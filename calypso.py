#!/usr/bin/env python3
"""
Calypso: Standalone taxonomy pipeline using inlined EMU alignment logic + ATLAS consensus

Combines:
- Phase 1: Outlier detection (inlined ATLAS score_blast logic) using EMU-derived CIGAR parsing
- Phase 2: Co-occurrence graph + Louvain clustering
- Phase 3: Read-to-partition assignment
- Phase 4: LCA consensus on outliers + partitions
- Phase 5: Final merge of both LCA outputs (includes all queries)

Usage:
  calypso.py \
    --query_dir fasta_dir \
    --ref ref_16S.fasta \
    --tax_file tax_map.tsv \
    --outdir results \
    [--qc 0.9] [--pid 99] [--raiseto 2.7] [--type map-ont] [--threads 8]

Requirements:
  pip install numpy scipy biopython networkx python-louvain pysam
"""
import os, sys, argparse, subprocess, math
from itertools import combinations
import glob
import numpy as np
from scipy import special
from Bio import SeqIO
import networkx as nx
import community as community_louvain
import pysam

# Edgar et al.'s identity thresholds
THRESHOLDS = {
    'species': 99,
    'genus': 97,
    'family': 92,
    'order': 90,
    'class': 88,
    'phylum': 80,
}
RANKS_FULL = ['domain','phylum','class','order','family','genus','species']

# --- EMU get_align_stats + get_align_len (inlined) --------------------
def get_align_stats(aln):
    ops = aln.cigartuples or []
    ins = sum(length for code, length in ops if code == 1)
    dels = sum(length for code, length in ops if code == 2)
    soft = sum(length for code, length in ops if code == 4)
    mm = aln.get_tag('NM') - ins - dels if aln.has_tag('NM') else 0
    return ins, dels, soft, mm

def get_align_len(aln):
    return sum(length for code, length in (aln.cigartuples or [])
               if code in (0,1,2,7,8))
# --- ATLAS MSA scoring (inlined) --------------------------------------
def execute_msa(seqs, width, raiseto, gh, g):
    arr = {b: np.cumsum([[1 if x==b else 0 for x in row] for row in seqs], axis=0)
           for b in ['A','C','T','G']}
    tot = [sum(arr[b][-1][i] for b in arr) for i in range(width)]
    entropy = [ -sum((arr[b][-1][i]/tot[i]) * math.log(arr[b][-1][i]/tot[i],4)
                      for b in arr if tot[i]>0 and arr[b][-1][i]>0)
                if tot[i]>0 else 0
                for i in range(width)]
    base = sum((entropy[i]**raiseto)*(
                   sum(special.gammaln(arr[b][-1][i]+0.5) for b in arr)
                   - len(arr)*special.gammaln(0.5)
                   - special.gammaln(tot[i]+2)
               ) for i in range(width) if tot[i]>0)
    scores = [base]
    for k in range(1,len(seqs)):
        pre = {b: arr[b][-1] - arr[b][k-1] for b in arr}
        post = {b: arr[b][k-1] for b in arr}
        tp = sum(pre[b][i] for b in arr for i in range(width))
        to = sum(post[b][i] for b in arr for i in range(width))
        sp = sum((entropy[i]**raiseto)*(
                      sum(special.gammaln(pre[b][i]+0.5) for b in arr)
                      - len(arr)*special.gammaln(0.5)
                      - special.gammaln(tp+2)
                  ) for i in range(width) if tp>0)
        so = sum((entropy[i]**raiseto)*(
                      sum(special.gammaln(post[b][i]+0.5) for b in arr)
                      - len(arr)*special.gammaln(0.5)
                      - special.gammaln(to+2)
                  ) for i in range(width) if to>0)
        scores.append(sp+so-base)
    if len(scores)>=3 and scores[-2]>scores[-1] and scores[-2]>scores[-3] and scores[-2]>0:
        return len(scores)-2, scores
    return -1, scores

# --- select outliers --------------------------------------------------
def write_outliers(name, recs, fout, fsum, sub, raiseto, gh, g):
    """Select outliers, record best‐hit PID/ref, then run ATLAS MSA logic."""
    # extract all PIDs & pick the best
    pid_vals = [ float(info.split('\t')[2]) for *_, info in recs ]
    best_idx = pid_vals.index(max(pid_vals)) if recs else None
    best_pid = pid_vals[best_idx] if recs else None
    best_ref = recs[best_idx][1]     if recs else None

    ids = [rid for _, rid, _ in recs]

    # no hits
    if not recs:
        fsum.write(f"{name}\tNA\tNA\tNA\tNA\n")
        fout.write(f"{name}\tNA\n")
        sub.write("\n")
        return

    # few hits (<3) -> just list them
    if len(recs) < 3:
        cand = ids
        fsum.write(
            f"{name}\t{best_pid:.2f}\t{best_ref}"
            f"\t{';'.join(cand)}\tNA\n"
        )
        fout.write(f"{name}\t{';'.join(cand)}\n")
        sub.write("\n".join(info for _,_,info in recs) + "\n")
        return

    # ATLAS MSA‐based outlier detection
    seqs = [alist for alist,_,_ in recs]
    idx, scores = execute_msa(seqs, len(seqs[0]), raiseto, gh, g)

    if idx <= 0:
        # fallback to all
        cand = ids
        fsum.write(
            f"{name}\t{best_pid:.2f}\t{best_ref}"
            f"\t{';'.join(cand)}\tNA\n"
        )
        fout.write(f"{name}\t{';'.join(cand)}\n")
        sub.write("\n".join(info for _,_,info in recs) + "\n")
    else:
        cand = ids[:idx]
        fsum.write(
            f"{name}\t{best_pid:.2f}\t{best_ref}"
            f"\t{';'.join(cand)}\t{scores[idx]:.3f}\n"
        )
        fout.write(f"{name}\t{';'.join(cand)}\n")
        sub.write("\n".join(recs[i][2] for i in range(idx)) + "\n")

# --- Phase 1: outlier detection ---------------------------------------
def phase1(query_fasta, sam_file, outdir, qc, pid_thr, raiseto):

    seqs = {r.id: str(r.seq) for r in SeqIO.parse(query_fasta, 'fasta')}
    qlen = {k: len(v) for k, v in seqs.items()}
    os.makedirs(outdir, exist_ok=True)

    fout = open(f"{outdir}/outliers.txt", 'w')
    fsum = open(f"{outdir}/outlier_summary.txt", 'w')
    sub  = open(f"{outdir}/subset.blast",    'w')

    # <-- new header with best_pid & best_ref columns -->
    fsum.write("SequenceID\tbest_pid\tbest_ref\tcandidates\tatlas_score\n")
    
    maxh = 500
    gh = [special.gammaln(i + 0.5) for i in range(maxh + 10)]
    g  = [special.gammaln(i + 2)   for i in range(maxh + 10)]

    recs_map = {}
    sam = pysam.AlignmentFile(sam_file, 'r')
    for aln in sam.fetch():
        if aln.is_secondary or aln.is_supplementary or aln.reference_name is None:
            continue
        ins, dels, soft, mm = get_align_stats(aln)
        alen = get_align_len(aln)
        pid  = (alen - mm) / alen * 100
        cov  = alen / qlen.get(aln.query_name, alen)
        if cov < qc or pid < pid_thr:
            continue
        recs_map.setdefault(aln.query_name, []).append(
            (list(aln.query_sequence.replace('-', '')),
             aln.reference_name,
             f"{aln.query_name}\t{aln.reference_name}\t{pid:.2f}\t{cov:.3f}")
        )
    sam.close()

    for q, recs in recs_map.items():
        write_outliers(q, recs, fout, fsum, sub, raiseto, gh, g)

    fout.close()
    fsum.close()
    sub.close()
    return f"{outdir}/outliers.txt"

# --- Phase 2: graph & clusters ----------------------------------------
def make_edge_list(in_f, out_f):
    os.makedirs(os.path.dirname(out_f), exist_ok=True)
    G=nx.Graph()
    for ln in open(in_f):
        _,c=ln.strip().split("\t")
        hits=c.split(';') if c else []
        for a,b in combinations(hits,2): G.add_edge(a,b)
    nx.write_edgelist(G,out_f,data=False)
    return out_f

def make_partitions(edges, out_f):
    os.makedirs(os.path.dirname(out_f), exist_ok=True)
    G=nx.read_edgelist(edges)
    part=community_louvain.best_partition(G)
    with open(out_f,'w') as o:
        for n,c in part.items(): o.write(f"{n}\t{c}\n")
    return out_f

# --- Phase 3: read->partition -----------------------------------------
def read_partition_assignment(in_f, part_f, out_f):
    os.makedirs(os.path.dirname(out_f), exist_ok=True)
    pm={ln.split()[0]:ln.split()[1] for ln in open(part_f)}
    with open(out_f,'w') as o:
        for ln in open(in_f):
            q,c=ln.strip().split("\t")
            ps=sorted({pm.get(h) for h in c.split(';') if h in pm})
            o.write(f"{q}\t{';'.join(ps)}\n")
    return out_f

# --- Phase 4: LCA ------------------------------------------------------
def load_taxonomy(path):
    """
    Read lines of “tid \\t domain;phylum;class;order;family;genus;species”
    into a dict tid -> [domain,phylum,class,order,family,genus,species].
    """
    tax = {}
    with open(path) as f:
        for ln in f:
            tid, lin = ln.rstrip().split('\t', 1)
            ranks = lin.split(';')
            # pad to exactly 7 ranks
            if len(ranks) < 7:
                ranks += ['Unclassified'] * (7 - len(ranks))
            tax[tid] = ranks[:7]
    return tax

def assign_lca(in_f, tax, out_f):
    os.makedirs(os.path.dirname(out_f), exist_ok=True)
    with open(out_f,'w') as o:
        o.write('SequenceID\t' + ';'.join(RANKS_FULL) + '\n')
        for ln in open(in_f):
            line = ln.strip()
            # skip empty or malformed lines
            if not line or '\t' not in line:
                continue
            seq, ids = line.split('\t', 1)
            # build list of taxonomic arrays (or empty if 'NA')
            vals = [tax.get(i, ['Unassigned']*7) for i in ids.split(';')] if ids != 'NA' else []
            # compute deepest common prefix
            lca = []
            for tpl in zip(*vals):
                if len(set(tpl)) == 1:
                    lca.append(tpl[0])
                else:
                    break
            o.write(f"{seq}\t{';'.join(lca or ['Unassigned'])}\n")
    return out_f

# --- Phase 5: merge (includes all queries) ------------------------------
def generate_output(out_lca, part_lca, query_fasta, final_f, taxmap):
    import os
    from Bio import SeqIO

    # load best‐hit summary
    summary = {}
    summary_fp = os.path.join(
        os.path.dirname(out_lca), '..', 'phase1', 'outlier_summary.txt'
    )
    for ln in open(summary_fp):
        if ln.startswith('SequenceID'): continue
        seq, pid, ref, *_ = ln.strip().split('\t')
        summary[seq] = (pid, ref)

    # load the two LCA fallbacks (not used if thresholds override)
    ot = {l.split('\t')[0]:l.split('\t')[1]
          for l in open(out_lca)   if not l.startswith('SequenceID')}
    pt = {l.split('\t')[0]:l.split('\t')[1]
          for l in open(part_lca) if not l.startswith('SequenceID')}

    # order thresholds from most‐stringent (species) to least
    TH_ORDER = ['species','genus','family','order','class','phylum']

    # prepare output
    query_ids = [r.id for r in SeqIO.parse(query_fasta, 'fasta')]
    os.makedirs(os.path.dirname(final_f), exist_ok=True)
    with open(final_f, 'w') as o:
        o.write('SequenceID\tbest_pid\tbest_ref\t' + ';'.join(RANKS_FULL) + '\n')

        for seq in query_ids:
            pid, ref = summary.get(seq, ('NA','NA'))
            # figure deepest rank we pass
            try:
                pidf = float(pid)
                rank_idx = None
                for rank in TH_ORDER:
                    if pidf >= THRESHOLDS[rank]:
                        rank_idx = RANKS_FULL.index(rank)
                        break
            except ValueError:
                rank_idx = None

            if rank_idx is not None and ref != 'NA':
                tid = ref.split(':',1)[0]
                full_lin = taxmap.get(tid, ['Unclassified'] * 7)
                # keep 0..rank_idx, blank below
                lineage = full_lin[:rank_idx+1] + ['Unclassified']*(6-rank_idx)
            else:
                lineage = ['Unassigned'] * 7

            o.write(f"{seq}\t{pid}\t{ref}\t{';'.join(lineage)}\n")

    return final_f

# --- Main -------------------------------------------------------------
if __name__=='__main__':
    p = argparse.ArgumentParser()
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument('--query', help="Single query FASTA")
    grp.add_argument('--query_dir', help="Directory of query FASTA files")
    p.add_argument('--ref',      required=True)
    p.add_argument('--tax_file', required=True)
    p.add_argument('--outdir',   default='results')
    p.add_argument('--qc',       type=float, default=0.9)
    p.add_argument('--pid',      type=float, default=99)
    p.add_argument('--raiseto',  type=float, default=2.7)
    p.add_argument('--type',     dest='mm2type', default='map-ont',
                   help='Minimap2 preset: map-ont,map-pb,map-hifi,etc.')
    p.add_argument('--threads',  type=int,   default=8)
    args = p.parse_args()

    # build minimap2 index once
    idx = f"{args.ref}.mmi"
    if not os.path.exists(idx):
        subprocess.run(['minimap2','-d', idx, args.ref], check=True)

    # collect all FASTA inputs
    if args.query_dir:
        fasta_files = sorted(glob.glob(os.path.join(args.query_dir, '*')))
        fasta_files = [f for f in fasta_files if f.endswith(('.fa','.fasta','.fna'))]
    else:
        fasta_files = [args.query]

    # make sure outdir exists
    os.makedirs(args.outdir, exist_ok=True)

    # prepare combined master output
    combined_fp = os.path.join(args.outdir, 'combined_taxonomic_annotation.txt')
    with open(combined_fp, 'w') as combo:
        combo.write('SequenceID\tbest_pid\tbest_ref\t' + ';'.join(RANKS_FULL) + '\n')

        # load taxonomy map once
        taxmap = load_taxonomy(args.tax_file)

        for query in fasta_files:
            sample = os.path.splitext(os.path.basename(query))[0]
            out_base = os.path.join(args.outdir, sample)
            # make phase directories
            for phase in ['phase1','phase2','phase3','phase4','phase5']:
                os.makedirs(os.path.join(out_base, phase), exist_ok=True)

            # 1) align
            sam = os.path.join(out_base, 'align.sam')
            cmd = ['minimap2','-ax', args.mm2type, idx, query, '-t', str(args.threads)]
            with open(sam, 'w') as fh:
                subprocess.run(cmd, stdout=fh, check=True)

            # 2) phase 1 → 5
            out1      = phase1(        query, sam, out_base+'/phase1', args.qc, args.pid, args.raiseto)
            edges     = make_edge_list(out1,      out_base+'/phase2/edges.txt')
            parts     = make_partitions(edges,    out_base+'/phase3/partitions.txt')
            r2p       = read_partition_assignment(out1, parts, out_base+'/phase4/read2part.txt')
            out_lca   = assign_lca(     out1, taxmap, out_base+'/phase5/tax_outliers.txt')
            part_lca  = assign_lca(      r2p, taxmap, out_base+'/phase5/tax_partitions.txt')
            final_tf  = generate_output( out_lca, part_lca, query,
                                         out_base+'/taxonomic_annotation.txt',
                                         taxmap)

            # 3) append this sample’s results (skip its header)
            with open(final_tf) as sf:
                next(sf)
                for ln in sf:
                    combo.write(ln)

            print(f"[{sample}] done.", file=sys.stderr)

    print(f"\nCombined taxonomy at: {combined_fp}", file=sys.stderr)

