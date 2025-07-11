#!/usr/bin/env python3
"""
Calypso: Standalone taxonomy pipeline using Parasail Smith–Waterman (BLASTn defaults) + ATLAS consensus

Combines:
- Phase 1: Outlier detection using inlined ATLAS score_blast logic
- Phase 2: Co-occurrence graph + Louvain clustering
- Phase 3: Read-to-partition assignment
- Phase 4: LCA consensus on outliers + partitions
- Phase 5: Final merge of both LCA outputs (includes all queries)

Usage:
  calypso.py \
    --query FILE.fasta    # single query file
    --query_dir DIR       # or directory of query FASTAs
    --ref ref_16S.fasta \
    --tax_file tax_map.tsv \
    --outdir results \
    [--qc 0.9] [--pid 99] [--raiseto 2.7]

Requirements:
  pip install numpy scipy parasail networkx python-louvain
"""
import os, sys, argparse, math, glob
from itertools import combinations
import numpy as np
from scipy import special
import parasail
from Bio import SeqIO
import networkx as nx
import community as community_louvain
from multiprocessing import Pool, cpu_count

# Edgar et al.'s identity thresholds
THRESHOLDS = {'species':99,'genus':97,'family':92,'order':90,'class':88,'phylum':80}
RANKS_FULL = ['domain','phylum','class','order','family','genus','species']

# --- ATLAS MSA scoring --------------------------------------
def execute_msa(seqs, width, raiseto, gh, g):
    arr = {b: np.cumsum([[1 if x==b else 0 for x in row] for row in seqs], axis=0)
           for b in ['A','C','T','G']}
    tot = [sum(arr[b][-1][i] for b in arr) for i in range(width)]
    entropy = [ -sum((arr[b][-1][i]/tot[i]) * math.log(arr[b][-1][i]/tot[i],4)
                      for b in arr if tot[i]>0 and arr[b][-1][i]>0)
                if tot[i]>0 else 0 for i in range(width)]
    base = sum((entropy[i]**raiseto)*(
                   sum(special.gammaln(arr[b][-1][i]+0.5) for b in arr)
                   - len(arr)*special.gammaln(0.5)
                   - special.gammaln(tot[i]+2))
                 for i in range(width) if tot[i]>0)
    scores = [base]
    for k in range(1, len(seqs)):
        pre  = {b: arr[b][-1] - arr[b][k-1] for b in arr}
        post = {b: arr[b][k-1] for b in arr}
        tp = sum(pre[b][i] for b in arr for i in range(width))
        to = sum(post[b][i] for b in arr for i in range(width))
        sp = sum((entropy[i]**raiseto)*(
                   sum(special.gammaln(pre[b][i]+0.5) for b in arr)
                   - len(arr)*special.gammaln(0.5)
                   - special.gammaln(tp+2))
                for i in range(width) if tp>0)
        so = sum((entropy[i]**raiseto)*(
                   sum(special.gammaln(post[b][i]+0.5) for b in arr)
                   - len(arr)*special.gammaln(0.5)
                   - special.gammaln(to+2))
                for i in range(width) if to>0)
        scores.append(sp + so - base)
    if len(scores) >= 3 and scores[-2] > scores[-1] and scores[-2] > scores[-3] and scores[-2] > 0:
        return len(scores) - 2, scores
    return -1, scores

# --- Select outliers and write summaries ------------------------------
def write_outliers(name, recs, fout, fsum, sub, raiseto, gh, g):
    pid_vals = [float(info.split('\t')[2]) for *_, info in recs]
    if pid_vals:
        min_pid, max_pid = min(pid_vals), max(pid_vals)
        mean_pid = sum(pid_vals) / len(pid_vals)
        median_pid = float(np.median(pid_vals))
        std_pid = float(np.std(pid_vals))
        best_idx = pid_vals.index(max_pid)
        best_ref = recs[best_idx][1]
    else:
        min_pid = max_pid = mean_pid = median_pid = std_pid = None
        best_ref = None

    if not recs:
        fsum.write(f"{name}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
        fout.write(f"{name}\tNA\n")
        sub.write("\n")
        return

    ids = [rid for _, rid, _ in recs]
    if len(recs) < 3:
        fsum.write(
            f"{name}\t{min_pid:.3f}\t{max_pid:.3f}\t"
            f"{mean_pid:.3f}\t{median_pid:.3f}\t{std_pid:.3f}\t"
            f"{best_ref}\t{';'.join(ids)}\tNA\n"
        )
        fout.write(f"{name}\t{';'.join(ids)}\n")
        sub.write("\n".join(info for _,_,info in recs) + "\n")
        return

    seqs = [alist for alist, _, _ in recs]
    width = max(len(row) for row in seqs)
    padded = [row + ['-'] * (width - len(row)) for row in seqs]
    idx, scores = execute_msa(padded, width, raiseto, gh, g)
    if idx <= 0:
        cand, atlas_score = ids, "NA"
    else:
        cand, atlas_score = ids[:idx], f"{scores[idx]:.3f}"

    fsum.write(
        f"{name}\t{min_pid:.3f}\t{max_pid:.3f}\t"
        f"{mean_pid:.3f}\t{median_pid:.3f}\t{std_pid:.3f}\t"
        f"{best_ref}\t{';'.join(cand)}\t{atlas_score}\n"
    )
    fout.write(f"{name}\t{';'.join(cand)}\n")
    sub.write("\n".join(recs[i][2] for i in range(len(cand))) + "\n")

# --- Parallel alignment helper using Parasail ----------------------------
def _align_query_parasail(args):
    qid, qseq, refs, qc, pid_thr = args
    mat = parasail.matrix_create("ACGT", 1, -3)
    recs = []
    for rid, rseq in refs.items():
        result = parasail.sw_trace_striped_16(qseq, rseq, 5, 2, mat)
        aligned_q = result.traceback.query
        aligned_r = result.traceback.ref
        matches = sum(1 for a, b in zip(aligned_q, aligned_r) if a == b and a != '-')
        aln_len = sum(1 for a, b in zip(aligned_q, aligned_r) if a != '-' and b != '-')
        if aln_len == 0:
            continue
        pid = matches / aln_len * 100
        cov = aln_len / len(qseq)
        if cov < qc or pid < pid_thr:
            continue
        info = f"{qid}\t{rid}\t{pid:.2f}\t{cov:.3f}"
        recs.append((list(aligned_r), rid, info))
    return qid, recs

# --- Phase 1: Outlier detection (parallel, Parasail) --------------------
def phase1(query_fasta, outdir, qc, pid_thr, raiseto, ref_fasta):
    os.makedirs(outdir, exist_ok=True)
    fout = open(f"{outdir}/outliers.txt", 'w')
    fsum = open(f"{outdir}/outlier_summary.txt", 'w')
    sub = open(f"{outdir}/subset.blast", 'w')
    fsum.write("SequenceID\tmin_pid\tmax_pid\tmean_pid\tmedian_pid\tstd_pid\tbest_ref\tcandidates\tatlas_score\n")
    gh = [special.gammaln(i + 0.5) for i in range(510)]
    g  = [special.gammaln(i + 2)   for i in range(510)]
    refs = {r.id: str(r.seq) for r in SeqIO.parse(ref_fasta, 'fasta')}
    queries = {r.id: str(r.seq) for r in SeqIO.parse(query_fasta, 'fasta')}
    jobs = [(qid, qseq, refs, qc, pid_thr) for qid, qseq in queries.items()]
    with Pool(processes=cpu_count()) as pool:
        for qid, recs in pool.imap_unordered(_align_query_parasail, jobs):
            write_outliers(qid, recs, fout, fsum, sub, raiseto, gh, g)
    fout.close(); fsum.close(); sub.close()
    return f"{outdir}/outliers.txt"

# --- Co-occurrence graph -----------------------------------------------
def make_edge_list(in_f, out_f):
    os.makedirs(os.path.dirname(out_f), exist_ok=True)
    G = nx.Graph()
    for ln in open(in_f):
        _, c = ln.strip().split("\t")
        hits = c.split(';') if c else []
        for a, b in combinations(hits, 2):
            G.add_edge(a, b)
    nx.write_edgelist(G, out_f, data=False)
    return out_f

# --- Louvain partitions ------------------------------------------------
def make_partitions(edges, out_f):
    os.makedirs(os.path.dirname(out_f), exist_ok=True)
    G = nx.read_edgelist(edges)
    part = community_louvain.best_partition(G)
    with open(out_f, 'w') as o:
        for n, c in part.items():
            o.write(f"{n}\t{c}\n")
    return out_f

# --- Read -> partition assignment --------------------------------------
def read_partition_assignment(in_f, part_f, out_f):
    os.makedirs(os.path.dirname(out_f), exist_ok=True)
    pm = {ln.split()[0]: ln.split()[1] for ln in open(part_f)}
    with open(out_f, 'w') as o:
        for ln in open(in_f):
            line = ln.strip()
            if not line or '\t' not in line:
                continue
            q, c = line.split("\t", 1)
            ps = sorted({pm.get(h) for h in c.split(';') if h in pm})
            o.write(f"{q}\t{';'.join(ps)}\n")
    return out_f

# --- LCA consensus -----------------------------------------------------
def load_taxonomy(path):
    tax = {}
    for ln in open(path):
        tid, lin = ln.rstrip().split('\t', 1)
        ranks = lin.split(';')
        if len(ranks) < 7:
            ranks += ['Unclassified'] * (7 - len(ranks))
        tax[tid] = ranks[:7]
    return tax

def assign_lca(in_f, tax, out_f):
    os.makedirs(os.path.dirname(out_f), exist_ok=True)
    with open(out_f, 'w') as o:
        o.write('SequenceID\t' + ';'.join(RANKS_FULL) + '\n')
        for ln in open(in_f):
            line = ln.strip()
            if not line or '\t' not in line:
                continue
            seq, ids = line.split('\t', 1)
            vals = [tax.get(i, ['Unassigned'] * 7) for i in ids.split(';')] if ids != 'NA' else []
            lca = []
            for tpl in zip(*vals):
                if len(set(tpl)) == 1:
                    lca.append(tpl[0])
                else:
                    break
            o.write(f"{seq}\t{';'.join(lca or ['Unassigned'])}\n")
    return out_f

# --- Final merge and annotation ----------------------------------------
def generate_output(out_lca, part_lca, query_fasta, final_f, taxmap):
    summary = {}
    summary_fp = os.path.join(os.path.dirname(out_lca), '..', 'phase1', 'outlier_summary.txt')
    for ln in open(summary_fp):
        if ln.startswith('SequenceID'): continue
        cols = ln.rstrip().split('\t')
        seq, median_pid, best_ref, cand = cols[0], cols[4], cols[6], cols[7]
        n_hits = 0 if cand == 'NA' else len(cand.split(';'))
        summary[seq] = (median_pid, best_ref, n_hits)
    TH_ORDER = ['species','genus','family','order','class','phylum']
    query_ids = [r.id for r in SeqIO.parse(query_fasta, 'fasta')]
    os.makedirs(os.path.dirname(final_f), exist_ok=True)
    with open(final_f, 'w') as o:
        o.write('SequenceID\tmedian_pid\tn_hits\tbest_ref\t' + ';'.join(RANKS_FULL) + '\n')
        for seq in query_ids:
            pid, ref, n_hits = summary.get(seq, ('NA','NA',0))
            rank_idx = None
            try:
                pidf = float(pid)
                for rank in TH_ORDER:
                    if pidf >= THRESHOLDS[rank]:
                        rank_idx = RANKS_FULL.index(rank)
                        break
            except ValueError:
                pass
            if rank_idx is not None and ref != 'NA':
                tid = ref.split(':', 1)[0]
                full = taxmap.get(tid, ['Unclassified'] * 7)
                lineage = full[:rank_idx+1] + ['Unclassified'] * (6 - rank_idx)
            else:
                lineage = ['Unassigned'] * 7
            o.write(f"{seq}\t{pid}\t{n_hits}\t{ref}\t{';'.join(lineage)}\n")
    return final_f

# --- Main ----------------------------------------------
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument('--query', help="Single query FASTA")
    grp.add_argument('--query_dir', help="Directory of query FASTA files")
    p.add_argument('--ref', required=True)
    p.add_argument('--tax_file', required=True)
    p.add_argument('--outdir', default='results')
    p.add_argument('--qc', type=float, default=0.9)
    p.add_argument('--pid', type=float, default=99)
    p.add_argument('--raiseto', type=float, default=2.7)
    args = p.parse_args()

    # collect FASTA inputs
    if args.query:
        fasta_files = [args.query]
    else:
        fasta_files = sorted(glob.glob(os.path.join(args.query_dir, '*')))
        fasta_files = [f for f in fasta_files if f.endswith(('.fa','.fasta','.fna'))]

    # load taxonomy map
    taxmap = load_taxonomy(args.tax_file)

    combined_fp = os.path.join(args.outdir, 'combined_taxonomic_annotation.txt')
    os.makedirs(args.outdir, exist_ok=True)
    with open(combined_fp, 'w') as combo:
        combo.write('SequenceID\tmedian_pid\tn_hits\tbest_ref\t' + ';'.join(RANKS_FULL) + '\n')
        for query in fasta_files:
            sample = os.path.splitext(os.path.basename(query))[0]
            base = os.path.join(args.outdir, sample)
            p1 = phase1(query,         base + '/phase1', args.qc, args.pid, args.raiseto, args.ref)
            e  = make_edge_list(p1,    base + '/phase2/edges.txt')
            pt = make_partitions(e,    base + '/phase3/partitions.txt')
            rp = read_partition_assignment(p1, pt,      base + '/phase4/read2part.txt')
            oL = assign_lca(p1,        taxmap, base + '/phase5/tax_outliers.txt')
            pL = assign_lca(rp,        taxmap, base + '/phase5/tax_partitions.txt')
            final = generate_output(oL, pL, query,    base + '/phase5/final_annotation.txt', taxmap)
            with open(final) as sf:
                next(sf)
                for ln in sf:
                    combo.write(ln)
            print(f"[{sample}] done.", file=sys.stderr)
    print(f"Combined taxonomy at: {combined_fp}", file=sys.stderr)