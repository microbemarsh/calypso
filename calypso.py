#!/usr/bin/env python3
"""
Self-contained taxonomy assignment pipeline combining:
- Phase 1: Outlier detection (ATLAS `score_blast` logic, inlined)
- Phase 2: Graph clustering (ATLAS `make_edge_list` + `make_partitions`)
- Phase 3: Read-to-partition assignment (ATLAS `read_partition_assignment`)
- Phase 4: Classification via Minimap2 + aln_stats + LCA consensus
- Phase 5: Final merge of outlier-based and partition-based LCA

Requirements:
  pip install numpy scipy biopython networkx python-louvain aln_stats minimap2

Usage:
  calypso.py \
    --query queries.fasta \
    --ref ref_16S.fasta \
    --tax_file tax_map.tsv \
    --outdir results \
    [--qc 0.9] [--pid 99] [--raiseto 2.7] [--threads 8]
"""
import os
import sys
import argparse
import subprocess
import math
from itertools import groupby, combinations

import numpy as np
from scipy import special
from Bio import SeqIO
import networkx as nx
import community as community_louvain

# Must have aln_stats library available
from aln_stats import parse_sam, compute_identity, compute_coverage

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

# --- ATLAS Phase 1: outlier detection (inlined score_blast) --------------

def fasta_iter(fasta_name):
    """Return dict of seq_name->sequence and length."""
    queries, qlen = {}, {}
    with open(fasta_name) as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0]=='>'))
        for header in faiter:
            hdr = next(header)[1:].strip().split()[0]
            seq = ''.join(s.strip() for s in next(faiter))
            queries[hdr] = seq
            qlen[hdr] = len(seq)
    return queries, qlen


def calc_entropy(counts, total):
    if total == 0: return 0.0
    probs = [c/total for c in counts]
    return -sum(p*math.log(p,4) for p in probs if p>0)


def calc_col_score(partial, total, gh, g):
    return sum(gh[p] for p in partial) - len(partial)*gh[0] - g[total]


def execute_msa(seqs, width, raiseto, gh, g):
    # seqs: list of equal-length lists of bases
    arr = {b: np.cumsum([[1 if x==b else 0 for x in row] for row in seqs], axis=0)
           for b in ['A','C','T','G']}
    tot = sum(arr[b][-1] for b in arr)
    entropy = [calc_entropy([arr[b][-1][i] for b in arr], tot[i]) for i in range(width)]
    base = sum((entropy[i]**raiseto)*calc_col_score([arr[b][-1][i] for b in arr], tot[i], gh, g)
               if tot[i]>0 else 0 for i in range(width))
    scores = [base]
    for k in range(1, len(seqs)):
        pre = {b: arr[b][-1] - arr[b][k-1] for b in arr}
        post = {b: arr[b][k-1] for b in arr}
        tp = sum(pre[b] for b in arr)
        to = sum(post[b] for b in arr)
        sp = sum((entropy[i]**raiseto)*calc_col_score([pre[b][i] for b in arr], pre[b][i] and tp, gh, g)
                 if tp>0 else 0 for i in range(width))
        so = sum((entropy[i]**raiseto)*calc_col_score([post[b][i] for b in arr], post[b][i] and to, gh, g)
                 if to>0 else 0 for i in range(width))
        scores.append(sp+so-base)
    if len(scores)>=3 and scores[-2]>max(scores[-1],scores[-3]) and scores[-2]>0:
        return len(scores)-2, scores
    return -1, scores


def write_outliers(name, recs, out_out, out_sum, subset):
    # recs: list of (aligned_query_list, ref_id)
    pid_names = [rid for _,rid in recs]
    seq_list = [aln for aln,_ in recs]
    idx, scores = execute_msa(seq_list, len(seq_list[0]), raiseto, gh, g) if recs else (-1,[])
    if idx<=0:
        out_sum.write(f"{name}\tNA\t{';'.join(pid_names)}\tNA\n")
        subset.write("\n".join(map(lambda x:x[2], recs))+"\n")
        return
    candidates = pid_names[:idx]
    out_sum.write(f"{name}\t{';'.join(candidates)}\t{scores[idx]:.3f}\n")
    out_out.write(f"{name}\t{';'.join(candidates)}\n")
    subset.write("\n".join([recs[i][2] for i in range(idx)])+"\n")


def phase1_minimap(query_fasta, sam_file, outdir, qc, pid_thr, raiseto):
    """Runs ATLAS outlier detection using precomputed SAM + aln_stats"""
    queries, qlen = fasta_iter(query_fasta)
    os.makedirs(outdir, exist_ok=True)
    out_out = open(os.path.join(outdir,'outliers.txt'),'w')
    out_sum = open(os.path.join(outdir,'outlier_summary.txt'),'w')
    subset = open(os.path.join(outdir,'subset.blast'),'w')
    # precompute gamma
    maxh=500; gh=[special.gammaln(i+0.5) for i in range(maxh+10)]; g=[special.gammaln(i+2) for i in range(maxh+10)]
    # parse alignments
    recs_map={}
    for rec in parse_sam(sam_file):
        pid=compute_identity(rec); cov=compute_coverage(rec)
        if cov<qc or pid<pid_thr: continue
        # aligned query as list
        aln_q = list(rec.aligned_query)
        recs_map.setdefault(rec.qname,[]).append((aln_q, rec.rname, f"{rec.qname}\t{rec.rname}	{pid:.2f}\t{cov:.3f}"))
    # write per-query
    for qname,recs in recs_map.items():
        write_outliers(qname, recs, out_out, out_sum, subset)
    out_out.close(); out_sum.close(); subset.close()
    return os.path.join(outdir,'outliers.txt')

# --- ATLAS Phase 2: co-occurrence edges ------------------------------
def make_edge_list(out_f, edges_f):
    G=nx.Graph()
    for ln in open(out_f):
        qid, cand=ln.strip().split("\t")
        hits=cand.split(';') if cand else []
        for a,b in combinations(hits,2): G.add_edge(a,b)
    nx.write_edgelist(G, edges_f, data=False)
    return edges_f

# --- ATLAS Phase 3: Louvain partitions ------------------------------
def make_partitions(edges_f, part_f):
    G=nx.read_edgelist(edges_f)
    part=community_louvain.best_partition(G)
    with open(part_f,'w') as out:
        for node,com in part.items(): out.write(f"{node}\t{com}\n")
    return part_f

# --- ATLAS Phase 4: query->partition assignment ---------------------
def read_partition_assignment(out_f, part_f, r2p_f):
    pmap={ln.split()[0]:ln.split()[1] for ln in open(part_f)}
    with open(r2p_f,'w') as out:
        for ln in open(out_f):
            qid,c=ln.strip().split("\t")
            parts=sorted({pmap.get(h) for h in c.split(';') if h in pmap})
            out.write(f"{qid}\t{';'.join(parts)}\n")
    return r2p_f

# --- Phase 5: LCA taxonomy -----------------------------------------
def load_taxonomy(path):
    idx=[0,2,5,8,10,12,13]
    tax={}
    for ln in open(path):
        tid,lin=ln.strip().split('\t',1)
        lv=lin.split(';'); lv+=['Unclassified']*(max(idx)+1-len(lv))
        tax[tid]=[lv[i] for i in idx]
    return tax


def assign_lca(in_f, tax, out_f):
    with open(out_f,'w') as out:
        out.write('SequenceID\t'+";".join(RANKS_FULL)+"\n")
        for ln in open(in_f):
            seq,ids=ln.strip().split('\t')
            lineages=[tax.get(i, ['Unassigned']*7) for i in ids.split(';')]
            lca=[]
            for vals in zip(*lineages):
                if len(set(vals))==1: lca.append(vals[0])
                else: break
            out.write(seq+'\t'+';'.join(lca or ['Unassigned'])+'\n')
    return out_f

# --- ATLAS Phase 6: merge -------------------------------------------
def generate_output(out_t, part_t, final_f):
    ot={ln.split()[0]:ln.split()[1] for ln in open(out_t) if not ln.startswith('SequenceID')}
    pt={ln.split()[0]:ln.split()[1] for ln in open(part_t) if not ln.startswith('SequenceID')}
    with open(final_f,'w') as out:
        out.write('SequenceID\t'+";".join(RANKS_FULL)+"\n")
        for seq in sorted(set(ot)|set(pt)):
            tax=ot.get(seq) or pt.get(seq) or 'Unassigned'
            out.write(seq+'\t'+tax+'\n')
    return final_f

# --- Main -----------------------------------------------------------
if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--query', required=True)
    p.add_argument('--ref', required=True)
    p.add_argument('--tax_file', required=True)
    p.add_argument('--outdir', default='results')
    p.add_argument('--qc', type=float, default=0.9)
    p.add_argument('--pid', type=float, default=99)
    p.add_argument('--raiseto', type=float, default=2.7)
    p.add_argument('--threads', type=int, default=8)
    args=p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    # 1) build index & align
    sam=os.path.join(args.outdir,'align.sam')
    idxfile=args.ref+'.mmi'
    if not os.path.exists(idxfile):
        subprocess.run(['minimap2','-d',idxfile,args.ref], check=True)
    with open(sam,'w') as outf:
        subprocess.run(['minimap2','-ax','sr',idxfile,args.query,'-t',str(args.threads)], stdout=outf, check=True)
    # 2) Phase 1 outliers
    out1=phase1_minimap(args.query, sam, os.path.join(args.outdir,'phase1'), args.qc, args.pid, args.raiseto)
    # 3) Phase 2-4
    edges=make_edge_list(out1, os.path.join(args.outdir,'phase2','edges.txt'))
    parts=make_partitions(edges, os.path.join(args.outdir,'phase3','partitions.txt'))
    r2p=read_partition_assignment(out1, parts, os.path.join(args.outdir,'phase4','read2part.txt'))
    # 4) Phases 5-6 classification
    taxmap=load_taxonomy(args.tax_file)
    out_lca=assign_lca(out1, taxmap, os.path.join(args.outdir,'phase5','tax_outliers.txt'))
    part_lca=assign_lca(r2p, taxmap, os.path.join(args.outdir,'phase5','tax_partitions.txt'))
    final=generate_output(out_lca, part_lca, os.path.join(args.outdir,'taxonomic_annotation.txt'))
    print(f"Done. Final taxonomy: {final}", file=sys.stderr)

