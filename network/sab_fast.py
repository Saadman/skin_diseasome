#!/usr/bin/env python3
"""Fast reimplementation of the Menche et al. (2015) network-separation statistic.

Computes the same quantities as separation.py -- d_A, d_B, d_AB and
S_AB = d_AB - (d_A + d_B)/2 -- but replaces the all-pairs shortest-path
computation with a breadth-first search that stops at the first member of the
target set, which makes all fifteen disease pairs tractable in one pass.

Verified against separation.py and against the values in
networksim/Similarity_scores.xlsx; see README.md in this folder.

Usage:
  python3 sab_fast.py interactome.tsv PS gene_lists/PS.txt AD gene_lists/AD.txt ...
"""

import sys, itertools
from collections import deque
def load_graph(path):
    adj={}
    for line in open(path):
        if line.startswith('#') or not line.strip(): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<2: continue
        a,b=f[0].strip(),f[1].strip()
        adj.setdefault(a,set()); adj.setdefault(b,set())
        if a!=b: adj[a].add(b); adj[b].add(a)
    return adj
def read_set(p):
    return {l.split('\t')[0].strip().strip('"') for l in open(p) if l.strip() and not l.startswith('#')}
def min_dist(adj, src, targets, exclude_self):
    if not exclude_self and src in targets: return 0
    seen={src}; q=deque([(src,0)])
    while q:
        v,d=q.popleft()
        if d>=10: return None
        for w in adj[v]:
            if w in seen: continue
            seen.add(w)
            if w in targets: return d+1
            q.append((w,d+1))
    return None
def mean_min(adj, A, B, self_excl):
    vals=[min_dist(adj,a,B,self_excl) for a in A]
    vals=[v for v in vals if v is not None]
    return sum(vals)/len(vals) if vals else float('nan')
def sab(adj, A, B):
    A={g for g in A if g in adj}; B={g for g in B if g in adj}
    dAA=mean_min(adj,A,A,True); dBB=mean_min(adj,B,B,True)
    vals=[min_dist(adj,a,B,False) for a in A]+[min_dist(adj,b,A,False) for b in B]
    vals=[v for v in vals if v is not None]
    dAB=sum(vals)/len(vals)
    return dAB-(dAA+dBB)/2.0, dAB, dAA, dBB, len(A), len(B)
if __name__=='__main__':
    adj=load_graph(sys.argv[1])
    labels=sys.argv[2::2]; paths=sys.argv[3::2]
    sets={l:read_set(p) for l,p in zip(labels,paths)}
    print("pair,S_AB,d_AB,d_A,d_B,n_A,n_B")
    for a,b in itertools.combinations(labels,2):
        s,dab,daa,dbb,na,nb=sab(adj,sets[a],sets[b])
        print("%s-%s,%.3f,%.3f,%.3f,%.3f,%d,%d"%(a,b,s,dab,daa,dbb,na,nb))
