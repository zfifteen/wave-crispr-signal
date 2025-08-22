#!/usr/bin/env python3
"""
Plotting script for bin-resonance CRISPR efficiency analysis.
Creates visualizations of phase-coherence vs efficiency across GC quartiles.
"""

import csv, math, argparse
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

def read_results(path):
    rows=[] 
    with open(path) as f:
        for d in csv.DictReader(f): rows.append(d)
    return rows

def read_data(path):
    rows=[]
    with open(path) as f:
        for d in csv.DictReader(f):
            s=(d.get("sequence","") or "").strip().upper()
            try: e=float(d.get("efficiency",""))
            except: continue
            if s: rows.append((s,e))
    return rows

def pc(seq):
    L=len(seq); mp={"A":1+0j,"T":-1+0j,"C":0+1j,"G":0-1j}
    s=0+0j
    for n,b in enumerate(seq):
        z=mp.get(b,0j); ang=2*math.pi*n/L
        s+= z*complex(math.cos(ang), math.sin(ang))
    return abs(s)/L

def gc(seq): return (seq.count("G")+seq.count("C"))/len(seq)

def quartiles(vals):
    v=sorted(vals); m=len(v)-1
    return v[int(.25*m)], v[int(.50*m)], v[int(.75*m)]

def main(data_csv, results_csv):
    data=read_data(data_csv)
    gcs=[gc(s) for s,_ in data]; q1,q2,q3=quartiles(gcs)
    bins=defaultdict(list)
    for (s,y),g in zip(data,gcs):
        b="Q1" if g<=q1 else "Q2" if g<=q2 else "Q3" if g<=q3 else "Q4"
        bins[b].append((pc(s), y, g))
    res_map={d["bin"]:d for d in read_results(results_csv)}
    
    # Fig 1: scatter + OLS line per quartile
    plt.figure(figsize=(10,6))
    for i,b in enumerate(["Q1","Q2","Q3","Q4"],1):
        xs=[x for x,_,_ in bins[b]]; ys=[y for _,y,_ in bins[b]]
        plt.subplot(2,2,i); plt.scatter(xs,ys,s=10,alpha=.6)
        if xs:
            mx,my=sum(xs)/len(xs),sum(ys)/len(ys)
            dx=[x-mx for x in xs]; dy=[y-my for y in ys]
            den=sum(d*d for d in dx) or 1e-9
            beta=sum(a*b for a,b in zip(dx,dy))/den; alpha=my-beta*mx
            x0,x1=min(xs),max(xs); plt.plot([x0,x1],[alpha+beta*x0,alpha+beta*x1])
        r=float(res_map.get(b,{}).get("r","nan"))
        p=res_map.get(b,{}).get("p_perm","nan")
        q=str(res_map.get(b,{}).get("q_pass",""))
        plt.title(f"{b}: r={r:.3f}, p={p}, FDR={q}")
        plt.xlabel("Phase coherence"); plt.ylabel("Efficiency")
    plt.tight_layout(); plt.savefig("fig1_scatter_quartiles.png",dpi=150)
    
    # Fig 2: r with 95% CI (error bars)
    plt.figure(figsize=(6,4))
    bs=["Q1","Q2","Q3","Q4"]
    r=[float(res_map.get(b,{}).get("r","nan")) for b in bs]
    lo=[float(res_map.get(b,{}).get("ci_low","nan")) for b in bs]
    hi=[float(res_map.get(b,{}).get("ci_high","nan")) for b in bs]
    x=np.arange(4); y=np.array(r); err=[y-np.array(lo), np.array(hi)-y]
    plt.errorbar(x,y,yerr=err,fmt='o'); plt.xticks(x,bs); plt.axhline(0,ls="--")
    plt.ylabel("Pearson r"); plt.title("Bootstrap 95% CI by GC quartile")
    plt.tight_layout(); plt.savefig("fig2_bootstrap_ci.png",dpi=150)
    
    # Fig 3: permutation p vs BH threshold (bars)
    plt.figure(figsize=(6,4))
    ps=[float(res_map.get(b,{}).get("p_perm","nan")) for b in bs]
    bh=[float(res_map.get("Q2",{}).get("bh_cutoff","nan"))]*4  # same cutoff for all bins
    x=np.arange(4); plt.bar(x-0.2,ps,width=0.4,label="p_perm")
    plt.bar(x+0.2,bh,width=0.4,label="BH cutoff"); plt.axhline(0.05,ls="--")
    plt.xticks(x,bs); plt.ylabel("p-value"); plt.title("Permutation p vs BH cutoff")
    plt.legend(); plt.tight_layout(); plt.savefig("fig3_permutation_bh.png",dpi=150)
    
    # Fig 4: 3D scatter (GC, PC, efficiency)
    fig=plt.figure(figsize=(6,5)); ax=fig.add_subplot(111,projection='3d')
    for b in ["Q1","Q2","Q3","Q4"]:
        xs=[x for x,_,g in bins[b]]; ys=[y for _,y,_ in bins[b]]; zs=[g for _,_,g in bins[b]]
        ax.scatter(zs,xs,ys,s=8,alpha=.6)
    ax.set_xlabel("GC fraction"); ax.set_ylabel("Phase coherence"); ax.set_zlabel("Efficiency")
    ax.set_title("GC × PC × Efficiency")
    plt.tight_layout(); plt.savefig("fig4_3d_scatter.png",dpi=150)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--data",required=True); ap.add_argument("--results",required=True)
    a=ap.parse_args(); main(a.data,a.results)