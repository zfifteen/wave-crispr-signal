#!/usr/bin/env python3
"""
ORCS v1.1.17 Physical Z-Metrics Testing Script

Tests four physical Z-metrics (base-pair opening, base-stacking dissociation, 
helical-twist fluctuation, denaturation/melting proxy) on real human CRISPR 
screen outcomes from BioGRID-ORCS 1.1.17.

Requirements: Python 3.12, numpy, scipy, biopython, mpmath (dps=50), pandas
Inputs:
  ORCS_DIR = "data/BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens"
  FASTA    = "data/hg38.fa"          # or guide FASTA (validated human-only)
  GTF      = "data/gencode.v44.annotation.gtf"  # if extracting CDS/TSS windows

Usage:
  python scripts/test_orcs_v1_1_17.py
"""

import os
import glob
import gzip
import json
import math
import re
import numpy as np
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import mpmath as mp
from scipy import stats
from pathlib import Path
import sys

# Add project root to path for imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from z_framework import ZFrameworkCalculator
from applications.crispr_physical_z_metrics import (
    PhysicalZMetricsCalculator,
    validate_human_dna_sequence,
    DNAValidationError
)

# Configure high precision
mp.mp.dps = 50
C_E2 = mp.e**2  # c = e^2

# === Safety guards ===
def validate_seq(s):
    """Validate DNA sequence contains only ACGTN nucleotides."""
    if not re.fullmatch(r"[ACGTN]+", s):
        raise ValueError("Invalid nucleotides (non-ACGTN) in sequence")
    return s

def compute_z(A, B, c=C_E2):
    """Compute Z-metric with causality check."""
    B = mp.fabs(B)
    if B >= c:
        raise ValueError(f"Causality violation: |B| >= c (B={B})")
    return A * (B / c)

# Golden ratio and geodesic resolution
phi = (1 + mp.sqrt(5)) / 2
KSTAR = mp.mpf("0.3")

def theta_prime(n, k=KSTAR, phi=phi):
    """Geodesic resolution function θ′(n,k) = φ·((n mod φ)/φ)^k."""
    n = mp.mpf(n)
    return phi * ((n % phi) / phi) ** k

# === Physical per-base tables (deterministic; replace with your calibrated set) ===
HBONDS = {'A':2,'C':3,'G':3,'T':2,'N':2}
STACK  = {'A':-1.0,'C':-1.5,'G':-2.0,'T':-0.6,'N':-1.0}
TWIST  = {'A':34.5,'C':34.0,'G':36.0,'T':35.5,'N':35.0}
MELT   = {'A':0.03,'C':0.05,'G':0.05,'T':0.03,'N':0.04}

PROP_TABLES = {
    'opening': HBONDS,
    'stack'  : STACK,
    'twist'  : TWIST,
    'melt'   : MELT,
}

def z_of_seq(seq, table):
    """Calculate Z-metric for sequence using property table."""
    # A = mean property; B = std of adjacent differences
    props = [table[b] for b in seq]
    A = mp.mpf(np.mean(props))
    diffs = np.abs(np.diff(props))
    B = mp.mpf(np.std(diffs) if len(diffs) else 0.001)
    return compute_z(A, B)

def enumerate_spcas9_guides(dna):
    """Return list of (20nt, pos_index) for + and - strands (pos_index is 1-based)."""
    dna = str(dna).upper()
    validate_seq(dna)
    guides = []
    # + strand: N20-NGG
    for i in range(len(dna)-23):
        if dna[i+20:i+23].endswith("GG"):
            guides.append((dna[i:i+20], i+1))
    # - strand: CCN-N20  (revcomp)
    comp = dna.translate(str.maketrans("ACGT","TGCA"))[::-1]
    for i in range(len(comp)-23):
        if comp[i+20:i+23].endswith("GG"):
            guides.append((comp[i:i+20], i+1))
    return guides

def aggregate_gene(guides, metric_key):
    """Aggregate Z-metric for gene across all guides."""
    table = PROP_TABLES[metric_key]
    z_vals = []
    for g, pos in guides:
        z = z_of_seq(g, table)
        wt = theta_prime(pos)
        z_vals.append(float(z * wt))
    if not z_vals:
        return None
    z_arr = np.array(z_vals, dtype=float)
    return {
        "mean": float(np.mean(z_arr)),
        "median": float(np.median(z_arr)),
        "top_decile_mean": float(np.mean(np.sort(z_arr)[int(0.9*len(z_arr)):])),
        "var": float(np.var(z_arr))
    }

def bootstrap_r(x, y, n=1000):
    """Bootstrap confidence intervals for Pearson correlation."""
    rng = np.random.default_rng(42)
    idx = np.arange(len(x))
    corrs = []
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        if np.std(x[s])==0 or np.std(y[s])==0:
            corrs.append(0.0)
        else:
            corrs.append(stats.pearsonr(x[s], y[s])[0])
    lo, hi = np.percentile(corrs, [2.5, 97.5])
    return float(lo), float(hi)

def cohen_d(x, y):
    """Cohen's d effect size between two groups."""
    # x,y are arrays for Hit=YES vs NO on Z aggregate
    nx, ny = len(x), len(y)
    sx2, sy2 = np.var(x, ddof=1), np.var(y, ddof=1)
    sp = math.sqrt(((nx-1)*sx2 + (ny-1)*sy2) / (nx+ny-2)) if (nx+ny-2)>0 else 0.0
    return float((np.mean(x)-np.mean(y)) / sp) if sp>0 else 0.0

def load_orcs_index(orcs_dir):
    """Load ORCS index files."""
    idx_files = glob.glob(os.path.join(orcs_dir, "*index.tab.txt"))
    frames = []
    for f in idx_files:
        frames.append(pd.read_csv(f, sep="\t", comment="#", header=None, dtype=str))
    if not frames:
        raise FileNotFoundError("No *.index.tab.txt files found")
    idx = pd.concat(frames, ignore_index=True)
    idx.columns = list(range(idx.shape[1]))  # keep numeric; we use positions below
    return idx

def select_screens(idx):
    """Select qualifying screens from ORCS index."""
    # Column positions based on actual index format (0-based):
    SCREEN_ID = 0; ORG_OFF = 34; ENZ = 23; LIBTYPE = 21; SCTYPE = 14; THR = 12; COUNT = 27; SCORE1_TYPE = 28
    
    # Convert columns to string and handle NaN values
    org_col = idx[ORG_OFF].fillna("").astype(str)
    enz_col = idx[ENZ].fillna("").astype(str)
    lib_col = idx[LIBTYPE].fillna("").astype(str)
    thr_col = idx[THR].fillna("").astype(str)
    cnt_col = pd.to_numeric(idx[COUNT], errors='coerce').fillna(0)
    
    idx2 = idx[
        (org_col.str.upper().str.contains("SAPIENS")) &
        (enz_col.str.upper()=="CAS9") &
        (lib_col.str.upper().isin(["CRISPRN","CRISPRKO"])) &
        (thr_col.str.contains("High", case=False, na=False)) &
        (cnt_col >= 1)
    ]
    
    if len(idx2) == 0:
        print("Debug: No screens found. Checking constraints...")
        print(f"Organisms: {org_col.unique()}")
        print(f"Enzymes: {enz_col.unique()}")
        print(f"Library types: {lib_col.unique()}")
        print(f"Throughput: {thr_col.unique()}")
    
    return idx2.assign(score1_type=idx2[SCORE1_TYPE], screen_id=idx2[SCREEN_ID]).loc[:, ["screen_id","score1_type"]]

def load_screen_rows(orcs_dir, screen_id):
    """Load screen data for specific screen ID."""
    # Each screen has a matching *.screen.tab.txt with first col = Screen ID
    rows = []
    for f in glob.glob(os.path.join(orcs_dir, "*.screen.tab.txt")):
        df = pd.read_csv(f, sep="\t", comment="#", header=None, dtype=str)
        df = df[df[1-1]==screen_id]
        if not df.empty:
            rows.append(df)
    if not rows:
        return None
    screen = pd.concat(rows, ignore_index=True)
    # Keep Official Symbol (4), Hit (13), Score.1 (8)
    screen = screen.loc[:, [4-1, 13-1, 8-1]].rename(columns={4-1:"symbol",13-1:"hit",8-1:"score1"})
    return screen

def seqs_by_symbol_from_fasta(fasta_path):
    """Extract sequences by gene symbol from FASTA file."""
    # Accepts FASTA with headers containing gene symbols (RefSeq/Ensembl names)
    # Heuristic: take longest transcript per symbol
    store = defaultdict(list)
    
    if not os.path.exists(fasta_path):
        print(f"Warning: FASTA file not found: {fasta_path}")
        return {}
    
    with (gzip.open(fasta_path, "rt") if fasta_path.endswith(".gz") else open(fasta_path)) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            header = rec.description.upper()
            seq = validate_seq(str(rec.seq).upper())
            # Try to extract gene symbol via common tokens
            m = re.search(r"GN[:=]([A-Z0-9\-\.]+)", header)
            if m:
                sym = m.group(1)
            else:
                # fallback: split on spaces and pick token without pipes
                tokens = [t for t in re.split(r"[ \|]", header) if t.isalnum()]
                sym = tokens[0] if tokens else rec.id
            store[sym].append(seq)
    out = {}
    for sym, seqs in store.items():
        out[sym] = max(seqs, key=len)
    return out

def main():
    """Main execution function."""
    ORCS_DIR = os.environ.get("ORCS_DIR","data/BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens")
    FASTA = os.environ.get("FASTA", "data/hg38_cdna.fasta.gz")  # must be human DNA
    out_csv = os.environ.get("OUT", "results/orcs_v1.1.17/summary.csv")
    max_screens = int(os.environ.get("MAX_SCREENS", "0"))  # 0 = no limit
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    
    print("Loading ORCS index...")
    idx = load_orcs_index(ORCS_DIR)
    todo = select_screens(idx)
    
    # Limit screens for testing
    if max_screens > 0:
        todo = todo.head(max_screens)
        print(f"Limited to first {max_screens} screens for testing")
    
    print(f"Found {len(todo)} qualifying screens")
    
    print("Loading sequences from FASTA...")
    sym2seq = seqs_by_symbol_from_fasta(FASTA)
    print(f"Loaded {len(sym2seq)} gene sequences")
    
    records = []
    for i, (_, row) in enumerate(todo.iterrows()):
        sid = row["screen_id"]
        stype = str(row["score1_type"]).lower()
        print(f"Processing screen {i+1}/{len(todo)}: {sid}")
        
        scr = load_screen_rows(ORCS_DIR, sid)
        if scr is None or scr.empty:
            print(f"  No data found for screen {sid}")
            continue
            
        # join sequences
        sym_u = [s.upper() for s in scr["symbol"].astype(str).tolist()]
        have = [(s, sym2seq.get(s)) for s in sym_u]
        scr["seq"] = [seq if seq is not None else "" for _, seq in have]
        scr = scr[scr["seq"]!=""]
        if scr.empty:
            print(f"  No sequences found for screen {sid}")
            continue

        print(f"  Found {len(scr)} genes with sequences")

        # build per-gene Z aggregates for each metric
        agg_by_metric = {k: {} for k in PROP_TABLES.keys()}
        for sym, seq in zip(scr["symbol"].astype(str).str.upper(), scr["seq"]):
            guides = enumerate_spcas9_guides(seq)
            if not guides:
                continue
            for mk in PROP_TABLES.keys():
                agg = aggregate_gene(guides, mk)
                if agg:
                    agg_by_metric[mk][sym] = agg

        for mk, mstore in agg_by_metric.items():
            if not mstore:
                continue
            dfm = pd.DataFrame.from_dict(mstore, orient="index")
            dfm["symbol"] = dfm.index
            base = scr.copy()
            # numeric score handling
            if "p-value" in stype:
                y = pd.to_numeric(base["score1"], errors="coerce")
                base["y"] = -np.log10(y.clip(lower=1e-300))
            else:
                base["y"] = pd.to_numeric(base["score1"], errors="coerce")
            base = base.merge(dfm, on="symbol", how="inner")
            
            for agg_key in ["mean","median","top_decile_mean"]:
                x = pd.to_numeric(base[agg_key], errors="coerce").to_numpy()
                y = pd.to_numeric(base["y"], errors="coerce").to_numpy()
                mask = np.isfinite(x) & np.isfinite(y)
                x, y = x[mask], y[mask]
                if len(x) < 10 or np.std(x)==0 or np.std(y)==0:
                    continue
                r, p = stats.pearsonr(x, y)
                lo, hi = bootstrap_r(x, y, n=1000)
                # Cohen's d on Z aggregate split by ORCS Hit
                hits = base.loc[mask, "hit"].astype(str).str.upper()=="YES"
                d = cohen_d(x[hits], x[~hits]) if hits.any() and (~hits).any() else np.nan
                records.append({
                    "screen_id": sid, "score1_type": stype,
                    "metric": mk, "aggregate": agg_key,
                    "N": int(len(x)),
                    "pearson_r": float(r), "p_value": float(p),
                    "ci95_lo": lo, "ci95_hi": hi,
                    "var_Z": float(np.var(pd.to_numeric(base["mean"], errors="coerce").to_numpy())),
                    "cohens_d_hit_vs_non": d
                })

    # Save results
    results_df = pd.DataFrame.from_records(records)
    results_df.to_csv(out_csv, index=False)
    print(f"\nWrote {len(records)} results to {out_csv}")
    
    # Generate top screens summary
    if not results_df.empty:
        top_screens = results_df[results_df["pearson_r"].abs() >= 0.5]
        if not top_screens.empty:
            top_file = os.path.join(os.path.dirname(out_csv), "top_screens.md")
            with open(top_file, 'w') as f:
                f.write("# Top ORCS Screens with |r| ≥ 0.5\n\n")
                f.write(f"Found {len(top_screens)} results from {top_screens['screen_id'].nunique()} screens:\n\n")
                for _, row in top_screens.sort_values('pearson_r', key=abs, ascending=False).head(20).iterrows():
                    f.write(f"- Screen {row['screen_id']}: {row['metric']} {row['aggregate']} "
                           f"r={row['pearson_r']:.3f} (p={row['p_value']:.3e}, N={row['N']})\n")
            print(f"Generated top screens summary: {top_file}")

if __name__ == "__main__":
    main()