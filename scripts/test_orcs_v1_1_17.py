#!/usr/bin/env python3
"""
ORCS v1.1.17 Discrete Biological Z-Metrics Testing Script

Tests four discrete biological Z-metrics (base-pair opening, base-stacking dissociation, 
helical-twist fluctuation, denaturation/melting proxy) on real human CRISPR 
screen outcomes from BioGRID-ORCS 1.1.17.

Uses discrete/biological Z-form: Z = A * (B / c) with c = e² ≈ 7.389
where A is sequence-dependent mean and B is rate of change between adjacent bases.

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
import argparse
import random
import hashlib
import subprocess

# Add project root to path for imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, os.path.join(project_root, 'scripts'))

from z_framework import ZFrameworkCalculator

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
    """Compute discrete biological Z-metric with domain validation."""
    B = mp.fabs(B)
    if B >= c:
        raise ValueError(f"Domain violation: |B| >= c in discrete domain (B={B})")
    return A * (B / c)

# Golden ratio and geodesic resolution
phi = (1 + mp.sqrt(5)) / 2
KSTAR = mp.mpf("0.3")

def theta_prime(n, k=KSTAR, phi=phi):
    """Geodesic resolution function θ′(n,k) = φ·((n mod φ)/φ)^k."""
    n = mp.mpf(n)
    return phi * ((n % phi) / phi) ** k

# === Discrete biological per-base tables (deterministic; replace with your calibrated set) ===
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
    """Calculate discrete biological Z-metric for sequence using property table."""
    # A = mean property; B = std of adjacent differences
    props = [table[b] for b in seq]
    A = mp.mpf(np.mean(props))
    diffs = np.abs(np.diff(props))
    B = mp.mpf(np.std(diffs) if len(diffs) else 0.001)
    return compute_z(A, B)

def enumerate_spcas9_guides(dna, max_guides=50):
    """Return list of (20nt, pos_index) for + and - strands (pos_index is 1-based).
    
    Args:
        dna: DNA sequence string
        max_guides: Maximum number of guides to return (for performance)
    """
    dna = str(dna).upper()
    validate_seq(dna)
    guides = []
    
    # + strand: N20-NGG
    for i in range(len(dna)-23):
        if len(guides) >= max_guides:
            break
        if dna[i+20:i+23].endswith("GG"):
            guides.append((dna[i:i+20], i+1))
    
    # - strand: CCN-N20  (revcomp) 
    if len(guides) < max_guides:
        comp = dna.translate(str.maketrans("ACGT","TGCA"))[::-1]
        for i in range(len(comp)-23):
            if len(guides) >= max_guides:
                break
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

def bootstrap_r(x, y, n=1000, seed=None):
    """Bootstrap confidence intervals for Pearson correlation."""
    if seed is not None:
        rng = np.random.default_rng(seed)
    else:
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

def permutation_test(x, y, n_permutations=1000, seed=None):
    """Permutation test for correlation significance."""
    if seed is not None:
        rng = np.random.default_rng(seed)
    else:
        rng = np.random.default_rng()
    
    # Observed correlation
    observed_r, _ = stats.pearsonr(x, y)
    
    # Permutation correlations
    null_corrs = []
    for _ in range(n_permutations):
        y_shuffled = rng.permutation(y)
        if np.std(x) > 0 and np.std(y_shuffled) > 0:
            r, _ = stats.pearsonr(x, y_shuffled)
            null_corrs.append(r)
    
    if not null_corrs:
        return 1.0  # No valid permutations
    
    # Two-tailed p-value
    null_corrs = np.array(null_corrs)
    p_value = np.mean(np.abs(null_corrs) >= np.abs(observed_r))
    return float(p_value)

def partial_correlation(x, y, z):
    """Calculate partial correlation between x and y controlling for z."""
    try:
        # Convert to numpy arrays
        x, y, z = np.array(x), np.array(y), np.array(z)
        
        # Check for valid data
        mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        if mask.sum() < 3:  # Need at least 3 points
            return np.nan, np.nan
        
        x, y, z = x[mask], y[mask], z[mask]
        
        # Check for zero variance
        if np.std(x) == 0 or np.std(y) == 0 or np.std(z) == 0:
            return np.nan, np.nan
        
        # Calculate partial correlation
        rxy, _ = stats.pearsonr(x, y)
        rxz, _ = stats.pearsonr(x, z)
        ryz, _ = stats.pearsonr(y, z)
        
        numerator = rxy - rxz * ryz
        denominator = np.sqrt((1 - rxz**2) * (1 - ryz**2))
        
        if denominator == 0:
            return np.nan, np.nan
        
        partial_r = numerator / denominator
        
        # Calculate p-value using t-distribution
        n = len(x)
        if n <= 3:
            return partial_r, np.nan
        
        t_stat = partial_r * np.sqrt((n - 3) / (1 - partial_r**2))
        p_value = 2 * (1 - stats.t.cdf(np.abs(t_stat), n - 3))
        
        return float(partial_r), float(p_value)
    except:
        return np.nan, np.nan

def cohen_d(x, y):
    """Cohen's d effect size between two groups."""
    # x,y are arrays for Hit=YES vs NO on Z aggregate
    nx, ny = len(x), len(y)
    if nx < 1 or ny < 1:
        return np.nan
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    pooled_std = np.sqrt(((nx-1)*vx + (ny-1)*vy) / (nx+ny-2))
    if pooled_std == 0:
        return np.nan
    return float((np.mean(x) - np.mean(y)) / pooled_std)

def calculate_gc_content(seq):
    """Calculate GC content of a DNA sequence."""
    if not seq:
        return 0.0
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq)

def get_git_commit_hash():
    """Get current git commit hash for reproducibility."""
    try:
        result = subprocess.run(['git', 'rev-parse', 'HEAD'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            return result.stdout.strip()[:7]  # Short hash
    except:
        pass
    return "unknown"
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
    SCREEN_ID = 0
    ORG_OFF = 34
    ENZ = 23
    LIBTYPE = 21
    SCTYPE = 14
    THR = 12
    COUNT = 27
    SCORE1_TYPE = 28
    
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
        df = df[df[0]==screen_id]
        if not df.empty:
            rows.append(df)
    if not rows:
        return None
    screen = pd.concat(rows, ignore_index=True)
    # Keep Official Symbol (4), Hit (13), Score.1 (8)
    screen = screen.loc[:, [3, 12, 7]].rename(columns={3:"symbol", 12:"hit", 7:"score1"})
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
    parser = argparse.ArgumentParser(description='ORCS v1.1.17 Discrete Biological Z-Metrics Testing')
    parser.add_argument('--orcs-dir', default=None, 
                       help='ORCS data directory (default: ORCS_DIR env var)')
    parser.add_argument('--fasta', default=None,
                       help='Human DNA sequences FASTA file (default: FASTA env var)')
    parser.add_argument('--output', default=None,
                       help='Output CSV file (default: OUT env var)')
    parser.add_argument('--max-screens', type=int, default=0,
                       help='Maximum screens to process (0=no limit)')
    parser.add_argument('--min-pairs', type=int, default=10,
                       help='Minimum gene pairs for correlation analysis')
    parser.add_argument('--seed', type=int, default=None,
                       help='Random seed for reproducibility')
    parser.add_argument('--bootstrap', type=int, default=1000,
                       help='Number of bootstrap samples for confidence intervals')
    parser.add_argument('--domain', choices=['discrete', 'biological'], default='discrete',
                       help='Z-metrics domain (default: discrete)')
    
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)
        print(f"Random seed set to: {args.seed}")
    
    # Get configuration from args or environment variables
    ORCS_DIR = args.orcs_dir or os.environ.get("ORCS_DIR","data/BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens")
    FASTA = args.fasta or os.environ.get("FASTA", "data/hg38_cdna.fasta.gz")  # must be human DNA
    out_csv = args.output or os.environ.get("OUT", "results/orcs_v1.1.17/summary.csv")
    max_screens = args.max_screens or int(os.environ.get("MAX_SCREENS", "0"))  # 0 = no limit
    min_pairs = args.min_pairs or int(os.environ.get("MIN_PAIRS", "10"))
    
    # Print configuration header
    commit_hash = get_git_commit_hash()
    print(f"ORCS Discrete Biological Z-Metrics Testing")
    print(f"==========================================")
    print(f"Version: {commit_hash}")
    print(f"Seed: {args.seed}")
    print(f"Bootstrap samples: {args.bootstrap}")
    print(f"Domain: {args.domain}")
    print(f"ORCS directory: {ORCS_DIR}")
    print(f"FASTA file: {FASTA}")
    print(f"Output: {out_csv}")
    print(f"Max screens: {max_screens if max_screens > 0 else 'unlimited'}")
    print(f"Min pairs: {min_pairs}")
    print()
    
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
        processed_genes = 0
        for sym, seq in zip(scr["symbol"].astype(str).str.upper(), scr["seq"]):
            guides = enumerate_spcas9_guides(seq, max_guides=20)  # Limit for performance
            if not guides:
                continue
            processed_genes += 1
            for mk in PROP_TABLES.keys():
                agg = aggregate_gene(guides, mk)
                if agg:
                    agg_by_metric[mk][sym] = agg

        print(f"  Processed {processed_genes} genes with guides")

        for mk, mstore in agg_by_metric.items():
            if not mstore:
                print(f"  No {mk} metrics computed")
                continue
            print(f"  Computed {mk} metrics for {len(mstore)} genes")
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
                
                # Calculate GC content and guide length for control variables
                sequences = base.loc[mask, "seq"].tolist()
                gc_content = np.array([calculate_gc_content(seq) for seq in sequences])
                guide_lengths = np.array([len(seq) for seq in sequences])
                
                print(f"    {mk} {agg_key}: {len(x)} valid pairs, std(x)={np.std(x):.4f}, std(y)={np.std(y):.4f}")
                if len(x) < min_pairs or np.std(x)==0 or np.std(y)==0:
                    print(f"    Skipping {mk} {agg_key}: insufficient data or zero variance (need >{min_pairs-1} pairs)")
                    continue
                    
                # Basic correlation
                r, p = stats.pearsonr(x, y)
                lo, hi = bootstrap_r(x, y, n=args.bootstrap, seed=args.seed)
                
                # Permutation test
                perm_p = permutation_test(x, y, n_permutations=1000, seed=args.seed)
                
                # Partial correlations controlling for GC content and guide length
                partial_r_gc, partial_p_gc = partial_correlation(x, y, gc_content)
                partial_r_len, partial_p_len = partial_correlation(x, y, guide_lengths)
                
                # Cohen's d on Z aggregate split by ORCS Hit
                hits = base.loc[mask, "hit"].astype(str).str.upper()=="YES"
                d = cohen_d(x[hits], x[~hits]) if hits.any() and (~hits).any() else np.nan
                
                record = {
                    "screen_id": sid, "score1_type": stype,
                    "metric": mk, "aggregate": agg_key,
                    "N": int(len(x)),
                    "pearson_r": float(r), "p_value": float(p),
                    "ci95_lo": lo, "ci95_hi": hi,
                    "permutation_p": perm_p,
                    "partial_r_gc": partial_r_gc, "partial_p_gc": partial_p_gc,
                    "partial_r_length": partial_r_len, "partial_p_length": partial_p_len,
                    "var_Z": float(np.var(x)),
                    "cohens_d_hit_vs_non": d,
                    "mean_gc_content": float(np.mean(gc_content)),
                    "mean_guide_length": float(np.mean(guide_lengths)),
                    "seed": args.seed,
                    "bootstrap_n": args.bootstrap,
                    "commit_hash": commit_hash
                }
                records.append(record)

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