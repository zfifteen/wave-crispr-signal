#!/usr/bin/env python3
import csv, math, random, argparse, statistics as stats
from typing import List, Tuple, Dict

# ---------- Phase-coherence on DNA (A/T/C/G) ----------
MAP = {"A": 1+0j, "T": -1+0j, "C": 0+1j, "G": 0-1j}

def phase_coherence(seq: str) -> float:
    s = 0+0j
    L = len(seq)
    if L == 0: return float("nan")
    for n, b in enumerate(seq.upper()):
        z = MAP.get(b, 0j)
        ang = 2*math.pi*n/L
        s += z * complex(math.cos(ang), math.sin(ang))
    return abs(s) / L

def gc_fraction(seq: str) -> float:
    seq = seq.upper()
    if not seq: return float("nan")
    g = seq.count("G"); c = seq.count("C")
    return (g + c) / len(seq)

# ---------- Core statistics (stdlib) ----------
def pearson_r(x: List[float], y: List[float]) -> float:
    mx, my = stats.fmean(x), stats.fmean(y)
    dx = [xi - mx for xi in x]
    dy = [yi - my for yi in y]
    num = sum(a*b for a, b in zip(dx, dy))
    den = math.sqrt(sum(a*a for a in dx) * sum(b*b for b in dy))
    return num/den if den else float("nan")

def bootstrap_ci_r(x: List[float], y: List[float], n_boot=2000, seed=42, alpha=0.05) -> Tuple[float,float]:
    rnd = random.Random(seed)
    n = len(x)
    rs = []
    for _ in range(n_boot):
        idx = [rnd.randrange(n) for __ in range(n)]
        xb = [x[i] for i in idx]; yb = [y[i] for i in idx]
        rs.append(pearson_r(xb, yb))
    rs.sort()
    lo = rs[int((alpha/2)*n_boot)]
    hi = rs[int((1-alpha/2)*n_boot)]
    return lo, hi

def permutation_pvalue(x: List[float], y: List[float], obs_r: float, n_perm=10000, seed=43, tail="two") -> float:
    rnd = random.Random(seed)
    n = len(y)
    y_perm = y[:]  # reused buffer
    hits = 0
    for _ in range(n_perm):
        rnd.shuffle(y_perm)
        r = pearson_r(x, y_perm)
        if tail == "greater":
            hits += (r >= obs_r)
        elif tail == "less":
            hits += (r <= obs_r)
        else:
            hits += (abs(r) >= abs(obs_r))
    # add-one smoothing (unbiased for finite perms)
    return (hits + 1) / (n_perm + 1)

def bh_fdr(pvals: List[float], alpha=0.05) -> Tuple[List[bool], float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    thresh = [((k+1)/m)*alpha for k in range(m)]
    cutoff_rank = -1
    for rank, i in enumerate(order):
        if pvals[i] <= thresh[rank]:
            cutoff_rank = rank
    passed = [False]*m
    cutoff = None
    if cutoff_rank >= 0:
        cutoff = thresh[cutoff_rank]
        for rank, i in enumerate(order):
            if pvals[i] <= cutoff:
                passed[i] = True
    return passed, (cutoff if cutoff is not None else float("nan"))

# ---------- Data, quartiles, and run ----------
def read_data(path: str):
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            seq = (row.get("sequence") or row.get("seq") or "").strip()
            eff = row.get("efficiency") or row.get("eff")
            if not seq or eff is None: continue
            try:
                rows.append((seq, float(eff)))
            except ValueError:
                continue
    return rows

def quartile_edges(values: List[float]) -> Tuple[float,float,float]:
    vs = sorted(v for v in values if not math.isnan(v))
    q1 = vs[int(0.25*(len(vs)-1))]
    q2 = vs[int(0.50*(len(vs)-1))]
    q3 = vs[int(0.75*(len(vs)-1))]
    return q1, q2, q3

def run(path: str, n_boot=2000, n_perm=10000, seed=42, tail="two", out=None):
    data = read_data(path)
    if not data: 
        raise SystemExit("No valid rows found. Need columns: sequence, efficiency")

    pcs = [phase_coherence(s) for s,_ in data]
    gcs = [gc_fraction(s) for s,_ in data]
    ys  = [e for _,e in data]

    q1,q2,q3 = quartile_edges(gcs)

    # Bin assignment
    grouped = {f"Q{i}":[[],[]] for i in range(1,5)}
    for pc,gc,y in zip(pcs,gcs,ys):
        if math.isnan(pc) or math.isnan(gc) or math.isnan(y): 
            continue
        if gc <= q1:        grouped["Q1"][0].append(pc); grouped["Q1"][1].append(y)
        elif gc <= q2:      grouped["Q2"][0].append(pc); grouped["Q2"][1].append(y)
        elif gc <= q3:      grouped["Q3"][0].append(pc); grouped["Q3"][1].append(y)
        else:               grouped["Q4"][0].append(pc); grouped["Q4"][1].append(y)

    # Deterministic per-bin seeds
    bin_seed_offset = {"Q1": 0, "Q2": 1, "Q3": 2, "Q4": 3}

    results = []
    pvals = []
    for bn in ["Q1","Q2","Q3","Q4"]:
        x, y = grouped[bn]
        n = len(x)
        if n < 10:
            results.append((bn, n, float("nan"), float("nan"), float("nan"), float("nan")))
            pvals.append(1.0)
            continue
        r = pearson_r(x,y)
        ci_lo, ci_hi = bootstrap_ci_r(x, y, n_boot=n_boot, seed=seed + 100 + bin_seed_offset[bn])
        p = permutation_pvalue(x, y, r, n_perm=n_perm, seed=seed + 200 + bin_seed_offset[bn], tail=tail)
        results.append((bn, n, r, ci_lo, ci_hi, p))
        pvals.append(p)

    passed, cutoff = bh_fdr(pvals, alpha=0.05)
    header = ["bin","n","r","ci_low","ci_high","p_perm","q_pass","bh_cutoff"]

    # Print table
    print("\nGC-quartile resonance test (bootstrap CI + permutation p, BH-FDR @ α=0.05)")
    print(f"Input: {path} | n_boot={n_boot} n_perm={n_perm} seed={seed} tail={tail}")
    print("{:>3} {:>5} {:>8} {:>10} {:>10} {:>10} {:>7}".format("bin","n","r","ci_low","ci_high","p_perm","q_pass"))
    for (bn,n,r,lo,hi,p), ok in zip(results, passed):
        print("{:>3} {:>5} {:>8.4f} {:>10.4f} {:>10.4f} {:>10.4g} {:>7}".format(bn,n,r,lo,hi,p,str(ok)))
    print(f"BH cutoff (p*): {cutoff:.4g}")

    # Optional CSV
    if out:
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            for (bn,n,r,lo,hi,p), ok in zip(results, passed):
                w.writerow([bn,n,r,lo,hi,p, ok, cutoff])

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Q1–Q4 resonance test with bootstrap CI, permutation p, BH-FDR.")
    ap.add_argument("--input", required=True, help="CSV with columns: sequence, efficiency")
    ap.add_argument("--output", default="", help="Optional output CSV path")
    ap.add_argument("--n_boot", type=int, default=2000)
    ap.add_argument("--n_perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tail", choices=["two","greater","less"], default="two")
    args = ap.parse_args()
    run(args.input, n_boot=args.n_boot, n_perm=args.n_perm, seed=args.seed, tail=args.tail, out=(args.output or None))