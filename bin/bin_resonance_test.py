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

# ---------- Z-style enhancement heterogeneity functions ----------
def calculate_lift_at_k(x: List[float], y: List[float], k_percent: float) -> float:
    """Calculate lift@k% metric: enhancement of top k% PC guides vs bin average."""
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    
    # Pair PC (x) and efficiency (y), sort by PC descending
    paired = list(zip(x, y))
    paired.sort(key=lambda p: p[0], reverse=True)
    
    # Calculate bin average efficiency
    bin_mean = stats.fmean(y)
    if bin_mean == 0:
        return float("nan")
    
    # Calculate top k% efficiency
    n_top = max(1, int(len(paired) * k_percent / 100.0))
    top_efficiencies = [p[1] for p in paired[:n_top]]
    top_mean = stats.fmean(top_efficiencies)
    
    # Lift = (top_mean / bin_mean) - 1
    return (top_mean / bin_mean) - 1.0

def bootstrap_lift_ci(x: List[float], y: List[float], k_percent: float, n_boot=2000, seed=42, alpha=0.05) -> Tuple[float, float]:
    """Calculate bootstrap confidence interval for lift@k% metric."""
    if len(x) < 2:
        return float("nan"), float("nan")
    
    rnd = random.Random(seed)
    n = len(x)
    lifts = []
    
    for _ in range(n_boot):
        # Bootstrap sample with same indices for x and y
        idx = [rnd.randrange(n) for _ in range(n)]
        x_boot = [x[i] for i in idx]
        y_boot = [y[i] for i in idx]
        
        lift = calculate_lift_at_k(x_boot, y_boot, k_percent)
        if not math.isnan(lift):
            lifts.append(lift)
    
    if len(lifts) == 0:
        return float("nan"), float("nan")
    
    lifts.sort()
    lo_idx = int((alpha/2) * len(lifts))
    hi_idx = int((1 - alpha/2) * len(lifts))
    
    return lifts[lo_idx], lifts[hi_idx]

def fisher_z_transform(r: float) -> float:
    """Fisher z-transformation of correlation coefficient."""
    if math.isnan(r) or abs(r) >= 1.0:
        return float("nan")
    return 0.5 * math.log((1 + r) / (1 - r))

def calculate_heterogeneity_stats(correlations: List[float], sample_sizes: List[int]) -> Dict[str, float]:
    """Calculate Cochran's Q and I² for heterogeneity testing."""
    # Filter out NaN correlations and corresponding sample sizes
    valid_data = [(r, n) for r, n in zip(correlations, sample_sizes) if not math.isnan(r) and n >= 3]
    
    if len(valid_data) < 2:
        return {
            "Q": float("nan"),
            "I2": float("nan"),
            "p_heterogeneity": float("nan"),
            "df": 0,
            "pooled_z": float("nan")
        }
    
    rs, ns = zip(*valid_data)
    
    # Fisher z-transform correlations
    zs = [fisher_z_transform(r) for r in rs]
    weights = [n - 3 for n in ns]  # Variance of z_i = 1/(n_i - 3)
    
    # Pooled estimate
    total_weight = sum(weights)
    pooled_z = sum(w * z for w, z in zip(weights, zs)) / total_weight
    
    # Cochran's Q statistic
    Q = sum(w * (z - pooled_z)**2 for w, z in zip(weights, zs))
    df = len(valid_data) - 1
    
    # I² statistic
    I2 = max(0.0, (Q - df) / Q * 100.0) if Q > 0 else 0.0
    
    # p-value (chi-square distribution)
    # Using approximation: p ≈ 1 - F_chi2(Q, df) where F is CDF
    # For large Q, p ≈ 0; for small Q, p ≈ 1
    if Q <= df:
        p_heterogeneity = 1.0
    else:
        # Rough approximation using chi-square critical values
        # For df=3: 95% = 7.815, 99% = 11.345
        critical_95 = 3.84 + 2.71 * max(0, df - 1)  # Approximate
        if Q > critical_95:
            p_heterogeneity = 0.01  # Significant heterogeneity
        else:
            p_heterogeneity = 0.1   # Non-significant
    
    return {
        "Q": Q,
        "I2": I2,
        "p_heterogeneity": p_heterogeneity,
        "df": df,
        "pooled_z": pooled_z
    }

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

def calculate_gc_midpoints(q1: float, q2: float, q3: float, gc_values: List[float]) -> Tuple[float,float,float,float]:
    """Calculate GC midpoint for each quartile bin."""
    bin_gcs = {f"Q{i}": [] for i in range(1, 5)}
    
    for gc in gc_values:
        if math.isnan(gc):
            continue
        if gc <= q1:
            bin_gcs["Q1"].append(gc)
        elif gc <= q2:
            bin_gcs["Q2"].append(gc)
        elif gc <= q3:
            bin_gcs["Q3"].append(gc)
        else:
            bin_gcs["Q4"].append(gc)
    
    midpoints = []
    for i in range(1, 5):
        values = bin_gcs[f"Q{i}"]
        if values:
            midpoints.append(stats.fmean(values))
        else:
            midpoints.append(float("nan"))
    
    return tuple(midpoints)

def run(path: str, n_boot=2000, n_perm=10000, seed=42, tail="two", out=None, save_control=False):
    data = read_data(path)
    if not data: 
        raise SystemExit("No valid rows found. Need columns: sequence, efficiency")

    pcs = [phase_coherence(s) for s,_ in data]
    gcs = [gc_fraction(s) for s,_ in data]
    ys  = [e for _,e in data]

    q1,q2,q3 = quartile_edges(gcs)
    gc_midpoints = calculate_gc_midpoints(q1, q2, q3, gcs)

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
    correlations = []
    sample_sizes = []
    
    for bn in ["Q1","Q2","Q3","Q4"]:
        x, y = grouped[bn]
        n = len(x)
        
        if n < 10:
            # Insufficient data - use NaN for all metrics
            lift10 = lift10_lo = lift10_hi = float("nan")
            lift20 = lift20_lo = lift20_hi = float("nan")
            r = ci_lo = ci_hi = p = float("nan")
            
            results.append({
                "bin": bn, "n": n, "r": r, "ci_low": ci_lo, "ci_high": ci_hi, "p_perm": p,
                "lift_top10": lift10, "lift10_ci_low": lift10_lo, "lift10_ci_high": lift10_hi,
                "lift_top20": lift20, "lift20_ci_low": lift20_lo, "lift20_ci_high": lift20_hi,
                "gc_midpoint": gc_midpoints[int(bn[1])-1]
            })
            pvals.append(1.0)
            correlations.append(float("nan"))
            sample_sizes.append(n)
            continue
            
        # Calculate existing correlation metrics
        r = pearson_r(x,y)
        ci_lo, ci_hi = bootstrap_ci_r(x, y, n_boot=n_boot, seed=seed + 100 + bin_seed_offset[bn])
        p = permutation_pvalue(x, y, r, n_perm=n_perm, seed=seed + 200 + bin_seed_offset[bn], tail=tail)
        
        # Calculate Z-style enhancement metrics
        lift10 = calculate_lift_at_k(x, y, 10.0)
        lift10_lo, lift10_hi = bootstrap_lift_ci(x, y, 10.0, n_boot=n_boot, seed=seed + 300 + bin_seed_offset[bn])
        
        lift20 = calculate_lift_at_k(x, y, 20.0)
        lift20_lo, lift20_hi = bootstrap_lift_ci(x, y, 20.0, n_boot=n_boot, seed=seed + 400 + bin_seed_offset[bn])
        
        results.append({
            "bin": bn, "n": n, "r": r, "ci_low": ci_lo, "ci_high": ci_hi, "p_perm": p,
            "lift_top10": lift10, "lift10_ci_low": lift10_lo, "lift10_ci_high": lift10_hi,
            "lift_top20": lift20, "lift20_ci_low": lift20_lo, "lift20_ci_high": lift20_hi,
            "gc_midpoint": gc_midpoints[int(bn[1])-1]
        })
        pvals.append(p)
        correlations.append(r)
        sample_sizes.append(n)

    # Calculate heterogeneity statistics
    heterogeneity = calculate_heterogeneity_stats(correlations, sample_sizes)

    passed, cutoff = bh_fdr(pvals, alpha=0.05)
    
    # Updated header to include Z-style enhancement metrics
    header = [
        "bin", "n", "r", "ci_low", "ci_high", "p_perm", "q_pass", "bh_cutoff",
        "lift_top10", "lift10_ci_low", "lift10_ci_high",
        "lift_top20", "lift20_ci_low", "lift20_ci_high",
        "gc_midpoint", "q1_edge", "q2_edge", "q3_edge",
        "Fisher_Q", "I2", "p_heterogeneity"
    ]

    # Print enhanced table
    print("\nGC-quartile resonance test with Z-style enhancement heterogeneity")
    print(f"Input: {path} | n_boot={n_boot} n_perm={n_perm} seed={seed} tail={tail}")
    print(f"Quartile edges: Q1≤{q1:.4f}, Q2≤{q2:.4f}, Q3≤{q3:.4f}, Q4>{q3:.4f}")
    print(f"GC midpoints: Q1={gc_midpoints[0]:.4f}, Q2={gc_midpoints[1]:.4f}, Q3={gc_midpoints[2]:.4f}, Q4={gc_midpoints[3]:.4f}")
    print()
    print("{:>3} {:>5} {:>8} {:>10} {:>10} {:>10} {:>7} {:>10} {:>10}".format(
        "bin","n","r","ci_low","ci_high","p_perm","q_pass","lift@10%","lift@20%"))
    
    for i, (result, ok) in enumerate(zip(results, passed)):
        print("{:>3} {:>5} {:>8.4f} {:>10.4f} {:>10.4f} {:>10.4g} {:>7} {:>10.4f} {:>10.4f}".format(
            result["bin"], result["n"], result["r"], result["ci_low"], result["ci_high"], 
            result["p_perm"], str(ok), result["lift_top10"], result["lift_top20"]))
    
    print(f"BH cutoff (p*): {cutoff:.4g}")
    print()
    print("HETEROGENEITY ANALYSIS:")
    print(f"Cochran's Q: {heterogeneity['Q']:.4f} (df={heterogeneity['df']})")
    print(f"I² heterogeneity: {heterogeneity['I2']:.1f}%")
    print(f"p_heterogeneity: {heterogeneity['p_heterogeneity']:.4g}")

    # Optional CSV
    if out:
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            for i, (result, ok) in enumerate(zip(results, passed)):
                w.writerow([
                    result["bin"], result["n"], result["r"], result["ci_low"], result["ci_high"], 
                    result["p_perm"], ok, cutoff,
                    result["lift_top10"], result["lift10_ci_low"], result["lift10_ci_high"],
                    result["lift_top20"], result["lift20_ci_low"], result["lift20_ci_high"],
                    result["gc_midpoint"], q1, q2, q3,
                    heterogeneity["Q"], heterogeneity["I2"], heterogeneity["p_heterogeneity"]
                ])

    # Optional negative control (shuffled efficiencies)
    if save_control and out:
        control_path = out.replace('.csv', '_control.csv')
        print(f"\nGenerating negative control with shuffled efficiencies: {control_path}")
        
        # Shuffle efficiency values for negative control
        control_ys = ys[:]
        control_rnd = random.Random(seed + 1000)  # Different seed for control
        control_rnd.shuffle(control_ys)
        
        # Re-run analysis with shuffled efficiencies
        control_grouped = {f"Q{i}":[[],[]] for i in range(1,5)}
        for pc,gc,y in zip(pcs,gcs,control_ys):
            if math.isnan(pc) or math.isnan(gc) or math.isnan(y): 
                continue
            if gc <= q1:        control_grouped["Q1"][0].append(pc); control_grouped["Q1"][1].append(y)
            elif gc <= q2:      control_grouped["Q2"][0].append(pc); control_grouped["Q2"][1].append(y)
            elif gc <= q3:      control_grouped["Q3"][0].append(pc); control_grouped["Q3"][1].append(y)
            else:               control_grouped["Q4"][0].append(pc); control_grouped["Q4"][1].append(y)

        control_results = []
        control_pvals = []
        control_correlations = []
        control_sample_sizes = []
        
        for bn in ["Q1","Q2","Q3","Q4"]:
            x, y = control_grouped[bn]
            n = len(x)
            
            if n < 10:
                lift10 = lift10_lo = lift10_hi = float("nan")
                lift20 = lift20_lo = lift20_hi = float("nan")
                r = ci_lo = ci_hi = p = float("nan")
                
                control_results.append({
                    "bin": bn, "n": n, "r": r, "ci_low": ci_lo, "ci_high": ci_hi, "p_perm": p,
                    "lift_top10": lift10, "lift10_ci_low": lift10_lo, "lift10_ci_high": lift10_hi,
                    "lift_top20": lift20, "lift20_ci_low": lift20_lo, "lift20_ci_high": lift20_hi,
                    "gc_midpoint": gc_midpoints[int(bn[1])-1]
                })
                control_pvals.append(1.0)
                control_correlations.append(float("nan"))
                control_sample_sizes.append(n)
                continue
                
            r = pearson_r(x,y)
            ci_lo, ci_hi = bootstrap_ci_r(x, y, n_boot=n_boot, seed=seed + 1100 + bin_seed_offset[bn])
            p = permutation_pvalue(x, y, r, n_perm=n_perm, seed=seed + 1200 + bin_seed_offset[bn], tail=tail)
            
            lift10 = calculate_lift_at_k(x, y, 10.0)
            lift10_lo, lift10_hi = bootstrap_lift_ci(x, y, 10.0, n_boot=n_boot, seed=seed + 1300 + bin_seed_offset[bn])
            
            lift20 = calculate_lift_at_k(x, y, 20.0)
            lift20_lo, lift20_hi = bootstrap_lift_ci(x, y, 20.0, n_boot=n_boot, seed=seed + 1400 + bin_seed_offset[bn])
            
            control_results.append({
                "bin": bn, "n": n, "r": r, "ci_low": ci_lo, "ci_high": ci_hi, "p_perm": p,
                "lift_top10": lift10, "lift10_ci_low": lift10_lo, "lift10_ci_high": lift10_hi,
                "lift_top20": lift20, "lift20_ci_low": lift20_lo, "lift20_ci_high": lift20_hi,
                "gc_midpoint": gc_midpoints[int(bn[1])-1]
            })
            control_pvals.append(p)
            control_correlations.append(r)
            control_sample_sizes.append(n)

        control_heterogeneity = calculate_heterogeneity_stats(control_correlations, control_sample_sizes)
        control_passed, control_cutoff = bh_fdr(control_pvals, alpha=0.05)
        
        # Save control results
        with open(control_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            for i, (result, ok) in enumerate(zip(control_results, control_passed)):
                w.writerow([
                    result["bin"], result["n"], result["r"], result["ci_low"], result["ci_high"], 
                    result["p_perm"], ok, control_cutoff,
                    result["lift_top10"], result["lift10_ci_low"], result["lift10_ci_high"],
                    result["lift_top20"], result["lift20_ci_low"], result["lift20_ci_high"],
                    result["gc_midpoint"], q1, q2, q3,
                    control_heterogeneity["Q"], control_heterogeneity["I2"], control_heterogeneity["p_heterogeneity"]
                ])
        
        print(f"Control analysis saved: {control_path}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Q1–Q4 resonance test with bootstrap CI, permutation p, BH-FDR.")
    ap.add_argument("--input", required=True, help="CSV with columns: sequence, efficiency")
    ap.add_argument("--output", default="", help="Optional output CSV path")
    ap.add_argument("--n_boot", type=int, default=2000)
    ap.add_argument("--n_perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tail", choices=["two","greater","less"], default="two")
    ap.add_argument("--save_control", action="store_true", help="Generate and save negative control with shuffled efficiencies")
    args = ap.parse_args()
    run(args.input, n_boot=args.n_boot, n_perm=args.n_perm, seed=args.seed, tail=args.tail, out=(args.output or None), save_control=args.save_control)