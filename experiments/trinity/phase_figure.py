"""
Rotational-phase figure for CRISPR guides with cut positions.
Plots average encoded magnitude vs. helical phase (0-2π).
CSV schema: guide,target,y,cut_pos
"""
from __future__ import annotations
import argparse, numpy as np, pandas as pd
from wave_crispr_signal.spectral import rotational_phase_curve
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="CSV with guide,target,y,cut_pos")
    ap.add_argument("--out", required=True, help="Output PNG file")
    ap.add_argument("--r", type=float, default=20.0)
    ap.add_argument("--bins", type=int, default=24)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    phases, mags = [], []
    for _, row in df.iterrows():
        p, m = rotational_phase_curve(row['target'], r=args.r, bins=args.bins)
        phases.append(p)
        mags.append(m)
    
    # Average across guides
    avg_mags = np.mean(mags, axis=0)
    
    plt.figure(figsize=(8,6))
    plt.plot(phases[0], avg_mags, 'o-')
    plt.xlabel('Helical Phase (radians)')
    plt.ylabel('Average Encoded Magnitude')
    plt.title('Rotational Phase Curve')
    plt.grid(True)
    plt.savefig(args.out)
    print(f"Saved to {args.out}")

if __name__ == "__main__":
    main()
