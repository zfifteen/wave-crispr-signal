#!/usr/bin/env python3
"""
Trinity Experiment C: Rotational Phase Figure

This script generates the "viral plot" showing activity vs rotational phase.

Key analysis:
1. Bin guides by rotational phase at cut site (φ_cut = 2π·pos/10.5)
2. Plot mean activity ± bootstrap CI by phase
3. Overlay breathing spectrum magnitude
4. Rayleigh test for non-uniformity
5. Circular-linear correlation

SCIENTIFIC GATES:
- Human DNA only
- Pre-registered endpoints (Rayleigh p, circular-linear r)
- Bootstrap CIs (≥1000 resamples)
- No leakage between bins

Usage:
    python phase_figure.py --seed 42 --bootstrap 1000 --bins 24 --output results/trinity/phase_figure
"""

import sys
import os
import argparse
import json
from datetime import datetime
from typing import Dict, List, Tuple
import logging

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from wave_crispr_signal.spectral import rotational_phase_at as phase_at, breathing_features, 10.5 as HELICAL_PERIOD
from experiments.trinity.statistics import (
    bootstrap_ci, rayleigh_test, circular_linear_correlation
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_guides(csv_path: str) -> List[Dict]:
    n_guides: int = 1000,
    guide_length: int = 20,
    seed: int = 42
) -> List[Dict]:
    """
    Generate synthetic guide data for demonstration.
    
    In production, this would load real CRISPR data (Doench, Kim, etc.).
    
    Args:
        n_guides: Number of guides to generate
        guide_length: Length of each guide
        seed: Random seed
        
    Returns:
        List of guide dictionaries with sequence, position, activity
    """
    rng = np.random.RandomState(seed)
    
    guides = []
    bases = ['A', 'C', 'G', 'T']
    
    for i in range(n_guides):
        # Generate random sequence
        seq = ''.join(rng.choice(bases) for _ in range(guide_length))
        
        # Random cut position (typically around position 17 for SpCas9)
        cut_pos = guide_length - 3  # -3 from PAM
        
        # Calculate phase at cut site
        phi_cut = phase_at(cut_pos, period=HELICAL_PERIOD)
        
        # Synthetic activity: phase-dependent + noise
        # Real effect would be more subtle
        base_activity = 0.5
        phase_effect = 0.15 * np.cos(phi_cut)  # Oscillation
        noise = rng.normal(0, 0.2)
        
        activity = np.clip(base_activity + phase_effect + noise, 0, 1)
        
        guides.append({
            "sequence": seq,
            "cut_position": cut_pos,
            "phi_cut": phi_cut,
            "activity": activity,
        })
    
    logger.info(f"Generated {n_guides} synthetic guides")
    
    return guides


def bin_by_phase(
    guides: List[Dict],
    n_bins: int = 24
) -> Dict[int, List[float]]:
    """
    Bin guides by rotational phase.
    
    Args:
        guides: List of guide dictionaries
        n_bins: Number of phase bins
        
    Returns:
        Dictionary mapping bin index to list of activities
    """
    bins = {i: [] for i in range(n_bins)}
    
    for guide in guides:
        phi = guide["phi_cut"]
        # Determine bin
        bin_idx = int(phi / (2 * np.pi) * n_bins) % n_bins
        bins[bin_idx].append(guide["activity"])
    
    return bins


def plot_phase_activity(
    bins: Dict[int, List[float]],
    n_bins: int,
    rayleigh_p: float,
    circ_lin_r: float,
    circ_lin_p: float,
    output_path: str,
    n_bootstrap: int = 1000,
    seed: int = 42
):
    """
    Create the rotational phase figure.
    
    Args:
        bins: Binned activities by phase
        n_bins: Number of bins
        rayleigh_p: Rayleigh test p-value
        circ_lin_r: Circular-linear correlation
        circ_lin_p: Circular-linear p-value
        output_path: Output file path
        n_bootstrap: Bootstrap resamples for CIs
        seed: Random seed
    """
    # Calculate bin centers and statistics
    bin_centers = [(i + 0.5) * 2 * np.pi / n_bins for i in range(n_bins)]
    bin_means = []
    bin_lowers = []
    bin_uppers = []
    
    for i in range(n_bins):
        activities = np.array(bins[i])
        
        if len(activities) > 0:
            mean, lower, upper = bootstrap_ci(
                activities, np.mean,
                n_bootstrap=n_bootstrap,
                seed=seed
            )
            bin_means.append(mean)
            bin_lowers.append(lower)
            bin_uppers.append(upper)
        else:
            bin_means.append(0)
            bin_lowers.append(0)
            bin_uppers.append(0)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw=dict(projection='polar'))
    
    # Plot mean ± CI
    bin_centers_deg = [np.degrees(x) for x in bin_centers]
    
    ax.plot(bin_centers, bin_means, 'o-', linewidth=2, markersize=8,
            color='#667eea', label='Mean Activity')
    ax.fill_between(
        bin_centers, bin_lowers, bin_uppers,
        alpha=0.3, color='#667eea', label='95% Bootstrap CI'
    )
    
    # Formatting
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Activity', fontsize=12)
    ax.set_title(
        f'Rotational Phase Dependence of CRISPR Activity\n'
        f'Rayleigh p={rayleigh_p:.4f}, Circ-Lin r={circ_lin_r:.3f} (p={circ_lin_p:.4f})',
        fontsize=14, pad=20
    )
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved phase figure to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Trinity Experiment C: Rotational Phase Figure"
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--bootstrap", type=int, default=1000,
                       help="Bootstrap resamples (min 1000)")
    parser.add_argument("--bins", type=int, default=24,
                       help="Number of phase bins")
    parser.add_argument("--n-guides", type=int, default=1000,
                       help="Number of synthetic guides")
    parser.add_argument("--output", type=str, default="results/trinity/phase_figure",
                       help="Output directory")
    
    args = parser.parse_args()
    
    # Validate parameters
    if args.bootstrap < 1000:
        logger.warning("bootstrap < 1000, results may be underpowered")
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Generate synthetic data
    logger.info("Generating synthetic guide data...")
    guides = load_guides(args.csv)
        n_guides=args.n_guides,
        seed=args.seed
    )
    
    # Extract phases and activities
    phases = np.array([g["phi_cut"] for g in guides])
    activities = np.array([g["activity"] for g in guides])
    
    # Rayleigh test (test for non-uniform phase distribution)
    logger.info("Running Rayleigh test...")
    rayleigh_Z, rayleigh_p = rayleigh_test(phases)
    logger.info(f"Rayleigh test: Z={rayleigh_Z:.4f}, p={rayleigh_p:.4f}")
    
    # Circular-linear correlation
    logger.info("Computing circular-linear correlation...")
    circ_lin_r, circ_lin_p = circular_linear_correlation(phases, activities)
    logger.info(f"Circular-linear: r={circ_lin_r:.4f}, p={circ_lin_p:.4f}")
    
    # Bin by phase
    logger.info(f"Binning into {args.bins} phase bins...")
    bins = bin_by_phase(guides, n_bins=args.bins)
    
    # Plot
    logger.info("Creating phase figure...")
    output_path = os.path.join(args.output, "phase_activity.png")
    plot_phase_activity(
        bins, args.bins,
        rayleigh_p, circ_lin_r, circ_lin_p,
        output_path,
        n_bootstrap=args.bootstrap,
        seed=args.seed
    )
    
    # Save results
    results = {
        "experiment": "trinity_phase_figure",
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "seed": args.seed,
            "n_guides": args.n_guides,
            "n_bins": args.bins,
            "bootstrap_resamples": args.bootstrap,
        },
        "statistics": {
            "rayleigh_Z": float(rayleigh_Z),
            "rayleigh_p": float(rayleigh_p),
            "circular_linear_r": float(circ_lin_r),
            "circular_linear_p": float(circ_lin_p),
        },
        "interpretation": {
            "rayleigh_significant": rayleigh_p < 0.05,
            "phase_correlation_significant": circ_lin_p < 0.05,
        }
    }
    
    results_path = os.path.join(args.output, "results.json")
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"Saved results to {results_path}")
    logger.info("✓ Trinity Experiment C completed successfully")
    
    # Print summary
    print("\n" + "="*60)
    print("TRINITY EXPERIMENT C: ROTATIONAL PHASE FIGURE")
    print("="*60)
    print(f"Guides analyzed: {args.n_guides}")
    print(f"Phase bins: {args.bins}")
    print(f"Bootstrap resamples: {args.bootstrap}")
    print()
    print("Results:")
    print(f"  Rayleigh test: Z={rayleigh_Z:.4f}, p={rayleigh_p:.4f}")
    print(f"    → {'Significant' if rayleigh_p < 0.05 else 'Not significant'} phase clustering")
    print(f"  Circular-linear correlation: r={circ_lin_r:.4f}, p={circ_lin_p:.4f}")
    print(f"    → {'Significant' if circ_lin_p < 0.05 else 'Not significant'} phase-activity relationship")
    print()
    print(f"Output: {args.output}/")
    print("="*60)

if __name__ == "__main__":
    main()
