"""
φ-Geometry Hypothesis Validation Experiment.

This script tests the hypothesis that φ (golden ratio) geometry
in B-DNA structure improves CRISPR guide prediction metrics.

HYPOTHESIS:
    φ-based features (phase alignment, spectral curvature) derived from the
    physical B-DNA geometry improve prediction of CRISPR guide efficiency
    compared to non-φ baselines (uniform-phase, random-phase, GC-only).

EMPIRICAL SUPPORT (Larsen, Symmetry 2021):
    - B-DNA helix length:diameter ≈ 1.6088 (close to φ ≈ 1.618)
    - Major/minor backbone spacing ≈ 1.64
    - Axial 10-fold symmetry with "golden diamonds"

FALSIFICATION CRITERIA:
    1. φ-features show NO significant correlation with efficiency (|r| < 0.1)
    2. φ-features do NOT improve over baselines (uniform, random, GC-only)
    3. Permutation test p-values are NOT significant (p > 0.05)

If any criterion fails (φ-features outperform), hypothesis is NOT falsified.
If all criteria hold, φ-features provide no predictive advantage.

Usage:
    python experiments/phi_geometry_validation_20251226/validate_phi_geometry.py \
        --seed 42 \
        --bootstrap 1000 \
        --permutation 1000 \
        --output experiments/phi_geometry_validation_20251226/results.json

Author: GitHub Copilot
Date: 2025-12-26
"""

import argparse
import json
import sys
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Any
import hashlib

import numpy as np
from scipy import stats

# Add repository root to path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from wave_crispr_signal.features.phi_geometry import (
    PHI,
    DNA_LENGTH_DIAMETER_RATIO,
    DNA_MAJOR_MINOR_SEP_RATIO,
    DNA_HELICAL_PERIOD,
    compute_phi_features,
    phi_phase_score,
    phi_curvature_score,
    uniform_phase_score,
    random_phase_score,
    simple_gc_content,
)


# ============================================================================
# CONFIGURATION
# ============================================================================

DEFAULT_DATASET = REPO_ROOT / "data" / "doench2016.csv"
BOOTSTRAP_SAMPLES = 1000
PERMUTATION_SAMPLES = 1000


# ============================================================================
# DATA LOADING
# ============================================================================

def load_crispr_data(filepath: Path) -> Tuple[List[str], np.ndarray]:
    """
    Load CRISPR guide sequences and efficiency labels from CSV.
    
    Args:
        filepath: Path to CSV with 'sequence' and 'efficiency' columns
        
    Returns:
        Tuple of (sequences, efficiencies)
        
    Raises:
        FileNotFoundError: If the data file does not exist
    """
    if not filepath.exists():
        raise FileNotFoundError(
            f"Dataset file not found: {filepath}\n"
            f"Please ensure the data file exists at the specified path."
        )
    
    sequences = []
    efficiencies = []
    total_rows = 0
    skipped_rows = 0
    
    with open(filepath, 'r') as f:
        header = f.readline().strip().split(',')
        seq_idx = header.index('sequence')
        eff_idx = header.index('efficiency')
        
        for line in f:
            total_rows += 1
            parts = line.strip().split(',')
            if len(parts) >= max(seq_idx, eff_idx) + 1:
                seq = parts[seq_idx].strip().upper()
                eff = float(parts[eff_idx])
                # Validate sequence (skip invalid)
                if all(c in 'ACGTN' for c in seq):
                    sequences.append(seq)
                    efficiencies.append(eff)
                else:
                    skipped_rows += 1
            else:
                skipped_rows += 1
    
    if skipped_rows > 0:
        print(f"Warning: Skipped {skipped_rows}/{total_rows} rows due to invalid sequences or malformed data")
    
    if len(sequences) == 0:
        raise ValueError("No valid sequences loaded from dataset")
    
    return sequences, np.array(efficiencies)


def compute_file_sha256(filepath: Path) -> str:
    """Compute SHA256 hash of file for reproducibility tracking."""
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            sha256.update(chunk)
    return sha256.hexdigest()


# ============================================================================
# FEATURE COMPUTATION
# ============================================================================

def compute_all_features(sequences: List[str]) -> Dict[str, np.ndarray]:
    """
    Compute all feature sets for comparison.
    
    Returns dictionary with:
        - phi_phase: φ-phase alignment scores
        - phi_curvature: φ-curvature spectral scores
        - phi_combined: geometric mean of phase + curvature
        - uniform_phase: baseline uniform scores
        - random_phase: baseline random scores
        - gc_content: GC fraction
    """
    n = len(sequences)
    
    features = {
        'phi_phase': np.zeros(n),
        'phi_curvature': np.zeros(n),
        'phi_combined': np.zeros(n),
        'uniform_phase': np.zeros(n),
        'random_phase': np.zeros(n),
        'gc_content': np.zeros(n),
    }
    
    for i, seq in enumerate(sequences):
        phi_feats = compute_phi_features(seq)
        features['phi_phase'][i] = phi_feats['phi_phase_score']
        features['phi_curvature'][i] = phi_feats['phi_curvature_score']
        features['phi_combined'][i] = phi_feats['phi_combined_score']
        features['uniform_phase'][i] = uniform_phase_score(seq)
        features['random_phase'][i] = random_phase_score(seq, seed=i)
        features['gc_content'][i] = simple_gc_content(seq)
    
    return features


# ============================================================================
# STATISTICAL ANALYSIS
# ============================================================================

def pearson_with_ci(x: np.ndarray, y: np.ndarray, 
                    n_bootstrap: int = 1000,
                    seed: int = 42) -> Dict[str, float]:
    """
    Compute Pearson correlation with 95% bootstrap CI.
    
    Returns:
        - r: correlation coefficient
        - ci_lower: 95% CI lower bound
        - ci_upper: 95% CI upper bound
    """
    r, p = stats.pearsonr(x, y)
    
    rng = np.random.default_rng(seed)
    n = len(x)
    bootstrap_rs = []
    
    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        boot_r, _ = stats.pearsonr(x[idx], y[idx])
        if not np.isnan(boot_r):
            bootstrap_rs.append(boot_r)
    
    if len(bootstrap_rs) > 0:
        ci_lower = np.percentile(bootstrap_rs, 2.5)
        ci_upper = np.percentile(bootstrap_rs, 97.5)
    else:
        ci_lower = ci_upper = r
    
    return {
        'r': float(r),
        'p_value': float(p),
        'ci_lower': float(ci_lower),
        'ci_upper': float(ci_upper),
    }


def permutation_test(x: np.ndarray, y: np.ndarray,
                     n_permutations: int = 1000,
                     seed: int = 42) -> Dict[str, float]:
    """
    Permutation test for correlation significance.
    
    Returns:
        - observed_r: observed correlation
        - p_value: permutation p-value (two-sided)
    """
    observed_r, _ = stats.pearsonr(x, y)
    
    rng = np.random.default_rng(seed)
    perm_rs = np.zeros(n_permutations)
    
    for i in range(n_permutations):
        perm_y = rng.permutation(y)
        perm_r, _ = stats.pearsonr(x, perm_y)
        perm_rs[i] = perm_r
    
    # Two-sided p-value
    p_value = np.mean(np.abs(perm_rs) >= np.abs(observed_r))
    
    return {
        'observed_r': float(observed_r),
        'p_value': float(p_value),
        'null_mean': float(np.mean(perm_rs)),
        'null_std': float(np.std(perm_rs)),
    }


def partial_correlation(x: np.ndarray, y: np.ndarray, 
                        z: np.ndarray) -> float:
    """
    Partial correlation between x and y controlling for z.
    
    Uses residuals from linear regression.
    """
    # Residuals of x on z
    slope_x = np.polyfit(z, x, 1)
    resid_x = x - np.polyval(slope_x, z)
    
    # Residuals of y on z
    slope_y = np.polyfit(z, y, 1)
    resid_y = y - np.polyval(slope_y, z)
    
    r, _ = stats.pearsonr(resid_x, resid_y)
    return float(r)


# ============================================================================
# COMPARISON ANALYSIS
# ============================================================================

def compare_features(features: Dict[str, np.ndarray],
                     efficiencies: np.ndarray,
                     n_bootstrap: int = 1000,
                     n_permutations: int = 1000,
                     seed: int = 42) -> Dict[str, Any]:
    """
    Compare all features against efficiency labels.
    
    Returns comprehensive comparison metrics.
    """
    results = {}
    
    # Skip features with constant values (e.g., uniform_phase)
    skip_features = set()
    for name, values in features.items():
        if np.std(values) < 1e-10:
            skip_features.add(name)
    
    # Analyze each feature
    for name, values in features.items():
        # Skip constant features
        if name in skip_features:
            results[name] = {
                'pearson_r': 0.0,
                'pearson_p': 1.0,
                'ci_lower': 0.0,
                'ci_upper': 0.0,
                'permutation_p': 1.0,
                'null_mean': 0.0,
                'null_std': 0.0,
                'partial_r_gc': 0.0,
                'skipped': True,
                'reason': 'constant value (no variance)',
            }
            continue
        
        # Basic correlation
        corr = pearson_with_ci(values, efficiencies, n_bootstrap, seed)
        
        # Permutation test
        perm = permutation_test(values, efficiencies, n_permutations, seed)
        
        # Partial correlation controlling for GC content
        if name != 'gc_content' and 'gc_content' not in skip_features:
            partial_r_gc = partial_correlation(values, efficiencies, 
                                               features['gc_content'])
        else:
            partial_r_gc = corr['r']
        
        results[name] = {
            'pearson_r': corr['r'],
            'pearson_p': corr['p_value'],
            'ci_lower': corr['ci_lower'],
            'ci_upper': corr['ci_upper'],
            'permutation_p': perm['p_value'],
            'null_mean': perm['null_mean'],
            'null_std': perm['null_std'],
            'partial_r_gc': partial_r_gc,
        }
    
    # Pairwise comparisons: φ-features vs baselines
    phi_features = ['phi_phase', 'phi_curvature', 'phi_combined']
    baselines = ['uniform_phase', 'random_phase', 'gc_content']
    
    comparisons = {}
    for phi_feat in phi_features:
        for baseline in baselines:
            r_phi = results[phi_feat]['pearson_r']
            r_baseline = results[baseline]['pearson_r']
            
            if np.isnan(r_phi):
                r_phi = 0.0
            if np.isnan(r_baseline):
                r_baseline = 0.0
            
            key = f"{phi_feat}_vs_{baseline}"
            comparisons[key] = {
                'phi_r': r_phi,
                'baseline_r': r_baseline,
                'phi_abs_r': abs(r_phi),
                'baseline_abs_r': abs(r_baseline),
                'phi_advantage': abs(r_phi) - abs(r_baseline),
                'phi_better': abs(r_phi) > abs(r_baseline),
            }
    
    results['pairwise_comparisons'] = comparisons
    
    return results


# ============================================================================
# FALSIFICATION LOGIC
# ============================================================================

def evaluate_falsification(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Evaluate whether the φ-geometry hypothesis is falsified.
    
    Falsification criteria:
    1. φ-features show NO significant correlation (|r| < 0.1)
    2. φ-features do NOT outperform baselines
    3. Permutation p-values are NOT significant (p > 0.05)
    
    Returns falsification verdict and reasoning.
    """
    phi_features = ['phi_phase', 'phi_curvature', 'phi_combined']
    
    # Criterion 1: Correlation magnitude
    max_phi_r = max(abs(results[f]['pearson_r']) for f in phi_features)
    criterion_1_falsified = max_phi_r < 0.1
    
    # Criterion 2: Comparison to baselines
    comparisons = results['pairwise_comparisons']
    phi_better_count = sum(1 for v in comparisons.values() if v['phi_better'])
    total_comparisons = len(comparisons)
    criterion_2_falsified = phi_better_count <= total_comparisons // 2
    
    # Criterion 3: Statistical significance
    min_phi_p = min(results[f]['permutation_p'] for f in phi_features)
    criterion_3_falsified = min_phi_p > 0.05
    
    # Overall falsification
    all_falsified = criterion_1_falsified and criterion_2_falsified and criterion_3_falsified
    
    return {
        'hypothesis_falsified': all_falsified,
        'criterion_1': {
            'description': 'φ-features show weak correlation (|r| < 0.1)',
            'met': criterion_1_falsified,
            'max_phi_r': max_phi_r,
        },
        'criterion_2': {
            'description': 'φ-features do NOT outperform baselines',
            'met': criterion_2_falsified,
            'phi_wins': phi_better_count,
            'total_comparisons': total_comparisons,
        },
        'criterion_3': {
            'description': 'φ-correlations are NOT significant (p > 0.05)',
            'met': criterion_3_falsified,
            'min_phi_p': min_phi_p,
        },
        'conclusion': (
            "HYPOTHESIS FALSIFIED: φ-geometry features provide no predictive advantage."
            if all_falsified else
            "HYPOTHESIS NOT FALSIFIED: φ-geometry features show some predictive signal."
        ),
    }


# ============================================================================
# MAIN EXPERIMENT
# ============================================================================

def run_experiment(
    data_path: Path,
    seed: int = 42,
    n_bootstrap: int = 1000,
    n_permutations: int = 1000,
    output_path: Path = None,
) -> Dict[str, Any]:
    """
    Run the complete validation experiment.
    """
    np.random.seed(seed)
    start_time = datetime.now()
    
    # Metadata
    metadata = {
        'experiment': 'phi_geometry_validation_20251226',
        'date': start_time.isoformat(),
        'seed': seed,
        'n_bootstrap': n_bootstrap,
        'n_permutations': n_permutations,
        'data_file': str(data_path),
        'data_sha256': compute_file_sha256(data_path),
        'phi_constants': {
            'phi': PHI,
            'dna_length_diameter_ratio': DNA_LENGTH_DIAMETER_RATIO,
            'dna_major_minor_sep_ratio': DNA_MAJOR_MINOR_SEP_RATIO,
            'dna_helical_period': DNA_HELICAL_PERIOD,
        },
    }
    
    # Load data
    print(f"Loading data from {data_path}...")
    sequences, efficiencies = load_crispr_data(data_path)
    print(f"Loaded {len(sequences)} sequences")
    
    metadata['n_sequences'] = len(sequences)
    metadata['efficiency_range'] = [float(efficiencies.min()), float(efficiencies.max())]
    
    # Compute features
    print("Computing features...")
    features = compute_all_features(sequences)
    
    # Statistical analysis
    print("Running statistical analysis...")
    results = compare_features(
        features, efficiencies,
        n_bootstrap=n_bootstrap,
        n_permutations=n_permutations,
        seed=seed,
    )
    
    # Falsification evaluation
    print("Evaluating falsification criteria...")
    falsification = evaluate_falsification(results)
    
    # Compile final results
    end_time = datetime.now()
    runtime = (end_time - start_time).total_seconds()
    
    output = {
        'metadata': metadata,
        'feature_results': results,
        'falsification': falsification,
        'runtime_seconds': runtime,
    }
    
    # Save if output path provided
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(output, f, indent=2)
        print(f"Results saved to {output_path}")
    
    return output


def format_results_report(results: Dict[str, Any]) -> str:
    """Format results as a human-readable report."""
    lines = []
    
    lines.append("=" * 80)
    lines.append("φ-GEOMETRY HYPOTHESIS VALIDATION EXPERIMENT")
    lines.append("=" * 80)
    lines.append("")
    
    # Executive Summary
    falsification = results['falsification']
    lines.append("EXECUTIVE SUMMARY")
    lines.append("-" * 40)
    lines.append(falsification['conclusion'])
    lines.append("")
    
    # Criterion details
    lines.append("FALSIFICATION CRITERIA")
    lines.append("-" * 40)
    for i, (crit_key, crit) in enumerate([
        ('criterion_1', results['falsification']['criterion_1']),
        ('criterion_2', results['falsification']['criterion_2']),
        ('criterion_3', results['falsification']['criterion_3']),
    ], 1):
        status = "✓ MET" if crit['met'] else "✗ NOT MET"
        lines.append(f"{i}. {crit['description']}")
        lines.append(f"   Status: {status}")
        if 'max_phi_r' in crit:
            lines.append(f"   Max |r|: {crit['max_phi_r']:.4f}")
        if 'phi_wins' in crit:
            lines.append(f"   φ wins: {crit['phi_wins']}/{crit['total_comparisons']}")
        if 'min_phi_p' in crit:
            lines.append(f"   Min p-value: {crit['min_phi_p']:.4f}")
        lines.append("")
    
    # Feature correlations
    lines.append("FEATURE CORRELATIONS WITH EFFICIENCY")
    lines.append("-" * 40)
    lines.append(f"{'Feature':<20} {'r':>8} {'95% CI':>20} {'p-perm':>10}")
    lines.append("-" * 60)
    
    for name in ['phi_phase', 'phi_curvature', 'phi_combined', 
                 'uniform_phase', 'random_phase', 'gc_content']:
        feat = results['feature_results'][name]
        ci = f"[{feat['ci_lower']:.3f}, {feat['ci_upper']:.3f}]"
        lines.append(f"{name:<20} {feat['pearson_r']:>8.4f} {ci:>20} {feat['permutation_p']:>10.4f}")
    
    lines.append("")
    
    # Metadata
    meta = results['metadata']
    lines.append("EXPERIMENT METADATA")
    lines.append("-" * 40)
    lines.append(f"Date: {meta['date']}")
    lines.append(f"Seed: {meta['seed']}")
    lines.append(f"Sequences: {meta['n_sequences']}")
    lines.append(f"Bootstrap samples: {meta['n_bootstrap']}")
    lines.append(f"Permutations: {meta['n_permutations']}")
    lines.append(f"Data SHA256: {meta['data_sha256'][:16]}...")
    lines.append(f"Runtime: {results['runtime_seconds']:.2f}s")
    lines.append("")
    
    # Physical constants
    lines.append("φ-GEOMETRY CONSTANTS (Larsen, Symmetry 2021)")
    lines.append("-" * 40)
    phi_const = meta['phi_constants']
    lines.append(f"φ (golden ratio): {phi_const['phi']:.6f}")
    lines.append(f"DNA helix L/D: {phi_const['dna_length_diameter_ratio']:.4f}")
    lines.append(f"Major/minor sep: {phi_const['dna_major_minor_sep_ratio']:.4f}")
    lines.append(f"Helical period: {phi_const['dna_helical_period']} bp")
    
    lines.append("")
    lines.append("=" * 80)
    
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="φ-Geometry Hypothesis Validation Experiment"
    )
    parser.add_argument(
        '--data', type=Path, default=DEFAULT_DATASET,
        help='Path to CRISPR efficiency CSV (default: data/doench2016.csv)'
    )
    parser.add_argument(
        '--seed', type=int, default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    parser.add_argument(
        '--bootstrap', type=int, default=BOOTSTRAP_SAMPLES,
        help='Number of bootstrap samples (default: 1000)'
    )
    parser.add_argument(
        '--permutation', type=int, default=PERMUTATION_SAMPLES,
        help='Number of permutation samples (default: 1000)'
    )
    parser.add_argument(
        '--output', type=Path, default=None,
        help='Output JSON path (default: print to stdout)'
    )
    parser.add_argument(
        '--quiet', action='store_true',
        help='Suppress progress output'
    )
    
    args = parser.parse_args()
    
    # Run experiment
    results = run_experiment(
        data_path=args.data,
        seed=args.seed,
        n_bootstrap=args.bootstrap,
        n_permutations=args.permutation,
        output_path=args.output,
    )
    
    # Print report
    if not args.quiet:
        print("\n")
        print(format_results_report(results))
    
    # Exit code based on falsification
    sys.exit(0 if results['falsification']['hypothesis_falsified'] else 1)


if __name__ == '__main__':
    main()
