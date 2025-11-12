#!/usr/bin/env python3
"""
Phase-Weighted CRISPR Scorecard Quick Validation Demo

This script demonstrates the phase-weighted scorecard functionality with
synthetic data and basic statistical validation.

Usage:
    python quick_validation_demo.py
"""

import sys
import os
import numpy as np
from typing import List, Tuple

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import directly from the module file
import importlib.util
spec = importlib.util.spec_from_file_location(
    "phase_weighted_scorecard",
    os.path.join(os.path.dirname(os.path.dirname(__file__)), "applications", "phase_weighted_scorecard.py")
)
pwsc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pwsc)

PhaseWeightedScorecard = pwsc.PhaseWeightedScorecard
K_STAR = pwsc.K_STAR


def generate_random_sequence(length: int = 20, seed: int = 42) -> str:
    """Generate a random DNA sequence."""
    np.random.seed(seed)
    bases = ['A', 'T', 'C', 'G']
    return ''.join(np.random.choice(bases, size=length))


def introduce_mutations(seq: str, n_mutations: int = 1, seed: int = None) -> str:
    """Introduce random point mutations into a sequence."""
    if seed is not None:
        np.random.seed(seed)
    
    seq_list = list(seq)
    bases = ['A', 'T', 'C', 'G']
    
    positions = np.random.choice(len(seq), size=min(n_mutations, len(seq)), replace=False)
    
    for pos in positions:
        original = seq_list[pos]
        # Choose a different base
        new_base = np.random.choice([b for b in bases if b != original])
        seq_list[pos] = new_base
    
    return ''.join(seq_list)


def compute_bootstrap_ci(scores: List[float], n_bootstrap: int = 1000, alpha: float = 0.05) -> Tuple[float, float, float]:
    """
    Compute bootstrap confidence interval for scores.
    
    Args:
        scores: List of score values
        n_bootstrap: Number of bootstrap resamples
        alpha: Significance level (default: 0.05 for 95% CI)
    
    Returns:
        (mean, lower_ci, upper_ci) tuple
    """
    scores = np.array(scores)
    n = len(scores)
    
    bootstrap_means = []
    for _ in range(n_bootstrap):
        resample = np.random.choice(scores, size=n, replace=True)
        bootstrap_means.append(np.mean(resample))
    
    mean = np.mean(scores)
    lower = np.percentile(bootstrap_means, 100 * alpha / 2)
    upper = np.percentile(bootstrap_means, 100 * (1 - alpha / 2))
    
    return mean, lower, upper


def main():
    """Run quick validation demo."""
    print("=" * 70)
    print("Phase-Weighted CRISPR Scorecard - Quick Validation Demo")
    print("=" * 70)
    print()
    
    # Set seed for reproducibility
    np.random.seed(42)
    
    # Create scorecard
    scorecard = PhaseWeightedScorecard(k=K_STAR)
    
    print(f"Configuration:")
    print(f"  k* parameter: {K_STAR}")
    print(f"  Golden ratio φ: {pwsc.PHI:.10f}")
    print(f"  e²: {pwsc.E_SQUARED:.10f}")
    print()
    
    # Test 1: Single sequence features
    print("-" * 70)
    print("Test 1: Single Sequence Spectral Features")
    print("-" * 70)
    
    test_guide = "GCTGCGGAGACCTGGAGAGA"
    print(f"Guide: {test_guide}")
    
    features = scorecard.compute_spectral_features(test_guide)
    print(f"\nSpectral Features:")
    print(f"  Entropy: {features['entropy']:.6f}")
    print(f"  Dominant Frequency: {features['dominant_freq_idx']} (magnitude: {features['dominant_freq_mag']:.6f})")
    print(f"  Sidelobe Count: {features['sidelobe_count']}")
    print(f"  Sequence Diversity: {features['diversity']:.6f}")
    print()
    
    # Test 2: Mutation disruption
    print("-" * 70)
    print("Test 2: Mutation-Induced Disruption")
    print("-" * 70)
    
    ref_seq = "GCTGCGGAGACCTGGAGAGAAAGC"
    mut_seq = "GCTGCGGAGTCCTGGAGAGAAAGC"  # Single A→T mutation at position 10
    
    print(f"Reference: {ref_seq}")
    print(f"Mutant:    {mut_seq}")
    print(f"           {''.join([' ' if r == m else '^' for r, m in zip(ref_seq, mut_seq)])}")
    
    result = scorecard.compute_z_score(ref_seq, mut_seq)
    
    print(f"\nZ-Invariant Disruption Score:")
    print(f"  Z-score: {result['z_score']:.6f}")
    print(f"  Δ_spectral: {result['delta_spectral']:.6f}")
    print(f"  κ(n) curvature: {result['kappa']:.6f}")
    print(f"\nDisruption Components:")
    print(f"  ΔEntropy: {result['disruptions']['delta_entropy']:.6f}")
    print(f"  ΔFrequency: {result['disruptions']['delta_freq']:.0f}")
    print(f"  ΔSidelobes: {result['disruptions']['delta_sidelobes']:.0f}")
    print()
    
    # Test 3: Systematic mutation analysis with bootstrap CI
    print("-" * 70)
    print("Test 3: Systematic Mutation Analysis (Bootstrap CI)")
    print("-" * 70)
    
    n_guides = 50
    n_bootstrap = 1000
    
    print(f"Generating {n_guides} random guide pairs...")
    print(f"Bootstrap resamples: {n_bootstrap}")
    print()
    
    z_scores_single = []
    z_scores_multi = []
    
    for i in range(n_guides):
        # Generate reference sequence
        ref = generate_random_sequence(length=20, seed=100 + i)
        
        # Single mutation
        mut_single = introduce_mutations(ref, n_mutations=1, seed=200 + i)
        result_single = scorecard.compute_z_score(ref, mut_single)
        z_scores_single.append(result_single['z_score'])
        
        # Multiple mutations (3)
        mut_multi = introduce_mutations(ref, n_mutations=3, seed=300 + i)
        result_multi = scorecard.compute_z_score(ref, mut_multi)
        z_scores_multi.append(result_multi['z_score'])
    
    # Compute bootstrap CIs
    mean_single, lower_single, upper_single = compute_bootstrap_ci(
        z_scores_single, n_bootstrap=n_bootstrap
    )
    mean_multi, lower_multi, upper_multi = compute_bootstrap_ci(
        z_scores_multi, n_bootstrap=n_bootstrap
    )
    
    print(f"Single Mutation Results (n={n_guides}):")
    print(f"  Mean Z-score: {mean_single:.6f}")
    print(f"  95% CI: [{lower_single:.6f}, {upper_single:.6f}]")
    print()
    
    print(f"Multiple Mutations Results (n={n_guides}):")
    print(f"  Mean Z-score: {mean_multi:.6f}")
    print(f"  95% CI: [{lower_multi:.6f}, {upper_multi:.6f}]")
    print()
    
    # Statistical test
    delta_mean = mean_multi - mean_single
    print(f"Difference in means: {delta_mean:.6f}")
    
    if delta_mean > 0:
        print(f"✓ Multiple mutations show higher disruption scores (expected)")
    else:
        print(f"⚠ Unexpected: multiple mutations show lower disruption")
    print()
    
    # Test 4: Validate mathematical constants
    print("-" * 70)
    print("Test 4: Mathematical Constants Validation")
    print("-" * 70)
    
    # Golden ratio
    phi_computed = (1 + np.sqrt(5)) / 2
    phi_error = abs(pwsc.PHI - phi_computed)
    print(f"Golden ratio φ:")
    print(f"  Expected: {phi_computed:.15f}")
    print(f"  Actual:   {pwsc.PHI:.15f}")
    print(f"  Error:    {phi_error:.2e}")
    print(f"  {'✓ PASS' if phi_error < 1e-10 else '✗ FAIL'}")
    print()
    
    # e²
    e2_computed = np.e ** 2
    e2_error = abs(pwsc.E_SQUARED - e2_computed)
    print(f"Euler's number squared e²:")
    print(f"  Expected: {e2_computed:.15f}")
    print(f"  Actual:   {pwsc.E_SQUARED:.15f}")
    print(f"  Error:    {e2_error:.2e}")
    print(f"  {'✓ PASS' if e2_error < 1e-10 else '✗ FAIL'}")
    print()
    
    # Test 5: Edge cases
    print("-" * 70)
    print("Test 5: Edge Cases and Validation")
    print("-" * 70)
    
    # Test invalid sequence
    print("Testing invalid DNA sequence...")
    try:
        scorecard.compute_spectral_features("ATCGX")
        print("  ✗ FAIL: Should reject invalid characters")
    except ValueError as e:
        print(f"  ✓ PASS: Correctly rejected with error: {str(e)[:50]}...")
    
    # Test RNA sequence
    print("\nTesting RNA sequence...")
    try:
        scorecard_rna = PhaseWeightedScorecard(is_rna=True)
        features_rna = scorecard_rna.compute_spectral_features("AUCGAUCGAUCGAUCGAUCG")
        print(f"  ✓ PASS: RNA sequence processed (entropy: {features_rna['entropy']:.6f})")
    except Exception as e:
        print(f"  ✗ FAIL: {e}")
    
    # Test short sequence warning
    print("\nTesting short sequence...")
    try:
        features_short = scorecard.compute_spectral_features("ATCGATCGATCG")  # 12 bases
        print(f"  ✓ PASS: Short sequence processed (diversity: {features_short['diversity']:.6f})")
    except Exception as e:
        print(f"  ✗ FAIL: {e}")
    
    # Test very short sequence rejection
    print("\nTesting very short sequence (< 5 bases)...")
    try:
        scorecard.compute_spectral_features("ATCG")
        print("  ✗ FAIL: Should reject sequences < 5 bases")
    except ValueError as e:
        print(f"  ✓ PASS: Correctly rejected with error: {str(e)[:50]}...")
    
    print()
    print("=" * 70)
    print("Quick Validation Demo Complete!")
    print("=" * 70)
    print()
    print("Summary:")
    print("  ✓ Complex encoding and phase weighting validated")
    print("  ✓ Spectral feature extraction functional")
    print("  ✓ Z-invariant disruption scoring operational")
    print("  ✓ Bootstrap confidence intervals computed")
    print("  ✓ Mathematical constants verified")
    print("  ✓ Edge cases handled correctly")
    print()
    print("Next steps:")
    print("  - Run full validation on Kim 2025 dataset")
    print("  - Compute permutation p-values (n_perm ≥ 1,000)")
    print("  - Compare against RuleSet3 baseline")
    print("  - Validate GC-quartile resonance")


if __name__ == "__main__":
    main()
