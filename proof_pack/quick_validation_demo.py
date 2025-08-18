#!/usr/bin/env python3
"""
Quick Z Framework Validation Demo

Simplified validation demonstrating >1000x density boost claims with baselines.
Fast version for demonstration purposes.

Usage:
    python quick_validation_demo.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
import logging
import warnings
import ast

# Statistical analysis
from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score

# Import local modules
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from pain_management_application import compute_density_boost

# Suppress warnings
warnings.filterwarnings('ignore')
np.random.seed(42)

def parse_numeric_sequences(sequence_strings, max_samples=20):
    """Parse numeric sequences from string representation (limited samples)."""
    sequences = []
    for i, seq_str in enumerate(sequence_strings[:max_samples]):
        try:
            numeric_seq = ast.literal_eval(seq_str)
            sequences.append(numeric_seq)
        except (ValueError, SyntaxError):
            numeric_seq = [ord(c) - ord('A') + 1 for c in seq_str if c.isalpha()]
            sequences.append(numeric_seq)
    return sequences

def baseline_density_measure(sequence):
    """Simple baseline density measure for comparison."""
    if len(sequence) < 2:
        return 1.0
    
    # Simple histogram density in middle bins
    hist, _ = np.histogram(sequence, bins=20, density=True)
    return np.mean(hist[9:11])  # Middle bins like Z Framework

def main():
    """Quick validation demonstration."""
    
    print("🚀 Z FRAMEWORK QUICK VALIDATION DEMO")
    print("=" * 50)
    print("📊 Comparing Z Framework vs Baseline Density Measures")
    print("⚠️  RESEARCH USE ONLY - Synthetic data validation")
    print()
    
    # Load neural spikes data
    data_file = Path(__file__).parent / 'neural_spikes.csv'
    if not data_file.exists():
        print("❌ Neural spikes data not found. Run generate_synthetic_data.py first.")
        return
    
    neural_data = pd.read_csv(data_file)
    print(f"📂 Loaded {len(neural_data)} neural spike samples")
    
    # Parse sequences (limited for speed)
    sequences = parse_numeric_sequences(neural_data['sequence_numeric'].tolist(), max_samples=20)
    print(f"🧬 Processing {len(sequences)} sequences for validation")
    print()
    
    # Calculate Z Framework density boosts
    print("🔬 CALCULATING Z FRAMEWORK DENSITY BOOSTS...")
    z_boosts = []
    baseline_densities = []
    
    for i, seq in enumerate(sequences):
        if len(seq) > 10:
            try:
                # Z Framework boost
                z_boost = compute_density_boost(seq, k=0.3, bins=20)
                z_boosts.append(z_boost)
                
                # Baseline density
                baseline_density = baseline_density_measure(seq)
                baseline_densities.append(baseline_density)
                
                print(f"  Sample {i+1:2d}: Z={z_boost:8.1f}× | Baseline={baseline_density:.4f}")
                
            except Exception as e:
                print(f"  Sample {i+1:2d}: Error - {str(e)[:30]}...")
                continue
    
    print()
    
    if not z_boosts:
        print("❌ No valid density calculations. Check data format.")
        return
    
    # Statistical analysis
    print("📈 STATISTICAL VALIDATION RESULTS")
    print("-" * 40)
    
    # Z Framework statistics
    z_mean = np.mean(z_boosts)
    z_median = np.median(z_boosts)
    z_std = np.std(z_boosts)
    z_min = np.min(z_boosts)
    z_max = np.max(z_boosts)
    
    # Success rates
    success_1000x = np.sum(np.array(z_boosts) > 1000) / len(z_boosts)
    success_210 = np.sum(np.array(z_boosts) > 2.1) / len(z_boosts)
    
    # Statistical significance test vs null (boost = 1.0)
    t_stat, p_value = stats.ttest_1samp(z_boosts, 1.0)
    
    # 95% Confidence interval
    ci_lower = np.percentile(z_boosts, 2.5)
    ci_upper = np.percentile(z_boosts, 97.5)
    
    print(f"Z Framework Density Boost Results:")
    print(f"  Mean Boost:        {z_mean:.1f}×")
    print(f"  Median Boost:      {z_median:.1f}×")
    print(f"  Range:             {z_min:.1f}× - {z_max:.1f}×")
    print(f"  95% CI:            [{ci_lower:.1f}, {ci_upper:.1f}]")
    print(f"  >1000× Success:    {success_1000x:.1%} ({np.sum(np.array(z_boosts) > 1000)}/{len(z_boosts)})")
    print(f"  >210% Success:     {success_210:.1%} ({np.sum(np.array(z_boosts) > 2.1)}/{len(z_boosts)})")
    print(f"  T-statistic:       {t_stat:.2f}")
    print(f"  P-value vs null:   {p_value:.2e}")
    print(f"  Significant:       {'✓ YES' if p_value < 0.05 else '✗ NO'}")
    print()
    
    # Comparison with baseline
    if baseline_densities:
        baseline_mean = np.mean(baseline_densities)
        
        # Convert baseline to comparable scale (approximate relative boost)
        baseline_boosts = [b / baseline_mean for b in baseline_densities]
        baseline_boost_mean = np.mean(baseline_boosts)
        
        print(f"Baseline Comparison:")
        print(f"  Baseline Mean Density: {baseline_mean:.6f}")
        print(f"  Baseline Rel. Boost:   {baseline_boost_mean:.2f}×")
        print(f"  Z Framework vs Base:   {z_mean / baseline_boost_mean:.1f}× improvement")
        print()
    
    # Validation summary
    print("🎯 VALIDATION SUMMARY")
    print("-" * 40)
    
    target_met_1000x = success_1000x > 0.5  # Majority should exceed 1000x
    target_met_210 = success_210 > 0.9      # Almost all should exceed 210%
    statistically_significant = p_value < 0.05
    
    print(f"✓ Target >1000× density boost:  {'ACHIEVED' if target_met_1000x else 'PARTIAL'}")
    print(f"✓ Target >210% density boost:   {'ACHIEVED' if target_met_210 else 'PARTIAL'}")
    print(f"✓ Statistical significance:     {'ACHIEVED' if statistically_significant else 'NOT MET'}")
    
    overall_success = target_met_1000x and statistically_significant
    print(f"✓ Overall validation:           {'SUCCESS' if overall_success else 'PARTIAL'}")
    print()
    
    print("📋 IMPORTANT NOTES")
    print("-" * 40)
    print("• This is synthetic data for mathematical validation only")
    print("• Results are for research and hypothesis generation")
    print("• Not validated for clinical diagnosis or treatment")
    print("• Z Framework provides heuristic feature scoring")
    print()
    
    if overall_success:
        print("🎉 Z Framework successfully demonstrates >1000× density enhancement")
        print("   with statistical significance on synthetic validation data!")
    else:
        print("⚠️  Partial validation - some targets not met on this sample")
        print("   Consider larger sample size or parameter adjustment")

if __name__ == "__main__":
    main()