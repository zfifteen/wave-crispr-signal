#!/usr/bin/env python3
"""
CORRECTED VERSION: DNA Breathing Dynamics vs Arbitrary Encoding

This version fixes the statistical flaw in the original analysis by comparing
ALL individual sequence scores (not aggregated means).

CORRECTION:
- Original (WRONG): n_breath=100 vs n_arb=10 encoder means
- Corrected (RIGHT): n_breath=100 vs n_arb=1000 individual scores (100 seqs × 10 encoders)
"""

# Import the original module
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from test_breathing_dynamics_encoding import (
    BreathingDynamicsEncoder,
    ArbitraryEncoder,
    BreathingDynamicsValidator,
    DiscreteZetaShift
)

import numpy as np
import random
from scipy import stats
from typing import Dict

class CorrectedValidator(BreathingDynamicsValidator):
    """Corrected statistical analysis"""

    def run_comparative_test_corrected(self, n_sequences: int = 100,
                                      n_arbitrary_trials: int = 10) -> Dict:
        """Run CORRECTED comparison: all individual scores"""

        print("\n" + "="*70)
        print("CORRECTED BREATHING DYNAMICS ENCODING TEST")
        print("="*70)
        print("\nFIX: Comparing ALL individual scores (not encoder means)")
        print(f"  - Breathing: {n_sequences} sequences × 1 encoder = {n_sequences} scores")
        print(f"  - Arbitrary: {n_sequences} sequences × {n_arbitrary_trials} encoders = {n_sequences * n_arbitrary_trials} scores")

        sequences = self.generate_crispr_sequences(n_sequences)

        # Breathing encoding: collect individual scores
        print("\n--- Testing Breathing Dynamics Encoding ---")
        breathing_scores_gc = []
        breathing_scores_at = []
        breathing_scores_random = []

        for seq in sequences:
            z_gc = self.analyze_mutation_sensitivity(seq, self.breathing_encoder, 'gc_affecting')
            z_at = self.analyze_mutation_sensitivity(seq, self.breathing_encoder, 'at_affecting')
            z_rand = self.analyze_mutation_sensitivity(seq, self.breathing_encoder, 'random')

            breathing_scores_gc.append(z_gc)
            breathing_scores_at.append(z_at)
            breathing_scores_random.append(z_rand)

        print(f"  Collected {len(breathing_scores_gc)} individual scores")
        print(f"  GC-affecting: Mean Z = {np.mean(breathing_scores_gc):.4f} ± {np.std(breathing_scores_gc):.4f}")

        # Arbitrary encoding: collect ALL individual scores (not means!)
        print(f"\n--- Testing Arbitrary Encodings ({n_arbitrary_trials} trials) ---")
        all_arbitrary_scores_gc = []
        all_arbitrary_scores_at = []
        all_arbitrary_scores_random = []

        for trial in range(n_arbitrary_trials):
            arb_encoder = ArbitraryEncoder(seed=trial)

            for seq in sequences:
                z_gc = self.analyze_mutation_sensitivity(seq, arb_encoder, 'gc_affecting')
                z_at = self.analyze_mutation_sensitivity(seq, arb_encoder, 'at_affecting')
                z_rand = self.analyze_mutation_sensitivity(seq, arb_encoder, 'random')

                all_arbitrary_scores_gc.append(z_gc)
                all_arbitrary_scores_at.append(z_at)
                all_arbitrary_scores_random.append(z_rand)

        print(f"  Collected {len(all_arbitrary_scores_gc)} individual scores")
        print(f"  GC-affecting: Mean Z = {np.mean(all_arbitrary_scores_gc):.4f} ± {np.std(all_arbitrary_scores_gc):.4f}")

        # Statistical comparison (CORRECTED)
        print("\n" + "="*70)
        print("CORRECTED STATISTICAL ANALYSIS")
        print("="*70)

        results = {}

        for mut_type, breath_scores, arb_scores in [
            ('GC-affecting', breathing_scores_gc, all_arbitrary_scores_gc),
            ('AT-affecting', breathing_scores_at, all_arbitrary_scores_at),
            ('Random', breathing_scores_random, all_arbitrary_scores_random)
        ]:
            # t-test (CORRECTED: comparing same-level data)
            t_stat, p_value = stats.ttest_ind(breath_scores, arb_scores)

            # Effect size (CORRECTED: proper pooled std)
            pooled_std = np.sqrt(
                ((len(breath_scores) - 1) * np.var(breath_scores, ddof=1) +
                 (len(arb_scores) - 1) * np.var(arb_scores, ddof=1)) /
                (len(breath_scores) + len(arb_scores) - 2)
            )
            cohens_d = (np.mean(breath_scores) - np.mean(arb_scores)) / pooled_std

            # Additional tests
            # Levene's test for equal variance
            levene_stat, levene_p = stats.levene(breath_scores, arb_scores)

            # Shapiro-Wilk for normality (sample if needed)
            if len(breath_scores) <= 5000:
                shapiro_breath_stat, shapiro_breath_p = stats.shapiro(breath_scores)
            else:
                shapiro_breath_p = None

            if len(arb_scores) <= 5000:
                shapiro_arb_stat, shapiro_arb_p = stats.shapiro(arb_scores)
            else:
                shapiro_arb_p = None

            # Mann-Whitney U (non-parametric alternative)
            u_stat, u_p = stats.mannwhitneyu(breath_scores, arb_scores, alternative='two-sided')

            breath_better = np.mean(breath_scores) > np.mean(arb_scores)

            print(f"\n{mut_type} Mutations:")
            print(f"  Breathing: n={len(breath_scores)}, mean={np.mean(breath_scores):.4f}, std={np.std(breath_scores, ddof=1):.4f}")
            print(f"  Arbitrary: n={len(arb_scores)}, mean={np.mean(arb_scores):.4f}, std={np.std(arb_scores, ddof=1):.4f}")
            print(f"  Difference: {np.mean(breath_scores) - np.mean(arb_scores):+.4f}")
            print(f"\n  Parametric t-test:")
            print(f"    t-statistic: {t_stat:.4f}")
            print(f"    p-value: {p_value:.6f}")
            print(f"    Cohen's d: {cohens_d:+.4f}")
            print(f"\n  Non-parametric Mann-Whitney U:")
            print(f"    U-statistic: {u_stat:.0f}")
            print(f"    p-value: {u_p:.6f}")
            print(f"\n  Assumption checks:")
            print(f"    Equal variance (Levene): p={levene_p:.4f} {'✓ PASS' if levene_p > 0.05 else '✗ FAIL (use Welch t-test)'}")
            if shapiro_breath_p is not None:
                print(f"    Normality (breath): p={shapiro_breath_p:.4f} {'✓' if shapiro_breath_p > 0.05 else '✗'}")
            if shapiro_arb_p is not None:
                print(f"    Normality (arb): p={shapiro_arb_p:.4f} {'✓' if shapiro_arb_p > 0.05 else '✗'}")
            print(f"\n  Winner: {'BREATHING' if breath_better else 'ARBITRARY'} (significant: {'YES' if p_value < 0.05 else 'NO'})")

            results[mut_type] = {
                'breathing_n': len(breath_scores),
                'arbitrary_n': len(arb_scores),
                'breathing_mean': np.mean(breath_scores),
                'arbitrary_mean': np.mean(arb_scores),
                'breathing_std': np.std(breath_scores, ddof=1),
                'arbitrary_std': np.std(arb_scores, ddof=1),
                'difference': np.mean(breath_scores) - np.mean(arb_scores),
                't_statistic': t_stat,
                'p_value': p_value,
                'cohens_d': cohens_d,
                'u_statistic': u_stat,
                'u_p_value': u_p,
                'levene_p': levene_p,
                'significant': p_value < 0.05,
                'breathing_wins': breath_better
            }

        return results

    def compare_with_original(self, corrected_results: Dict) -> None:
        """Compare corrected vs original (flawed) results"""
        print("\n" + "="*70)
        print("COMPARISON: CORRECTED vs ORIGINAL (FLAWED) ANALYSIS")
        print("="*70)

        # Note: Original had n_arb=10 encoder means
        # Corrected has n_arb=1000 individual scores

        print("\nGC-affecting Mutations:")
        print("  ORIGINAL (WRONG):")
        print("    Comparison: 100 breath scores vs 10 encoder MEANS")
        print("    Cohen's d: ~+4.130 (INFLATED!)")
        print("    p-value: <0.000001")
        print("\n  CORRECTED (RIGHT):")
        gc_result = corrected_results['GC-affecting']
        print(f"    Comparison: {gc_result['breathing_n']} breath scores vs {gc_result['arbitrary_n']} arb scores")
        print(f"    Cohen's d: {gc_result['cohens_d']:+.4f} (ACCURATE)")
        print(f"    p-value: {gc_result['p_value']:.6f}")
        print(f"\n  IMPACT:")
        print(f"    Effect size {'REDUCED' if abs(gc_result['cohens_d']) < 4.0 else 'SIMILAR'}")
        print(f"    Still significant: {'YES' if gc_result['significant'] else 'NO'}")
        print(f"    Core finding: {'STANDS' if gc_result['breathing_wins'] and gc_result['significant'] else 'REJECTED'}")


def main():
    """Run corrected analysis"""

    # Set random seed for reproducibility
    random.seed(42)
    np.random.seed(42)

    print("\n" + "="*70)
    print("STATISTICAL CORRECTION: DNA BREATHING DYNAMICS")
    print("="*70)
    print("\nPROBLEM IDENTIFIED:")
    print("  Original analysis compared 100 individual scores against 10 encoder")
    print("  MEANS, which inflates effect size and violates t-test assumptions.")
    print("\nCORRECTION:")
    print("  Compare ALL individual scores at the same granularity level:")
    print("  - Breathing: 100 sequences × 1 encoder = 100 scores")
    print("  - Arbitrary: 100 sequences × 10 encoders = 1000 scores")

    # Run corrected analysis
    validator = CorrectedValidator()
    results = validator.run_comparative_test_corrected(n_sequences=100, n_arbitrary_trials=10)

    # Compare with original
    validator.compare_with_original(results)

    # Summary
    print("\n" + "="*70)
    print("CORRECTED CONCLUSIONS")
    print("="*70)

    gc_result = results['GC-affecting']

    if gc_result['breathing_wins'] and gc_result['significant']:
        print("\n✓ HYPOTHESIS STILL SUPPORTED (with corrected statistics)")
        print(f"  Cohen's d = {gc_result['cohens_d']:+.3f}")

        if abs(gc_result['cohens_d']) > 0.8:
            print("  Effect size: LARGE")
        elif abs(gc_result['cohens_d']) > 0.5:
            print("  Effect size: MEDIUM")
        elif abs(gc_result['cohens_d']) > 0.2:
            print("  Effect size: SMALL")
        else:
            print("  Effect size: NEGLIGIBLE")

        print(f"  p-value: {gc_result['p_value']:.6f}")
        print(f"  Non-parametric p-value: {gc_result['u_p_value']:.6f}")
        print("\n  The finding is REAL, though the original magnitude was overstated.")
    else:
        print("\n✗ HYPOTHESIS REJECTED (with corrected statistics)")
        print("  The apparent effect was a statistical artifact.")

    print("\n" + "="*70)
    print("CORRECTED ANALYSIS COMPLETE")
    print("="*70)

    return results


if __name__ == '__main__':
    results = main()
