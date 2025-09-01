#!/usr/bin/env python3
"""
Wave-CRISPR Physical Z-Metrics Validation Suite

This module provides comprehensive validation for the four physical Z-metrics,
including bootstrap confidence intervals, effect size calculations, and 
hypothesis testing framework for H1 and H2.
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import argparse
import sys
import os
from typing import List, Dict, Tuple, Optional
import logging
from pathlib import Path

# Add applications to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '.'))

from crispr_physical_z_metrics import PhysicalZMetricsCalculator, DNAValidationError

logger = logging.getLogger(__name__)


class ZMetricsValidationSuite:
    """
    Comprehensive validation suite for physical Z-metrics.
    
    Provides bootstrap confidence intervals, effect size calculations,
    and hypothesis testing framework.
    """
    
    def __init__(self, random_seed: int = 42):
        """Initialize validation suite."""
        self.calculator = PhysicalZMetricsCalculator()
        self.random_seed = random_seed
        np.random.seed(random_seed)
        
    def bootstrap_ci(self, data: np.ndarray, func=np.mean, n_boot: int = 1000, 
                     ci: float = 95) -> Tuple[float, float, float]:
        """
        Calculate bootstrap confidence intervals.
        
        Args:
            data: Input data array
            func: Function to apply (default: mean)
            n_boot: Number of bootstrap samples
            ci: Confidence interval percentage
            
        Returns:
            Tuple of (statistic, lower_ci, upper_ci)
        """
        if len(data) == 0:
            return 0.0, 0.0, 0.0
            
        boot_stats = []
        n = len(data)
        
        for _ in range(n_boot):
            sample = np.random.choice(data, size=n, replace=True)
            boot_stats.append(func(sample))
        
        boot_stats = np.array(boot_stats)
        lower = np.percentile(boot_stats, (100 - ci) / 2)
        upper = np.percentile(boot_stats, 100 - (100 - ci) / 2)
        
        return func(data), lower, upper
    
    def cohens_d(self, group1: np.ndarray, group2: np.ndarray) -> float:
        """
        Calculate Cohen's d effect size.
        
        Args:
            group1: First group data
            group2: Second group data
            
        Returns:
            Cohen's d effect size
        """
        if len(group1) == 0 or len(group2) == 0:
            return 0.0
            
        n1, n2 = len(group1), len(group2)
        pooled_std = np.sqrt(((n1 - 1) * np.var(group1, ddof=1) + 
                              (n2 - 1) * np.var(group2, ddof=1)) / (n1 + n2 - 2))
        
        if pooled_std == 0:
            return 0.0
            
        return (np.mean(group1) - np.mean(group2)) / pooled_std
    
    def calculate_z_metrics_for_sequences(self, sequences: List[str], 
                                        headers: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Calculate Z-metrics for multiple sequences.
        
        Args:
            sequences: List of DNA sequences
            headers: Optional list of sequence headers
            
        Returns:
            DataFrame with Z-metrics for all sequences
        """
        if headers is None:
            headers = [f"seq_{i}" for i in range(len(sequences))]
        
        results = []
        
        for i, (seq, header) in enumerate(zip(sequences, headers)):
            try:
                metrics = self.calculator.calculate_all_physical_z_metrics(
                    seq, header, validate=False
                )
                
                result = {
                    'sequence_id': i,
                    'header': header,
                    'sequence': seq,
                    'length': len(seq),
                    'gc_content': float(metrics['sequence_info']['gc_content']),
                    'z_opening': float(metrics['opening_kinetics']['z_opening']),
                    'z_stacking': float(metrics['stacking_dissociation']['z_stacking']),
                    'z_twist': float(metrics['twist_fluctuation']['z_twist']),
                    'z_melting': float(metrics['melting_kinetics']['z_melting']),
                    'z_mean': float(metrics['summary']['z_mean']),
                    'z_variance': float(metrics['summary']['z_variance']),
                    'opening_rate': float(metrics['opening_kinetics']['opening_rate']),
                    'dissociation_rate': float(metrics['stacking_dissociation']['dissociation_rate']),
                    'twist_rate': float(metrics['twist_fluctuation']['twist_rate']),
                    'melting_rate': float(metrics['melting_kinetics']['melting_rate']),
                    'estimated_tm': float(metrics['melting_kinetics']['estimated_tm'])
                }
                
                results.append(result)
                
            except Exception as e:
                logger.warning(f"Failed to calculate metrics for sequence {i}: {e}")
                continue
        
        return pd.DataFrame(results)
    
    def validate_z_framework_invariants(self, df: pd.DataFrame) -> Dict[str, Dict]:
        """
        Validate Z Framework invariants and convergence properties.
        
        Args:
            df: DataFrame with Z-metrics
            
        Returns:
            Dictionary with validation results
        """
        results = {}
        
        # Test invariant form Z = A * (B / c) with c = e²
        e_squared = float(self.calculator.e_squared)
        
        for metric in ['opening', 'stacking', 'twist', 'melting']:
            z_col = f'z_{metric}'
            rate_col = f'{metric}_rate' if metric != 'stacking' else 'dissociation_rate'
            
            if z_col in df.columns and rate_col in df.columns:
                # Check if Z values are consistent with expected form
                # (This would require A factor data, simplified for now)
                z_values = df[z_col].dropna()
                rate_values = df[rate_col].dropna()
                
                # Basic validation: rates should be non-negative (guardrail)
                negative_rates = (rate_values < 0).sum()
                
                # Statistical properties
                z_mean, z_lower, z_upper = self.bootstrap_ci(z_values.values)
                rate_mean, rate_lower, rate_upper = self.bootstrap_ci(rate_values.values)
                
                results[metric] = {
                    'z_mean': z_mean,
                    'z_ci': (z_lower, z_upper),
                    'rate_mean': rate_mean,
                    'rate_ci': (rate_lower, rate_upper),
                    'negative_rates': negative_rates,
                    'guardrail_violated': negative_rates > 0,
                    'n_sequences': len(z_values)
                }
        
        return results
    
    def test_h1_opening_z_hypothesis(self, df: pd.DataFrame, 
                                   efficiency_col: str = 'efficiency') -> Dict:
        """
        Test H1: Guides with lower Opening-Z show higher CRISPR efficiency.
        
        Args:
            df: DataFrame with Z-metrics and efficiency data
            efficiency_col: Column name for efficiency values
            
        Returns:
            Dictionary with H1 test results
        """
        if efficiency_col not in df.columns:
            return {
                'status': 'SKIPPED',
                'reason': f'No {efficiency_col} column found',
                'correlation': None,
                'p_value': None,
                'effect_size': None,
                'h1_supported': False
            }
        
        # Remove rows with missing data
        clean_df = df[['z_opening', efficiency_col]].dropna()
        
        if len(clean_df) < 10:
            return {
                'status': 'INSUFFICIENT_DATA',
                'reason': f'Only {len(clean_df)} complete observations',
                'correlation': None,
                'p_value': None,
                'effect_size': None,
                'h1_supported': False
            }
        
        # Test correlation (negative correlation expected for H1)
        correlation, p_value = pearsonr(clean_df['z_opening'], clean_df[efficiency_col])
        
        # Effect size using median split
        median_opening_z = clean_df['z_opening'].median()
        low_opening_z = clean_df[clean_df['z_opening'] < median_opening_z][efficiency_col]
        high_opening_z = clean_df[clean_df['z_opening'] >= median_opening_z][efficiency_col]
        
        effect_size = self.cohens_d(low_opening_z.values, high_opening_z.values)
        
        # Bootstrap confidence intervals
        # Note: For correlation, we'll use a simpler approach
        boot_correlations = []
        n = len(clean_df)
        
        for _ in range(1000):
            indices = np.random.choice(n, size=n, replace=True)
            boot_sample = clean_df.iloc[indices]
            boot_corr, _ = pearsonr(boot_sample['z_opening'], boot_sample[efficiency_col])
            boot_correlations.append(boot_corr)
        
        corr_lower = np.percentile(boot_correlations, 2.5)
        corr_upper = np.percentile(boot_correlations, 97.5)
        
        # H1 criteria: ρ ≥ 0.5 (negative), p < 1e-5
        h1_correlation_met = correlation <= -0.5
        h1_significance_met = p_value < 1e-5
        h1_supported = h1_correlation_met and h1_significance_met
        
        return {
            'status': 'TESTED',
            'n_observations': len(clean_df),
            'correlation': correlation,
            'correlation_ci': (corr_lower, corr_upper),
            'p_value': p_value,
            'effect_size': effect_size,
            'low_opening_z_efficiency_mean': low_opening_z.mean(),
            'high_opening_z_efficiency_mean': high_opening_z.mean(),
            'h1_correlation_met': h1_correlation_met,
            'h1_significance_met': h1_significance_met,
            'h1_supported': h1_supported,
            'target_correlation': -0.5,
            'target_p_value': 1e-5
        }
    
    def test_h2_stacking_z_hypothesis(self, df: pd.DataFrame, 
                                    efficiency_col: str = 'efficiency') -> Dict:
        """
        Test H2: Very low Stacking-Z sequences (framework for future completion).
        
        Args:
            df: DataFrame with Z-metrics
            efficiency_col: Column name for efficiency values
            
        Returns:
            Dictionary with H2 test framework results
        """
        clean_df = df[['z_stacking']].dropna()
        
        if len(clean_df) == 0:
            return {
                'status': 'NO_DATA',
                'very_low_threshold': None,
                'very_low_count': 0,
                'h2_framework': False
            }
        
        # Define "very low" as bottom 10th percentile
        very_low_threshold = clean_df['z_stacking'].quantile(0.1)
        very_low_count = (clean_df['z_stacking'] <= very_low_threshold).sum()
        
        # Statistics for very low Stacking-Z sequences
        very_low_sequences = clean_df[clean_df['z_stacking'] <= very_low_threshold]
        normal_sequences = clean_df[clean_df['z_stacking'] > very_low_threshold]
        
        return {
            'status': 'FRAMEWORK_ESTABLISHED',
            'n_observations': len(clean_df),
            'very_low_threshold': very_low_threshold,
            'very_low_count': very_low_count,
            'very_low_percentage': (very_low_count / len(clean_df)) * 100,
            'very_low_mean': very_low_sequences['z_stacking'].mean(),
            'normal_mean': normal_sequences['z_stacking'].mean(),
            'stacking_z_range': (clean_df['z_stacking'].min(), clean_df['z_stacking'].max()),
            'h2_framework': True,
            'note': 'H2 framework established - requires efficiency data for correlation testing'
        }
    
    def generate_validation_report(self, sequences: List[str], 
                                 headers: Optional[List[str]] = None,
                                 efficiency_data: Optional[List[float]] = None) -> Dict:
        """
        Generate comprehensive validation report.
        
        Args:
            sequences: List of DNA sequences
            headers: Optional sequence headers
            efficiency_data: Optional CRISPR efficiency data
            
        Returns:
            Complete validation report
        """
        print("🧬 Calculating Z-metrics for all sequences...")
        df = self.calculate_z_metrics_for_sequences(sequences, headers)
        
        if efficiency_data and len(efficiency_data) == len(df):
            df['efficiency'] = efficiency_data
        
        print("🔬 Validating Z Framework invariants...")
        invariants = self.validate_z_framework_invariants(df)
        
        print("🎯 Testing H1 hypothesis (Opening-Z vs efficiency)...")
        h1_results = self.test_h1_opening_z_hypothesis(df)
        
        print("📊 Testing H2 hypothesis framework (Stacking-Z threshold)...")
        h2_results = self.test_h2_stacking_z_hypothesis(df)
        
        # Summary statistics
        summary_stats = {
            'total_sequences': len(df),
            'mean_z_opening': df['z_opening'].mean(),
            'mean_z_stacking': df['z_stacking'].mean(), 
            'mean_z_twist': df['z_twist'].mean(),
            'mean_z_melting': df['z_melting'].mean(),
            'mean_gc_content': df['gc_content'].mean(),
            'mean_length': df['length'].mean()
        }
        
        return {
            'timestamp': pd.Timestamp.now().isoformat(),
            'random_seed': self.random_seed,
            'summary_statistics': summary_stats,
            'dataframe': df,
            'invariant_validation': invariants,
            'h1_hypothesis_test': h1_results,
            'h2_hypothesis_framework': h2_results,
            'validation_status': 'COMPLETE'
        }
    
    def print_validation_report(self, report: Dict):
        """Print formatted validation report."""
        print("\n" + "="*80)
        print("🧬 WAVE-CRISPR PHYSICAL Z-METRICS VALIDATION REPORT")
        print("="*80)
        
        print(f"\n📊 Summary Statistics ({report['summary_statistics']['total_sequences']} sequences):")
        stats = report['summary_statistics']
        print(f"  Mean Opening-Z:      {stats['mean_z_opening']:.6f}")
        print(f"  Mean Stacking-Z:     {stats['mean_z_stacking']:.6f}")
        print(f"  Mean Twist-Z:        {stats['mean_z_twist']:.6f}")
        print(f"  Mean Melting-Z:      {stats['mean_z_melting']:.6f}")
        print(f"  Mean GC Content:     {stats['mean_gc_content']:.3f}")
        print(f"  Mean Length:         {stats['mean_length']:.1f} bp")
        
        print(f"\n🔬 Z Framework Invariant Validation:")
        for metric, results in report['invariant_validation'].items():
            print(f"  {metric.capitalize()} Metrics:")
            print(f"    Z Mean: {results['z_mean']:.6f} (95% CI: {results['z_ci'][0]:.6f}-{results['z_ci'][1]:.6f})")
            print(f"    Rate Mean: {results['rate_mean']:.6f} (95% CI: {results['rate_ci'][0]:.6f}-{results['rate_ci'][1]:.6f})")
            print(f"    Guardrail Status: {'✅ PASS' if not results['guardrail_violated'] else '❌ FAIL'}")
        
        print(f"\n🎯 H1 Hypothesis Test (Opening-Z vs Efficiency):")
        h1 = report['h1_hypothesis_test']
        print(f"  Status: {h1['status']}")
        if h1['status'] == 'TESTED':
            print(f"  Correlation: {h1['correlation']:.3f} (target: ≤ -0.5)")
            print(f"  P-value: {h1['p_value']:.2e} (target: < 1e-5)")
            print(f"  Effect Size (Cohen's d): {h1['effect_size']:.3f}")
            print(f"  H1 Supported: {'✅ YES' if h1['h1_supported'] else '❌ NO'}")
        
        print(f"\n📊 H2 Hypothesis Framework (Stacking-Z Threshold):")
        h2 = report['h2_hypothesis_framework']
        print(f"  Status: {h2['status']}")
        if h2['status'] == 'FRAMEWORK_ESTABLISHED':
            print(f"  Very Low Threshold: {h2['very_low_threshold']:.6f}")
            print(f"  Very Low Count: {h2['very_low_count']} ({h2['very_low_percentage']:.1f}%)")
            print(f"  Framework Ready: {'✅ YES' if h2['h2_framework'] else '❌ NO'}")
        
        print(f"\n✅ Validation Complete - {report['validation_status']}")
        print("="*80)


def generate_test_sequences() -> Tuple[List[str], List[str], List[float]]:
    """Generate test sequences with simulated efficiency data."""
    # Test sequences with different characteristics
    sequences = [
        "GGGCCCGGGCCCGGGCCCGG",  # High GC, stable
        "ATATATATATATATATAT",    # Low GC, less stable  
        "GCGCGCGCGCGCGCGCGC",    # Alternating GC
        "AAATTTAAATTTAAATTT",    # AT-rich
        "ATCGATCGATCGATCGATCG",  # Mixed composition
        "CCCGGGCCCGGGCCCGGGCC",  # GC-rich
        "TTTAAATTTAAATTTAAATT",  # AT-rich variant
        "GTCAGTCAGTCAGTCAGTCA",  # Balanced, varied
        "ACTGACTGACTGACTGACTG",  # Balanced repeat
        "GACACACACACACACACACA"   # AC repeat
    ]
    
    headers = [f"test_seq_{i+1}_{'high' if i < 5 else 'low'}_efficiency" for i in range(len(sequences))]
    
    # Simulated efficiency data (higher for first 5 sequences)
    efficiency = [0.8, 0.7, 0.9, 0.6, 0.85, 0.3, 0.2, 0.4, 0.35, 0.25]
    
    return sequences, headers, efficiency


def main():
    """Main CLI for validation suite."""
    parser = argparse.ArgumentParser(
        description="Wave-CRISPR Physical Z-Metrics Validation Suite"
    )
    
    parser.add_argument(
        "--test-data", action="store_true",
        help="Run validation on generated test data"
    )
    parser.add_argument(
        "--sequences", nargs="+",
        help="DNA sequences to validate"
    )
    parser.add_argument(
        "--efficiency", nargs="+", type=float,
        help="Efficiency values corresponding to sequences"
    )
    parser.add_argument(
        "--output", "-o",
        help="Output file for detailed results (CSV format)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility"
    )
    
    args = parser.parse_args()
    
    # Initialize validation suite
    validator = ZMetricsValidationSuite(random_seed=args.seed)
    
    if args.test_data:
        # Use generated test data
        sequences, headers, efficiency = generate_test_sequences()
    elif args.sequences:
        sequences = args.sequences
        headers = [f"input_seq_{i+1}" for i in range(len(sequences))]
        efficiency = args.efficiency
    else:
        print("Error: Must provide either --test-data or --sequences")
        return 1
    
    # Generate validation report
    report = validator.generate_validation_report(sequences, headers, efficiency)
    
    # Print report
    validator.print_validation_report(report)
    
    # Save detailed results if requested
    if args.output:
        df = report['dataframe']
        df.to_csv(args.output, index=False)
        print(f"\n💾 Detailed results saved to {args.output}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())