#!/usr/bin/env python3
"""
Simple Validation Script for Z Framework Density Boost Claims

This script loads minimal synthetic datasets and computes fold-change between
z5d_* and baseline_* measurements to validate the >1000x density boost claims.

Usage:
    python validate.py

Expected Output:
    - Fold-change statistics for each dataset
    - Bootstrap confidence intervals  
    - Verification of >1000x boost claims
"""

import pandas as pd
import numpy as np
import argparse
import os
from pathlib import Path
import matplotlib.pyplot as plt

def bootstrap_ci(x, func=np.mean, n_boot=1000, ci=95):
    """Calculate bootstrap confidence intervals."""
    boot_stats = []
    n = len(x)
    for _ in range(n_boot):
        sample = np.random.choice(x, size=n, replace=True)
        boot_stats.append(func(sample))
    lower = np.percentile(boot_stats, (100-ci)/2)
    upper = np.percentile(boot_stats, 100 - (100-ci)/2)
    return func(x), lower, upper

def analyze_file(filepath):
    """Analyze a single dataset file."""
    try:
        df = pd.read_csv(filepath)
        print(f"\n📊 Analyzing {filepath.name}")
        print(f"   Loaded {len(df)} samples")
        
        # Determine column names based on file
        if 'nav1_8' in filepath.name:
            baseline_col = 'baseline_density'
            z5d_col = 'z5d_density'
            unit = 'density'
        elif 'bcl11a' in filepath.name:
            baseline_col = 'baseline_edit_eff'
            z5d_col = 'z5d_edit_eff'
            unit = 'edit efficiency'
        elif 'neural' in filepath.name:
            baseline_col = 'baseline_spike_rate'
            z5d_col = 'z5d_spike_rate'
            unit = 'spike rate'
        else:
            # Generic approach - find baseline and z5d columns
            baseline_cols = [col for col in df.columns if 'baseline' in col.lower()]
            z5d_cols = [col for col in df.columns if 'z5d' in col.lower()]
            
            if not baseline_cols or not z5d_cols:
                print(f"   ❌ Could not find baseline/z5d columns in {filepath.name}")
                return None
                
            baseline_col = baseline_cols[0]
            z5d_col = z5d_cols[0]
            unit = 'value'
        
        # Check if columns exist
        if baseline_col not in df.columns or z5d_col not in df.columns:
            print(f"   ❌ Missing columns: {baseline_col}, {z5d_col}")
            return None
            
        # Calculate fold changes (avoid division by zero)
        baseline_vals = df[baseline_col].values
        z5d_vals = df[z5d_col].values
        
        # Replace zeros with small value to avoid division by zero
        baseline_vals_safe = np.where(baseline_vals == 0, 0.001, baseline_vals)
        fold_changes = z5d_vals / baseline_vals_safe
        
        # Calculate statistics
        mean_fold, mean_ci_low, mean_ci_high = bootstrap_ci(fold_changes, np.mean)
        median_fold, med_ci_low, med_ci_high = bootstrap_ci(fold_changes, np.median)
        
        # Count successes
        success_1000 = np.sum(fold_changes >= 1000)
        success_100 = np.sum(fold_changes >= 100)
        success_10 = np.sum(fold_changes >= 10)
        
        results = {
            'dataset': filepath.stem,
            'n_samples': len(fold_changes),
            'unit': unit,
            'fold_changes': fold_changes,
            'mean_fold': mean_fold,
            'mean_ci': (mean_ci_low, mean_ci_high),
            'median_fold': median_fold,
            'median_ci': (med_ci_low, med_ci_high),
            'min_fold': np.min(fold_changes),
            'max_fold': np.max(fold_changes),
            'success_1000x': success_1000,
            'success_100x': success_100,
            'success_10x': success_10,
            'success_rate_1000x': success_1000 / len(fold_changes),
            'success_rate_100x': success_100 / len(fold_changes),
            'success_rate_10x': success_10 / len(fold_changes)
        }
        
        # Print results
        print(f"   📈 Mean fold-change: {mean_fold:.1f}x (95% CI: {mean_ci_low:.1f}-{mean_ci_high:.1f})")
        print(f"   📊 Median fold-change: {median_fold:.1f}x (95% CI: {med_ci_low:.1f}-{med_ci_high:.1f})")
        print(f"   📏 Range: {np.min(fold_changes):.1f}x to {np.max(fold_changes):.1f}x")
        print(f"   🚀 >1000x success: {success_1000}/{len(fold_changes)} ({success_1000/len(fold_changes):.1%})")
        print(f"   💯 >100x success: {success_100}/{len(fold_changes)} ({success_100/len(fold_changes):.1%})")
        print(f"   🔟 >10x success: {success_10}/{len(fold_changes)} ({success_10/len(fold_changes):.1%})")
        
        return results
        
    except Exception as e:
        print(f"   ❌ Error analyzing {filepath.name}: {e}")
        return None

def create_summary_plot(results_list, output_dir=None):
    """Create a summary bar chart of fold-changes."""
    if not results_list:
        return
        
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot 1: Mean fold-changes with error bars
        datasets = [r['dataset'] for r in results_list if r]
        means = [r['mean_fold'] for r in results_list if r]
        ci_lows = [r['mean_ci'][0] for r in results_list if r]
        ci_highs = [r['mean_ci'][1] for r in results_list if r]
        
        yerr = [np.array(means) - np.array(ci_lows), 
                np.array(ci_highs) - np.array(means)]
        
        bars1 = ax1.bar(datasets, means, yerr=yerr, capsize=5, alpha=0.7, color=['skyblue', 'lightcoral', 'lightgreen'])
        ax1.set_ylabel('Mean Fold-Change')
        ax1.set_title('Z Framework vs Baseline Performance')
        ax1.set_yscale('log')
        ax1.axhline(y=1000, color='red', linestyle='--', label='1000x target')
        ax1.axhline(y=100, color='orange', linestyle='--', label='100x')
        ax1.legend()
        ax1.tick_params(axis='x', rotation=45)
        
        # Plot 2: Success rates
        success_1000 = [r['success_rate_1000x'] * 100 for r in results_list if r]
        success_100 = [r['success_rate_100x'] * 100 for r in results_list if r]
        success_10 = [r['success_rate_10x'] * 100 for r in results_list if r]
        
        x = np.arange(len(datasets))
        width = 0.25
        
        ax2.bar(x - width, success_1000, width, label='>1000x', alpha=0.8, color='red')
        ax2.bar(x, success_100, width, label='>100x', alpha=0.8, color='orange') 
        ax2.bar(x + width, success_10, width, label='>10x', alpha=0.8, color='green')
        
        ax2.set_ylabel('Success Rate (%)')
        ax2.set_title('Achievement Rates by Threshold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(datasets, rotation=45)
        ax2.legend()
        ax2.set_ylim(0, 100)
        
        plt.tight_layout()
        
        if output_dir:
            plt.savefig(Path(output_dir) / 'validation_summary.png', dpi=150, bbox_inches='tight')
            print(f"📊 Plot saved to {Path(output_dir) / 'validation_summary.png'}")
        else:
            plt.savefig('validation_summary.png', dpi=150, bbox_inches='tight')
            print("📊 Plot saved to validation_summary.png")
            
        plt.show()
        
    except Exception as e:
        print(f"❌ Error creating plot: {e}")

def main():
    """Main validation function."""
    parser = argparse.ArgumentParser(description='Validate Z Framework density boost claims')
    parser.add_argument('--data-dir', default='data', help='Directory containing CSV files')
    parser.add_argument('--plot', action='store_true', help='Generate summary plots')
    parser.add_argument('--output-dir', default='.', help='Output directory for plots')
    
    args = parser.parse_args()
    
    # Set random seed for reproducible bootstrap
    np.random.seed(42)
    
    print("🔬 Z Framework Validation Script")
    print("=" * 50)
    print("🎯 Validating >1000x density boost claims")
    print("📂 Scanning for datasets...")
    
    # Find CSV files
    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        print(f"❌ Data directory {data_dir} does not exist")
        return
        
    csv_files = list(data_dir.glob('*.csv'))
    if not csv_files:
        print(f"❌ No CSV files found in {data_dir}")
        return
        
    print(f"📁 Found {len(csv_files)} dataset(s)")
    
    # Analyze each file
    results_list = []
    for filepath in sorted(csv_files):
        result = analyze_file(filepath)
        if result:
            results_list.append(result)
    
    # Print overall summary
    if results_list:
        print("\n" + "=" * 50)
        print("📋 OVERALL SUMMARY")
        print("=" * 50)
        
        total_samples = sum(r['n_samples'] for r in results_list)
        total_1000x = sum(r['success_1000x'] for r in results_list)
        overall_1000x_rate = total_1000x / total_samples if total_samples > 0 else 0
        
        print(f"📊 Total samples analyzed: {total_samples}")
        print(f"🚀 Overall >1000x success: {total_1000x}/{total_samples} ({overall_1000x_rate:.1%})")
        
        # Statistical significance test
        fold_changes_all = np.concatenate([r['fold_changes'] for r in results_list])
        
        # Test if mean fold-change is significantly > 1000
        from scipy import stats
        t_stat, p_value = stats.ttest_1samp(fold_changes_all, 1000)
        
        print(f"📈 Combined mean fold-change: {np.mean(fold_changes_all):.1f}x")
        print(f"📊 Combined median fold-change: {np.median(fold_changes_all):.1f}x")
        print(f"🧪 T-test vs 1000x: t={t_stat:.2f}, p={p_value:.6f}")
        print(f"✅ Significant (p<0.05): {'Yes' if p_value < 0.05 else 'No'}")
        
        # Determine overall conclusion
        if overall_1000x_rate >= 0.8 and p_value < 0.05:
            print("\n🎉 CONCLUSION: >1000x density boost claims VALIDATED")
            print("   ✅ High success rate (≥80%)")
            print("   ✅ Statistically significant (p<0.05)")
        elif overall_1000x_rate >= 0.5:
            print("\n⚠️  CONCLUSION: >1000x density boost claims PARTIALLY VALIDATED")
            print("   ⚠️  Moderate success rate (≥50%)")
        else:
            print("\n❌ CONCLUSION: >1000x density boost claims NOT VALIDATED")
            print("   ❌ Low success rate (<50%)")
            
        # Create plots if requested
        if args.plot:
            create_summary_plot(results_list, args.output_dir)
            
    else:
        print("❌ No valid results to summarize")
    
    print("\n" + "=" * 50)
    print("⚠️  RESEARCH USE ONLY")
    print("   Not for diagnosis, treatment, or clinical decisions")
    print("=" * 50)

if __name__ == "__main__":
    main()