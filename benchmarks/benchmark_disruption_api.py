#!/usr/bin/env python3
"""
Benchmarking Script for Genomic Disruption Analyzer

Tests throughput and latency metrics for the SaaS API to ensure:
- Single guide scoring: <500ms
- Batch processing: 1,000 guides/min
- Statistical validity: Bootstrap CI <1%

Scientific Gates:
- All measurements with controlled seed for reproducibility
- Bootstrap CI (≥1,000 resamples) for variance estimates
- KS test validation for score distributions

Usage:
    python benchmarks/benchmark_disruption_api.py --n-guides 100 --seed 42
    python benchmarks/benchmark_disruption_api.py --full --output results/benchmark.json
"""

import sys
import os
import time
import json
import argparse
from typing import List, Dict, Any
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy import stats

from applications.genomic_disruption_api import DisruptionAnalyzer, K_STAR


def generate_synthetic_guides(n_guides: int, length: int = 20, seed: int = 42) -> List[str]:
    """
    Generate synthetic guide sequences for benchmarking.
    
    Args:
        n_guides: Number of guides to generate
        length: Length of each guide (default: 20)
        seed: Random seed
    
    Returns:
        List of synthetic guide sequences
    """
    np.random.seed(seed)
    bases = ['A', 'C', 'G', 'T']
    
    guides = []
    for _ in range(n_guides):
        guide = ''.join(np.random.choice(bases, size=length))
        guides.append(guide)
    
    return guides


def benchmark_single_guide(
    analyzer: DisruptionAnalyzer,
    n_trials: int = 100,
) -> Dict[str, float]:
    """
    Benchmark single guide scoring latency.
    
    Args:
        analyzer: DisruptionAnalyzer instance
        n_trials: Number of trials to run
    
    Returns:
        Dictionary with timing statistics
    """
    guides = generate_synthetic_guides(n_trials)
    
    latencies = []
    
    for guide in guides:
        start = time.time()
        result = analyzer.score_guide(guide)
        end = time.time()
        
        latencies.append((end - start) * 1000)  # Convert to ms
    
    latencies = np.array(latencies)
    
    return {
        'n_trials': n_trials,
        'mean_ms': float(np.mean(latencies)),
        'median_ms': float(np.median(latencies)),
        'std_ms': float(np.std(latencies)),
        'min_ms': float(np.min(latencies)),
        'max_ms': float(np.max(latencies)),
        'p95_ms': float(np.percentile(latencies, 95)),
        'p99_ms': float(np.percentile(latencies, 99)),
    }


def benchmark_batch_processing(
    analyzer: DisruptionAnalyzer,
    batch_sizes: List[int] = [10, 50, 100, 500, 1000],
) -> Dict[str, Any]:
    """
    Benchmark batch processing throughput.
    
    Args:
        analyzer: DisruptionAnalyzer instance
        batch_sizes: List of batch sizes to test
    
    Returns:
        Dictionary with throughput statistics
    """
    results = {}
    
    for batch_size in batch_sizes:
        guides = generate_synthetic_guides(batch_size)
        
        start = time.time()
        batch_results = analyzer.batch_score(guides)
        end = time.time()
        
        elapsed_s = end - start
        throughput = batch_size / elapsed_s if elapsed_s > 0 else 0
        
        results[f'batch_{batch_size}'] = {
            'batch_size': batch_size,
            'elapsed_s': elapsed_s,
            'throughput_guides_per_min': throughput * 60,
            'avg_latency_ms': (elapsed_s / batch_size) * 1000,
        }
    
    return results


def benchmark_bootstrap_ci(
    analyzer: DisruptionAnalyzer,
    n_guides: int = 10,
    n_bootstrap: int = 1000,
) -> Dict[str, Any]:
    """
    Benchmark bootstrap CI computation time.
    
    Args:
        analyzer: DisruptionAnalyzer instance
        n_guides: Number of guides to test
        n_bootstrap: Number of bootstrap resamples
    
    Returns:
        Dictionary with CI timing statistics
    """
    guides = generate_synthetic_guides(n_guides)
    
    ci_times = []
    ci_widths = []
    
    for guide in guides:
        start = time.time()
        result = analyzer.score_guide(guide, compute_ci=True, n_bootstrap=n_bootstrap)
        end = time.time()
        
        ci_times.append((end - start) * 1000)
        
        if 'confidence_interval' in result:
            ci = result['confidence_interval']
            width = ci['upper_95'] - ci['lower_95']
            ci_widths.append(width)
    
    return {
        'n_guides': n_guides,
        'n_bootstrap': n_bootstrap,
        'mean_time_ms': float(np.mean(ci_times)),
        'median_time_ms': float(np.median(ci_times)),
        'mean_ci_width': float(np.mean(ci_widths)),
        'median_ci_width': float(np.median(ci_widths)),
    }


def validate_score_distribution(
    analyzer: DisruptionAnalyzer,
    n_guides: int = 1000,
) -> Dict[str, Any]:
    """
    Validate that score distributions are reasonable using KS test.
    
    Args:
        analyzer: DisruptionAnalyzer instance
        n_guides: Number of guides to sample
    
    Returns:
        Dictionary with distribution statistics and KS test results
    """
    guides = generate_synthetic_guides(n_guides)
    
    scores = []
    for guide in guides:
        result = analyzer.score_guide(guide)
        scores.append(result['disruption_score'])
    
    scores = np.array(scores)
    
    # KS test against uniform distribution (null hypothesis)
    # We expect scores NOT to be uniform (so p < 0.05 is good)
    ks_stat, ks_p = stats.kstest(scores, 'uniform')
    
    # Also test normality
    shapiro_stat, shapiro_p = stats.shapiro(scores)
    
    return {
        'n_guides': n_guides,
        'mean_score': float(np.mean(scores)),
        'median_score': float(np.median(scores)),
        'std_score': float(np.std(scores)),
        'min_score': float(np.min(scores)),
        'max_score': float(np.max(scores)),
        'ks_statistic': float(ks_stat),
        'ks_pvalue': float(ks_p),
        'shapiro_statistic': float(shapiro_stat),
        'shapiro_pvalue': float(shapiro_p),
    }


def run_full_benchmark(
    k: float = K_STAR,
    seed: int = 42,
    output_file: str = None,
) -> Dict[str, Any]:
    """
    Run complete benchmark suite.
    
    Args:
        k: Resolution exponent
        seed: Random seed
        output_file: Optional output JSON file
    
    Returns:
        Dictionary with all benchmark results
    """
    print("=" * 60)
    print("GENOMIC DISRUPTION ANALYZER - BENCHMARK SUITE")
    print("=" * 60)
    print()
    
    # Initialize analyzer
    print("Initializing DisruptionAnalyzer...")
    analyzer = DisruptionAnalyzer(k=k, seed=seed)
    print(f"  k = {k}")
    print(f"  seed = {seed}")
    print()
    
    results = {
        'metadata': {
            'k': k,
            'seed': seed,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        }
    }
    
    # Benchmark 1: Single guide latency
    print("Benchmark 1: Single Guide Latency")
    print("-" * 60)
    latency_results = benchmark_single_guide(analyzer, n_trials=100)
    results['single_guide_latency'] = latency_results
    
    print(f"  Trials: {latency_results['n_trials']}")
    print(f"  Mean latency: {latency_results['mean_ms']:.2f} ms")
    print(f"  Median latency: {latency_results['median_ms']:.2f} ms")
    print(f"  P95 latency: {latency_results['p95_ms']:.2f} ms")
    print(f"  P99 latency: {latency_results['p99_ms']:.2f} ms")
    
    # Check if meets target (<500ms)
    if latency_results['p95_ms'] < 500:
        print(f"  ✓ PASSED: P95 latency < 500ms target")
    else:
        print(f"  ✗ FAILED: P95 latency exceeds 500ms target")
    print()
    
    # Benchmark 2: Batch processing throughput
    print("Benchmark 2: Batch Processing Throughput")
    print("-" * 60)
    batch_results = benchmark_batch_processing(analyzer)
    results['batch_processing'] = batch_results
    
    for key, result in batch_results.items():
        print(f"  Batch size {result['batch_size']}:")
        print(f"    Throughput: {result['throughput_guides_per_min']:.1f} guides/min")
        print(f"    Avg latency: {result['avg_latency_ms']:.2f} ms/guide")
    
    # Check if largest batch meets target (1,000 guides/min)
    largest_batch = batch_results[f'batch_{max([10, 50, 100, 500, 1000])}']
    if largest_batch['throughput_guides_per_min'] >= 1000:
        print(f"  ✓ PASSED: Throughput ≥ 1,000 guides/min target")
    else:
        print(f"  ✗ FAILED: Throughput < 1,000 guides/min target")
    print()
    
    # Benchmark 3: Bootstrap CI computation
    print("Benchmark 3: Bootstrap CI Computation")
    print("-" * 60)
    ci_results = benchmark_bootstrap_ci(analyzer, n_guides=10, n_bootstrap=1000)
    results['bootstrap_ci'] = ci_results
    
    print(f"  Guides: {ci_results['n_guides']}")
    print(f"  Bootstrap resamples: {ci_results['n_bootstrap']}")
    print(f"  Mean time: {ci_results['mean_time_ms']:.2f} ms")
    print(f"  Median time: {ci_results['median_time_ms']:.2f} ms")
    print(f"  Mean CI width: {ci_results['mean_ci_width']:.4f}")
    
    # Check if CI width meets target (<1%)
    if ci_results['mean_ci_width'] < 0.01:
        print(f"  ✓ PASSED: CI width < 1% target")
    else:
        print(f"  ✗ FAILED: CI width ≥ 1% target")
    print()
    
    # Benchmark 4: Score distribution validation
    print("Benchmark 4: Score Distribution Validation")
    print("-" * 60)
    dist_results = validate_score_distribution(analyzer, n_guides=1000)
    results['score_distribution'] = dist_results
    
    print(f"  Guides: {dist_results['n_guides']}")
    print(f"  Mean score: {dist_results['mean_score']:.4f}")
    print(f"  Median score: {dist_results['median_score']:.4f}")
    print(f"  Std dev: {dist_results['std_score']:.4f}")
    print(f"  Range: [{dist_results['min_score']:.4f}, {dist_results['max_score']:.4f}]")
    print(f"  KS test p-value: {dist_results['ks_pvalue']:.4f}")
    print(f"  Shapiro test p-value: {dist_results['shapiro_pvalue']:.4f}")
    
    # Check if distribution is non-uniform (p < 0.05 for KS test)
    if dist_results['ks_pvalue'] < 0.05:
        print(f"  ✓ PASSED: Distribution is non-uniform (informative)")
    else:
        print(f"  ✗ FAILED: Distribution appears uniform (not informative)")
    print()
    
    # Summary
    print("=" * 60)
    print("BENCHMARK SUMMARY")
    print("=" * 60)
    
    all_passed = (
        latency_results['p95_ms'] < 500 and
        largest_batch['throughput_guides_per_min'] >= 1000 and
        ci_results['mean_ci_width'] < 0.01 and
        dist_results['ks_pvalue'] < 0.05
    )
    
    if all_passed:
        print("✓ ALL BENCHMARKS PASSED")
    else:
        print("✗ SOME BENCHMARKS FAILED")
    
    print()
    
    # Save results
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"Results saved to: {output_file}")
    
    return results


def main():
    """CLI interface."""
    parser = argparse.ArgumentParser(
        description='Benchmark Genomic Disruption Analyzer'
    )
    parser.add_argument(
        '--k',
        type=float,
        default=K_STAR,
        help=f'Resolution exponent (default: {K_STAR})'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed (default: 42)'
    )
    parser.add_argument(
        '--n-guides',
        type=int,
        default=100,
        help='Number of guides for quick benchmark (default: 100)'
    )
    parser.add_argument(
        '--full',
        action='store_true',
        help='Run full benchmark suite'
    )
    parser.add_argument(
        '--output',
        help='Output JSON file for results'
    )
    
    args = parser.parse_args()
    
    if args.full:
        run_full_benchmark(k=args.k, seed=args.seed, output_file=args.output)
    else:
        # Quick benchmark
        analyzer = DisruptionAnalyzer(k=args.k, seed=args.seed)
        print(f"Quick benchmark with {args.n_guides} guides...")
        
        start = time.time()
        guides = generate_synthetic_guides(args.n_guides, seed=args.seed)
        results = analyzer.batch_score(guides)
        end = time.time()
        
        elapsed = end - start
        throughput = args.n_guides / elapsed * 60  # guides/min
        
        print(f"Elapsed time: {elapsed:.2f} s")
        print(f"Throughput: {throughput:.1f} guides/min")
        print(f"Avg latency: {(elapsed / args.n_guides) * 1000:.2f} ms/guide")


if __name__ == '__main__':
    main()
