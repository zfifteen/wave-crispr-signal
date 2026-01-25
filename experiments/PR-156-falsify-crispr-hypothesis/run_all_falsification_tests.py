#!/usr/bin/env python3
"""
Combined Runner for PR-156 Falsification Experiments

This script runs both Hypothesis 1 and Hypothesis 2 falsification tests
and generates a combined report.

Usage:
  python run_all_falsification_tests.py [--smoke] [--output-dir DIR]
"""

import sys
import os
import argparse
import json
import subprocess
from datetime import datetime
from pathlib import Path

# Add parent directories to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))


def run_hypothesis1(args, smoke=False):
    """Run Hypothesis 1 falsification test."""
    print("=" * 70)
    print("RUNNING HYPOTHESIS 1 TEST")
    print("=" * 70)
    
    script_path = os.path.join(
        os.path.dirname(__file__), 
        "falsify_hypothesis1.py"
    )
    
    cmd = [
        sys.executable,
        script_path,
        "--seed", str(args.seed),
        "--k-parameter", str(args.k_parameter),
        "--output-dir", args.output_dir,
    ]
    
    if smoke:
        cmd.extend([
            "--n-samples", "20",
            "--n-bootstrap", "50",
            "--n-folds", "3",
        ])
    else:
        cmd.extend([
            "--n-samples", str(args.n_samples),
            "--n-bootstrap", str(args.n_bootstrap),
            "--n-folds", str(args.n_folds),
        ])
    
    result = subprocess.run(cmd, capture_output=False)
    
    return result.returncode


def run_hypothesis2(args, smoke=False):
    """Run Hypothesis 2 falsification test."""
    print("\n")
    print("=" * 70)
    print("RUNNING HYPOTHESIS 2 TEST")
    print("=" * 70)
    
    script_path = os.path.join(
        os.path.dirname(__file__), 
        "falsify_hypothesis2.py"
    )
    
    cmd = [
        sys.executable,
        script_path,
        "--seed", str(args.seed),
        "--k-parameter", str(args.k_parameter),
        "--output-dir", args.output_dir,
    ]
    
    if smoke:
        cmd.extend([
            "--min-length", "20",
            "--max-length", "40",
            "--length-step", "10",
            "--n-per-length", "5",
        ])
    else:
        cmd.extend([
            "--min-length", str(args.min_length),
            "--max-length", str(args.max_length),
            "--length-step", str(args.length_step),
            "--n-per-length", str(args.n_per_length),
        ])
    
    if args.test_motif:
        cmd.append("--test-motif")
    
    result = subprocess.run(cmd, capture_output=False)
    
    return result.returncode


def generate_combined_report(args):
    """Generate combined report from both hypothesis results."""
    print("\n")
    print("=" * 70)
    print("GENERATING COMBINED REPORT")
    print("=" * 70)
    
    # Read results
    output_dir = os.path.join(
        os.path.dirname(__file__), "..", "..", args.output_dir
    )
    
    h1_file = os.path.join(output_dir, "hypothesis1_results.json")
    h2_file = os.path.join(output_dir, "hypothesis2_results.json")
    
    try:
        with open(h1_file, 'r') as f:
            h1_results = json.load(f)
    except FileNotFoundError:
        print(f"Warning: Hypothesis 1 results not found at {h1_file}")
        h1_results = None
    
    try:
        with open(h2_file, 'r') as f:
            h2_results = json.load(f)
    except FileNotFoundError:
        print(f"Warning: Hypothesis 2 results not found at {h2_file}")
        h2_results = None
    
    # Generate combined report
    report = {
        "experiment_id": "PR-156-falsify-crispr-hypothesis",
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "seed": args.seed,
            "k_parameter": args.k_parameter,
        },
        "hypothesis_1": h1_results,
        "hypothesis_2": h2_results,
    }
    
    # Summary
    summary = {
        "hypothesis_1_falsified": h1_results["falsification"]["falsified"] if h1_results else None,
        "hypothesis_2_falsified": h2_results["falsification"]["falsified"] if h2_results else None,
    }
    
    # Overall conclusion
    if summary["hypothesis_1_falsified"] is not None and summary["hypothesis_2_falsified"] is not None:
        if summary["hypothesis_1_falsified"] or summary["hypothesis_2_falsified"]:
            conclusion = "At least one hypothesis was falsified. Review individual results."
        else:
            conclusion = "Both hypotheses are supported by the data (not falsified)."
    else:
        conclusion = "Incomplete results. Check individual hypothesis outputs."
    
    report["summary"] = summary
    report["conclusion"] = conclusion
    
    # Save combined report
    combined_file = os.path.join(output_dir, "combined_report.json")
    with open(combined_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\nCombined report saved to: {combined_file}")
    print()
    print("SUMMARY:")
    print(f"  Hypothesis 1 falsified: {summary['hypothesis_1_falsified']}")
    print(f"  Hypothesis 2 falsified: {summary['hypothesis_2_falsified']}")
    print(f"  Conclusion: {conclusion}")
    print()
    
    return report


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Run all PR-156 falsification tests"
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run smoke tests with reduced parameters"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--k-parameter",
        type=float,
        default=0.3,
        help="Resolution exponent k"
    )
    
    # Hypothesis 1 parameters
    parser.add_argument(
        "--n-samples",
        type=int,
        default=100,
        help="Number of samples (Hypothesis 1)"
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=1000,
        help="Number of bootstrap resamples (Hypothesis 1)"
    )
    parser.add_argument(
        "--n-folds",
        type=int,
        default=10,
        help="Number of CV folds (Hypothesis 1)"
    )
    
    # Hypothesis 2 parameters
    parser.add_argument(
        "--min-length",
        type=int,
        default=20,
        help="Minimum sequence length (Hypothesis 2)"
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=100,
        help="Maximum sequence length (Hypothesis 2)"
    )
    parser.add_argument(
        "--length-step",
        type=int,
        default=10,
        help="Length step size (Hypothesis 2)"
    )
    parser.add_argument(
        "--n-per-length",
        type=int,
        default=20,
        help="Sequences per length (Hypothesis 2)"
    )
    parser.add_argument(
        "--test-motif",
        action="store_true",
        help="Test with preserved motif (Hypothesis 2)"
    )
    
    # Output
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/PR-156-falsify-combined",
        help="Output directory for all results"
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("PR-156 Falsification Experiment Suite")
    print("=" * 70)
    print(f"Mode: {'SMOKE TEST' if args.smoke else 'FULL RUN'}")
    print(f"Seed: {args.seed}")
    print(f"k parameter: {args.k_parameter}")
    print()
    
    # Run both hypotheses
    exit_code_h1 = run_hypothesis1(args, smoke=args.smoke)
    exit_code_h2 = run_hypothesis2(args, smoke=args.smoke)
    
    # Generate combined report
    generate_combined_report(args)
    
    # Exit with non-zero if any test failed
    if exit_code_h1 != 0 or exit_code_h2 != 0:
        print("WARNING: One or more tests exited with non-zero code")
        return max(exit_code_h1, exit_code_h2)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
