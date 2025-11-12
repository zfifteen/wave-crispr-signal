#!/usr/bin/env python3
"""
Test script for Z Framework implementation

This script tests the core Z Framework functionality to ensure
proper mathematical calculations and convergence validation.
"""

import sys
import os

# Add the repository root to the path for importing modules
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, repo_root)

from scripts.z_framework import (
    ZFrameworkCalculator,
    format_mpmath_for_json,
    format_mpmath_for_display,
)
import json


def test_z_framework():
    """Test the Z Framework implementation with example sequences."""

    print("=" * 60)
    print("Z FRAMEWORK IMPLEMENTATION TEST")
    print("=" * 60)

    # Initialize calculator
    calc = ZFrameworkCalculator(
        precision_dps=30
    )  # Use lower precision for faster testing

    # Test sequences
    test_sequences = [
        ("Simple test", "ATCGATCGATCGATCG"),
        ("PCSK9 sample", "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGG"),
        ("High GC", "GGCCGGCCGGCCGGCCGGCC"),
        ("Low complexity", "AAATTTAAATTTAAAA"),
        ("Balanced", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
    ]

    for name, sequence in test_sequences:
        print(f"\n--- Testing: {name} ({len(sequence)} bp) ---")
        print(f"Sequence: {sequence}")

        try:
            # Basic Z Framework analysis
            z_results = calc.calculate_z_values(sequence)

            print(f"Z mean: {format_mpmath_for_display(z_results['z_mean'])}")
            print(f"Z variance: {format_mpmath_for_display(z_results['z_variance'])}")
            print(f"Z std: {format_mpmath_for_display(z_results['z_std'])}")
            print(f"Converges to φ-1: {z_results['converges_to_phi_conjugate']}")
            print(
                f"Converges to target variance: {z_results['converges_to_target_variance']}"
            )
            print(
                f"φ-1 deviation: {format_mpmath_for_display(z_results['phi_conjugate_convergence'])}"
            )
            print(
                f"Variance deviation: {format_mpmath_for_display(z_results['variance_convergence'])}"
            )

            # Quick falsification test
            print("\nPerforming quick falsification test...")
            falsification = calc.perform_falsification_test(
                sequence, num_perturbations=10, perturbation_rate=0.1
            )
            print(
                f"Convergence rate: {format_mpmath_for_display(falsification['convergence_rate'], 2)}"
            )
            print(f"Falsification passed: {falsification['falsification_passed']}")

        except Exception as e:
            print(f"Error testing {name}: {e}")

    print("\n" + "=" * 60)
    print("Z FRAMEWORK TEST COMPLETED")
    print("=" * 60)


def test_precision():
    """Test high-precision calculations."""
    print("\n--- High Precision Test ---")

    calc = ZFrameworkCalculator(precision_dps=50)
    sequence = "ATGCTGCGGAGACCTGGAGAGA"

    print("Testing with 50 decimal precision...")
    z_results = calc.calculate_z_values(sequence)

    print(f"φ value: {calc.phi}")
    print(f"φ-1 value: {calc.phi_conjugate}")
    print(f"Z mean (high precision): {format_mpmath_for_display(z_results['z_mean'])}")
    print(
        f"Z variance (high precision): {format_mpmath_for_display(z_results['z_variance'])}"
    )


def test_json_serialization():
    """Test JSON serialization of results."""
    print("\n--- JSON Serialization Test ---")

    calc = ZFrameworkCalculator(precision_dps=30)
    sequence = "ATCGATCGATCGATCG"

    z_results = calc.calculate_z_values(sequence)
    density_results = calc.calculate_density_enhancement(sequence)

    # Convert to JSON-serializable format
    results = {"z_framework": z_results, "density_enhancement": density_results}

    json_results = format_mpmath_for_json(results)

    try:
        json_str = json.dumps(json_results, indent=2)
        print("JSON serialization successful!")
        print(f"JSON size: {len(json_str)} characters")
    except Exception as e:
        print(f"JSON serialization failed: {e}")


if __name__ == "__main__":
    test_z_framework()
    test_precision()
    test_json_serialization()
