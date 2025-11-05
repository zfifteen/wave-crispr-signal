#!/usr/bin/env python3
"""
Simple integration test for Molecular Dynamics Z Framework

This script provides a quick validation that the MD framework
integrates properly with the existing Z Framework infrastructure.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "modules"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from molecular_dynamics_framework import (
    MolecularDynamicsZFramework,
    compare_md_vs_base_framework,
)
from z_framework import format_mpmath_for_display


def quick_md_integration_test():
    """Quick test to validate MD framework integration."""
    print("=" * 50)
    print("MOLECULAR DYNAMICS INTEGRATION TEST")
    print("=" * 50)

    # Test basic initialization
    print("\n1. Testing MD Framework Initialization...")
    try:
        md_calc = MolecularDynamicsZFramework(precision_dps=30, md_weight=0.3)
        print("✓ MD Framework initialized successfully")
    except Exception as e:
        print(f"✗ Initialization failed: {e}")
        return False

    # Test basic calculation
    print("\n2. Testing Basic MD Calculation...")
    test_sequence = "ATCGATCGATCGATCG"
    try:
        results = md_calc.calculate_md_z_values(test_sequence)
        print("✓ MD calculation successful")
        print(f"  MD Z mean: {format_mpmath_for_display(results['md_z_mean'])}")
        print(
            f"  Enhancement factor: {format_mpmath_for_display(results['md_enhancement_factor'])}"
        )
    except Exception as e:
        print(f"✗ MD calculation failed: {e}")
        return False

    # Test comparative analysis
    print("\n3. Testing Comparative Analysis...")
    try:
        comparison = compare_md_vs_base_framework([test_sequence], md_weight=0.3)
        if comparison["results"]:
            result = comparison["results"][0]
            print("✓ Comparative analysis successful")
            print(f"  Base Z mean: {result['base_z_mean']:.6f}")
            print(f"  MD Z mean: {result['md_z_mean']:.6f}")
            print(f"  Enhancement: {result['enhancement_factor']:.6f}")
        else:
            print("✗ No comparison results generated")
            return False
    except Exception as e:
        print(f"✗ Comparative analysis failed: {e}")
        return False

    # Test falsification
    print("\n4. Testing MD Parameter Perturbation...")
    try:
        falsification = md_calc.perform_md_parameter_perturbation_test(
            test_sequence, perturbation_range=0.1, num_tests=5
        )
        if "error" not in falsification:
            print("✓ Falsification test successful")
            print(f"  Convergence rate: {falsification['convergence_rate']:.3f}")
            print(f"  Successful tests: {falsification['num_successful_tests']}/5")
        else:
            print(f"✗ Falsification test failed: {falsification['error']}")
            return False
    except Exception as e:
        print(f"✗ Falsification test failed: {e}")
        return False

    print("\n" + "=" * 50)
    print("MD INTEGRATION TEST PASSED")
    print("=" * 50)
    return True


if __name__ == "__main__":
    success = quick_md_integration_test()
    if success:
        print("\n✓ All integration tests passed!")
        print("The MD framework is ready for use.")
    else:
        print("\n✗ Integration tests failed!")
        print("Please check the implementation.")

    sys.exit(0 if success else 1)
