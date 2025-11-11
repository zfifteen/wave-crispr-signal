#!/usr/bin/env python3
"""
Test Suite Runner for Z Framework Geodesic Bridge Implementation

This script runs all the tests to validate the Z Framework implementation
including the new geodesic bridge test.
"""

import subprocess
import sys
import os


def run_test(test_name, test_file):
    """Run a single test and return result."""
    print(f"\n{'='*60}")
    print(f"RUNNING: {test_name}")
    print(f"{'='*60}")

    try:
        result = subprocess.run(
            [sys.executable, test_file], capture_output=True, text=True, timeout=300
        )

        if result.returncode == 0:
            print(f"✓ {test_name} PASSED")
            return True
        else:
            print(f"❌ {test_name} FAILED")
            if result.stderr:
                print("STDERR:", result.stderr)
            return False
    except subprocess.TimeoutExpired:
        print(f"❌ {test_name} TIMED OUT")
        return False
    except Exception as e:
        print(f"❌ {test_name} ERROR: {e}")
        return False


def main():
    """Run all tests."""
    print("Z FRAMEWORK TEST SUITE")
    print("=" * 60)

    # Test definitions - updated paths for new structure
    tests = [
        ("Z Framework Core", "tests/test_z_framework.py"),
        ("Invariant Features", "tests/test_invariant_features.py"),
        ("Geodesic Bridge", "tests/test_geodesic_bridge.py"),
        ("Bin-Resonance Test", "tests/test_bin_resonance.py"),
        ("ORCS Physical Z-Metrics", "tests/test_orcs_physical_z_metrics.py"),
        ("MRI Z5D Analysis", "tests/test_mri_z5d_analysis.py"),
    ]

    results = []

    # Change to repository root to ensure proper imports
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    original_cwd = os.getcwd()
    os.chdir(repo_root)

    try:
        # Run each test
        for test_name, test_file in tests:
            if os.path.exists(test_file):
                success = run_test(test_name, test_file)
                results.append((test_name, success))
            else:
                print(f"⚠ {test_name} - File not found: {test_file}")
                results.append((test_name, False))
    finally:
        # Restore original working directory
        os.chdir(original_cwd)

    # Summary
    print(f"\n{'='*60}")
    print("TEST SUMMARY")
    print(f"{'='*60}")

    passed = sum(1 for _, success in results if success)
    total = len(results)

    for test_name, success in results:
        status = "✓ PASSED" if success else "❌ FAILED"
        print(f"{test_name}: {status}")

    print(f"\nOVERALL: {passed}/{total} tests passed")

    if passed == total:
        print("🎉 ALL TESTS PASSED!")
        return 0
    else:
        print("💥 SOME TESTS FAILED!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
