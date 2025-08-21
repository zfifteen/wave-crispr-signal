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
        # Set up environment for proper imports
        env = os.environ.copy()
        env['PYTHONPATH'] = os.getcwd()
        
        result = subprocess.run(
            [sys.executable, test_file], 
            capture_output=True, 
            text=True, 
            timeout=300,
            env=env
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

    # Test Setup 1: Core Framework Tests
    core_tests = [
        ("Z Framework Core", "tests/test_z_framework.py"),
        ("Invariant Features", "tests/test_invariant_features.py"),
        ("Geodesic Bridge", "tests/test_geodesic_bridge.py"),
    ]
    
    # Test Setup 2: Mathematical Analysis Tests
    analysis_tests = [
        ("Topological Analysis", "tests/test_topological_analysis.py"),
    ]
    
    # Test Setup 3: Extended Test Coverage (Application Tests)
    application_tests = [
        ("CRISPR Simple", "tests/test_crispr_simple.py"),
    ]
    
    # Combine all test configurations
    tests = core_tests + analysis_tests + application_tests

    results = []

    # Change to repository root to ensure proper imports
    repo_root = os.path.dirname(os.path.abspath(__file__))
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

    # Summary with test setup breakdown
    print(f"\n{'='*60}")
    print("TEST SUMMARY")
    print(f"{'='*60}")

    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    # Detailed breakdown by test setup
    setup_results = {
        "Core Framework Tests": results[:3],
        "Mathematical Analysis Tests": results[3:4] if len(results) > 3 else [],
        "Application Tests": results[4:] if len(results) > 4 else []
    }
    
    for setup_name, setup_tests in setup_results.items():
        if setup_tests:
            setup_passed = sum(1 for _, success in setup_tests if success)
            setup_total = len(setup_tests)
            print(f"\n{setup_name}: {setup_passed}/{setup_total}")
            for test_name, success in setup_tests:
                status = "✓ PASSED" if success else "❌ FAILED"
                print(f"  {test_name}: {status}")

    print(f"\nOVERALL: {passed}/{total} tests passed")

    if passed == total:
        print("🎉 ALL TESTS PASSED!")
        return 0
    else:
        print("💥 SOME TESTS FAILED!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
