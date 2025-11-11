#!/usr/bin/env python3
"""
Quick Z Framework Validation Demo

Single-command demonstration of >1000x density boost with statistical significance.
Loads minimal synthetic datasets and validates claims in under 2 minutes.

Usage:
    python quick_validation_demo.py

Expected Output:
    - Fast validation of >1000x density boost claims
    - Statistical significance confirmation (p < 0.05)
    - Summary of Z Framework vs baseline performance
"""

import sys
from pathlib import Path


def main():
    """Quick validation demonstration."""
    print("🚀 Z Framework Quick Validation Demo")
    print("=" * 45)
    print("⚡ Fast validation of >1000x density boost claims")
    print("📊 Statistical significance testing (p < 0.05)")
    print("⏱️  Estimated runtime: <2 minutes")
    print()

    # Check if data directory exists
    data_dir = Path(__file__).parent / "data"
    if not data_dir.exists():
        print(
            "❌ Data directory not found. Please run generate_synthetic_data.py first."
        )
        return 1

    # Check if validation script exists
    validate_script = Path(__file__).parent / "validate.py"
    if not validate_script.exists():
        print("❌ Validation script not found.")
        return 1

    print("🔍 Running validation on synthetic datasets...")
    print()

    # Run the main validation script
    import subprocess
    import time

    start_time = time.time()

    try:
        result = subprocess.run(
            [sys.executable, str(validate_script), "--data-dir", str(data_dir)],
            capture_output=True,
            text=True,
            timeout=120,
        )

        end_time = time.time()
        runtime = end_time - start_time

        # Print the output
        print(result.stdout)

        if result.stderr:
            print("⚠️  Warnings/Errors:")
            print(result.stderr)

        # Print runtime
        print(f"⏱️  Runtime: {runtime:.1f} seconds")

        # Check if validation was successful
        if ">1000x density boost claims VALIDATED" in result.stdout:
            print("\n✅ QUICK DEMO RESULT: Z Framework claims VALIDATED")
            print("   🚀 >1000x density boost achieved")
            print("   📊 Statistical significance confirmed (p < 0.05)")
            return_code = 0
        elif ">1000x density boost claims PARTIALLY VALIDATED" in result.stdout:
            print("\n⚠️  QUICK DEMO RESULT: Z Framework claims PARTIALLY VALIDATED")
            return_code = 0
        else:
            print("\n❌ QUICK DEMO RESULT: Z Framework claims NOT VALIDATED")
            return_code = 1

        print("\n" + "=" * 45)
        print("📋 Summary:")
        print("   • Synthetic datasets loaded and analyzed")
        print("   • Baseline vs Z Framework comparison completed")
        print("   • Statistical significance testing performed")
        print("   • >1000x density boost claims evaluated")
        print()
        print("⚠️  DISCLAIMER: Research use only")
        print("   Not for diagnosis, treatment, or clinical decisions")
        print("=" * 45)

        return return_code

    except subprocess.TimeoutExpired:
        print("❌ Validation timed out (>2 minutes)")
        return 1
    except Exception as e:
        print(f"❌ Error running validation: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
