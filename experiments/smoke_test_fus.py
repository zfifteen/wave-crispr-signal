#!/usr/bin/env python3
"""
Smoke Test: Focused Ultrasound MVE

Quick validation test that runs the FUS experiment with minimal parameters
to ensure it completes successfully in CI environment.

Target runtime: <5 seconds
"""

import subprocess
import sys
import time
from pathlib import Path

def run_smoke_test():
    """Run FUS experiment with minimal parameters for CI."""
    print("=" * 80)
    print("SMOKE TEST: Focused Ultrasound MVE")
    print("=" * 80)
    
    # Change to repository root
    repo_root = Path(__file__).parent.parent
    
    # Run experiment with minimal parameters
    import tempfile
    temp_dir = tempfile.mkdtemp(prefix='fus_smoke_test_')
    
    cmd = [
        sys.executable,
        str(repo_root / "experiments" / "focused_ultrasound_mve.py"),
        "--seed", "42",
        "--bootstrap", "100",  # Reduced for speed
        "--permutation", "100",  # Reduced for speed
        "--splits", "single",
        "--domain", "discrete",
        "--k-parameter", "0.3",
        "--grid-size", "50",  # Reduced grid size
        "--n-trials", "100",  # Reduced trials
        "--output-dir", temp_dir,
        "--smoke-test"  # Enable smoke test mode
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    print()
    
    start_time = time.time()
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=30  # 30 second timeout for CI
        )
        
        elapsed_time = time.time() - start_time
        
        # Print output
        print(result.stdout)
        
        # Check for success indicators
        if "HYPOTHESIS SUPPORTED" in result.stdout or "HYPOTHESIS NOT SUPPORTED" in result.stdout:
            print("\n" + "=" * 80)
            print("✓ SMOKE TEST PASSED")
            print(f"  Runtime: {elapsed_time:.2f} seconds")
            print("=" * 80)
            return 0
        else:
            print("\n" + "=" * 80)
            print("✗ SMOKE TEST FAILED: No hypothesis conclusion found")
            print("=" * 80)
            return 1
            
    except subprocess.CalledProcessError as e:
        elapsed_time = time.time() - start_time
        print(f"✗ SMOKE TEST FAILED: Command returned non-zero exit code")
        print(f"  Exit code: {e.returncode}")
        print(f"  Runtime: {elapsed_time:.2f} seconds")
        print("\nSTDOUT:")
        print(e.stdout)
        print("\nSTDERR:")
        print(e.stderr)
        return 1
        
    except subprocess.TimeoutExpired:
        print("✗ SMOKE TEST FAILED: Timeout exceeded (>30 seconds)")
        return 1
        
    except Exception as e:
        print(f"✗ SMOKE TEST FAILED: Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(run_smoke_test())
