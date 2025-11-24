#!/usr/bin/env python3
"""
Smoke test for prime approximation falsification experiment.

This test runs quickly (<5s) to validate basic functionality for CI.
Uses smaller N and fewer bootstrap samples.
"""

import sys
import time
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from falsify_prime_density import (
    mu, riemann_R, count_primes_exact, 
    calculate_density_boost
)
import mpmath as mp

def test_moebius_function():
    """Test Möbius function μ(n)."""
    assert mu(1) == 1, "μ(1) should be 1"
    assert mu(2) == -1, "μ(2) should be -1 (one prime)"
    assert mu(4) == 0, "μ(4) should be 0 (squared factor)"
    assert mu(6) == 1, "μ(6) should be 1 (two distinct primes)"
    assert mu(30) == -1, "μ(30) should be -1 (three distinct primes)"
    print("✓ Möbius function tests passed")


def test_prime_counting():
    """Test exact prime counting."""
    assert count_primes_exact(10) == 4, "π(10) = 4"
    assert count_primes_exact(100) == 25, "π(100) = 25"
    assert count_primes_exact(1000) == 168, "π(1000) = 168"
    print("✓ Prime counting tests passed")


def test_riemann_R():
    """Test Riemann R(x) calculation."""
    # R(100) should be close to π(100) = 25
    R_100 = riemann_R(mp.mpf(100), max_terms=50)
    assert 20 < float(R_100) < 30, f"R(100) = {R_100} should be near 25"
    
    # R(1000) should be close to π(1000) = 168
    R_1000 = riemann_R(mp.mpf(1000), max_terms=50)
    assert 160 < float(R_1000) < 175, f"R(1000) = {R_1000} should be near 168"
    
    print("✓ Riemann R(x) tests passed")


def test_density_calculation_smoke():
    """Smoke test for density boost calculation (small N, few bootstrap)."""
    print("\nRunning smoke test with N=10000, bootstrap=10...")
    start = time.time()
    
    results = calculate_density_boost(
        N=10000,
        num_bootstrap=10,
        seed=42
    )
    
    elapsed = time.time() - start
    print(f"  Elapsed time: {elapsed:.2f}s")
    
    # Check basic result structure
    assert 'pi_N_exact' in results, "Missing pi_N_exact"
    assert 'R_N' in results, "Missing R_N"
    assert 'actual_boost_pct' in results, "Missing actual_boost_pct"
    assert 'bootstrap_ci_lower' in results, "Missing bootstrap_ci_lower"
    assert 'hypothesis_falsified' in results, "Missing hypothesis_falsified"
    
    # Sanity checks
    assert results['pi_N_exact'] > 0, "π(N) should be positive"
    assert results['R_N'] > 0, "R(N) should be positive"
    assert 50 < results['actual_boost_pct'] < 200, "Boost should be reasonable"
    
    # The hypothesis should be falsified (actual boost << 210%)
    if results['actual_boost_pct'] < 150:
        print(f"  ✓ Actual boost ({results['actual_boost_pct']:.1f}%) << claimed (210%)")
    
    print("✓ Density calculation smoke test passed")


def main():
    """Run all smoke tests."""
    print("=" * 60)
    print("Prime Approximation Falsification - Smoke Test")
    print("=" * 60)
    
    start_total = time.time()
    
    try:
        test_moebius_function()
        test_prime_counting()
        test_riemann_R()
        test_density_calculation_smoke()
        
        total_elapsed = time.time() - start_total
        
        print("\n" + "=" * 60)
        print(f"All smoke tests passed in {total_elapsed:.2f}s")
        
        if total_elapsed < 5.0:
            print("✓ CI compliance: runtime < 5s")
        else:
            print(f"⚠ Warning: runtime {total_elapsed:.2f}s exceeds 5s target")
        
        print("=" * 60)
        return 0
        
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
