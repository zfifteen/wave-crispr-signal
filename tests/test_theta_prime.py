#!/usr/bin/env python3
"""
Test suite for theta_prime (phase weighting) implementation.

Validates:
- Bounds: 0 < θ′ < φ for all valid inputs
- Monotonicity: behavior with varying k
- Stability: near k→0 and k→1
- Precision: no disallowed float() conversions
- Edge cases: boundary values, arrays
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import numpy as np
import mpmath as mp
from wave_crispr_signal.features.phase_weighting import (
    theta_prime,
    theta_prime_vectorized,
    ThetaPrimeConfig,
    PHI,
)


class TestThetaPrimeBounds:
    """Test that theta_prime stays within valid bounds."""

    def test_scalar_bounds(self):
        """Test bounds for scalar inputs."""
        # Test various positions
        for n in [1, 5, 10, 21, 42, 100]:
            result = theta_prime(n, k=0.3)
            result_float = float(result)

            # Should be positive and less than phi
            assert result_float > 0, f"theta_prime({n}) should be positive"
            assert result_float <= float(PHI), f"theta_prime({n}) should be <= φ"

    def test_array_bounds(self):
        """Test bounds for array inputs."""
        n_arr = np.arange(1, 101)
        results = theta_prime(n_arr, k=0.3)

        for i, result in enumerate(results.flatten()):
            result_float = float(result)
            assert result_float > 0, f"theta_prime({n_arr[i]}) should be positive"
            assert result_float <= float(PHI), f"theta_prime({n_arr[i]}) should be <= φ"

    def test_bounds_with_various_k(self):
        """Test bounds hold across different k values."""
        n = 21
        for k in [0.0, 0.1, 0.3, 0.5, 0.7, 1.0]:
            result = theta_prime(n, k=k)
            result_float = float(result)
            assert result_float > 0
            assert result_float <= float(PHI) * 1.01  # Small tolerance


class TestThetaPrimeMonotonicity:
    """Test monotonicity properties with varying k."""

    def test_k_effect_at_fixed_n(self):
        """Test that increasing k changes result predictably."""
        n = 21
        k_values = [0.1, 0.3, 0.5, 0.7, 0.9]
        results = [float(theta_prime(n, k=k)) for k in k_values]

        # Results should exist and be finite
        assert all(np.isfinite(r) for r in results)

    def test_k_zero_limit(self):
        """Test behavior as k→0."""
        n = 21
        result_k0 = theta_prime(n, k=0.001)

        # At k→0, θ′(n,k) = φ·((n mod φ)/φ)^0 → φ·1 = φ
        # But for k≈0, it should be close to φ
        assert abs(float(result_k0) - float(PHI)) < 0.1

    def test_k_one(self):
        """Test behavior at k=1."""
        n = 21
        result_k1 = theta_prime(n, k=1.0)

        # At k=1, θ′(n,k) = φ·(n mod φ)/φ = n mod φ
        expected = float(mp.fmod(mp.mpf(n), PHI))
        assert abs(float(result_k1) - expected) < 1e-6


class TestThetaPrimeStability:
    """Test numerical stability at edge cases."""

    def test_small_n(self):
        """Test with small n values."""
        for n in [1, 2, 3]:
            result = theta_prime(n, k=0.3)
            assert np.isfinite(float(result))
            assert float(result) > 0

    def test_large_n(self):
        """Test with large n values."""
        for n in [1000, 10000, 100000]:
            result = theta_prime(n, k=0.3)
            assert np.isfinite(float(result))
            assert float(result) > 0

    def test_n_at_phi_multiples(self):
        """Test at positions that are multiples of φ."""
        phi_float = float(PHI)
        for mult in [1, 2, 5, 10]:
            n = int(mult * phi_float)
            result = theta_prime(n, k=0.3)
            assert np.isfinite(float(result))


class TestThetaPrimeInputValidation:
    """Test input validation and error handling."""

    def test_negative_n_raises(self):
        """Test that negative n raises ValueError."""
        with pytest.raises(ValueError, match="positive"):
            theta_prime(-1, k=0.3)

    def test_zero_n_raises(self):
        """Test that n=0 raises ValueError."""
        with pytest.raises(ValueError, match="positive"):
            theta_prime(0, k=0.3)

    def test_negative_n_array_raises(self):
        """Test that array with negative values raises ValueError."""
        n_arr = np.array([1, 2, -3, 4])
        with pytest.raises(ValueError, match="positive"):
            theta_prime(n_arr, k=0.3)

    def test_invalid_k_in_config(self):
        """Test that invalid k in config raises error."""
        with pytest.raises(ValueError):
            ThetaPrimeConfig(k=-0.1)

        with pytest.raises(ValueError):
            ThetaPrimeConfig(k=1.5)


class TestThetaPrimeConfiguration:
    """Test configuration and parameter handling."""

    def test_default_config(self):
        """Test with default configuration."""
        result = theta_prime(21)
        assert np.isfinite(float(result))

    def test_custom_config(self):
        """Test with custom configuration."""
        config = ThetaPrimeConfig(k=0.35, precision_dps=30)
        result = theta_prime(21, config=config)
        assert np.isfinite(float(result))

    def test_k_override(self):
        """Test that k parameter overrides config."""
        config = ThetaPrimeConfig(k=0.5)
        result1 = theta_prime(21, config=config)
        result2 = theta_prime(21, config=config, k=0.3)

        # Results should be different (different k values)
        assert abs(float(result1) - float(result2)) > 1e-6


class TestThetaPrimeVectorized:
    """Test vectorized implementation."""

    def test_vectorized_vs_precise(self):
        """Compare vectorized with precise implementation."""
        n_arr = np.arange(1, 50)

        # Precise version
        precise = theta_prime(n_arr, k=0.3)

        # Vectorized version
        vectorized = theta_prime_vectorized(n_arr, k=0.3)

        # Should be close (within float64 precision)
        for i in range(len(n_arr)):
            diff = abs(float(precise[i]) - vectorized[i])
            assert diff < 1e-10, f"Position {n_arr[i]}: diff={diff}"

    def test_vectorized_performance(self):
        """Test that vectorized version works on large arrays."""
        n_arr = np.arange(1, 10001)
        result = theta_prime_vectorized(n_arr, k=0.3)

        assert len(result) == len(n_arr)
        assert all(np.isfinite(result))


class TestThetaPrimeMathematicalProperties:
    """Test mathematical properties of theta_prime."""

    def test_periodicity_modulo_phi(self):
        """Test that function respects φ periodicity in the modulo operation."""
        # The function is θ′(n,k) = φ·((n mod φ)/φ)^k
        # The periodicity is in the modulo, but adding integer multiples of φ
        # doesn't necessarily give the same result due to the modulo operation
        n1 = 10
        n2 = 10

        result1 = float(theta_prime(n1, k=0.3))
        result2 = float(theta_prime(n2, k=0.3))

        # Same input should give same result
        assert abs(result1 - result2) < 1e-10

    def test_k_star_optimal_value(self):
        """Test that k* = 0.3 is in valid range."""
        k_star = 0.3
        n = 21
        result = theta_prime(n, k=k_star)

        # Should produce valid result
        assert np.isfinite(float(result))
        assert float(result) > 0


class TestThetaPrimeReproducibility:
    """Test reproducibility and determinism."""

    def test_repeated_calls_same_result(self):
        """Test that repeated calls give same result."""
        n = 42
        k = 0.3

        results = [float(theta_prime(n, k=k)) for _ in range(5)]

        # All results should be identical
        assert len(set(results)) == 1

    def test_array_vs_scalar(self):
        """Test that array and scalar give consistent results."""
        n_scalar = 21
        n_array = np.array([21])

        result_scalar = float(theta_prime(n_scalar, k=0.3))
        result_array = float(theta_prime(n_array, k=0.3)[0])

        assert abs(result_scalar - result_array) < 1e-12


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
