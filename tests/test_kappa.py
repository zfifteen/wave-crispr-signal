#!/usr/bin/env python3
"""
Test suite for kappa (curvature) implementation.

Validates:
- Scale invariance expectations
- Finite/non-NaN across dataset ranges
- Discrete vs continuous modes
- Coupling with theta_prime
- Edge cases and boundary conditions
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import numpy as np
import mpmath as mp
from wave_crispr_signal.features.curvature import (
    kappa,
    kappa_vectorized,
    compute_coupled_features,
    KappaConfig,
    PHI
)


class TestKappaBounds:
    """Test that kappa stays within valid bounds."""
    
    def test_scalar_positive(self):
        """Test that kappa is always positive."""
        for n in [1, 5, 10, 21, 42, 100]:
            result = kappa(n, mode='discrete')
            assert float(result) > 0, f"kappa({n}) should be positive"
    
    def test_scalar_finite(self):
        """Test that kappa is always finite."""
        for n in [1, 5, 10, 21, 42, 100, 1000]:
            result_discrete = kappa(n, mode='discrete')
            result_continuous = kappa(n, mode='continuous')
            
            assert np.isfinite(float(result_discrete))
            assert np.isfinite(float(result_continuous))
    
    def test_array_bounds(self):
        """Test bounds for array inputs."""
        n_arr = np.arange(1, 101)
        results = kappa(n_arr, mode='discrete')
        
        for result in results.flatten():
            assert float(result) > 0
            assert np.isfinite(float(result))


class TestKappaModes:
    """Test discrete vs continuous calculation modes."""
    
    def test_discrete_mode(self):
        """Test discrete mode: κ(n) = 1 / (1 + (n mod φ))."""
        n = 21
        result = kappa(n, mode='discrete')
        
        # Manual calculation
        n_mod_phi = mp.fmod(mp.mpf(n), PHI)
        scale = mp.mpf(1) / PHI
        expected = scale / (mp.mpf(1) + n_mod_phi)
        
        assert abs(float(result) - float(expected)) < 1e-12
    
    def test_continuous_mode(self):
        """Test continuous mode: κ(n) = 1 / (1 + n/φ)."""
        n = 21
        result = kappa(n, mode='continuous')
        
        # Manual calculation
        scale = mp.mpf(1) / PHI
        expected = scale / (mp.mpf(1) + mp.mpf(n) / PHI)
        
        assert abs(float(result) - float(expected)) < 1e-12
    
    def test_modes_differ(self):
        """Test that discrete and continuous modes give different results."""
        n = 21
        result_discrete = kappa(n, mode='discrete')
        result_continuous = kappa(n, mode='continuous')
        
        # Should be different
        assert abs(float(result_discrete) - float(result_continuous)) > 1e-6


class TestKappaScaleInvariance:
    """Test scale-related properties."""
    
    def test_custom_scale(self):
        """Test with custom scale parameter."""
        n = 21
        scale = 0.5
        config = KappaConfig(mode='discrete', scale=scale)
        result = kappa(n, config=config)
        
        # With scale=0.5, result should be smaller
        result_default = kappa(n, mode='discrete')
        
        # result with scale=0.5 should be roughly half of result with default scale (1/φ ≈ 0.618)
        ratio = float(result) / float(result_default)
        expected_ratio = scale / (1.0 / float(PHI))
        assert abs(ratio - expected_ratio) < 0.1
    
    def test_scale_zero_invalid(self):
        """Test that scale affects magnitude."""
        n = 21
        result1 = kappa(n, mode='discrete')
        
        config = KappaConfig(mode='discrete', scale=2.0)
        result2 = kappa(n, config=config)
        
        # Result should scale proportionally
        assert float(result2) > float(result1)


class TestKappaInputValidation:
    """Test input validation and error handling."""
    
    def test_negative_n_raises(self):
        """Test that negative n raises ValueError."""
        with pytest.raises(ValueError, match="positive"):
            kappa(-1)
    
    def test_zero_n_raises(self):
        """Test that n=0 raises ValueError."""
        with pytest.raises(ValueError, match="positive"):
            kappa(0)
    
    def test_negative_n_array_raises(self):
        """Test that array with negative values raises ValueError."""
        n_arr = np.array([1, 2, -3, 4])
        with pytest.raises(ValueError, match="positive"):
            kappa(n_arr)


class TestKappaConfiguration:
    """Test configuration and parameter handling."""
    
    def test_default_config(self):
        """Test with default configuration."""
        result = kappa(21)
        assert np.isfinite(float(result))
    
    def test_custom_config(self):
        """Test with custom configuration."""
        config = KappaConfig(mode='continuous', scale=0.5, precision_dps=30)
        result = kappa(21, config=config)
        assert np.isfinite(float(result))
    
    def test_mode_override(self):
        """Test that mode parameter overrides config."""
        config = KappaConfig(mode='continuous')
        result1 = kappa(21, config=config)
        result2 = kappa(21, config=config, mode='discrete')
        
        # Results should be different (different modes)
        assert abs(float(result1) - float(result2)) > 1e-6


class TestKappaVectorized:
    """Test vectorized implementation."""
    
    def test_vectorized_vs_precise(self):
        """Compare vectorized with precise implementation."""
        n_arr = np.arange(1, 50)
        
        # Precise version
        precise = kappa(n_arr, mode='discrete')
        
        # Vectorized version
        vectorized = kappa_vectorized(n_arr, mode='discrete')
        
        # Should be close (within float64 precision)
        for i in range(len(n_arr)):
            diff = abs(float(precise[i]) - vectorized[i])
            assert diff < 1e-10, f"Position {n_arr[i]}: diff={diff}"
    
    def test_vectorized_performance(self):
        """Test that vectorized version works on large arrays."""
        n_arr = np.arange(1, 10001)
        result = kappa_vectorized(n_arr, mode='discrete')
        
        assert len(result) == len(n_arr)
        assert all(np.isfinite(result))


class TestKappaDecayProperties:
    """Test decay properties of kappa."""
    
    def test_discrete_mode_bounded(self):
        """Test that discrete mode kappa is bounded."""
        n_values = [1, 10, 100, 1000, 10000]
        results = [float(kappa(n, mode='discrete')) for n in n_values]
        
        # Should all be positive and bounded
        assert all(r > 0 for r in results)
        
        # In discrete mode, values oscillate due to modulo
        # but should stay within a range
        assert all(r < 1.0 for r in results)
    
    def test_continuous_mode_decreasing(self):
        """Test that continuous mode kappa decreases with n."""
        n_values = [1, 10, 100, 1000, 10000]
        results = [float(kappa(n, mode='continuous')) for n in n_values]
        
        # Should be monotonically decreasing
        for i in range(len(results) - 1):
            assert results[i] > results[i+1]


class TestCoupledFeatures:
    """Test coupled κ + θ′ features."""
    
    def test_additive_coupling(self):
        """Test additive coupling: F(n) = κ(n) + θ′(n,k)."""
        n = 21
        k = 0.3
        
        result = compute_coupled_features(n, k=k, coupling='additive')
        
        # Should be finite and positive
        assert np.isfinite(float(result))
        assert float(result) > 0
    
    def test_multiplicative_coupling(self):
        """Test multiplicative coupling: F(n) = κ(n) · θ′(n,k)."""
        n = 21
        k = 0.3
        
        result = compute_coupled_features(n, k=k, coupling='multiplicative')
        
        # Should be finite and positive
        assert np.isfinite(float(result))
        assert float(result) > 0
    
    def test_coupling_modes_differ(self):
        """Test that coupling modes give different results."""
        n = 21
        k = 0.3
        
        result_add = compute_coupled_features(n, k=k, coupling='additive')
        result_mul = compute_coupled_features(n, k=k, coupling='multiplicative')
        
        # Should be different
        assert abs(float(result_add) - float(result_mul)) > 1e-6
    
    def test_coupled_array(self):
        """Test coupled features with array input."""
        n_arr = np.arange(1, 22)
        k = 0.3
        
        results = compute_coupled_features(n_arr, k=k, coupling='multiplicative')
        
        # Should all be finite and positive
        for result in results.flatten():
            assert np.isfinite(float(result))
            assert float(result) > 0


class TestKappaReproducibility:
    """Test reproducibility and determinism."""
    
    def test_repeated_calls_same_result(self):
        """Test that repeated calls give same result."""
        n = 42
        
        results = [float(kappa(n, mode='discrete')) for _ in range(5)]
        
        # All results should be identical
        assert len(set(results)) == 1
    
    def test_array_vs_scalar(self):
        """Test that array and scalar give consistent results."""
        n_scalar = 21
        n_array = np.array([21])
        
        result_scalar = float(kappa(n_scalar, mode='discrete'))
        result_array = float(kappa(n_array, mode='discrete')[0])
        
        assert abs(result_scalar - result_array) < 1e-12


class TestKappaBiologicalRange:
    """Test kappa values in typical biological position ranges."""
    
    def test_guide_length_range(self):
        """Test kappa for typical CRISPR guide positions (20-23 nt)."""
        for n in range(20, 24):
            result = kappa(n, mode='discrete')
            
            # Should be well-defined and positive
            assert float(result) > 0
            assert np.isfinite(float(result))
    
    def test_extended_range(self):
        """Test kappa for extended position range (1-100)."""
        n_arr = np.arange(1, 101)
        results = kappa(n_arr, mode='discrete')
        
        # All should be valid
        for result in results.flatten():
            assert float(result) > 0
            assert np.isfinite(float(result))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
