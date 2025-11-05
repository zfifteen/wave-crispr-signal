#!/usr/bin/env python3
"""
Smoke test for bridge feature pipeline.

Fast (<30s) end-to-end test with synthetic data to validate:
- Feature extraction works
- No runtime errors
- Output has expected structure
- Basic sanity checks pass
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import numpy as np
import pandas as pd
from wave_crispr_signal.features import theta_prime, kappa
from wave_crispr_signal.features.curvature import compute_coupled_features


class TestPipelineSmoke:
    """Fast smoke tests for the bridge feature pipeline."""

    def test_single_position_pipeline(self):
        """Test complete pipeline for single position."""
        n = 21  # Typical guide position

        # Extract features
        theta = theta_prime(n, k=0.3)
        kappa_val = kappa(n, mode="discrete")
        coupled = compute_coupled_features(n, k=0.3, coupling="multiplicative")

        # Convert to float for assertions
        theta_f = float(theta)
        kappa_f = float(kappa_val)
        coupled_f = float(coupled)

        # Sanity checks
        assert np.isfinite(theta_f), "theta_prime should be finite"
        assert np.isfinite(kappa_f), "kappa should be finite"
        assert np.isfinite(coupled_f), "coupled feature should be finite"

        assert theta_f > 0, "theta_prime should be positive"
        assert kappa_f > 0, "kappa should be positive"
        assert coupled_f > 0, "coupled feature should be positive"

        # Check coupling is correct
        expected_coupled = kappa_f * theta_f
        assert abs(coupled_f - expected_coupled) < 1e-10

    def test_guide_sequence_pipeline(self):
        """Test pipeline for typical CRISPR guide (20-23 nt)."""
        # Positions for a 20nt guide
        positions = np.arange(1, 21)

        # Extract features for all positions
        theta_vals = theta_prime(positions, k=0.3)
        kappa_vals = kappa(positions, mode="discrete")
        coupled_vals = compute_coupled_features(
            positions, k=0.3, coupling="multiplicative"
        )

        # Convert to float arrays
        theta_arr = np.array([float(t) for t in theta_vals])
        kappa_arr = np.array([float(k) for k in kappa_vals])
        coupled_arr = np.array([float(c) for c in coupled_vals])

        # All values should be finite
        assert np.all(np.isfinite(theta_arr))
        assert np.all(np.isfinite(kappa_arr))
        assert np.all(np.isfinite(coupled_arr))

        # All values should be positive
        assert np.all(theta_arr > 0)
        assert np.all(kappa_arr > 0)
        assert np.all(coupled_arr > 0)

        # Check shapes match
        assert len(theta_arr) == len(positions)
        assert len(kappa_arr) == len(positions)
        assert len(coupled_arr) == len(positions)

    def test_dataframe_integration(self):
        """Test integration with pandas DataFrame."""
        # Create synthetic guide data
        n_guides = 100

        data = {
            "guide_id": [f"guide_{i:03d}" for i in range(n_guides)],
            "sequence": ["ACGT" * 5 for _ in range(n_guides)],  # 20nt guides
            "position": np.arange(1, n_guides + 1),
            "gc_percent": np.random.uniform(0.3, 0.7, n_guides),
        }

        df = pd.DataFrame(data)

        # Extract bridge features
        theta_features = []
        kappa_features = []
        coupled_features = []

        for pos in df["position"]:
            theta = float(theta_prime(pos, k=0.3))
            kappa_val = float(kappa(pos, mode="discrete"))
            coupled = float(
                compute_coupled_features(pos, k=0.3, coupling="multiplicative")
            )

            theta_features.append(theta)
            kappa_features.append(kappa_val)
            coupled_features.append(coupled)

        df["theta_prime"] = theta_features
        df["kappa"] = kappa_features
        df["coupled"] = coupled_features

        # Verify columns exist and have correct types
        assert "theta_prime" in df.columns
        assert "kappa" in df.columns
        assert "coupled" in df.columns

        assert df["theta_prime"].dtype == np.float64
        assert df["kappa"].dtype == np.float64
        assert df["coupled"].dtype == np.float64

        # Check no NaNs
        assert not df["theta_prime"].isna().any()
        assert not df["kappa"].isna().any()
        assert not df["coupled"].isna().any()

        # Check all positive
        assert (df["theta_prime"] > 0).all()
        assert (df["kappa"] > 0).all()
        assert (df["coupled"] > 0).all()

    def test_coupling_modes(self):
        """Test both coupling modes work."""
        positions = np.arange(1, 21)

        # Additive coupling
        coupled_add = compute_coupled_features(positions, k=0.3, coupling="additive")

        # Multiplicative coupling
        coupled_mul = compute_coupled_features(
            positions, k=0.3, coupling="multiplicative"
        )

        # Convert to float arrays
        coupled_add_arr = np.array([float(c) for c in coupled_add])
        coupled_mul_arr = np.array([float(c) for c in coupled_mul])

        # Both should be finite and positive
        assert np.all(np.isfinite(coupled_add_arr))
        assert np.all(np.isfinite(coupled_mul_arr))
        assert np.all(coupled_add_arr > 0)
        assert np.all(coupled_mul_arr > 0)

        # Results should differ
        assert not np.allclose(coupled_add_arr, coupled_mul_arr)

    def test_k_parameter_sweep(self):
        """Test that varying k produces different results."""
        n = 21
        k_values = [0.1, 0.2, 0.3, 0.4, 0.5]

        results = []
        for k in k_values:
            theta = float(theta_prime(n, k=k))
            results.append(theta)

        # All should be finite and positive
        assert all(np.isfinite(r) for r in results)
        assert all(r > 0 for r in results)

        # Results should vary with k
        assert len(set(results)) == len(
            k_values
        ), "Different k values should give different results"

    def test_performance_1000_guides(self):
        """Test that processing 1000 guides completes quickly."""
        import time

        positions = np.arange(1, 1001)

        start = time.time()

        # Extract all features
        theta_vals = theta_prime(positions, k=0.3)
        kappa_vals = kappa(positions, mode="discrete")
        coupled_vals = compute_coupled_features(
            positions, k=0.3, coupling="multiplicative"
        )

        elapsed = time.time() - start

        # Should complete in reasonable time (<10s for 1000 positions)
        assert (
            elapsed < 10.0
        ), f"Processing 1000 guides took {elapsed:.2f}s, should be <10s"

        # Verify results are valid
        assert len(theta_vals) == 1000
        assert len(kappa_vals) == 1000
        assert len(coupled_vals) == 1000


class TestFeatureStatistics:
    """Test statistical properties of extracted features."""

    def test_theta_prime_distribution(self):
        """Test distribution of theta_prime over typical guide range."""
        positions = np.arange(1, 101)
        theta_vals = theta_prime(positions, k=0.3)
        theta_arr = np.array([float(t) for t in theta_vals])

        # Check range
        phi = 1.618033988749895
        assert np.all(theta_arr > 0)
        assert np.all(theta_arr <= phi * 1.01)  # Small tolerance

        # Check mean and std make sense
        mean = np.mean(theta_arr)
        std = np.std(theta_arr)

        assert mean > 0.5  # Should be reasonably sized
        assert std > 0.1  # Should have variation

    def test_kappa_distribution(self):
        """Test distribution of kappa over typical guide range."""
        positions = np.arange(1, 101)
        kappa_vals = kappa(positions, mode="discrete")
        kappa_arr = np.array([float(k) for k in kappa_vals])

        # Check range
        assert np.all(kappa_arr > 0)
        assert np.all(kappa_arr < 1.0)  # Scale is 1/φ ≈ 0.618

        # Check variation
        std = np.std(kappa_arr)
        assert std > 0.01  # Should have some variation

    def test_coupled_feature_correlation(self):
        """Test that coupled features are related to components."""
        positions = np.arange(1, 101)

        theta_vals = np.array([float(theta_prime(p, k=0.3)) for p in positions])
        kappa_vals = np.array([float(kappa(p, mode="discrete")) for p in positions])
        coupled_vals = np.array(
            [
                float(compute_coupled_features(p, k=0.3, coupling="multiplicative"))
                for p in positions
            ]
        )

        # For multiplicative coupling, coupled = theta * kappa
        # Check this relationship holds
        expected_coupled = theta_vals * kappa_vals

        # Should be very close (within numerical precision)
        assert np.allclose(
            coupled_vals, expected_coupled, rtol=1e-10
        ), "Multiplicative coupling should equal theta * kappa"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
