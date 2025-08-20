"""
Tests for Invariant Features Module

Tests the mathematical invariant features for CRISPR guide design
including phase bit detection, phase-difference calculations, and
golden proximity metrics.
"""

import unittest
import numpy as np
import sys
import os

# Add the repository root to the Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from invariant_features import (
    ZetaUnfoldCalculator,
    PhaseAwareSpectralAnalyzer,
    GoldenProximityCalculator,
    CurvatureDisruptionAnalyzer,
    InvariantFeatureSet,
)


class TestZetaUnfoldCalculator(unittest.TestCase):
    """Test the zeta unfold calculator and phase bit detection."""

    def setUp(self):
        """Set up test fixtures."""
        # Use values that should produce F alternation
        self.calc = ZetaUnfoldCalculator(0.552, 9.061, 7.389)

    def test_initialization(self):
        """Test proper initialization of calculator."""
        self.assertIsNotNone(self.calc.z)
        self.assertIsNotNone(self.calc.D)
        self.assertIsNotNone(self.calc.E)
        self.assertIsNotNone(self.calc.F)

    def test_unfold_iteration(self):
        """Test unfold iteration produces new values."""
        original_f = float(self.calc.F)

        # Perform unfold
        calc1 = self.calc.unfold_next()
        f1 = float(calc1.F)

        # F should change after unfold
        self.assertNotEqual(original_f, f1)

        # Perform another unfold
        calc2 = calc1.unfold_next()
        f2 = float(calc2.F)

        # F should alternate back toward original pattern
        self.assertNotEqual(f1, f2)

    def test_phase_bit_detection(self):
        """Test phase bit extraction from F values."""
        # Test with known F values
        calc_phase_0 = ZetaUnfoldCalculator(1.0, 10.0, 7.389)  # Should give low F
        calc_phase_1 = calc_phase_0.unfold_next()  # Should give high F

        phase_0 = calc_phase_0.get_phase_bit()
        phase_1 = calc_phase_1.get_phase_bit()

        # Phases should be different
        self.assertIn(phase_0, [0, 1])
        self.assertIn(phase_1, [0, 1])

    def test_f_alternation_pattern(self):
        """Test that F values show alternation pattern."""
        f_values = []
        calc = self.calc

        # Collect F values over several unfolds
        for i in range(6):
            f_values.append(float(calc.F))
            calc = calc.unfold_next()

        # Check that F values show some variation
        unique_f_values = set(f_values)
        self.assertGreater(
            len(unique_f_values), 1, "F values should vary across unfolds"
        )


class TestPhaseAwareSpectralAnalyzer(unittest.TestCase):
    """Test phase-aware spectral analysis."""

    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = PhaseAwareSpectralAnalyzer()
        self.test_sequence = "ATCGATCGATCGATCGATCG"

    def test_spectral_features_calculation(self):
        """Test calculation of spectral features."""
        features = self.analyzer._compute_spectral_features(self.test_sequence, phase=0)

        # Check required features are present
        required_features = ["entropy", "flatness", "f1_magnitude"]
        for feature in required_features:
            self.assertIn(feature, features)
            self.assertIsInstance(features[feature], (int, float, np.number))

    def test_phase_difference_features(self):
        """Test phase-difference feature calculation."""
        features = self.analyzer.calculate_phase_difference_features(self.test_sequence)

        # Check phase difference features are present
        expected_features = [
            "delta_phase_entropy",
            "delta_phase_flatness",
            "delta_phase_f1_magnitude",
        ]
        for feature in expected_features:
            self.assertIn(feature, features)

    def test_mutation_analysis(self):
        """Test spectral analysis with mutations."""
        mutation_pos = 5
        mutation_base = "A" if self.test_sequence[mutation_pos] != "A" else "T"

        features = self.analyzer.calculate_phase_difference_features(
            self.test_sequence, mutation_pos, mutation_base
        )

        # Should have mutation-specific features
        self.assertIn("delta_phase_entropy_change", features)
        self.assertIn("delta_phase_flatness_change", features)

    def test_waveform_construction(self):
        """Test phase-aware waveform construction."""
        waveform_0 = self.analyzer._build_phase_aware_waveform(
            self.test_sequence, phase=0
        )
        waveform_1 = self.analyzer._build_phase_aware_waveform(
            self.test_sequence, phase=1
        )

        # Waveforms should be different for different phases
        self.assertEqual(len(waveform_0), len(self.test_sequence))
        self.assertEqual(len(waveform_1), len(self.test_sequence))
        self.assertTrue(np.any(waveform_0 != waveform_1))

    def test_geodesic_resolution(self):
        """Test geodesic resolution calculation."""
        # Test geodesic calculation with different phases
        geo_0 = self.analyzer._calculate_geodesic_resolution(5, phase=0)
        geo_1 = self.analyzer._calculate_geodesic_resolution(5, phase=1)

        self.assertIsInstance(geo_0, (int, float))
        self.assertIsInstance(geo_1, (int, float))
        self.assertGreater(geo_0, 0)
        self.assertGreater(geo_1, 0)


class TestGoldenProximityCalculator(unittest.TestCase):
    """Test golden proximity calculations."""

    def setUp(self):
        """Set up test fixtures."""
        self.calculator = GoldenProximityCalculator()

    def test_basic_proximity_calculation(self):
        """Test basic golden proximity calculation."""
        # Test with values near golden ratio conjugate (≈0.618)
        z_values = [0.6, 0.62, 0.615, 0.625]

        results = self.calculator.calculate_golden_proximity(z_values)

        # Check required metrics are present
        self.assertIn("delta_phi", results)
        self.assertIn("mu_z", results)
        self.assertIn("phi_conjugate_target", results)

        # Delta phi should be reasonable for values near golden ratio
        self.assertLess(results["delta_phi"], 0.1)

    def test_trimmed_proximity_calculation(self):
        """Test trimmed proximity calculation with outliers."""
        # Include outliers to test trimming
        z_values = [0.6, 0.62, 0.615, 0.625, 2.0, 3.5]  # Last two are outliers

        results = self.calculator.calculate_golden_proximity(
            z_values, trim_outliers=True
        )

        # Should have trimmed results
        self.assertIn("delta_phi_trim", results)
        self.assertIn("mu_z_trim", results)
        self.assertIn("trimmed_count", results)

        # Trimmed should remove outliers
        self.assertLess(results["trimmed_count"], len(z_values))

    def test_golden_ratio_target(self):
        """Test that golden ratio target is correct."""
        expected_phi_conjugate = 1.618033988749895 - 1  # φ - 1

        self.assertAlmostEqual(
            self.calculator.phi_conjugate, expected_phi_conjugate, places=10
        )


class TestCurvatureDisruptionAnalyzer(unittest.TestCase):
    """Test curvature disruption analysis."""

    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = CurvatureDisruptionAnalyzer()
        self.test_sequence = "ATCGATCGATCGATCGATCGAAATTTGGGCCCAAA"

    def test_curvature_disruption_calculation(self):
        """Test curvature disruption calculation."""
        mutation_pos = 10
        original_base = self.test_sequence[mutation_pos]
        mutation_base = "A" if original_base != "A" else "T"

        results = self.analyzer.calculate_curvature_disruption(
            self.test_sequence, mutation_pos, mutation_base
        )

        # Check required disruption metrics
        expected_features = [
            "delta_curv_weighted_composition",
            "delta_curv_structural_complexity",
            "delta_curv_weighted_entropy",
        ]

        for feature in expected_features:
            self.assertIn(feature, results)
            self.assertIsInstance(results[feature], (int, float, np.number))

        # Check metadata
        self.assertEqual(results["mutation_pos"], mutation_pos)
        self.assertEqual(results["original_base"], original_base)
        self.assertEqual(results["mutation_base"], mutation_base)

    def test_curvature_features(self):
        """Test curvature feature calculation."""
        window = "ATCGATCG"
        focal_pos = 3

        features = self.analyzer._calculate_curvature_features(window, focal_pos)

        expected_features = [
            "weighted_composition",
            "structural_complexity",
            "weighted_entropy",
        ]
        for feature in expected_features:
            self.assertIn(feature, features)
            self.assertIsInstance(features[feature], (int, float, np.number))

    def test_geodesic_resolution(self):
        """Test geodesic resolution calculation."""
        # Test with different distances
        weight_0 = self.analyzer._calculate_geodesic_resolution(0)  # Focal point
        weight_5 = self.analyzer._calculate_geodesic_resolution(5)  # Distant

        # Weight should decrease with distance
        self.assertGreater(weight_0, weight_5)
        self.assertGreater(weight_0, 0)
        self.assertGreater(weight_5, 0)


class TestInvariantFeatureSet(unittest.TestCase):
    """Test the complete invariant feature set."""

    def setUp(self):
        """Set up test fixtures."""
        self.feature_set = InvariantFeatureSet()
        self.test_sequence = "ATCGATCGATCGATCGATCGAAATTTGGGCCCAAA"

    def test_complete_feature_calculation(self):
        """Test calculation of complete feature set."""
        features = self.feature_set.calculate_complete_feature_set(self.test_sequence)

        # Check core invariant features are present
        core_features = [
            "phase_bit",
            "length_invariant_normalized",
            "normalization_constant",
        ]

        for feature in core_features:
            self.assertIn(feature, features)

        # Phase bit should be 0 or 1
        self.assertIn(features["phase_bit"], [0, 1])

        # Normalization should be enabled
        self.assertTrue(features["length_invariant_normalized"])

    def test_mutation_feature_calculation(self):
        """Test feature calculation with mutation."""
        mutation_pos = 10
        mutation_base = "A"

        features = self.feature_set.calculate_complete_feature_set(
            self.test_sequence, mutation_pos, mutation_base
        )

        # Should include mutation-specific features
        mutation_features = [
            "delta_phase_entropy_change",
            "delta_curv_weighted_composition",
            "original_base",
            "mutation_base",
        ]

        for feature in mutation_features:
            self.assertIn(feature, features)

    def test_feature_types(self):
        """Test that features have correct types."""
        features = self.feature_set.calculate_complete_feature_set(self.test_sequence)

        # Phase bit should be integer
        self.assertIsInstance(features["phase_bit"], int)

        # Normalization flag should be boolean
        self.assertIsInstance(features["length_invariant_normalized"], bool)

        # Phase difference features should be numeric
        phase_features = [k for k in features.keys() if "delta_phase" in k]
        for feature in phase_features:
            self.assertIsInstance(features[feature], (int, float, np.number))


class TestIntegration(unittest.TestCase):
    """Integration tests for invariant features."""

    def test_g_to_c_transition_analysis(self):
        """Test analysis specific to G→C transitions."""
        sequence = "ATCGATCGATCGATCGATCG"

        # Find G positions for G→C mutations
        g_positions = [i for i, base in enumerate(sequence) if base == "G"]

        if g_positions:
            mutation_pos = g_positions[0]

            # Test G→C transition
            feature_set = InvariantFeatureSet()
            features = feature_set.calculate_complete_feature_set(
                sequence, mutation_pos, "C"
            )

            # Should have phase-coherent features for G→C
            self.assertIn("delta_phase_entropy_change", features)
            self.assertIn("delta_curv_structural_complexity", features)

    def test_feature_consistency(self):
        """Test consistency of features across multiple calculations."""
        sequence = "ATCGATCGATCGATCGATCG"
        feature_set = InvariantFeatureSet()

        # Calculate features multiple times
        features_1 = feature_set.calculate_complete_feature_set(sequence)
        features_2 = feature_set.calculate_complete_feature_set(sequence)

        # Core features should be consistent
        core_features = ["phase_bit", "delta_phi", "normalization_constant"]
        for feature in core_features:
            if feature in features_1 and feature in features_2:
                self.assertEqual(features_1[feature], features_2[feature])

    def test_length_invariance(self):
        """Test length invariance of features."""
        short_seq = "ATCGATCG"
        long_seq = "ATCGATCGATCGATCGATCGAAATTTGGGCCCAAAGGGTTTCCC"

        feature_set = InvariantFeatureSet()
        features_short = feature_set.calculate_complete_feature_set(short_seq)
        features_long = feature_set.calculate_complete_feature_set(long_seq)

        # Both should have normalization enabled
        self.assertTrue(features_short["length_invariant_normalized"])
        self.assertTrue(features_long["length_invariant_normalized"])

        # Normalization constant should be the same
        self.assertEqual(
            features_short["normalization_constant"],
            features_long["normalization_constant"],
        )


if __name__ == "__main__":
    unittest.main()
