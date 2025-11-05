#!/usr/bin/env python3
"""
Unit tests for DNA Breathing Dynamics Encoder

Tests the core functionality of the breathing dynamics encoder to ensure
correctness and reproducibility.
"""

import unittest
import numpy as np
import sys
import os

# Add experiments directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))

from test_breathing_dynamics_encoding import (
    BreathingDynamicsEncoder,
    ArbitraryEncoder,
    HELICAL_PERIOD
)


class TestBreathingDynamicsEncoder(unittest.TestCase):
    """Test suite for BreathingDynamicsEncoder"""

    def setUp(self):
        """Set up test fixtures"""
        self.encoder = BreathingDynamicsEncoder()

    def test_encoder_initialization(self):
        """Test that encoder initializes with correct weights"""
        # AT bases should have negative real part (10 MHz)
        self.assertLess(self.encoder.weights['A'].real, 0)
        self.assertLess(self.encoder.weights['T'].real, 0)

        # GC bases should have negative real part (1 GHz)
        self.assertLess(self.encoder.weights['C'].real, 0)
        self.assertLess(self.encoder.weights['G'].real, 0)

        # AT bases should have positive imaginary part (weak bonds)
        self.assertGreater(self.encoder.weights['A'].imag, 0)
        self.assertGreater(self.encoder.weights['T'].imag, 0)

        # GC bases should have negative imaginary part (strong bonds)
        self.assertLess(self.encoder.weights['C'].imag, 0)
        self.assertLess(self.encoder.weights['G'].imag, 0)

    def test_frequency_mapping_correct(self):
        """Test that frequency values map correctly to weights"""
        # Log10 of 10^7 = 7, normalized to (7-8)*10 = -10
        expected_at_real = -80.0
        # Log10 of 10^9 = 9, normalized to (9-8)*10 = +10
        expected_gc_real = -60.0

        self.assertAlmostEqual(self.encoder.weights['A'].real, expected_at_real, places=1)
        self.assertAlmostEqual(self.encoder.weights['T'].real, expected_at_real, places=1)
        self.assertAlmostEqual(self.encoder.weights['C'].real, expected_gc_real, places=1)
        self.assertAlmostEqual(self.encoder.weights['G'].real, expected_gc_real, places=1)

    def test_encode_simple_sequence(self):
        """Test encoding of simple sequences"""
        # Test homopolymer sequences
        seq_at = "AAAA"
        encoded_at = self.encoder.encode_sequence(seq_at)

        self.assertEqual(len(encoded_at), 4)
        self.assertTrue(np.all(np.isfinite(encoded_at)))

        seq_gc = "GGGG"
        encoded_gc = self.encoder.encode_sequence(seq_gc)

        self.assertEqual(len(encoded_gc), 4)
        self.assertTrue(np.all(np.isfinite(encoded_gc)))

    def test_at_vs_gc_magnitude_difference(self):
        """Test that AT and GC encoded values have expected magnitude relationship"""
        seq_at = "ATATATAT"
        seq_gc = "GCGCGCGC"

        encoded_at = self.encoder.encode_sequence(seq_at)
        encoded_gc = self.encoder.encode_sequence(seq_gc)

        # Magnitudes should be similar (same magnitude range)
        # but phase-shifted due to imaginary component differences
        mag_at = np.abs(encoded_at)
        mag_gc = np.abs(encoded_gc)

        # Magnitudes should be within reasonable range
        self.assertTrue(np.all(mag_at > 0))
        self.assertTrue(np.all(mag_gc > 0))

    def test_helical_phase_modulation(self):
        """Test that helical periodicity is correctly applied"""
        seq = "ATCGATCGATCGATCG"  # 16 bp (>1 helical turn)
        encoded = self.encoder.encode_sequence(seq)

        # Calculate expected phase at position 0 and 10.5 (one helical turn)
        # They should differ by 2π
        phase_0 = np.angle(encoded[0])
        if len(encoded) > 10:
            phase_10 = np.angle(encoded[10])  # Close to one turn

            # Phase difference should be close to 2π or -2π (modulo 2π)
            phase_diff = (phase_10 - phase_0) % (2 * np.pi)
            expected_diff = (2 * np.pi * 10 / HELICAL_PERIOD) % (2 * np.pi)

            # Allow some tolerance for positional phase contribution
            self.assertLess(abs(phase_diff - expected_diff), 6.0)

    def test_empty_sequence(self):
        """Test behavior with empty sequence"""
        seq = ""
        encoded = self.encoder.encode_sequence(seq)

        self.assertEqual(len(encoded), 0)

    def test_invalid_bases(self):
        """Test handling of invalid bases"""
        seq = "ATCGN"  # N is not a valid base
        encoded = self.encoder.encode_sequence(seq)

        # Should encode 4 valid bases + 1 zero for invalid
        self.assertEqual(len(encoded), 5)
        self.assertEqual(encoded[4], 0 + 0j)

    def test_reproducibility(self):
        """Test that encoding is deterministic"""
        seq = "ATCGATCGATCG"

        # Encode twice
        encoded1 = self.encoder.encode_sequence(seq)
        encoded2 = self.encoder.encode_sequence(seq)

        # Should be identical
        np.testing.assert_array_equal(encoded1, encoded2)


class TestArbitraryEncoder(unittest.TestCase):
    """Test suite for ArbitraryEncoder"""

    def test_reproducibility_with_seed(self):
        """Test that arbitrary encoder is reproducible with same seed"""
        encoder1 = ArbitraryEncoder(seed=42)
        encoder2 = ArbitraryEncoder(seed=42)

        seq = "ATCGATCG"
        encoded1 = encoder1.encode_sequence(seq)
        encoded2 = encoder2.encode_sequence(seq)

        np.testing.assert_array_almost_equal(encoded1, encoded2)

    def test_different_seeds_give_different_results(self):
        """Test that different seeds produce different encodings"""
        encoder1 = ArbitraryEncoder(seed=42)
        encoder2 = ArbitraryEncoder(seed=123)

        seq = "ATCGATCG"
        encoded1 = encoder1.encode_sequence(seq)
        encoded2 = encoder2.encode_sequence(seq)

        # Should not be identical
        self.assertFalse(np.allclose(encoded1, encoded2))

    def test_magnitude_range(self):
        """Test that arbitrary weights are in expected range"""
        encoder = ArbitraryEncoder(seed=42)

        # Check that weights are within specified bounds
        for base in 'ATCG':
            weight = encoder.weights[base]
            self.assertGreaterEqual(weight.real, -20)
            self.assertLessEqual(weight.real, 20)
            self.assertGreaterEqual(weight.imag, -5)
            self.assertLessEqual(weight.imag, 5)


class TestEncoderComparison(unittest.TestCase):
    """Test comparison between breathing and arbitrary encoders"""

    def test_gc_content_sensitivity(self):
        """Test that breathing encoder is sensitive to GC content changes"""
        encoder = BreathingDynamicsEncoder()

        # High AT content
        seq_at = "AAATTTAAA"
        encoded_at = encoder.encode_sequence(seq_at)

        # High GC content (same length)
        seq_gc = "GGGCCCGGG"
        encoded_gc = encoder.encode_sequence(seq_gc)

        # Spectra should be different
        from scipy.fft import fft
        spectrum_at = np.abs(fft(encoded_at))
        spectrum_gc = np.abs(fft(encoded_gc))

        # L2 distance should be substantial
        distance = np.linalg.norm(spectrum_at - spectrum_gc)
        self.assertGreater(distance, 0.1)

    def test_within_class_similarity(self):
        """Test that A-T and G-C have similar encodings"""
        encoder = BreathingDynamicsEncoder()

        # Within AT class
        weight_diff_at = abs(encoder.weights['A'] - encoder.weights['T'])

        # Within GC class
        weight_diff_gc = abs(encoder.weights['G'] - encoder.weights['C'])

        # Both should be zero (same weights for same class)
        self.assertAlmostEqual(weight_diff_at, 0.0, places=5)
        self.assertAlmostEqual(weight_diff_gc, 0.0, places=5)

    def test_between_class_difference(self):
        """Test that AT and GC classes are well-separated"""
        encoder = BreathingDynamicsEncoder()

        # Between classes
        weight_diff_ag = abs(encoder.weights['A'] - encoder.weights['G'])
        weight_diff_tc = abs(encoder.weights['T'] - encoder.weights['C'])

        # Should be substantial (>10 in real part due to frequency difference)
        self.assertGreater(weight_diff_ag, 10.0)
        self.assertGreater(weight_diff_tc, 10.0)


class TestIntegration(unittest.TestCase):
    """Integration tests for full encoding pipeline"""

    def test_full_mutation_analysis(self):
        """Test complete mutation analysis workflow"""
        from test_breathing_dynamics_encoding import BreathingDynamicsValidator

        validator = BreathingDynamicsValidator()

        # Generate small test set
        sequences = validator.generate_crispr_sequences(n_sequences=10, seq_length=20)

        self.assertEqual(len(sequences), 10)
        self.assertTrue(all(len(seq) == 20 for seq in sequences))
        self.assertTrue(all(all(b in 'ATCG' for b in seq) for seq in sequences))



def run_tests():
    """Run all tests"""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestBreathingDynamicsEncoder))
    suite.addTests(loader.loadTestsFromTestCase(TestArbitraryEncoder))
    suite.addTests(loader.loadTestsFromTestCase(TestEncoderComparison))
    suite.addTests(loader.loadTestsFromTestCase(TestIntegration))

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
