"""
Test suite for Phase-Weighted CRISPR Scorecard

This module tests the phase-weighted spectral analysis implementation
according to Z Framework principles and scientific gates.
"""

import unittest
import numpy as np
import sys
import os

# Add applications directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "applications"))

from phase_weighted_scorecard import (
    validate_dna_sequence,
    encode_complex,
    theta_prime,
    apply_phase_weighting,
    compute_spectrum,
    compute_spectral_entropy,
    find_dominant_frequency,
    count_sidelobes,
    compute_sequence_diversity,
    kappa_curvature,
    sigmoid_aggregator,
    PhaseWeightedScorecard,
    score_guide_batch,
    PHI,
    K_STAR,
    E_SQUARED,
)


class TestSequenceValidation(unittest.TestCase):
    """Test DNA/RNA sequence validation."""
    
    def test_valid_dna_sequence(self):
        """Test valid DNA sequences pass validation."""
        valid_seqs = [
            "ATCG",
            "ATCGATCGATCG",
            "AAATTTCCCGGG",
            "NNNAAATTT",  # N is allowed
        ]
        for seq in valid_seqs:
            try:
                validate_dna_sequence(seq, allow_rna=False)
            except ValueError:
                self.fail(f"Valid DNA sequence rejected: {seq}")
    
    def test_valid_rna_sequence(self):
        """Test valid RNA sequences pass validation."""
        valid_seqs = [
            "AUCG",
            "AUCGAUCGAUCG",
            "AAAUUUCCCGGG",
            "NNNAAAUUU",
        ]
        for seq in valid_seqs:
            try:
                validate_dna_sequence(seq, allow_rna=True)
            except ValueError:
                self.fail(f"Valid RNA sequence rejected: {seq}")
    
    def test_invalid_dna_with_u(self):
        """Test DNA sequences with U are rejected."""
        with self.assertRaises(ValueError) as cm:
            validate_dna_sequence("AUCG", allow_rna=False)
        self.assertIn("'U'", str(cm.exception))
    
    def test_invalid_rna_with_t(self):
        """Test RNA sequences with T are rejected."""
        with self.assertRaises(ValueError) as cm:
            validate_dna_sequence("ATCG", allow_rna=True)
        self.assertIn("'T'", str(cm.exception))
    
    def test_invalid_characters(self):
        """Test sequences with invalid characters are rejected."""
        invalid_seqs = ["ATCGX", "AT-CG", "ATCG!", "ATCGR"]
        for seq in invalid_seqs:
            with self.assertRaises(ValueError):
                validate_dna_sequence(seq, allow_rna=False)
    
    def test_iupac_codes_rejected(self):
        """Test IUPAC ambiguity codes (beyond N) are rejected."""
        iupac_codes = ["ATCGR", "ATCGY", "ATCGW", "ATCGS"]
        for seq in iupac_codes:
            with self.assertRaises(ValueError):
                validate_dna_sequence(seq, allow_rna=False)


class TestComplexEncoding(unittest.TestCase):
    """Test complex encoding of DNA/RNA sequences."""
    
    def test_dna_encoding(self):
        """Test DNA complex encoding."""
        seq = "ATCG"
        encoded = encode_complex(seq, is_rna=False)
        
        expected = np.array([
            1.0 + 0.0j,   # A
            -1.0 + 0.0j,  # T
            0.0 + 1.0j,   # C
            0.0 - 1.0j,   # G
        ], dtype=np.complex128)
        
        np.testing.assert_array_almost_equal(encoded, expected)
    
    def test_rna_encoding(self):
        """Test RNA complex encoding."""
        seq = "AUCG"
        encoded = encode_complex(seq, is_rna=True)
        
        expected = np.array([
            1.0 + 0.0j,   # A
            -1.0 + 0.0j,  # U
            0.0 + 1.0j,   # C
            0.0 - 1.0j,   # G
        ], dtype=np.complex128)
        
        np.testing.assert_array_almost_equal(encoded, expected)
    
    def test_ambiguous_base_n(self):
        """Test N encodes to zero."""
        seq = "ANTN"
        encoded = encode_complex(seq, is_rna=False)
        
        # N should be 0+0j
        self.assertEqual(encoded[1], 0.0 + 0.0j)
        self.assertEqual(encoded[3], 0.0 + 0.0j)
    
    def test_case_insensitive(self):
        """Test encoding is case-insensitive."""
        seq_upper = "ATCG"
        seq_lower = "atcg"
        
        encoded_upper = encode_complex(seq_upper, is_rna=False)
        encoded_lower = encode_complex(seq_lower, is_rna=False)
        
        np.testing.assert_array_almost_equal(encoded_upper, encoded_lower)


class TestThetaPrime(unittest.TestCase):
    """Test geometric resolution function θ′(n,k)."""
    
    def test_theta_prime_scalar(self):
        """Test θ′ with scalar input."""
        n = 10
        k = 0.3
        result = theta_prime(n, k)
        
        # Should return float
        self.assertIsInstance(result, (float, np.floating))
        
        # Should be positive
        self.assertGreater(result, 0)
        
        # Should be close to phi for typical values
        self.assertLess(result, 2 * PHI)
    
    def test_theta_prime_array(self):
        """Test θ′ with array input."""
        n_array = np.array([0, 1, 2, 3, 4, 5])
        k = K_STAR
        result = theta_prime(n_array, k)
        
        # Should return array
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(len(result), len(n_array))
        
        # All values should be non-negative (0 at n=0 is valid)
        self.assertTrue(np.all(result >= 0))
    
    def test_theta_prime_k_parameter(self):
        """Test θ′ with different k values."""
        n = 20
        k_values = [0.1, 0.3, 0.5, 0.7]
        results = [theta_prime(n, k) for k in k_values]
        
        # Results should be different for different k
        self.assertEqual(len(set(results)), len(k_values))
    
    def test_theta_prime_periodicity(self):
        """Test θ′ has quasi-periodic structure with phi."""
        # Values should show modular pattern with phi
        n_values = np.arange(0, 100)
        results = theta_prime(n_values, K_STAR)
        
        # Should have reasonable range
        self.assertLess(np.max(results), 3 * PHI)
        self.assertGreaterEqual(np.min(results), 0)  # Allow 0 at n=0


class TestPhaseWeighting(unittest.TestCase):
    """Test phase-weighted transform."""
    
    def test_phase_weighting_length(self):
        """Test phase weighting preserves array length."""
        seq = "ATCGATCGATCGATCG"
        encoded = encode_complex(seq)
        weighted = apply_phase_weighting(encoded, k=K_STAR)
        
        self.assertEqual(len(weighted), len(encoded))
    
    def test_phase_weighting_complex(self):
        """Test phase weighting produces complex output."""
        seq = "ATCG"
        encoded = encode_complex(seq)
        weighted = apply_phase_weighting(encoded, k=K_STAR)
        
        self.assertEqual(weighted.dtype, np.complex128)
    
    def test_phase_weighting_modifies_values(self):
        """Test phase weighting modifies encoded values."""
        seq = "ATCGATCGATCGATCG"
        encoded = encode_complex(seq)
        weighted = apply_phase_weighting(encoded, k=K_STAR)
        
        # Should modify values (not identical)
        self.assertFalse(np.allclose(encoded, weighted))


class TestSpectralFeatures(unittest.TestCase):
    """Test spectral feature extraction."""
    
    def test_compute_spectrum(self):
        """Test FFT spectrum computation."""
        seq = "ATCGATCGATCGATCG"
        encoded = encode_complex(seq)
        weighted = apply_phase_weighting(encoded, k=K_STAR)
        spectrum = compute_spectrum(weighted)
        
        # Should have same length
        self.assertEqual(len(spectrum), len(encoded))
        
        # Should be real and non-negative
        self.assertTrue(np.all(spectrum >= 0))
    
    def test_spectral_entropy(self):
        """Test spectral entropy calculation."""
        # Create simple spectrum
        spectrum = np.array([1.0, 2.0, 3.0, 4.0])
        entropy = compute_spectral_entropy(spectrum)
        
        # Should be positive
        self.assertGreater(entropy, 0)
        
        # Should be finite
        self.assertTrue(np.isfinite(entropy))
    
    def test_entropy_uniform_distribution(self):
        """Test entropy for uniform distribution."""
        # Uniform distribution should have maximum entropy
        spectrum_uniform = np.ones(10)
        entropy_uniform = compute_spectral_entropy(spectrum_uniform)
        
        # Compare with peaked distribution
        spectrum_peaked = np.array([10.0] + [1.0] * 9)
        entropy_peaked = compute_spectral_entropy(spectrum_peaked)
        
        # Uniform should have higher entropy
        self.assertGreater(entropy_uniform, entropy_peaked)
    
    def test_dominant_frequency(self):
        """Test dominant frequency detection."""
        # Create spectrum with clear peak at index 3
        spectrum = np.array([1.0, 2.0, 3.0, 10.0, 2.0, 1.0])
        freq_idx, freq_mag = find_dominant_frequency(spectrum)
        
        # Should find peak at index 3
        self.assertEqual(freq_idx, 3)
        self.assertAlmostEqual(freq_mag, 10.0)
    
    def test_sidelobe_count(self):
        """Test sidelobe counting."""
        # Create spectrum with multiple peaks
        spectrum = np.array([1.0, 5.0, 1.0, 8.0, 1.0, 4.0, 1.0])
        count = count_sidelobes(spectrum, threshold_ratio=0.3)
        
        # Should detect peaks above threshold
        self.assertGreater(count, 0)
        self.assertIsInstance(count, (int, np.integer))


class TestSequenceDiversity(unittest.TestCase):
    """Test sequence diversity calculation."""
    
    def test_diversity_balanced(self):
        """Test diversity for balanced sequence."""
        # Equal distribution of all bases
        seq = "AAATTTCCCGGG"
        diversity = compute_sequence_diversity(seq)
        
        # Should be close to 1 (maximum diversity)
        self.assertGreater(diversity, 0.9)
    
    def test_diversity_homogeneous(self):
        """Test diversity for homogeneous sequence."""
        # All same base
        seq = "AAAAAAAAAA"
        diversity = compute_sequence_diversity(seq)
        
        # Should be 0 (no diversity)
        self.assertAlmostEqual(diversity, 0.0)
    
    def test_diversity_range(self):
        """Test diversity is in valid range."""
        seqs = [
            "ATCGATCG",
            "AAAATTTT",
            "ATATATATATAT",
            "GGGGCCCC",
        ]
        for seq in seqs:
            diversity = compute_sequence_diversity(seq)
            self.assertGreaterEqual(diversity, 0.0)
            self.assertLessEqual(diversity, 1.0)


class TestKappaCurvature(unittest.TestCase):
    """Test curvature weight calculation."""
    
    def test_kappa_short_sequence(self):
        """Test κ(n) for short sequences defaults to 1."""
        for n in range(1, 10):
            kappa = kappa_curvature(n, d_n=0.5)
            self.assertEqual(kappa, 1.0)
    
    def test_kappa_long_sequence(self):
        """Test κ(n) for long sequences."""
        n = 20
        diversity = 0.8
        kappa = kappa_curvature(n, d_n=diversity)
        
        # Should not be 1.0 for n >= 10
        self.assertNotEqual(kappa, 1.0)
        
        # Should be positive
        self.assertGreater(kappa, 0)
    
    def test_kappa_increases_with_length(self):
        """Test κ(n) increases with sequence length."""
        diversity = 0.5
        kappa_10 = kappa_curvature(10, d_n=diversity)
        kappa_20 = kappa_curvature(20, d_n=diversity)
        kappa_50 = kappa_curvature(50, d_n=diversity)
        
        # Should increase with length
        self.assertLess(kappa_10, kappa_20)
        self.assertLess(kappa_20, kappa_50)
    
    def test_kappa_increases_with_diversity(self):
        """Test κ(n) increases with diversity."""
        n = 20
        kappa_low = kappa_curvature(n, d_n=0.2)
        kappa_high = kappa_curvature(n, d_n=0.9)
        
        # Should increase with diversity
        self.assertLess(kappa_low, kappa_high)


class TestSigmoidAggregator(unittest.TestCase):
    """Test sigmoid aggregator function."""
    
    def test_sigmoid_range(self):
        """Test sigmoid output is in (0, 1)."""
        values = [-10, -1, 0, 1, 10]
        for x in values:
            result = sigmoid_aggregator(x, kappa=1.0)
            self.assertGreater(result, 0)
            self.assertLess(result, 1)
    
    def test_sigmoid_zero(self):
        """Test sigmoid at zero is 0.5."""
        result = sigmoid_aggregator(0, kappa=1.0)
        self.assertAlmostEqual(result, 0.5)
    
    def test_sigmoid_monotonic(self):
        """Test sigmoid is monotonically increasing."""
        x_values = np.linspace(-5, 5, 100)
        results = [sigmoid_aggregator(x, kappa=1.0) for x in x_values]
        
        # Should be monotonically increasing
        for i in range(len(results) - 1):
            self.assertLessEqual(results[i], results[i + 1])


class TestPhaseWeightedScorecard(unittest.TestCase):
    """Test PhaseWeightedScorecard class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.scorecard = PhaseWeightedScorecard()
        self.test_guide = "GCTGCGGAGACCTGGAGAGA"
        self.test_target = "GCTGCGGAGACCTGGAGAGA"
        self.test_mutant = "GCTGCGGAGTCCTGGAGAGA"  # Single mutation at position 10
    
    def test_initialization(self):
        """Test scorecard initialization."""
        # Default initialization
        sc = PhaseWeightedScorecard()
        self.assertEqual(sc.k, K_STAR)
        self.assertFalse(sc.is_rna)
        
        # Custom initialization
        sc_custom = PhaseWeightedScorecard(k=0.5, is_rna=True)
        self.assertEqual(sc_custom.k, 0.5)
        self.assertTrue(sc_custom.is_rna)
    
    def test_compute_spectral_features(self):
        """Test spectral feature computation."""
        features = self.scorecard.compute_spectral_features(self.test_guide)
        
        # Should contain expected keys
        self.assertIn("entropy", features)
        self.assertIn("dominant_freq_idx", features)
        self.assertIn("dominant_freq_mag", features)
        self.assertIn("sidelobe_count", features)
        self.assertIn("diversity", features)
        
        # Values should be valid
        self.assertGreater(features["entropy"], 0)
        self.assertGreaterEqual(features["dominant_freq_idx"], 0)
        self.assertGreater(features["dominant_freq_mag"], 0)
        self.assertGreaterEqual(features["sidelobe_count"], 0)
        self.assertGreaterEqual(features["diversity"], 0)
        self.assertLessEqual(features["diversity"], 1)
    
    def test_compute_spectral_features_short_sequence(self):
        """Test error handling for short sequences."""
        short_seq = "ATCG"
        with self.assertRaises(ValueError) as cm:
            self.scorecard.compute_spectral_features(short_seq)
        self.assertIn("too short", str(cm.exception).lower())
    
    def test_compute_disruption_features(self):
        """Test disruption feature computation."""
        disruptions = self.scorecard.compute_disruption_features(
            self.test_target, self.test_mutant
        )
        
        # Should contain expected keys
        self.assertIn("delta_entropy", disruptions)
        self.assertIn("delta_freq", disruptions)
        self.assertIn("delta_sidelobes", disruptions)
        self.assertIn("ref_features", disruptions)
        self.assertIn("mut_features", disruptions)
        
        # Deltas should be non-negative (absolute values)
        self.assertGreaterEqual(disruptions["delta_freq"], 0)
        self.assertGreaterEqual(disruptions["delta_sidelobes"], 0)
    
    def test_compute_disruption_features_length_mismatch(self):
        """Test error handling for length mismatch."""
        with self.assertRaises(ValueError) as cm:
            self.scorecard.compute_disruption_features(
                "ATCGATCG", "ATCG"
            )
        self.assertIn("same length", str(cm.exception).lower())
    
    def test_compute_z_score(self):
        """Test Z-score computation."""
        result = self.scorecard.compute_z_score(
            self.test_target, self.test_mutant
        )
        
        # Should contain expected keys
        self.assertIn("z_score", result)
        self.assertIn("delta_spectral", result)
        self.assertIn("kappa", result)
        self.assertIn("disruptions", result)
        
        # Z-score should be in (0, 1)
        self.assertGreater(result["z_score"], 0)
        self.assertLess(result["z_score"], 1)
        
        # Kappa should be positive
        self.assertGreater(result["kappa"], 0)
    
    def test_score_guide_without_target(self):
        """Test guide scoring without target."""
        result = self.scorecard.score_guide(self.test_guide)
        
        # Should contain features
        self.assertIn("guide_features", result)
        
        # Should not have z_score
        self.assertIsNone(result["z_score"])
        self.assertIsNone(result["delta_spectral"])
    
    def test_score_guide_with_target(self):
        """Test guide scoring with target."""
        result = self.scorecard.score_guide(
            self.test_mutant, target_seq=self.test_target
        )
        
        # Should contain z_score
        self.assertIn("z_score", result)
        self.assertIsNotNone(result["z_score"])
        
        # Should contain features
        self.assertIn("guide_features", result)
        self.assertIn("target_features", result)


class TestBatchScoring(unittest.TestCase):
    """Test batch scoring functionality."""
    
    def test_score_guide_batch_without_targets(self):
        """Test batch scoring without targets."""
        guides = [
            "GCTGCGGAGACCTGGAGAGA",
            "ATCGATCGATCGATCGATCG",
            "AAAATTTTCCCCGGGGAAAA",
        ]
        results = score_guide_batch(guides)
        
        self.assertEqual(len(results), len(guides))
        
        for result in results:
            self.assertIn("guide_features", result)
    
    def test_score_guide_batch_with_targets(self):
        """Test batch scoring with targets."""
        guides = [
            "GCTGCGGAGACCTGGAGAGA",
            "ATCGATCGATCGATCGATCG",
        ]
        targets = [
            "GCTGCGGAGACCTGGAGAGA",
            "ATCGATCGATCGATCGATCG",
        ]
        results = score_guide_batch(guides, targets)
        
        self.assertEqual(len(results), len(guides))
        
        for result in results:
            if "error" not in result:
                self.assertIn("z_score", result)
    
    def test_score_guide_batch_length_mismatch(self):
        """Test error handling for mismatched lengths."""
        guides = ["ATCGATCGATCGATCGATCG"]
        targets = ["ATCGATCGATCGATCGATCG", "AAAATTTTCCCCGGGGAAAA"]
        
        with self.assertRaises(ValueError):
            score_guide_batch(guides, targets)


class TestIntegration(unittest.TestCase):
    """Integration tests for end-to-end workflows."""
    
    def test_complete_workflow(self):
        """Test complete scoring workflow."""
        # Create scorecard
        scorecard = PhaseWeightedScorecard(k=K_STAR)
        
        # Define sequences
        ref_seq = "GCTGCGGAGACCTGGAGAGAAAGC"
        mut_seq = "GCTGCGGAGTCCTGGAGAGAAAGC"  # A→T at position 10
        
        # Compute Z-score
        result = scorecard.compute_z_score(ref_seq, mut_seq)
        
        # Verify result structure
        self.assertIn("z_score", result)
        self.assertIn("delta_spectral", result)
        self.assertIn("kappa", result)
        
        # Verify values are reasonable
        self.assertGreater(result["z_score"], 0)
        self.assertLess(result["z_score"], 1)
        self.assertIsInstance(result["delta_spectral"], (float, np.floating))
        self.assertIsInstance(result["kappa"], (float, np.floating))
    
    def test_rna_workflow(self):
        """Test workflow with RNA sequences."""
        scorecard = PhaseWeightedScorecard(k=K_STAR, is_rna=True)
        
        # RNA sequences
        ref_rna = "GCUGCGGAGACCUGGAGAGAAAGC"
        mut_rna = "GCUGCGGAGUCCUGGAGAGAAAGC"  # A→U
        
        # Compute Z-score
        result = scorecard.compute_z_score(ref_rna, mut_rna)
        
        # Should complete successfully
        self.assertIn("z_score", result)
        self.assertGreater(result["z_score"], 0)


class TestConstants(unittest.TestCase):
    """Test mathematical constants."""
    
    def test_phi_value(self):
        """Test golden ratio value."""
        # φ = (1 + √5) / 2 ≈ 1.618
        expected_phi = (1 + np.sqrt(5)) / 2
        self.assertAlmostEqual(PHI, expected_phi, places=10)
    
    def test_e_squared_value(self):
        """Test e² value."""
        expected_e2 = np.e ** 2
        self.assertAlmostEqual(E_SQUARED, expected_e2, places=10)
    
    def test_k_star_value(self):
        """Test K_STAR is 0.3."""
        self.assertEqual(K_STAR, 0.3)


if __name__ == "__main__":
    unittest.main()
