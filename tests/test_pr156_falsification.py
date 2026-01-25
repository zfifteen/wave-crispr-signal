"""
Unit tests for PR-156 Falsification Experiments

This module tests the falsification experiment implementations for
Hypothesis 1 and Hypothesis 2.

Tests validate:
- Spectral feature computation
- Z-score calculation
- Statistical test execution
- Result serialization
"""

import unittest
import sys
import os
import numpy as np
from pathlib import Path

# Add parent directories to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments", "PR-156-falsify-crispr-hypothesis"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "applications"))

from falsify_hypothesis1 import Hypothesis1Falsifier, generate_synthetic_crispr_data
from falsify_hypothesis2 import Hypothesis2Falsifier, generate_variable_length_sequences


class TestHypothesis1Falsifier(unittest.TestCase):
    """Test Hypothesis 1 falsification implementation."""
    
    def setUp(self):
        """Set up test fixture."""
        self.falsifier = Hypothesis1Falsifier(seed=42, k_parameter=0.3)
    
    def test_initialization(self):
        """Test falsifier initialization."""
        self.assertEqual(self.falsifier.seed, 42)
        self.assertEqual(self.falsifier.k, 0.3)
    
    def test_encode_and_phase(self):
        """Test DNA encoding and phase weighting."""
        sequence = "ATCG"
        result = self.falsifier.encode_and_phase(sequence)
        
        # Should be complex array
        self.assertEqual(result.dtype, np.complex128)
        self.assertEqual(len(result), 4)
    
    def test_spectral_features(self):
        """Test spectral feature computation."""
        wt_seq = "ATCGATCGATCGATCGATCG"
        mut_seq = "ATCGATCGAGCGATCGATCG"  # One mutation
        
        features = self.falsifier.compute_spectral_features(wt_seq, mut_seq)
        
        # Check expected keys
        expected_keys = [
            "delta_entropy", "delta_sidelobes", "freq_shift",
            "wt_entropy", "mut_entropy", "wt_sidelobes", "mut_sidelobes"
        ]
        for key in expected_keys:
            self.assertIn(key, features)
        
        # All values should be numeric
        for key, value in features.items():
            self.assertIsInstance(value, (int, float, np.number))
    
    def test_baseline_features(self):
        """Test baseline feature computation."""
        sequence = "ATCGATCGATCGATCGATCG"
        features = self.falsifier.compute_baseline_features(sequence)
        
        # Check expected keys
        expected_keys = ["gc_content", "length", "gc_first_5", "gc_last_5"]
        for key in expected_keys:
            self.assertIn(key, features)
        
        # GC content should be 0.5 for ATCG repeat
        self.assertAlmostEqual(features["gc_content"], 0.5, places=2)
        self.assertEqual(features["length"], 20)
    
    def test_synthetic_data_generation(self):
        """Test synthetic CRISPR data generation."""
        wt_seqs, mut_seqs, labels = generate_synthetic_crispr_data(
            n_samples=10, seed=42
        )
        
        self.assertEqual(len(wt_seqs), 10)
        self.assertEqual(len(mut_seqs), 10)
        self.assertEqual(len(labels), 10)
        
        # All sequences should be 20 nt
        for seq in wt_seqs + mut_seqs:
            self.assertEqual(len(seq), 20)
            # All bases should be valid
            self.assertTrue(all(b in 'ATCG' for b in seq))
        
        # Labels should be binary
        self.assertTrue(all(label in [0, 1] for label in labels))


class TestHypothesis2Falsifier(unittest.TestCase):
    """Test Hypothesis 2 falsification implementation."""
    
    def setUp(self):
        """Set up test fixture."""
        self.falsifier = Hypothesis2Falsifier(seed=42, k_parameter=0.3)
    
    def test_initialization(self):
        """Test falsifier initialization."""
        self.assertEqual(self.falsifier.seed, 42)
        self.assertEqual(self.falsifier.k, 0.3)
    
    def test_z_score_computation(self):
        """Test Z-score computation."""
        sequence = "ATCGATCGATCGATCGATCG"
        z_score = self.falsifier.compute_z_score(sequence)
        
        # Z-score should be in [0, 1] range (sigmoid output)
        self.assertGreaterEqual(z_score, 0.0)
        self.assertLessEqual(z_score, 1.0)
    
    def test_z_score_different_lengths(self):
        """Test Z-score computation for different lengths."""
        sequences = {
            20: "ATCGATCGATCGATCGATCG",
            30: "ATCGATCGATCGATCGATCGATCGATCGATCG",
            40: "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        }
        
        z_scores = {}
        for length, seq in sequences.items():
            z_scores[length] = self.falsifier.compute_z_score(seq)
        
        # All should be valid
        for z in z_scores.values():
            self.assertGreaterEqual(z, 0.0)
            self.assertLessEqual(z, 1.0)
    
    def test_variable_length_sequences(self):
        """Test variable-length sequence generation."""
        lengths = [20, 30, 40]
        sequences_by_length = generate_variable_length_sequences(
            lengths, n_per_length=5, seed=42
        )
        
        self.assertEqual(len(sequences_by_length), 3)
        
        for length, seqs in sequences_by_length.items():
            self.assertEqual(len(seqs), 5)
            for seq in seqs:
                self.assertEqual(len(seq), length)
                self.assertTrue(all(b in 'ATCG' for b in seq))
    
    def test_length_invariance_test(self):
        """Test length invariance testing."""
        lengths = [20, 30, 40]
        sequences_by_length = generate_variable_length_sequences(
            lengths, n_per_length=5, seed=42
        )
        
        results = self.falsifier.test_length_invariance(sequences_by_length)
        
        # Check expected keys
        expected_keys = [
            "f_statistic", "p_value", "invariant", 
            "overall_variance", "mean_z_by_length"
        ]
        for key in expected_keys:
            self.assertIn(key, results)
        
        # Invariance should be boolean
        self.assertIsInstance(results["invariant"], bool)
    
    def test_periodicity_test(self):
        """Test periodicity detection."""
        # Create simple sequence list
        sequences = ["ATCGATCGATCGATCGATCG"] * 10
        
        results = self.falsifier.test_periodicity(sequences)
        
        # Check expected keys
        expected_keys = [
            "max_autocorrelation", "periodicity_detected", 
            "autocorrelation_values", "z_scores"
        ]
        for key in expected_keys:
            self.assertIn(key, results)
        
        # Periodicity detected should be boolean
        self.assertIsInstance(results["periodicity_detected"], bool)


class TestIntegration(unittest.TestCase):
    """Integration tests for complete workflows."""
    
    def test_hypothesis1_smoke_run(self):
        """Test Hypothesis 1 smoke run (minimal parameters)."""
        falsifier = Hypothesis1Falsifier(seed=42, k_parameter=0.3)
        
        # Generate small dataset
        wt_seqs, mut_seqs, labels = generate_synthetic_crispr_data(
            n_samples=10, seed=42
        )
        
        # Extract features
        spectral_features = []
        baseline_features = []
        
        for wt_seq, mut_seq in zip(wt_seqs, mut_seqs):
            spec_feat = falsifier.compute_spectral_features(wt_seq, mut_seq)
            spectral_features.append([
                spec_feat["delta_entropy"],
                spec_feat["delta_sidelobes"],
                spec_feat["freq_shift"],
            ])
            
            base_feat = falsifier.compute_baseline_features(wt_seq)
            baseline_features.append([
                base_feat["gc_content"],
                base_feat["length"],
            ])
        
        X_spectral = np.array(spectral_features)
        X_baseline = np.array(baseline_features)
        
        # Should have correct shapes
        self.assertEqual(X_spectral.shape, (10, 3))
        self.assertEqual(X_baseline.shape, (10, 2))
    
    def test_hypothesis2_smoke_run(self):
        """Test Hypothesis 2 smoke run (minimal parameters)."""
        falsifier = Hypothesis2Falsifier(seed=42, k_parameter=0.3)
        
        # Generate sequences
        lengths = [20, 30, 40]
        sequences_by_length = generate_variable_length_sequences(
            lengths, n_per_length=5, seed=42
        )
        
        # Test invariance
        results = falsifier.test_length_invariance(sequences_by_length)
        
        # Should have results
        self.assertIsNotNone(results["f_statistic"])
        self.assertIsNotNone(results["p_value"])


if __name__ == "__main__":
    unittest.main()
