#!/usr/bin/env python3
"""
Test Suite for Genomic Disruption Analyzer API

This module tests the core DisruptionAnalyzer class and GenomicDisruptionAPI wrapper,
ensuring scientific gates are enforced and results are reproducible.

Scientific Gates Tested:
- Human DNA/RNA only (A/C/G/T or A/C/G/U, N allowed)
- Fail-fast validation with clear errors
- Reproducibility with seed control
- Bootstrap CI validity
- Z-invariant normalization

Test Coverage:
- Single guide scoring
- Batch processing
- Guide design
- Off-target analysis
- Edge cases (short sequences, high GC, low GC)
- Invalid sequences (reject U in DNA, T in RNA, other bases)
"""

import sys
import os
import unittest
from typing import List, Dict

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from applications.genomic_disruption_api import (
    DisruptionAnalyzer,
    GenomicDisruptionAPI,
    K_STAR,
)


class TestDisruptionAnalyzer(unittest.TestCase):
    """Test DisruptionAnalyzer core functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = DisruptionAnalyzer(k=K_STAR, seed=42)
    
    def test_score_single_guide_dna(self):
        """Test scoring a single DNA guide."""
        guide = "GCTGCGGAGACCTGGAGAGA"
        result = self.analyzer.score_guide(guide)
        
        self.assertIn('disruption_score', result)
        self.assertIn('spectral_features', result)
        self.assertIn('metadata', result)
        
        # Check score is in valid range
        self.assertGreaterEqual(result['disruption_score'], 0.0)
        self.assertLessEqual(result['disruption_score'], 1.0)
        
        # Check metadata
        self.assertEqual(result['metadata']['guide_length'], 20)
        self.assertEqual(result['metadata']['guide_type'], 'DNA')
    
    def test_score_single_guide_rna(self):
        """Test scoring a single RNA guide (with U)."""
        guide = "GCUGCGGAGACCUGGAGAGA"  # RNA with U
        result = self.analyzer.score_guide(guide)
        
        self.assertIn('disruption_score', result)
        self.assertEqual(result['metadata']['guide_type'], 'RNA')
    
    def test_score_with_target(self):
        """Test scoring guide with target sequence."""
        guide = "GCTGCGGAGACCTGGAGAGA"
        target = "GCTGCGGAGACCTGGAGAGA"  # Same as guide
        
        result = self.analyzer.score_guide(guide, target)
        
        self.assertIn('disruption_score', result)
        self.assertIn('on_target_metrics', result)
        self.assertIn('delta_entropy', result['on_target_metrics'])
    
    def test_score_with_bootstrap_ci(self):
        """Test bootstrap confidence interval computation."""
        guide = "GCTGCGGAGACCTGGAGAGA"
        result = self.analyzer.score_guide(guide, compute_ci=True, n_bootstrap=100)
        
        self.assertIn('confidence_interval', result)
        ci = result['confidence_interval']
        
        self.assertIn('median', ci)
        self.assertIn('lower_95', ci)
        self.assertIn('upper_95', ci)
        
        # CI bounds should be ordered
        self.assertLessEqual(ci['lower_95'], ci['median'])
        self.assertLessEqual(ci['median'], ci['upper_95'])
    
    def test_batch_score(self):
        """Test batch scoring of multiple guides."""
        guides = [
            "GCTGCGGAGACCTGGAGAGA",
            "ATCGATCGATCGATCGATCG",
            "GGGGGGGGGGGGGGGGGGGG",
        ]
        
        results = self.analyzer.batch_score(guides)
        
        self.assertEqual(len(results), 3)
        
        for i, result in enumerate(results):
            self.assertEqual(result['index'], i)
            if 'error' not in result:
                self.assertIn('disruption_score', result)
    
    def test_design_guides(self):
        """Test guide design from target sequence."""
        target = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGGATCGATC"
        results = self.analyzer.design_guides(target, n_guides=3)
        
        self.assertGreater(len(results), 0)
        self.assertLessEqual(len(results), 3)
        
        for result in results:
            self.assertIn('guide', result)
            self.assertIn('position', result)
            self.assertIn('pam_site', result)
            self.assertIn('disruption_score', result)
            
            # Guide should be 20 bp
            self.assertEqual(len(result['guide']), 20)
    
    def test_offtarget_analysis(self):
        """Test off-target analysis."""
        guide = "GCTGCGGAGACCTGGAGAGA"
        offtargets = [
            "GCTGCGGAGACCTGGAGAGA",  # Perfect match
            "ACTGCGGAGACCTGGAGAGA",  # 1 mismatch
            "ACTGCGGAGACCTGGAGACA",  # 2 mismatches
        ]
        
        results = self.analyzer.analyze_offtarget(guide, offtargets, mismatch_threshold=2)
        
        self.assertGreater(len(results), 0)
        
        for result in results:
            self.assertIn('sequence', result)
            self.assertIn('mismatches', result)
            self.assertIn('spectral_distance', result)
            self.assertLessEqual(result['mismatches'], 2)
    
    def test_invalid_dna_sequence(self):
        """Test rejection of invalid DNA sequences."""
        # Mixed U and T should fail (ambiguous DNA/RNA)
        with self.assertRaises(ValueError) as cm:
            self.analyzer.score_guide("GCUGCGTACGTACGTA")  # Mixed U (RNA) and T (DNA)
        # Should detect U first and then complain about T in RNA context
        self.assertIn("contains", str(cm.exception).lower())
        
        # Invalid bases should fail
        with self.assertRaises(ValueError):
            self.analyzer.score_guide("GCTXCGGAGACCTGGAGAGA")
    
    def test_high_gc_guide(self):
        """Test scoring high GC content guide."""
        guide = "GCGCGCGCGCGCGCGCGCGC"  # 100% GC
        result = self.analyzer.score_guide(guide)
        
        self.assertIn('disruption_score', result)
        # High GC should still produce valid score
        self.assertGreaterEqual(result['disruption_score'], 0.0)
    
    def test_low_gc_guide(self):
        """Test scoring low GC content guide."""
        guide = "ATATATATATATATATATAT"  # 0% GC
        result = self.analyzer.score_guide(guide)
        
        self.assertIn('disruption_score', result)
        # Low GC should still produce valid score
        self.assertGreaterEqual(result['disruption_score'], 0.0)
    
    def test_reproducibility_with_seed(self):
        """Test that results are reproducible with same seed."""
        guide = "GCTGCGGAGACCTGGAGAGA"
        
        analyzer1 = DisruptionAnalyzer(seed=42)
        result1 = analyzer1.score_guide(guide, compute_ci=True, n_bootstrap=100)
        
        analyzer2 = DisruptionAnalyzer(seed=42)
        result2 = analyzer2.score_guide(guide, compute_ci=True, n_bootstrap=100)
        
        # Scores should be identical with same seed
        self.assertAlmostEqual(
            result1['disruption_score'],
            result2['disruption_score'],
            places=10
        )


class TestGenomicDisruptionAPI(unittest.TestCase):
    """Test GenomicDisruptionAPI wrapper functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.api = GenomicDisruptionAPI(k=K_STAR, seed=42)
    
    def test_handle_score(self):
        """Test API score handler."""
        request = {
            'guide': 'GCTGCGGAGACCTGGAGAGA',
        }
        
        response = self.api.handle_score(request)
        
        self.assertTrue(response['success'])
        self.assertIn('data', response)
        self.assertIn('disruption_score', response['data'])
    
    def test_handle_score_with_target(self):
        """Test API score handler with target."""
        request = {
            'guide': 'GCTGCGGAGACCTGGAGAGA',
            'target': 'GCTGCGGAGACCTGGAGAGA',
        }
        
        response = self.api.handle_score(request)
        
        self.assertTrue(response['success'])
        self.assertIn('on_target_metrics', response['data'])
    
    def test_handle_batch(self):
        """Test API batch handler."""
        request = {
            'guides': [
                'GCTGCGGAGACCTGGAGAGA',
                'ATCGATCGATCGATCGATCG',
            ]
        }
        
        response = self.api.handle_batch(request)
        
        self.assertTrue(response['success'])
        self.assertEqual(response['count'], 2)
        self.assertIn('data', response)
    
    def test_handle_design(self):
        """Test API design handler."""
        request = {
            'target': 'ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGGATCGATC',
            'n_guides': 3,
        }
        
        response = self.api.handle_design(request)
        
        self.assertTrue(response['success'])
        self.assertGreater(response['count'], 0)
        self.assertIn('data', response)
    
    def test_handle_offtarget(self):
        """Test API off-target handler."""
        request = {
            'guide': 'GCTGCGGAGACCTGGAGAGA',
            'offtargets': [
                'GCTGCGGAGACCTGGAGAGA',
                'ACTGCGGAGACCTGGAGAGA',
            ],
        }
        
        response = self.api.handle_offtarget(request)
        
        self.assertTrue(response['success'])
        self.assertGreater(response['count'], 0)
    
    def test_error_handling_missing_guide(self):
        """Test error handling for missing guide."""
        request = {}  # Missing 'guide'
        
        response = self.api.handle_score(request)
        
        self.assertIn('error', response)
    
    def test_error_handling_invalid_sequence(self):
        """Test error handling for invalid sequences."""
        request = {
            'guide': 'GCTXCGGAGACCTGGAGAGA',  # Invalid base X
        }
        
        response = self.api.handle_score(request)
        
        self.assertFalse(response['success'])
        self.assertIn('error', response)


class TestScientificGates(unittest.TestCase):
    """Test enforcement of scientific gates."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = DisruptionAnalyzer(seed=42)
    
    def test_gate_human_dna_only(self):
        """Test that only valid DNA/RNA bases are accepted."""
        # Valid DNA
        valid_dna = "ACGTACGTACGT"
        result = self.analyzer.score_guide(valid_dna)
        self.assertIn('disruption_score', result)
        
        # Valid RNA
        valid_rna = "ACGUACGUACGU"
        result = self.analyzer.score_guide(valid_rna)
        self.assertIn('disruption_score', result)
        
        # Invalid: U in DNA context
        with self.assertRaises(ValueError):
            self.analyzer.score_guide("ACGUACGTACGT")  # Mixed U and T
    
    def test_gate_no_fabrication(self):
        """Test that fabricated/invalid bases are rejected."""
        invalid_sequences = [
            "ACGTXYZ",  # Invalid bases
            "ACGT123",  # Numbers
            "ACGT---",  # Dashes
        ]
        
        for seq in invalid_sequences:
            with self.assertRaises(ValueError):
                self.analyzer.score_guide(seq)
    
    def test_gate_z_invariant_range(self):
        """Test that Z-invariant scores are in valid range."""
        guides = [
            "GCTGCGGAGACCTGGAGAGA",
            "ATATATATATATATATATAT",
            "GCGCGCGCGCGCGCGCGCGC",
        ]
        
        for guide in guides:
            result = self.analyzer.score_guide(guide)
            score = result['disruption_score']
            
            # Z-score should be in [0, 1] due to sigmoid
            self.assertGreaterEqual(score, 0.0)
            self.assertLessEqual(score, 1.0)
    
    def test_gate_reproducibility(self):
        """Test that results are reproducible with seed control."""
        guide = "GCTGCGGAGACCTGGAGAGA"
        
        # Same seed should give identical results
        analyzer1 = DisruptionAnalyzer(seed=123)
        result1 = analyzer1.score_guide(guide)
        
        analyzer2 = DisruptionAnalyzer(seed=123)
        result2 = analyzer2.score_guide(guide)
        
        self.assertEqual(result1['disruption_score'], result2['disruption_score'])


def run_tests():
    """Run all tests and return results."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestDisruptionAnalyzer))
    suite.addTests(loader.loadTestsFromTestCase(TestGenomicDisruptionAPI))
    suite.addTests(loader.loadTestsFromTestCase(TestScientificGates))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result


if __name__ == '__main__':
    result = run_tests()
    
    # Exit with non-zero status if tests failed
    if not result.wasSuccessful():
        sys.exit(1)
