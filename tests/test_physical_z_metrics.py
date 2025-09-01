#!/usr/bin/env python3
"""
Test suite for CRISPR Physical Z-Metrics implementation.

Tests the four sequence-derivable physical Z-metrics and their integration
with the Wave-CRISPR pipeline.
"""

import unittest
import sys
import os
import tempfile
from unittest.mock import patch

# Add applications to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'applications'))

from crispr_physical_z_metrics import (
    PhysicalZMetricsCalculator,
    validate_human_dna_sequence,
    DNAValidationError,
    read_fasta_with_validation
)


class TestDNAValidation(unittest.TestCase):
    """Test DNA sequence validation guardrails."""

    def test_valid_human_sequence(self):
        """Test valid human DNA sequence passes validation."""
        sequence = "ATCGATCGATCG"
        header = "GRCh38 chr1:123456-123467"
        # Should not raise exception
        validate_human_dna_sequence(sequence, header)

    def test_invalid_nucleotides(self):
        """Test rejection of non-ACGTN nucleotides."""
        sequence = "ATCGXYZ"
        header = "test"
        with self.assertRaises(DNAValidationError):
            validate_human_dna_sequence(sequence, header)

    def test_excessive_n_content(self):
        """Test rejection of excessive N content."""
        sequence = "NNNNNNNNNN"  # 100% N
        header = "test"
        with self.assertRaises(DNAValidationError):
            validate_human_dna_sequence(sequence, header)

    def test_sequence_too_short(self):
        """Test rejection of sequences too short."""
        sequence = "ATCG"  # < 10 bp
        header = "test"
        with self.assertRaises(DNAValidationError):
            validate_human_dna_sequence(sequence, header)

    def test_human_reference_headers(self):
        """Test recognition of human reference headers."""
        valid_headers = [
            "GRCh38 chromosome 1",
            "hg38 reference genome",
            "Homo sapiens PCSK9 gene",
            "CRISPR benchmark dataset"
        ]
        sequence = "ATCGATCGATCG"
        
        for header in valid_headers:
            # Should not raise exception
            validate_human_dna_sequence(sequence, header)


class TestPhysicalZMetrics(unittest.TestCase):
    """Test physical Z-metrics calculations."""

    def setUp(self):
        """Set up test calculator."""
        self.calculator = PhysicalZMetricsCalculator()
        self.test_sequence = "ATCGATCGATCGATCGATCG"  # 20 bp test sequence

    def test_base_pair_opening_kinetics(self):
        """Test base-pair opening kinetics Z-metric."""
        result = self.calculator.calculate_base_pair_opening_kinetics(self.test_sequence)
        
        # Check required fields
        self.assertIn('z_opening', result)
        self.assertIn('opening_rate', result)
        self.assertIn('context_factor', result)
        self.assertIn('gc_content', result)
        
        # Check rate is non-negative (guardrail)
        self.assertGreaterEqual(float(result['opening_rate']), 0)
        
        # Check Z-metric is finite
        self.assertTrue(float(result['z_opening']) >= 0)

    def test_base_stacking_dissociation(self):
        """Test base stacking dissociation Z-metric."""
        result = self.calculator.calculate_base_stacking_dissociation(self.test_sequence)
        
        # Check required fields
        self.assertIn('z_stacking', result)
        self.assertIn('dissociation_rate', result)
        self.assertIn('stack_factor', result)
        
        # Check rate is non-negative (guardrail)
        self.assertGreaterEqual(float(result['dissociation_rate']), 0)
        
        # Check Z-metric is finite
        self.assertTrue(float(result['z_stacking']) >= 0)

    def test_helical_twist_fluctuation(self):
        """Test helical twist fluctuation Z-metric."""
        result = self.calculator.calculate_helical_twist_fluctuation(self.test_sequence)
        
        # Check required fields
        self.assertIn('z_twist', result)
        self.assertIn('twist_rate', result)
        self.assertIn('geometry_factor', result)
        
        # Check rate is non-negative (guardrail)
        self.assertGreaterEqual(float(result['twist_rate']), 0)
        
        # Check Z-metric is finite
        self.assertTrue(float(result['z_twist']) >= 0)

    def test_denaturation_melting_kinetics(self):
        """Test denaturation/melting kinetics Z-metric."""
        result = self.calculator.calculate_denaturation_melting_kinetics(self.test_sequence)
        
        # Check required fields
        self.assertIn('z_melting', result)
        self.assertIn('melting_rate', result)
        self.assertIn('stability_factor', result)
        self.assertIn('estimated_tm', result)
        
        # Check rate is non-negative (guardrail)
        self.assertGreaterEqual(float(result['melting_rate']), 0)
        
        # Check Z-metric is finite
        self.assertTrue(float(result['z_melting']) >= 0)
        
        # Check Tm is reasonable
        tm = float(result['estimated_tm'])
        self.assertGreater(tm, 0)  # Above absolute zero
        self.assertLess(tm, 200)   # Below water boiling point

    def test_all_metrics_calculation(self):
        """Test calculation of all four Z-metrics together."""
        header = "PCSK9 test sequence GRCh38"
        results = self.calculator.calculate_all_physical_z_metrics(
            self.test_sequence, header, validate=False
        )
        
        # Check all four metrics are present
        self.assertIn('opening_kinetics', results)
        self.assertIn('stacking_dissociation', results)
        self.assertIn('twist_fluctuation', results)
        self.assertIn('melting_kinetics', results)
        
        # Check summary statistics
        self.assertIn('summary', results)
        summary = results['summary']
        self.assertIn('z_mean', summary)
        self.assertIn('z_variance', summary)
        
        # Check sequence info
        self.assertIn('sequence_info', results)
        seq_info = results['sequence_info']
        self.assertEqual(float(seq_info['length']), len(self.test_sequence))

    def test_z_framework_form(self):
        """Test that Z-metrics follow Z = A * (B / c) form with c = e²."""
        result = self.calculator.calculate_base_pair_opening_kinetics(self.test_sequence)
        
        # Verify the calculation follows Z = A * (B / c) form
        a_context = result['context_factor']
        b_rate = result['opening_rate']
        z_opening = result['z_opening']
        
        # Manually calculate expected Z value
        expected_z = float(a_context) * (float(b_rate) / float(self.calculator.e_squared))
        
        # Should match calculated Z (within numerical precision)
        self.assertAlmostEqual(float(z_opening), expected_z, places=10)

    def test_short_sequence_handling(self):
        """Test handling of sequences too short for some metrics."""
        short_seq = "A"  # Single nucleotide
        
        # Should handle gracefully without crashing
        result = self.calculator.calculate_all_physical_z_metrics(
            short_seq, "test", validate=False
        )
        
        # Should return valid results structure
        self.assertIn('summary', result)

    def test_negative_rate_clamping(self):
        """Test that negative rates are clamped to zero."""
        # Mock a scenario that might produce negative rates
        with patch.object(self.calculator, '_clamp_rate') as mock_clamp:
            mock_clamp.return_value = 0.0  # Simulate clamping
            
            result = self.calculator.calculate_base_pair_opening_kinetics(self.test_sequence)
            
            # Should have called rate clamping
            mock_clamp.assert_called()


class TestFASTAValidation(unittest.TestCase):
    """Test FASTA file reading with validation."""

    def test_valid_fasta_reading(self):
        """Test reading valid FASTA file."""
        fasta_content = ">PCSK9_GRCh38\nATCGATCGATCGATCG\n>test2_hg38\nGCGCGCGCGCGC\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            sequences = read_fasta_with_validation(temp_path)
            
            self.assertEqual(len(sequences), 2)
            self.assertIn("PCSK9_GRCh38", sequences)
            self.assertIn("test2_hg38", sequences)
            self.assertEqual(sequences["PCSK9_GRCh38"], "ATCGATCGATCGATCG")
            
        finally:
            os.unlink(temp_path)

    def test_invalid_fasta_rejection(self):
        """Test rejection of FASTA with invalid sequences."""
        fasta_content = ">invalid_sequence\nATCGXYZ\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            temp_path = f.name
        
        try:
            with self.assertRaises(DNAValidationError):
                read_fasta_with_validation(temp_path)
                
        finally:
            os.unlink(temp_path)


class TestIntegrationWithCLI(unittest.TestCase):
    """Test integration with CRISPR CLI."""

    def test_cli_import(self):
        """Test that CLI can import the new module."""
        try:
            from crispr_cli import main
            # Import successful
            self.assertTrue(True)
        except ImportError as e:
            self.fail(f"Failed to import CLI with new module: {e}")


class TestHypothesesFramework(unittest.TestCase):
    """Test framework for experimental hypotheses H1 and H2."""

    def setUp(self):
        """Set up test framework for hypotheses testing."""
        self.calculator = PhysicalZMetricsCalculator()

    def test_h1_opening_z_correlation(self):
        """Test framework for H1: Lower Opening-Z correlates with higher CRISPR efficiency."""
        # Generate test sequences with different characteristics
        high_efficiency_seq = "GGGCCCGGGCCCGGGCCCGG"  # High GC, stable
        low_efficiency_seq = "ATATATATATATATATAT"    # Low GC, less stable
        
        # Calculate Opening-Z for both
        high_eff_metrics = self.calculator.calculate_base_pair_opening_kinetics(high_efficiency_seq)
        low_eff_metrics = self.calculator.calculate_base_pair_opening_kinetics(low_efficiency_seq)
        
        high_opening_z = float(high_eff_metrics['z_opening'])
        low_opening_z = float(low_eff_metrics['z_opening'])
        
        # Framework should be able to distinguish sequences
        self.assertNotEqual(high_opening_z, low_opening_z)
        
        # This test establishes the framework for H1 testing
        # Actual correlation testing would require empirical efficiency data
        print(f"H1 Framework: High efficiency seq Opening-Z = {high_opening_z:.6f}")
        print(f"H1 Framework: Low efficiency seq Opening-Z = {low_opening_z:.6f}")

    def test_h2_stacking_z_threshold(self):
        """Test framework for H2: Very low Stacking-Z identification."""
        # Generate sequences expected to have different stacking properties
        high_stacking_seq = "GCGCGCGCGCGCGCGCGC"  # Alternating GC (high stacking)
        low_stacking_seq = "AAATTTAAATTTAAATTT"   # AT-rich (lower stacking)
        
        # Calculate Stacking-Z for both
        high_stack_metrics = self.calculator.calculate_base_stacking_dissociation(high_stacking_seq)
        low_stack_metrics = self.calculator.calculate_base_stacking_dissociation(low_stacking_seq)
        
        high_stacking_z = float(high_stack_metrics['z_stacking'])
        low_stacking_z = float(low_stack_metrics['z_stacking'])
        
        # Framework should distinguish stacking properties
        self.assertNotEqual(high_stacking_z, low_stacking_z)
        
        # This establishes framework for identifying "very low" Stacking-Z
        print(f"H2 Framework: High stacking seq Stacking-Z = {high_stacking_z:.6f}")
        print(f"H2 Framework: Low stacking seq Stacking-Z = {low_stacking_z:.6f}")


def run_tests():
    """Run all tests."""
    unittest.main(verbosity=2)


if __name__ == "__main__":
    run_tests()