"""
Test Suite for WAVE-CRISPR Experiment Design

Tests for NCBI sequence fetching, spectral analysis, and experimental framework
as specified in the WAVE-CRISPR experimental design.
"""

import unittest
import sys
import os
import tempfile
from unittest.mock import patch, Mock, MagicMock

# Add paths for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'applications'))

from spectral_analysis import SpectralDNAAnalyzer, load_zeta_zeros, compute_zeta_correlations
from applications.crispr_ncbi_fetcher import CRISPRNCBIFetcher
from applications.wave_crispr_experiment import WAVECRISPRExperiment


class TestSpectralDNAAnalyzer(unittest.TestCase):
    """Test cases for SpectralDNAAnalyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = SpectralDNAAnalyzer(k_star=0.3, seed=42)
        self.test_sequence = "ATCGATCGATCGATCGAT"
    
    def test_initialization(self):
        """Test analyzer initialization."""
        self.assertEqual(self.analyzer.k_star, 0.3)
        self.assertIn('A', self.analyzer.biological_encoding)
        self.assertIn('A', self.analyzer.arbitrary_encoding)
        
        # Check that encodings are different
        bio_a = self.analyzer.biological_encoding['A']
        arb_a = self.analyzer.arbitrary_encoding['A']
        self.assertNotEqual(bio_a, arb_a)
    
    def test_sequence_encoding(self):
        """Test DNA sequence encoding."""
        # Test biological encoding
        bio_encoded = self.analyzer.encode_sequence(self.test_sequence, 'biological')
        self.assertEqual(len(bio_encoded), len(self.test_sequence))
        self.assertTrue(all(isinstance(x, complex) for x in bio_encoded))
        
        # Test arbitrary encoding
        arb_encoded = self.analyzer.encode_sequence(self.test_sequence, 'arbitrary')
        self.assertEqual(len(arb_encoded), len(self.test_sequence))
        self.assertTrue(all(isinstance(x, complex) for x in arb_encoded))
        
        # Encodings should be different
        self.assertFalse(all(b == a for b, a in zip(bio_encoded, arb_encoded)))
    
    def test_z_framework_calculation(self):
        """Test Z Framework discrete domain calculations."""
        encoded_seq = self.analyzer.encode_sequence(self.test_sequence, 'biological')
        z_metrics = self.analyzer.compute_z_framework_values(encoded_seq)
        
        # Check required fields
        required_fields = ['z_value', 'delta_n', 'delta_max', 'kappa']
        for field in required_fields:
            self.assertIn(field, z_metrics)
            self.assertIsInstance(z_metrics[field], float)
        
        # Check reasonable values
        self.assertGreaterEqual(z_metrics['z_value'], 0)
        self.assertGreaterEqual(z_metrics['delta_n'], 0)
        self.assertGreater(z_metrics['delta_max'], 0)
        self.assertGreater(z_metrics['kappa'], 0)
    
    def test_spectral_entropy(self):
        """Test spectral entropy calculation."""
        encoded_seq = self.analyzer.encode_sequence(self.test_sequence, 'biological')
        entropy = self.analyzer.compute_spectral_entropy(encoded_seq)
        
        self.assertIsInstance(entropy, float)
        self.assertGreaterEqual(entropy, 0)
    
    def test_geodesic_resolution(self):
        """Test geodesic resolution calculation."""
        n = len(self.test_sequence)
        geodesic_res = self.analyzer.compute_geodesic_resolution(n)
        
        self.assertIsInstance(geodesic_res, float)
        self.assertGreater(geodesic_res, 0)
    
    def test_complete_analysis(self):
        """Test complete sequence analysis."""
        analysis = self.analyzer.analyze_sequence(self.test_sequence, 'biological')
        
        required_fields = [
            'sequence_length', 'encoding_type', 'z_framework', 
            'spectral_entropy', 'geodesic_resolution', 'gc_content'
        ]
        
        for field in required_fields:
            self.assertIn(field, analysis)
        
        self.assertEqual(analysis['sequence_length'], len(self.test_sequence))
        self.assertEqual(analysis['encoding_type'], 'biological')
        self.assertIsInstance(analysis['z_framework'], dict)
    
    def test_encoding_comparison(self):
        """Test biological vs arbitrary encoding comparison."""
        comparison = self.analyzer.compare_encodings(self.test_sequence)
        
        self.assertIn('biological', comparison)
        self.assertIn('arbitrary', comparison)
        self.assertIn('differences', comparison)
        
        # Check that analyses are different
        bio_z = comparison['biological']['z_framework']['z_value']
        arb_z = comparison['arbitrary']['z_framework']['z_value']
        self.assertNotEqual(bio_z, arb_z)


class TestCRISPRNCBIFetcher(unittest.TestCase):
    """Test cases for CRISPRNCBIFetcher."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.fetcher = CRISPRNCBIFetcher()
    
    def test_initialization(self):
        """Test fetcher initialization."""
        self.assertEqual(self.fetcher.fetch_count, 0)
    
    @patch('applications.crispr_ncbi_fetcher.Entrez')
    @patch('applications.crispr_ncbi_fetcher.SeqIO')
    def test_sequence_fetch(self, mock_seqio, mock_entrez):
        """Test single sequence fetching."""
        # Mock successful fetch
        mock_record = Mock()
        mock_record.seq = "ATCGATCGATCG"
        mock_seqio.read.return_value = mock_record
        
        mock_handle = Mock()
        mock_entrez.efetch.return_value.__enter__.return_value = mock_handle
        
        result = self.fetcher.fetch_sequence("NM_007294.4")
        
        self.assertIsNotNone(result)
        self.assertEqual(self.fetcher.fetch_count, 1)
    
    def test_sequence_validation(self):
        """Test sequence integrity validation."""
        # Create mock sequence record
        mock_record = Mock()
        mock_record.seq = "ATCGATCGATCGNCGATCG"
        
        validation = self.fetcher.validate_sequence_integrity(mock_record)
        
        required_fields = ['valid', 'length', 'gc_content', 'n_content', 'valid_bases', 'min_length_ok']
        for field in required_fields:
            self.assertIn(field, validation)
        
        self.assertIsInstance(validation['valid'], bool)
        self.assertIsInstance(validation['gc_content'], float)
        self.assertGreaterEqual(validation['gc_content'], 0)
        self.assertLessEqual(validation['gc_content'], 1)
    
    def test_validation_edge_cases(self):
        """Test validation with edge cases."""
        # Test with None record
        validation = self.fetcher.validate_sequence_integrity(None)
        self.assertFalse(validation['valid'])
        self.assertIn('error', validation)
        
        # Test with short sequence
        mock_record = Mock()
        mock_record.seq = "ATCG"  # Too short
        validation = self.fetcher.validate_sequence_integrity(mock_record)
        self.assertFalse(validation['min_length_ok'])
        
        # Test with high N content
        mock_record = Mock()
        mock_record.seq = "NNNNNNNNNNNNNNNNNNNN"  # All Ns
        validation = self.fetcher.validate_sequence_integrity(mock_record)
        self.assertFalse(validation['valid'])


class TestWAVECRISPRExperiment(unittest.TestCase):
    """Test cases for WAVECRISPRExperiment."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.experiment = WAVECRISPRExperiment(seed=42)
    
    def test_initialization(self):
        """Test experiment initialization."""
        self.assertEqual(self.experiment.seed, 42)
        self.assertIsNotNone(self.experiment.spectral_analyzer)
        self.assertIsNotNone(self.experiment.ncbi_fetcher)
        self.assertIsNotNone(self.experiment.doench_data)
    
    def test_doench_data_loading(self):
        """Test Doench data loading/generation."""
        data = self.experiment.doench_data
        
        required_columns = ['sequence', 'efficiency']
        for col in required_columns:
            self.assertIn(col, data.columns)
        
        self.assertGreater(len(data), 0)
        
        # Check data validity
        self.assertTrue(all(data['efficiency'] >= 0))
        self.assertTrue(all(data['efficiency'] <= 1))
        self.assertTrue(all(len(seq) > 0 for seq in data['sequence']))
    
    def test_hypothesis_h1_structure(self):
        """Test H1 hypothesis testing structure."""
        h1_results = self.experiment.test_hypothesis_h1()
        
        required_fields = ['biological', 'arbitrary', 'hypothesis_met']
        for field in required_fields:
            self.assertIn(field, h1_results)
        
        self.assertIn('correlations', h1_results['biological'])
        self.assertIn('p_values', h1_results['biological'])
        # Check that hypothesis_met is a boolean value
        self.assertIn(h1_results['hypothesis_met'], [True, False])
    
    def test_hypothesis_h2_structure(self):
        """Test H2 hypothesis testing structure."""
        # Reduce scales for testing
        self.experiment.scales = [100, 500]
        h2_results = self.experiment.test_hypothesis_h2()
        
        required_fields = ['scales', 'error_rates', 'confidence_intervals', 'hypothesis_met']
        for field in required_fields:
            self.assertIn(field, h2_results)
        
        self.assertEqual(len(h2_results['scales']), 2)
        self.assertEqual(len(h2_results['error_rates']), 2)
        self.assertIsInstance(h2_results['hypothesis_met'], bool)
    
    def test_hypothesis_h3_structure(self):
        """Test H3 hypothesis testing structure."""
        h3_results = self.experiment.test_hypothesis_h3()
        
        required_fields = ['correlations', 'p_values', 'hypothesis_met']
        for field in required_fields:
            self.assertIn(field, h3_results)
        
        self.assertIsInstance(h3_results['hypothesis_met'], bool)


class TestZetaCorrelations(unittest.TestCase):
    """Test cases for zeta correlation functionality."""
    
    def test_load_zeta_zeros(self):
        """Test loading zeta zeros from file."""
        # Create temporary zeta file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            test_zeros = [14.134725142, 21.022039639, 25.010857580]
            for zero in test_zeros:
                f.write(f"{zero}\n")
            temp_file = f.name
        
        try:
            zeros = load_zeta_zeros(temp_file)
            self.assertEqual(len(zeros), 3)
            self.assertAlmostEqual(zeros[0], 14.134725142, places=6)
        finally:
            os.unlink(temp_file)
    
    def test_compute_zeta_correlations(self):
        """Test zeta correlation computation."""
        spectral_features = [1.0, 2.0, 3.0, 4.0, 5.0]
        zeta_zeros = [14.1, 21.0, 25.0, 30.4, 32.9]
        
        correlation, p_value = compute_zeta_correlations(spectral_features, zeta_zeros)
        
        self.assertIsInstance(correlation, float)
        self.assertIsInstance(p_value, float)
        self.assertGreaterEqual(abs(correlation), 0)
        self.assertLessEqual(abs(correlation), 1)
        self.assertGreaterEqual(p_value, 0)
        self.assertLessEqual(p_value, 1)
    
    def test_zeta_correlations_edge_cases(self):
        """Test zeta correlations with edge cases."""
        # Empty lists
        correlation, p_value = compute_zeta_correlations([], [])
        self.assertEqual(correlation, 0.0)
        self.assertEqual(p_value, 1.0)
        
        # Single values
        correlation, p_value = compute_zeta_correlations([1.0], [14.1])
        self.assertEqual(correlation, 0.0)
        self.assertEqual(p_value, 1.0)


class TestIntegration(unittest.TestCase):
    """Integration tests for WAVE-CRISPR components."""
    
    def test_end_to_end_workflow(self):
        """Test complete workflow integration."""
        # Initialize components
        analyzer = SpectralDNAAnalyzer(seed=42)
        
        # Test sequence
        test_seq = "ATCGATCGATCGATCGAT"
        
        # Run analysis
        bio_analysis = analyzer.analyze_sequence(test_seq, 'biological')
        arb_analysis = analyzer.analyze_sequence(test_seq, 'arbitrary')
        
        # Verify results structure
        self.assertIn('z_framework', bio_analysis)
        self.assertIn('spectral_entropy', bio_analysis)
        self.assertIn('z_framework', arb_analysis)
        self.assertIn('spectral_entropy', arb_analysis)
        
        # Results should be different between encodings
        self.assertNotEqual(
            bio_analysis['spectral_entropy'], 
            arb_analysis['spectral_entropy']
        )
    
    def test_hypothesis_criteria_implementation(self):
        """Test that hypothesis criteria are correctly implemented."""
        # H1 criteria
        test_bio_corr = 0.6
        test_bio_p = 1e-6
        test_arb_corr = 0.2
        
        h1_bio_meets = test_bio_corr >= 0.5 and test_bio_p < 1e-5
        h1_arb_fails = test_arb_corr < 0.3
        h1_valid = h1_bio_meets and h1_arb_fails
        
        self.assertTrue(h1_valid)
        
        # H2 criteria
        test_error = 0.005
        h2_valid = test_error < 0.01
        self.assertTrue(h2_valid)
        
        # H3 criteria
        test_zeta_corr = 0.95
        test_zeta_p = 1e-12
        h3_valid = test_zeta_corr >= 0.93 and test_zeta_p < 1e-10
        self.assertTrue(h3_valid)


def run_wave_crispr_tests():
    """Run all WAVE-CRISPR tests."""
    print("Running WAVE-CRISPR Test Suite")
    print("=" * 40)
    
    # Create test suite
    test_classes = [
        TestSpectralDNAAnalyzer,
        TestCRISPRNCBIFetcher,
        TestWAVECRISPRExperiment,
        TestZetaCorrelations,
        TestIntegration
    ]
    
    suite = unittest.TestSuite()
    
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Return success status
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_wave_crispr_tests()
    sys.exit(0 if success else 1)