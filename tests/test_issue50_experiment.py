#!/usr/bin/env python3
"""
Test suite for Issue #50 Experiment Framework

Validates the experimental infrastructure and hypothesis testing framework
for the WAVE-CRISPR Z-Framework claims.
"""

import unittest
import sys
import os
import numpy as np
import pandas as pd
import json
from pathlib import Path

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from experiments.issue50_experiment import Issue50ExperimentFramework


class TestIssue50Experiment(unittest.TestCase):
    """Test cases for Issue #50 experiment framework."""
    
    def setUp(self):
        """Set up test environment."""
        self.experiment = Issue50ExperimentFramework(random_seed=42)
    
    def test_experiment_initialization(self):
        """Test that the experiment initializes correctly."""
        self.assertEqual(self.experiment.random_seed, 42)
        self.assertIsNotNone(self.experiment.bio_encoder)
        self.assertIsNotNone(self.experiment.arbitrary_encoder)
        self.assertIsNotNone(self.experiment.z_calculator)
        self.assertIn('experiment_metadata', self.experiment.results)
    
    def test_synthetic_dataset_creation(self):
        """Test synthetic dataset creation."""
        df = self.experiment._create_synthetic_dataset(n_sequences=100)
        
        self.assertEqual(len(df), 100)
        self.assertIn('sequence', df.columns)
        self.assertIn('efficiency', df.columns)
        self.assertIn('gene_group', df.columns)
        
        # Check sequence format
        self.assertTrue(all(len(seq) == 20 for seq in df['sequence']))
        self.assertTrue(all(set(seq).issubset(set('ATCG')) for seq in df['sequence']))
        
        # Check efficiency range
        self.assertTrue(all(0 <= eff <= 1 for eff in df['efficiency']))
    
    def test_feature_extraction(self):
        """Test Z-Framework feature extraction."""
        test_sequence = "ATCGATCGATCGATCGATCG"
        
        # Test bio-anchored features
        bio_features = self.experiment.extract_features(test_sequence, 'bio', k=0.3)
        self.assertIn('spectral_entropy', bio_features)
        self.assertIn('z_mean', bio_features)
        self.assertIn('delta_phi', bio_features)
        self.assertIn('gc_content', bio_features)
        self.assertIn('phase_bit', bio_features)
        
        # Test arbitrary features
        arb_features = self.experiment.extract_features(test_sequence, 'arbitrary', k=0.3)
        self.assertIn('spectral_entropy', arb_features)
        self.assertIn('z_mean', arb_features)
        
        # Features should be different between encodings
        self.assertNotEqual(bio_features['spectral_entropy'], arb_features['spectral_entropy'])
    
    def test_baseline_features(self):
        """Test baseline feature creation."""
        test_sequence = "ATCGATCGATCGATCGATCG"
        baseline_features = self.experiment.create_baseline_features(test_sequence)
        
        self.assertIn('gc_content', baseline_features)
        self.assertIn('at_content', baseline_features)
        self.assertIn('sequence_length', baseline_features)
        
        # Check GC content calculation
        expected_gc = 0.5  # Half A/T, half G/C
        self.assertAlmostEqual(baseline_features['gc_content'], expected_gc, places=2)
    
    def test_h1_hypothesis_test(self):
        """Test H1 hypothesis testing framework."""
        # Create small test dataset
        df = self.experiment._create_synthetic_dataset(n_sequences=50)
        
        # Run H1 test
        h1_results = self.experiment.run_hypothesis_h1_lift_vs_baseline(df, scale=50)
        
        # Check result structure
        self.assertIn('bio_band', h1_results)
        self.assertIn('arbitrary_band', h1_results)
        
        for band in ['bio_band', 'arbitrary_band']:
            band_result = h1_results[band]
            self.assertIn('rmr_mean', band_result)
            self.assertIn('rmr_ci_lower', band_result)
            self.assertIn('rmr_ci_upper', band_result)
            self.assertIn('p_value', band_result)
            self.assertIn('passes_h1', band_result)
            
            # Check that confidence intervals are ordered correctly
            self.assertLessEqual(band_result['rmr_ci_lower'], band_result['rmr_mean'])
            self.assertLessEqual(band_result['rmr_mean'], band_result['rmr_ci_upper'])
    
    def test_results_saving(self):
        """Test that results are saved correctly."""
        # Run a minimal experiment
        df = self.experiment._create_synthetic_dataset(n_sequences=30)
        h1_results = self.experiment.run_hypothesis_h1_lift_vs_baseline(df, scale=30)
        self.experiment.results['hypothesis_tests']['h1_scale_30'] = h1_results
        
        # Save results
        self.experiment.save_results()
        
        # Check files exist
        self.assertTrue(Path('results/issue50_summary.json').exists())
        self.assertTrue(Path('results/issue50_table.csv').exists())
        
        # Check JSON content
        with open('results/issue50_summary.json', 'r') as f:
            summary = json.load(f)
        
        self.assertIn('experiment_metadata', summary)
        self.assertIn('hypothesis_tests', summary)
        self.assertEqual(summary['experiment_metadata']['random_seed'], 42)
        
        # Check CSV content
        df_results = pd.read_csv('results/issue50_table.csv')
        self.assertIn('scale', df_results.columns)
        self.assertIn('hypothesis', df_results.columns)
        self.assertIn('band', df_results.columns)
        self.assertIn('rmr_mean', df_results.columns)
    
    def test_reproducibility(self):
        """Test that experiments are reproducible with same seed."""
        # Create two identical experiments
        exp1 = Issue50ExperimentFramework(random_seed=123)
        exp2 = Issue50ExperimentFramework(random_seed=123)
        
        # Generate datasets
        df1 = exp1._create_synthetic_dataset(n_sequences=20)
        df2 = exp2._create_synthetic_dataset(n_sequences=20)
        
        # Should be identical
        pd.testing.assert_frame_equal(df1, df2)
        
        # Test feature extraction reproducibility
        test_seq = "ATCGATCGATCGATCGATCG"
        features1 = exp1.extract_features(test_seq, 'bio')
        features2 = exp2.extract_features(test_seq, 'bio')
        
        for key in features1:
            if isinstance(features1[key], (int, float)):
                self.assertAlmostEqual(features1[key], features2[key], places=6)
    
    def test_data_loading(self):
        """Test CRISPR data loading functionality."""
        # Test loading real dataset
        df = self.experiment.load_crispr_data()
        
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn('sequence', df.columns)
        self.assertIn('efficiency', df.columns)
        self.assertIn('gene_group', df.columns)
        self.assertGreater(len(df), 0)
    
    def test_experiment_metadata(self):
        """Test experiment metadata is properly recorded."""
        metadata = self.experiment.results['experiment_metadata']
        
        self.assertIn('timestamp', metadata)
        self.assertIn('random_seed', metadata)
        self.assertIn('scales', metadata)
        self.assertIn('k_values', metadata)
        self.assertIn('n_bootstrap', metadata)
        self.assertIn('n_permutations', metadata)
        
        # Check values
        self.assertEqual(metadata['random_seed'], 42)
        self.assertEqual(metadata['scales'], [100, 1000, 10000])
        self.assertEqual(metadata['n_bootstrap'], 1000)
    
    def test_visualization_creation(self):
        """Test that visualizations are created correctly."""
        # Set up minimal results
        self.experiment.results['hypothesis_tests']['h1_scale_100'] = {
            'bio_band': {'rmr_mean': 0.20},
            'arbitrary_band': {'rmr_mean': 0.10}
        }
        
        # Create visualization
        self.experiment._create_simple_visualization()
        
        # Check file exists
        self.assertTrue(Path('results/issue50_demo.png').exists())


class TestExperimentIntegration(unittest.TestCase):
    """Integration tests for the complete experiment workflow."""
    
    def test_complete_workflow(self):
        """Test the complete experimental workflow."""
        experiment = Issue50ExperimentFramework(random_seed=42)
        
        # Run complete experiment (simplified version)
        results = experiment.run_complete_experiment()
        
        # Check results structure
        self.assertIsInstance(results, dict)
        self.assertIn('experiment_metadata', results)
        self.assertIn('hypothesis_tests', results)
        
        # Check that H1 was tested
        h1_key = 'h1_scale_100'
        self.assertIn(h1_key, results['hypothesis_tests'])
        
        h1_results = results['hypothesis_tests'][h1_key]
        self.assertIn('bio_band', h1_results)
        self.assertIn('arbitrary_band', h1_results)
        
        # Check files were created
        self.assertTrue(Path('results/issue50_summary.json').exists())
        self.assertTrue(Path('results/issue50_table.csv').exists())
        self.assertTrue(Path('results/issue50_demo.png').exists())


def main():
    """Run all tests."""
    # Create results directory if it doesn't exist
    Path('results').mkdir(exist_ok=True)
    Path('logs').mkdir(exist_ok=True)
    
    # Run tests
    unittest.main(verbosity=2)


if __name__ == '__main__':
    main()