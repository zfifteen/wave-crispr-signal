#!/usr/bin/env python3
"""
Test suite for Focused Ultrasound MVE experiment

Tests the minimal viable experiment for Z Framework spatial targeting precision
in simulated focused ultrasound applications.
"""

import sys
import tempfile
import unittest
from pathlib import Path
import subprocess
import json

# Add parent directory for imports
sys.path.append('.')
sys.path.append('..')

from experiments.focused_ultrasound_mve import (
    ExperimentConfig, 
    AcousticGrid,
    BaselineTargetingModel,
    ZFrameworkTargetingModel,
    FocusedUltrasoundExperiment,
    StatisticalAnalysis,
    get_git_commit_sha
)


class TestFocusedUltrasoundMVE(unittest.TestCase):
    """Test suite for Focused Ultrasound MVE components."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.config = ExperimentConfig(
            grid_size=20,  # Small for fast testing
            n_trials=10,   # Few trials for speed
            n_bootstrap=100,  # Reduced for speed but still > threshold in actual run
            n_permutation=100,  # Reduced for speed
            seed=42,
            output_dir="test_output"
        )
        
    def test_acoustic_grid_initialization(self):
        """Test acoustic grid initialization and properties."""
        grid = AcousticGrid(
            size=20,
            base_velocity=1540.0,
            variance=0.1,
            seed=42
        )
        
        # Check basic properties
        self.assertEqual(grid.size, 20)
        self.assertEqual(grid.base_velocity, 1540.0)
        self.assertEqual(grid.variance, 0.1)
        
        # Check velocity field properties
        self.assertEqual(grid.velocity_field.shape, (20, 20))
        self.assertGreater(grid.mean_velocity, 1400)  # Should be around base velocity
        self.assertLess(grid.mean_velocity, 1700)
        self.assertGreater(grid.std_velocity, 0)  # Should have some variation
        
    def test_baseline_targeting_model(self):
        """Test baseline targeting model calculations."""
        grid = AcousticGrid(20, 1540.0, 0.1, 42)
        model = BaselineTargetingModel(grid)
        
        # Test targeting calculation
        error, time_to_target = model.calculate_targeting_error(5, 5, 15, 15)
        
        # Check that results are reasonable
        self.assertGreater(error, 0)
        self.assertGreater(time_to_target, 0)
        self.assertIsInstance(error, float)
        self.assertIsInstance(time_to_target, float)
        
    def test_z_framework_targeting_model(self):
        """Test Z Framework targeting model calculations."""
        grid = AcousticGrid(20, 1540.0, 0.1, 42)
        model = ZFrameworkTargetingModel(grid, k_parameter=0.3)
        
        # Test targeting calculation
        error, time_to_target = model.calculate_targeting_error(5, 5, 15, 15)
        
        # Check that results are reasonable
        self.assertGreater(error, 0)
        self.assertGreater(time_to_target, 0)
        self.assertIsInstance(error, float)
        self.assertIsInstance(time_to_target, float)
        
        # Test theta_prime calculation
        theta_val = model.theta_prime(10, 0.3)
        self.assertGreater(theta_val, 0)
        self.assertLess(theta_val, 5)  # Should be reasonable value
        
    def test_focused_ultrasound_experiment(self):
        """Test complete experiment execution."""
        experiment = FocusedUltrasoundExperiment(self.config)
        
        # Run single trial
        result = experiment.run_single_trial(0)
        
        # Check result structure
        self.assertEqual(result.trial_id, 0)
        self.assertGreater(result.target_x, 0)
        self.assertGreater(result.target_y, 0)
        self.assertGreater(result.baseline_error, 0)
        self.assertGreater(result.z_framework_error, 0)
        self.assertGreater(result.baseline_time, 0)
        self.assertGreater(result.z_framework_time, 0)
        
    def test_statistical_analysis(self):
        """Test statistical analysis functionality."""
        # Run minimal experiment
        experiment = FocusedUltrasoundExperiment(self.config)
        results = experiment.run_experiment()
        
        # Test analysis
        analysis = StatisticalAnalysis(results, self.config)
        stats = analysis.analyze_results()
        
        # Check required statistical outputs
        required_keys = [
            'baseline_mean_error',
            'z_framework_mean_error', 
            'improvement_percentage',
            'cohens_d',
            'correlation_r',
            'permutation_p_value',
            'statistically_significant'
        ]
        
        for key in required_keys:
            self.assertIn(key, stats)
            
        # Check types
        self.assertIsInstance(stats['statistically_significant'], bool)
        self.assertIsInstance(stats['improvement_percentage'], float)
        
    def test_cli_interface_smoke(self):
        """Test CLI interface with minimal parameters."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test that CLI can be called without errors
            cmd = [
                'python', 'experiments/focused_ultrasound_mve.py',
                '--seed', '42',
                '--bootstrap', '1000', 
                '--permutation', '1000',
                '--splits', 'single',
                '--domain', 'discrete',
                '--k-parameter', '0.3',
                '--n-trials', '10',
                '--grid-size', '20',
                '--output-dir', temp_dir
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd='.')
            
            # Check that it completes successfully
            self.assertEqual(result.returncode, 0, 
                           f"CLI failed with stderr: {result.stderr}")
            
            # Check that output files are created
            results_dir = Path(temp_dir) / "focused_ultrasound_mve"
            self.assertTrue(results_dir.exists())
            
            # Find the run directory
            run_dirs = list(results_dir.glob("run-*"))
            self.assertGreater(len(run_dirs), 0)
            
            run_dir = run_dirs[0]
            
            # Check required files exist
            required_files = ['results.csv', 'analysis.json', 'metadata.json', 'experiment.log']
            for filename in required_files:
                file_path = run_dir / filename
                self.assertTrue(file_path.exists(), f"Missing file: {filename}")
                
            # Check analysis.json is valid JSON
            with open(run_dir / 'analysis.json') as f:
                analysis_data = json.load(f)
                self.assertIn('statistically_significant', analysis_data)
                
    def test_git_commit_format(self):
        """Test that git commit SHA is captured in correct format."""
        sha = get_git_commit_sha()
        
        # Should be either 'unknown' or a 40-character hex string
        if sha != 'unknown':
            self.assertEqual(len(sha), 40, f"Git SHA should be 40 characters, got {len(sha)}")
            self.assertTrue(all(c in '0123456789abcdef' for c in sha.lower()), 
                          f"Git SHA should be hex, got: {sha}")
        else:
            # In environments where git is not available, 'unknown' is acceptable
            self.assertEqual(sha, 'unknown')
            
    def test_z_framework_import_compatibility(self):
        """Test that z_framework module can be imported and instantiated."""
        try:
            from z_framework import ZFrameworkCalculator
            calc = ZFrameworkCalculator()
            self.assertIsNotNone(calc)
            # Test that basic methods exist that are used by the MVE
            self.assertTrue(hasattr(calc, 'calculate_z_values'))
            self.assertTrue(hasattr(calc, 'calculate_geodesic_resolution'))
        except ImportError as e:
            self.fail(f"Could not import ZFrameworkCalculator: {e}")
        except Exception as e:
            self.fail(f"Could not instantiate ZFrameworkCalculator: {e}")


class TestScientificGatesCompliance(unittest.TestCase):
    """Test compliance with scientific gates and repository policy."""
    
    def test_statistical_rigor_requirements(self):
        """Test that statistical rigor requirements are enforced."""
        from experiments.focused_ultrasound_mve import main
        import sys
        
        # Test that bootstrap < 1000 is rejected
        original_argv = sys.argv
        try:
            sys.argv = [
                'focused_ultrasound_mve.py',
                '--seed', '42',
                '--bootstrap', '999',  # Below threshold
                '--permutation', '1000',
                '--splits', 'single',
                '--domain', 'discrete'
            ]
            
            with self.assertRaises(ValueError) as context:
                main()
                
            self.assertIn("Bootstrap samples must be ≥1,000", str(context.exception))
            
        finally:
            sys.argv = original_argv
            
    def test_reproducibility_with_seed(self):
        """Test that results are reproducible with same seed."""
        config1 = ExperimentConfig(
            grid_size=20, n_trials=5, n_bootstrap=100, 
            n_permutation=100, seed=42
        )
        config2 = ExperimentConfig(
            grid_size=20, n_trials=5, n_bootstrap=100,
            n_permutation=100, seed=42
        )
        
        # Run two experiments with same seed
        exp1 = FocusedUltrasoundExperiment(config1)
        results1 = exp1.run_experiment()
        
        exp2 = FocusedUltrasoundExperiment(config2)
        results2 = exp2.run_experiment()
        
        # Compare first trial results (should be identical with same seed)
        self.assertAlmostEqual(results1[0].target_x, results2[0].target_x, places=6)
        self.assertAlmostEqual(results1[0].target_y, results2[0].target_y, places=6)
        self.assertAlmostEqual(results1[0].baseline_error, results2[0].baseline_error, places=6)


def run_smoke_test():
    """Run minimal smoke test for CI."""
    print("Running Focused Ultrasound MVE Smoke Test...")
    
    # Create minimal test
    config = ExperimentConfig(
        grid_size=10,
        n_trials=5, 
        n_bootstrap=1000,  # Meet scientific gates requirement
        n_permutation=1000,
        seed=42
    )
    
    # Run experiment
    experiment = FocusedUltrasoundExperiment(config)
    results = experiment.run_experiment()
    
    # Run analysis
    analysis = StatisticalAnalysis(results, config).analyze_results()
    
    print(f"✓ Experiment completed successfully")
    print(f"✓ Trials: {len(results)}")
    print(f"✓ Improvement: {analysis['improvement_percentage']:.1f}%")
    print(f"✓ Significant: {analysis['statistically_significant']}")
    print(f"✓ Runtime: {analysis.get('runtime_seconds', 0):.2f}s")
    
    # For smoke test, just check that experiment ran without errors
    # and shows some improvement (significance not required with few trials)
    runtime = analysis.get('runtime_seconds', 0)
    success = (
        len(results) > 0 and 
        analysis['improvement_percentage'] > 0 and
        runtime < 60
    )
    
    return success


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke', action='store_true', help='Run smoke test only')
    args = parser.parse_args()
    
    if args.smoke:
        success = run_smoke_test()
        exit(0 if success else 1)
    else:
        unittest.main()