#!/usr/bin/env python3
"""
Test Suite for FUS Enhancer: Vectorized Focused Ultrasound Targeting

This module provides comprehensive tests for the optimized fus_enhancer.py 
implementation, ensuring correctness, performance, and scientific validity.
"""

import unittest
import sys
import time
import tempfile
import shutil
from pathlib import Path

import torch
import numpy as np

# Add current directory to path for imports
sys.path.append('.')

from scripts.fus_enhancer import (
    VectorizedConfig, VectorizedAcousticGrid, VectorizedZFramework,
    VectorizedTargetingModels, VectorizedStatistics, VectorizedFUSExperiment
)


class TestVectorizedAcousticGrid(unittest.TestCase):
    """Test vectorized acoustic grid functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.device = 'cpu'  # Use CPU for consistent testing
        self.grid_size = 50
        self.base_velocity = 1540.0
        self.variance = 0.1
        self.seed = 42
        
    def test_grid_initialization(self):
        """Test acoustic grid initialization."""
        grid = VectorizedAcousticGrid(
            size=self.grid_size,
            base_velocity=self.base_velocity,
            variance=self.variance,
            seed=self.seed,
            device=self.device
        )
        
        self.assertEqual(grid.size, self.grid_size)
        self.assertEqual(grid.base_velocity, self.base_velocity)
        self.assertEqual(grid.variance, self.variance)
        self.assertEqual(grid.velocity_field.shape, (self.grid_size, self.grid_size))
        self.assertTrue(isinstance(grid.velocity_field, torch.Tensor))
        
    def test_velocity_sampling_batch(self):
        """Test batch velocity sampling along paths."""
        grid = VectorizedAcousticGrid(
            size=self.grid_size,
            base_velocity=self.base_velocity,
            variance=self.variance,
            seed=self.seed,
            device=self.device
        )
        
        batch_size = 10
        x_coords = torch.rand(batch_size, 2, device=self.device) * (self.grid_size - 10) + 5
        y_coords = torch.rand(batch_size, 2, device=self.device) * (self.grid_size - 10) + 5
        
        velocity_heterogeneity, mean_path_velocities = grid.sample_velocities_batch(x_coords, y_coords)
        
        self.assertEqual(velocity_heterogeneity.shape, (batch_size,))
        self.assertEqual(mean_path_velocities.shape, (batch_size,))
        self.assertTrue(torch.all(velocity_heterogeneity >= 0))
        self.assertTrue(torch.all(mean_path_velocities > 0))


class TestVectorizedZFramework(unittest.TestCase):
    """Test vectorized Z Framework calculations."""
    
    def setUp(self):
        """Set up test environment."""
        self.device = 'cpu'
        self.k_parameter = 0.3
        
    def test_z_framework_initialization(self):
        """Test Z Framework calculator initialization."""
        z_framework = VectorizedZFramework(
            k_parameter=self.k_parameter,
            device=self.device
        )
        
        self.assertEqual(z_framework.k_parameter, self.k_parameter)
        self.assertTrue(z_framework.phi > 1.6)  # Golden ratio check
        self.assertTrue(z_framework.e_squared > 7.0)  # e² check
        
    def test_theta_prime_vectorized(self):
        """Test vectorized geometric resolution calculation."""
        z_framework = VectorizedZFramework(
            k_parameter=self.k_parameter,
            device=self.device
        )
        
        batch_size = 100
        n_values = torch.rand(batch_size, device=self.device) * 50 + 1
        
        theta_prime = z_framework.theta_prime_vectorized(n_values)
        
        self.assertEqual(theta_prime.shape, (batch_size,))
        self.assertTrue(torch.all(theta_prime > 0))
        self.assertTrue(torch.all(theta_prime <= z_framework.phi))
        
    def test_z_enhancement_vectorized(self):
        """Test vectorized Z Framework enhancement calculation."""
        z_framework = VectorizedZFramework(
            k_parameter=self.k_parameter,
            device=self.device
        )
        
        batch_size = 100
        distances = torch.rand(batch_size, device=self.device) * 50 + 1
        velocity_heterogeneity = torch.rand(batch_size, device=self.device) * 0.2
        
        enhancements = z_framework.calculate_z_enhancement_vectorized(distances, velocity_heterogeneity)
        
        self.assertEqual(enhancements.shape, (batch_size,))
        self.assertTrue(torch.all(enhancements >= 0.05))  # Minimum enhancement
        self.assertTrue(torch.all(enhancements <= 0.45))  # Maximum enhancement


class TestVectorizedTargetingModels(unittest.TestCase):
    """Test vectorized targeting models."""
    
    def setUp(self):
        """Set up test environment."""
        self.device = 'cpu'
        self.grid = VectorizedAcousticGrid(
            size=50, base_velocity=1540.0, variance=0.1, seed=42, device=self.device
        )
        self.z_framework = VectorizedZFramework(k_parameter=0.3, device=self.device)
        
    def test_targeting_models_initialization(self):
        """Test targeting models initialization."""
        models = VectorizedTargetingModels(self.grid, self.z_framework)
        
        self.assertEqual(models.grid, self.grid)
        self.assertEqual(models.z_framework, self.z_framework)
        self.assertEqual(models.device, self.grid.device)
        
    def test_vectorized_error_calculation(self):
        """Test vectorized targeting error calculation."""
        models = VectorizedTargetingModels(self.grid, self.z_framework)
        
        batch_size = 50
        source_coords = torch.tensor([[5.0, 5.0]], device=self.device).expand(batch_size, 2)
        target_coords = torch.rand(batch_size, 2, device=self.device) * 40 + 10
        
        baseline_errors, z_framework_errors, baseline_times, z_framework_times = \
            models.calculate_targeting_errors_vectorized(source_coords, target_coords)
        
        # Check shapes
        self.assertEqual(baseline_errors.shape, (batch_size,))
        self.assertEqual(z_framework_errors.shape, (batch_size,))
        self.assertEqual(baseline_times.shape, (batch_size,))
        self.assertEqual(z_framework_times.shape, (batch_size,))
        
        # Check that Z Framework performs better (lower errors)
        self.assertTrue(torch.all(z_framework_errors <= baseline_errors))
        self.assertTrue(torch.all(baseline_errors > 0))
        self.assertTrue(torch.all(baseline_times > 0))


class TestVectorizedStatistics(unittest.TestCase):
    """Test vectorized statistical analysis."""
    
    def setUp(self):
        """Set up test environment."""
        self.device = 'cpu'
        self.statistics = VectorizedStatistics(device=self.device)
        
    def test_bootstrap_correlation(self):
        """Test vectorized bootstrap correlation analysis."""
        # Generate correlated test data
        n_samples = 1000
        x = torch.randn(n_samples, device=self.device)
        y = 0.7 * x + 0.3 * torch.randn(n_samples, device=self.device)  # Known correlation ~0.7
        
        correlation, ci_low, ci_high, p_value = self.statistics.bootstrap_correlation_vectorized(
            x, y, n_boot=100
        )
        
        # Check correlation is reasonable (theoretical ~0.92 for y=0.7x+0.3e with Var=1)
        self.assertTrue(0.85 < correlation < 0.97)
        self.assertTrue(ci_low < correlation < ci_high)
        self.assertTrue(0 <= p_value <= 1)
        
    def test_permutation_test(self):
        """Test vectorized permutation test."""
        # Generate test data with known difference
        n_samples = 500
        x = torch.randn(n_samples, device=self.device) + 0.5  # Mean ~0.5
        y = torch.randn(n_samples, device=self.device)        # Mean ~0.0
        
        p_value = self.statistics.permutation_test_vectorized(x, y, n_perm=100)
        
        self.assertTrue(0 <= p_value <= 1)
        # With this setup, we expect a significant difference (low p-value)
        self.assertTrue(p_value < 0.1)


class TestVectorizedFUSExperiment(unittest.TestCase):
    """Test complete vectorized FUS experiment."""
    
    def setUp(self):
        """Set up test environment."""
        self.config = VectorizedConfig(
            batch_size=100,
            n_trials=500,  # Small for testing
            grid_size=30,
            n_bootstrap=50,  # Reduced for testing
            n_permutation=50,  # Reduced for testing
            seed=42,
            device='cpu'
        )
        
    def test_experiment_initialization(self):
        """Test experiment initialization."""
        experiment = VectorizedFUSExperiment(self.config)
        
        self.assertEqual(experiment.config, self.config)
        self.assertTrue(hasattr(experiment, 'acoustic_grid'))
        self.assertTrue(hasattr(experiment, 'z_framework'))
        self.assertTrue(hasattr(experiment, 'targeting_models'))
        self.assertTrue(hasattr(experiment, 'statistics'))
        
    def test_coordinate_generation(self):
        """Test batch coordinate generation."""
        experiment = VectorizedFUSExperiment(self.config)
        
        batch_size = 100
        source_coords, target_coords = experiment.generate_trial_coordinates_batch(batch_size)
        
        self.assertEqual(source_coords.shape, (batch_size, 2))
        self.assertEqual(target_coords.shape, (batch_size, 2))
        
        # Check source coordinates are fixed
        self.assertTrue(torch.allclose(source_coords, torch.tensor([5.0, 5.0])))
        
        # Check target coordinates are in valid range
        self.assertTrue(torch.all(target_coords >= 10))
        self.assertTrue(torch.all(target_coords <= self.config.grid_size - 10))
        
    def test_full_experiment_run(self):
        """Test complete experiment execution."""
        experiment = VectorizedFUSExperiment(self.config)
        
        start_time = time.time()
        results = experiment.run_vectorized_experiment()
        elapsed_time = time.time() - start_time
        
        # Check result structure
        self.assertIn('experiment_config', results)
        self.assertIn('performance_metrics', results)
        self.assertIn('statistical_analysis', results)
        self.assertIn('raw_results', results)
        
        # Check experiment config
        exp_config = results['experiment_config']
        self.assertEqual(exp_config['n_trials'], self.config.n_trials)
        self.assertEqual(exp_config['batch_size'], self.config.batch_size)
        
        # Check performance metrics
        perf = results['performance_metrics']
        self.assertTrue(perf['elapsed_time_seconds'] > 0)
        self.assertTrue(perf['trials_per_second'] > 0)
        
        # Check statistical analysis
        stats = results['statistical_analysis']
        self.assertIn('improvement_percentage', stats)
        self.assertIn('cohens_d', stats)
        self.assertIn('correlation_r', stats)
        self.assertIn('permutation_p_value', stats)
        self.assertIn('statistically_significant', stats)
        
        # Check that Z Framework shows improvement
        self.assertTrue(stats['improvement_percentage'] > 0)
        
        print(f"✓ Full experiment completed in {elapsed_time:.2f}s")
        print(f"✓ Processing rate: {perf['trials_per_second']:.0f} trials/second")
        print(f"✓ Improvement: {stats['improvement_percentage']:.1f}%")


class TestPerformanceBenchmark(unittest.TestCase):
    """Performance benchmark tests."""
    
    def test_performance_scaling(self):
        """Test performance with different batch sizes."""
        test_configs = [
            (1000, 100),   # 1K trials, 100 batch
            (10000, 1000), # 10K trials, 1K batch
        ]
        
        for n_trials, batch_size in test_configs:
            config = VectorizedConfig(
                batch_size=batch_size,
                n_trials=n_trials,
                grid_size=50,
                n_bootstrap=10,  # Minimal for performance test
                n_permutation=10,
                seed=42,
                device='cpu'
            )
            
            experiment = VectorizedFUSExperiment(config)
            
            start_time = time.time()
            results = experiment.run_vectorized_experiment()
            elapsed_time = time.time() - start_time
            
            trials_per_second = n_trials / elapsed_time
            
            print(f"✓ {n_trials:,} trials: {elapsed_time:.2f}s ({trials_per_second:.0f} trials/sec)")
            
            # Performance requirement: should process at least 1000 trials/second
            self.assertGreater(trials_per_second, 1000, 
                             f"Performance too slow: {trials_per_second:.0f} trials/sec")


def run_smoke_test():
    """Run a quick smoke test for CI."""
    print("Running FUS Enhancer smoke test...")
    
    config = VectorizedConfig(
        batch_size=50,
        n_trials=200,
        grid_size=20,
        n_bootstrap=20,
        n_permutation=20,
        seed=42,
        device='cpu'
    )
    
    experiment = VectorizedFUSExperiment(config)
    results = experiment.run_vectorized_experiment()
    
    # Basic checks
    assert results['experiment_config']['n_trials'] == 200
    assert results['statistical_analysis']['improvement_percentage'] > 0
    assert results['performance_metrics']['trials_per_second'] > 100
    
    print("✓ FUS Enhancer smoke test passed")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke', action='store_true', help='Run smoke test only')
    args = parser.parse_args()
    
    if args.smoke:
        run_smoke_test()
    else:
        # Run full test suite
        unittest.main(verbosity=2)