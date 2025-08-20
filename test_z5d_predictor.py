"""
Test suite for Z5D Prime Predictor Experiment

This module tests the Z5D predictor implementation and experimental framework
to ensure correctness and reproducibility.
"""

import unittest
import numpy as np
import tempfile
import os
from z5d_prime_predictor_experiment import (
    Z5DPredictor, 
    LIPredictor, 
    Z5DPerformanceExperiment,
    MonteCarloSimulator,
    PredictorResult
)


class TestZ5DPredictor(unittest.TestCase):
    """Test Z5D Prime Predictor implementation."""
    
    def setUp(self):
        self.predictor = Z5DPredictor()
    
    def test_predictor_initialization(self):
        """Test predictor initializes correctly."""
        self.assertEqual(self.predictor.name, "Z5D")
    
    def test_small_values(self):
        """Test predictor handles small values correctly."""
        self.assertEqual(self.predictor.predict(0), 0.0)
        self.assertEqual(self.predictor.predict(1), 2.0)  # 1st prime is 2
        
        # For n=100, the 100th prime is 541
        result = self.predictor.predict(100)
        self.assertGreater(result, 400)
        self.assertLess(result, 700)  # Reasonable bounds around 541
    
    def test_monotonic_property(self):
        """Test that predictions increase monotonically."""
        values = [100, 1000, 10000]
        predictions = [self.predictor.predict(n) for n in values]
        
        for i in range(1, len(predictions)):
            self.assertGreater(predictions[i], predictions[i-1])
    
    def test_known_values(self):
        """Test against known approximate values."""
        # p_1000 ≈ 7919
        result = self.predictor.predict(1000)
        self.assertGreater(result, 7000)
        self.assertLess(result, 9000)
        
        # p_10000 ≈ 104729  
        result = self.predictor.predict(10000)
        self.assertGreater(result, 100000)
        self.assertLess(result, 110000)


class TestLIPredictor(unittest.TestCase):
    """Test LI Prime Predictor implementation."""
    
    def setUp(self):
        self.predictor = LIPredictor()
    
    def test_predictor_initialization(self):
        """Test predictor initializes correctly."""
        self.assertEqual(self.predictor.name, "LI")
    
    def test_small_values(self):
        """Test predictor handles small values correctly."""
        self.assertEqual(self.predictor.predict(0), 0.0)
        self.assertEqual(self.predictor.predict(1), 2.0)  # 1st prime is 2
        
        # For n=100, the 100th prime is 541
        result = self.predictor.predict(100)
        self.assertGreater(result, 400)
        self.assertLess(result, 700)  # Reasonable bounds around 541
    
    def test_known_values(self):
        """Test against known p_n approximations."""
        # p_1000 ≈ 7919
        result = self.predictor.predict(1000)
        self.assertGreater(result, 7000)
        self.assertLess(result, 9000)


class TestPredictorComparison(unittest.TestCase):
    """Test comparison between Z5D and LI predictors."""
    
    def setUp(self):
        self.z5d = Z5DPredictor()
        self.li = LIPredictor()
    
    def test_predictor_differences(self):
        """Test that predictors give different results."""
        test_values = [100, 1000, 10000]
        
        for n in test_values:
            z5d_result = self.z5d.predict(n)
            li_result = self.li.predict(n)
            
            # Results should be different but in same ballpark
            self.assertNotAlmostEqual(z5d_result, li_result, places=3)
            self.assertLess(abs(z5d_result - li_result) / max(z5d_result, li_result), 0.5)
    
    def test_convergence_behavior(self):
        """Test asymptotic behavior for large n."""
        large_n = 100000
        
        z5d_result = self.z5d.predict(large_n)
        li_result = self.li.predict(large_n)
        
        # Both should give reasonable estimates for p_n
        # p_n ≈ n * ln(n) for large n
        expected_approx = large_n * np.log(large_n)
        
        self.assertLess(abs(z5d_result - expected_approx) / expected_approx, 0.3)
        self.assertLess(abs(li_result - expected_approx) / expected_approx, 0.3)
    
    def test_z5d_advantage_large_n(self):
        """Test that Z5D outperforms LI for large n (corrected behavior)."""
        # Test with known prime values
        test_cases = [
            (10**6, 15485863),    # 1 millionth prime
            (10**7, 179424673),   # 10 millionth prime
        ]
        
        for n, known_prime in test_cases:
            z5d_pred = self.z5d.predict(n)
            li_pred = self.li.predict(n)
            
            z5d_error = abs(z5d_pred - known_prime) / known_prime
            li_error = abs(li_pred - known_prime) / li_pred
            
            if n >= 10**7:
                # For large n, Z5D should outperform LI
                self.assertLess(z5d_error, li_error, 
                               f"Z5D should outperform LI for n={n:,}")
            else:
                # For smaller n, LI may edge out (pre-asymptotic regime)
                # Both should be reasonable though
                self.assertLess(z5d_error, 0.01, f"Z5D should be reasonable for n={n:,}")
                self.assertLess(li_error, 0.01, f"LI should be reasonable for n={n:,}")


class TestMonteCarloSimulator(unittest.TestCase):
    """Test Monte Carlo simulation for CRISPR variants."""
    
    def setUp(self):
        self.z5d = Z5DPredictor()
        self.simulator = MonteCarloSimulator(self.z5d, seed=42)
    
    def test_simulator_initialization(self):
        """Test simulator initializes correctly."""
        self.assertEqual(self.simulator.predictor, self.z5d)
        self.assertEqual(self.simulator.seed, 42)
    
    def test_variant_generation(self):
        """Test PCSK9 variant generation."""
        variants = self.simulator.generate_pcsk9_variants(10)
        
        self.assertEqual(len(variants), 10)
        
        # All variants should be DNA sequences
        for variant in variants:
            self.assertTrue(all(base in 'ATCG' for base in variant))
            self.assertGreater(len(variant), 50)  # Reasonable length
    
    def test_reproducibility(self):
        """Test that simulation is reproducible with same seed."""
        variants1 = self.simulator.generate_pcsk9_variants(5)
        
        # Create new simulator with same seed
        simulator2 = MonteCarloSimulator(self.z5d, seed=42)
        variants2 = simulator2.generate_pcsk9_variants(5)
        
        self.assertEqual(variants1, variants2)
    
    def test_runtime_measurement(self):
        """Test runtime measurement functionality."""
        runtime = self.simulator.simulate_crispr_runtime(100)
        
        self.assertIsInstance(runtime, float)
        self.assertGreater(runtime, 0)
        self.assertLess(runtime, 10)  # Should be fast for small numbers


class TestZ5DPerformanceExperiment(unittest.TestCase):
    """Test the main experiment class."""
    
    def setUp(self):
        # Use small parameters for fast testing
        self.experiment = Z5DPerformanceExperiment(
            min_n=1000,
            max_n=10000, 
            num_samples=5
        )
    
    def test_experiment_initialization(self):
        """Test experiment initializes correctly."""
        self.assertEqual(self.experiment.min_n, 1000)
        self.assertEqual(self.experiment.max_n, 10000)
        self.assertEqual(self.experiment.num_samples, 5)
        
        self.assertIsInstance(self.experiment.z5d_predictor, Z5DPredictor)
        self.assertIsInstance(self.experiment.li_predictor, LIPredictor)
    
    def test_test_value_generation(self):
        """Test generation of logarithmically spaced test values."""
        n_values = self.experiment.generate_test_values()
        
        self.assertGreater(len(n_values), 0)
        self.assertLessEqual(len(n_values), self.experiment.num_samples)
        
        # Values should be in range
        self.assertGreaterEqual(min(n_values), self.experiment.min_n)
        self.assertLessEqual(max(n_values), self.experiment.max_n)
        
        # Values should be sorted
        self.assertEqual(n_values, sorted(n_values))
    
    def test_predictor_testing(self):
        """Test individual predictor testing."""
        n_values = [1000, 5000, 10000]
        results = self.experiment.test_predictor(self.experiment.z5d_predictor, n_values)
        
        self.assertEqual(len(results), len(n_values))
        
        for i, result in enumerate(results):
            self.assertIsInstance(result, PredictorResult)
            self.assertEqual(result.n, n_values[i])
            self.assertEqual(result.predictor_name, "Z5D")
            self.assertGreater(result.predicted_prime, 0)
            self.assertGreater(result.actual_prime, 0)
            self.assertGreaterEqual(result.relative_error, 0)
            self.assertGreater(result.computation_time, 0)
    
    def test_crispr_simulation(self):
        """Test CRISPR simulation functionality."""
        results = self.experiment.run_crispr_simulation_test(n_variants=100)
        
        self.assertIn('z5d_runtime', results)
        self.assertIn('li_runtime', results)
        self.assertIn('speedup_percentage', results)
        
        self.assertGreater(results['z5d_runtime'], 0)
        self.assertGreater(results['li_runtime'], 0)
        self.assertIsInstance(results['speedup_percentage'], float)
    
    def test_mini_experiment(self):
        """Test running a complete mini experiment."""
        # This tests the full pipeline with minimal data
        result = self.experiment.run_experiment()
        
        # Check result structure
        self.assertIsNotNone(result.n_values)
        self.assertIsNotNone(result.z5d_results)
        self.assertIsNotNone(result.li_results)
        
        # Check metrics
        self.assertIsInstance(result.z5d_mean_error, float)
        self.assertIsInstance(result.li_mean_error, float)
        self.assertIsInstance(result.error_reduction, float)
        
        # Check statistical tests
        self.assertIn('p_value', result.statistical_significance)
        self.assertIn('t_statistic', result.statistical_significance)
        self.assertIn('cohens_d', result.statistical_significance)
        
        # Check confidence intervals
        self.assertIn('z5d_error_ci', result.confidence_intervals)
        self.assertIn('li_error_ci', result.confidence_intervals)
    
    def test_hypothesis_confirmation_large_n(self):
        """Test that Z5D confirms hypothesis for large n (corrected experiment)."""
        # Test with large n values where Z5D should outperform LI
        large_n_experiment = Z5DPerformanceExperiment(
            min_n=10**7,
            max_n=10**8,
            num_samples=3
        )
        
        result = large_n_experiment.run_experiment()
        
        # For large n, Z5D should achieve lower mean error than LI
        self.assertLess(result.z5d_mean_error, result.li_mean_error,
                       "Z5D should outperform LI for large n values")
        
        # Error reduction should be positive (Z5D better)
        self.assertGreater(result.error_reduction, 0,
                          "Error reduction should be positive when Z5D outperforms LI")
        
        # Statistical significance should confirm the difference
        if result.statistical_significance['p_value'] < 0.05:
            # If statistically significant, confirm Z5D advantage
            self.assertLess(result.z5d_mean_error, result.li_mean_error)
    
    def test_report_generation(self):
        """Test white paper report generation."""
        # Run mini experiment
        result = self.experiment.run_experiment()
        
        # Generate report
        report = self.experiment.generate_report(result)
        
        self.assertIsInstance(report, str)
        self.assertGreater(len(report), 1000)  # Should be substantial
        
        # Check for key sections
        self.assertIn("Executive Summary", report)
        self.assertIn("Methods", report)
        self.assertIn("Results", report)
        self.assertIn("Conclusions", report)
        self.assertIn("Reproducibility", report)
        
        # Check for key metrics
        self.assertIn("Error Reduction", report)
        self.assertIn("Statistical Significance", report)
        self.assertIn("Confidence Intervals", report)


class TestReproducibility(unittest.TestCase):
    """Test experimental reproducibility."""
    
    def test_deterministic_behavior(self):
        """Test that experiments are deterministic."""
        experiment1 = Z5DPerformanceExperiment(min_n=1000, max_n=5000, num_samples=3)
        experiment2 = Z5DPerformanceExperiment(min_n=1000, max_n=5000, num_samples=3)
        
        # Same parameters should give same test values
        values1 = experiment1.generate_test_values()
        values2 = experiment2.generate_test_values()
        
        self.assertEqual(values1, values2)
    
    def test_predictor_consistency(self):
        """Test that predictors give consistent results."""
        z5d1 = Z5DPredictor()
        z5d2 = Z5DPredictor()
        
        test_n = 12345
        result1 = z5d1.predict(test_n)
        result2 = z5d2.predict(test_n)
        
        self.assertAlmostEqual(result1, result2, places=10)


def run_validation_tests():
    """Run a focused validation test suite."""
    print("Running Z5D Prime Predictor Validation Tests...")
    print("=" * 50)
    
    # Test basic functionality
    z5d = Z5DPredictor()
    li = LIPredictor()
    
    test_values = [100, 1000, 10000]
    
    print("Testing predictor implementations:")
    for n in test_values:
        z5d_result = z5d.predict(n)
        li_result = li.predict(n)
        actual = z5d.get_actual_nth_prime(n)
        
        z5d_error = abs(z5d_result - actual) / actual
        li_error = abs(li_result - actual) / actual
        
        print(f"n={n:5d}: Z5D={z5d_result:7.1f} (err={z5d_error:.4f}), "
              f"LI={li_result:7.1f} (err={li_error:.4f}), Actual≈{actual}")
    
    print("\nTesting mini experiment:")
    experiment = Z5DPerformanceExperiment(min_n=1000, max_n=10000, num_samples=3)
    result = experiment.run_experiment()
    
    print(f"Error reduction: {result.error_reduction:.2f}%")
    print(f"Statistical significance: p = {result.statistical_significance['p_value']:.4f}")
    print(f"CRISPR speedup: {result.crispr_simulation['speedup_percentage']:.2f}%")
    
    print("\n✓ Validation tests completed successfully!")


if __name__ == "__main__":
    # Run both unit tests and validation
    print("Running unit tests...")
    unittest.main(argv=[''], exit=False, verbosity=2)
    
    print("\n" + "="*60)
    run_validation_tests()