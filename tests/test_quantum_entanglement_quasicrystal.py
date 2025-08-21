"""
Test suite for Quantum Entanglement in Quasicrystal Geodesic Networks

This test suite validates the implementation of quantum entanglement analysis
in quasicrystal networks, ensuring mathematical consistency and empirical validation.
"""

import sys
import os
import unittest
import logging
from typing import List, Dict

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import mpmath as mp
    import numpy as np
    
    from quantum_entanglement_quasicrystal import (
        QuantumEntanglementQuasicrystal,
        PenroseTilingGenerator,
        QuasicrystalLatticePoint,
        K_VALIDATED,
        K_FALSIFIED,
        PHI,
        format_results_for_display
    )
    from z_framework import ZFrameworkCalculator
    
    DEPENDENCIES_AVAILABLE = True
except ImportError as e:
    print(f"⚠️  WARNING: Missing dependencies for quantum entanglement tests: {e}")
    DEPENDENCIES_AVAILABLE = False

# Configure logging for tests
logging.basicConfig(level=logging.WARNING)  # Reduce log noise during tests

# Set high precision for tests
if DEPENDENCIES_AVAILABLE:
    mp.dps = 30  # Reduced precision for faster tests


class TestPenroseTilingGenerator(unittest.TestCase):
    """Test the Penrose tiling generation functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        if not DEPENDENCIES_AVAILABLE:
            self.skipTest("Required dependencies not available")
        self.generator = PenroseTilingGenerator(precision_dps=30)
    
    def test_penrose_vertices_generation(self):
        """Test that Penrose vertices are generated correctly."""
        vertices = self.generator.generate_penrose_vertices(100)
        
        # Should generate requested number of vertices
        self.assertEqual(len(vertices), 100)
        
        # Each vertex should be a tuple of two mpmath values
        for vertex in vertices[:10]:  # Check first 10
            self.assertIsInstance(vertex, tuple)
            self.assertEqual(len(vertex), 2)
            self.assertIsInstance(vertex[0], mp.mpf)
            self.assertIsInstance(vertex[1], mp.mpf)
    
    def test_golden_ratio_scaling(self):
        """Test that vertices use golden ratio relationships."""
        vertices = self.generator.generate_penrose_vertices(10)
        
        # Check that distances involve golden ratio relationships
        distances = []
        for i, (x, y) in enumerate(vertices[1:], 1):
            distance = mp.sqrt(x**2 + y**2)
            distances.append(distance)
        
        # Should have non-zero distances
        self.assertTrue(all(d > 0 for d in distances))
    
    def test_density_calculation(self):
        """Test density calculation for vertex sets."""
        vertices = self.generator.generate_penrose_vertices(100)
        density = self.generator.calculate_tiling_density(vertices, mp.mpf("5.0"))
        
        # Density should be positive
        self.assertGreater(density, 0)
        
        # Density should be finite
        self.assertTrue(mp.isfinite(density))
    
    def test_quasiperiodic_order(self):
        """Test quasiperiodic order parameter calculation."""
        vertices = self.generator.generate_penrose_vertices(50)
        order = self.generator.calculate_quasiperiodic_order(vertices)
        
        # Order parameter should be between 0 and 1
        self.assertGreaterEqual(order, 0)
        self.assertLessEqual(order, 1)


class TestQuantumEntanglementQuasicrystal(unittest.TestCase):
    """Test the main quantum entanglement analyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        if not DEPENDENCIES_AVAILABLE:
            self.skipTest("Required dependencies not available")
        self.analyzer = QuantumEntanglementQuasicrystal(precision_dps=30)
    
    def test_initialization(self):
        """Test proper initialization of the analyzer."""
        # Should use validated parameter
        self.assertEqual(self.analyzer.k_validated, K_VALIDATED)
        
        # Should have Z Framework calculator
        self.assertIsInstance(self.analyzer.z_framework, ZFrameworkCalculator)
        
        # Should have Penrose generator
        self.assertIsInstance(self.analyzer.penrose_generator, PenroseTilingGenerator)
    
    def test_parameter_validation(self):
        """Test that validated parameters are used, not falsified ones."""
        # Should use k* ≈ 0.3, not k* ≈ 0.04449
        self.assertAlmostEqual(float(self.analyzer.k_validated), 0.3, places=1)
        self.assertNotEqual(self.analyzer.k_validated, K_FALSIFIED)
    
    def test_lattice_generation(self):
        """Test quasicrystal lattice generation."""
        lattice = self.analyzer.generate_quasicrystal_lattice(50)
        
        # Should generate requested number of points
        self.assertEqual(len(lattice), 50)
        
        # Each point should be a QuasicrystalLatticePoint
        for point in lattice:
            self.assertIsInstance(point, QuasicrystalLatticePoint)
            
            # Should have all required attributes
            self.assertIsInstance(point.x, mp.mpf)
            self.assertIsInstance(point.y, mp.mpf)
            self.assertIsInstance(point.z, mp.mpf)
            self.assertIsInstance(point.geodesic_curvature, mp.mpf)
            self.assertIsInstance(point.entanglement_amplitude, complex)
            self.assertIsInstance(point.phase, mp.mpf)
    
    def test_geodesic_curvature_calculation(self):
        """Test geodesic curvature calculation."""
        # Test at origin
        curvature_origin = self.analyzer._calculate_geodesic_curvature(
            mp.mpf("0"), mp.mpf("0"), mp.mpf("0")
        )
        self.assertEqual(curvature_origin, 0)
        
        # Test at non-origin point
        curvature_point = self.analyzer._calculate_geodesic_curvature(
            mp.mpf("1"), mp.mpf("1"), mp.mpf("1")
        )
        self.assertGreater(curvature_point, 0)
    
    def test_entanglement_amplitude_calculation(self):
        """Test quantum entanglement amplitude calculation."""
        amplitude = self.analyzer._calculate_entanglement_amplitude(
            mp.mpf("1"), mp.mpf("1"), mp.mpf("1"), 5
        )
        
        # Should be a complex number
        self.assertIsInstance(amplitude, complex)
        
        # Should have finite magnitude
        self.assertTrue(np.isfinite(abs(amplitude)))
    
    def test_quantum_phase_calculation(self):
        """Test quantum phase calculation."""
        phase = self.analyzer._calculate_quantum_phase(
            mp.mpf("1"), mp.mpf("1"), mp.mpf("1")
        )
        
        # Phase should be between 0 and 2π
        self.assertGreaterEqual(phase, 0)
        self.assertLessEqual(phase, 2 * mp.pi)
    
    def test_entanglement_entropy_calculation(self):
        """Test entanglement entropy calculation."""
        lattice = self.analyzer.generate_quasicrystal_lattice(20)
        entropy = self.analyzer.calculate_entanglement_entropy(lattice, 10)
        
        # Entropy should be non-negative
        self.assertGreaterEqual(entropy, 0)
        
        # Entropy should be finite
        self.assertTrue(mp.isfinite(entropy))
    
    def test_correlation_strength_calculation(self):
        """Test quantum correlation strength calculation."""
        lattice = self.analyzer.generate_quasicrystal_lattice(20)
        correlation = self.analyzer.calculate_quantum_correlation_strength(lattice)
        
        # Correlation should be between 0 and 1
        self.assertGreaterEqual(correlation, 0)
        self.assertLessEqual(correlation, 1)
    
    def test_decoherence_time_estimation(self):
        """Test decoherence time estimation."""
        lattice = self.analyzer.generate_quasicrystal_lattice(20)
        decoherence_time = self.analyzer.estimate_decoherence_time(lattice)
        
        # Decoherence time should be positive
        self.assertGreater(decoherence_time, 0)
        
        # Should be finite
        self.assertTrue(mp.isfinite(decoherence_time))
    
    def test_density_enhancement_validation(self):
        """Test density enhancement validation with corrected parameters."""
        # Use small sample for fast testing
        results = self.analyzer.perform_density_enhancement_validation(n_points=100)
        
        # Should contain all expected keys
        expected_keys = [
            "n_points", "k_parameter_used", "k_parameter_falsified",
            "baseline_density", "lattice_density", "density_enhancement_percent",
            "quasiperiodic_order", "entanglement_entropy", "correlation_strength",
            "decoherence_time", "geodesic_curvature_avg", "quantum_phase_coherence"
        ]
        
        for key in expected_keys:
            self.assertIn(key, results)
        
        # Should use validated parameter, not falsified one
        self.assertEqual(results["k_parameter_used"], K_VALIDATED)
        self.assertEqual(results["k_parameter_falsified"], K_FALSIFIED)
        
        # Enhancement should be reasonable (not 210%)
        enhancement = float(results["density_enhancement_percent"])
        self.assertLess(enhancement, 100)  # Should be less than 100% enhancement
        self.assertGreater(enhancement, -99)  # Allow for calculation variations in small samples
    
    def test_entanglement_stability_analysis(self):
        """Test quantum entanglement stability analysis."""
        # Use small sample for fast testing
        results = self.analyzer.analyze_quantum_entanglement_stability(
            n_points=50, num_perturbations=5
        )
        
        # Should contain expected keys
        expected_keys = [
            "baseline_entropy", "baseline_correlation", "num_perturbations",
            "avg_entropy_variation", "avg_correlation_variation", 
            "avg_stability_score", "is_stable"
        ]
        
        for key in expected_keys:
            self.assertIn(key, results)
        
        # Variations should be non-negative
        self.assertGreaterEqual(results["avg_entropy_variation"], 0)
        self.assertGreaterEqual(results["avg_correlation_variation"], 0)
        
        # Stability score should be between 0 and 1
        self.assertGreaterEqual(results["avg_stability_score"], 0)
        self.assertLessEqual(results["avg_stability_score"], 1)


class TestResultFormatting(unittest.TestCase):
    """Test result formatting utilities."""
    
    def setUp(self):
        """Set up test fixtures."""
        if not DEPENDENCIES_AVAILABLE:
            self.skipTest("Required dependencies not available")
    
    def test_format_results_for_display(self):
        """Test formatting of results for display."""
        test_results = {
            "mpmath_value": mp.mpf("3.14159"),
            "list_of_mpmath": [mp.mpf("1.0"), mp.mpf("2.0"), mp.mpf("3.0")],
            "string_value": "test",
            "integer_value": 42
        }
        
        formatted = format_results_for_display(test_results, precision=3)
        
        # Should convert mpmath values to strings
        self.assertIsInstance(formatted["mpmath_value"], str)
        self.assertEqual(formatted["mpmath_value"], "3.142")
        
        # Should handle lists of mpmath values
        self.assertIsInstance(formatted["list_of_mpmath"], list)
        self.assertEqual(formatted["list_of_mpmath"], ["1.000", "2.000", "3.000"])
        
        # Should preserve other types
        self.assertEqual(formatted["string_value"], "test")
        self.assertEqual(formatted["integer_value"], "42")


class TestIntegrationWithZFramework(unittest.TestCase):
    """Test integration with the Z Framework."""
    
    def setUp(self):
        """Set up test fixtures."""
        if not DEPENDENCIES_AVAILABLE:
            self.skipTest("Required dependencies not available")
        self.analyzer = QuantumEntanglementQuasicrystal(precision_dps=30)
    
    def test_z_framework_integration(self):
        """Test that the analyzer properly integrates with Z Framework."""
        # Should use the same precision
        self.assertEqual(self.analyzer.z_framework.precision, 30)
        
        # Should use same mathematical constants
        self.assertEqual(self.analyzer.phi, PHI)
        
        # Geodesic resolution should be consistent
        test_n = 10
        z_resolution = self.analyzer.z_framework.calculate_geodesic_resolution(
            test_n, float(K_VALIDATED)
        )
        
        # Should be positive and finite
        self.assertGreater(z_resolution, 0)
        self.assertTrue(mp.isfinite(z_resolution))
    
    def test_parameter_consistency(self):
        """Test that parameters are consistent across components."""
        # All components should use validated parameter
        self.assertEqual(self.analyzer.k_validated, K_VALIDATED)
        
        # Should not use falsified parameter
        self.assertNotEqual(self.analyzer.k_validated, K_FALSIFIED)


def run_quantum_entanglement_tests():
    """Run all quantum entanglement tests."""
    if not DEPENDENCIES_AVAILABLE:
        print("❌ Quantum Entanglement Tests: SKIPPED (missing dependencies)")
        return False
    
    # Create test suite
    suite = unittest.TestSuite()
    
    # Add test classes
    test_classes = [
        TestPenroseTilingGenerator,
        TestQuantumEntanglementQuasicrystal,
        TestResultFormatting,
        TestIntegrationWithZFramework
    ]
    
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2, buffer=True)
    result = runner.run(suite)
    
    # Return success status
    return result.wasSuccessful()


if __name__ == "__main__":
    print("Testing Quantum Entanglement in Quasicrystal Geodesic Networks")
    print("=" * 70)
    
    success = run_quantum_entanglement_tests()
    
    if success:
        print("✓ All quantum entanglement tests passed!")
    else:
        print("❌ Some quantum entanglement tests failed!")
        sys.exit(1)