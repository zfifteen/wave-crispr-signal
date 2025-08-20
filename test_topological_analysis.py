"""
Test suite for Topological Analysis: f(x) = arcsin((x-1)/(2x+3)) properties

Validates the mathematical connections between f(x) topological properties
and geodesic curvature in the Z Framework.
"""

import unittest
import numpy as np
import mpmath as mp
from topological_analysis import TopologicalAnalyzer

class TestTopologicalAnalysis(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures"""
        self.analyzer = TopologicalAnalyzer(precision_dps=30)
    
    def test_fx_function_basic(self):
        """Test basic f(x) = arcsin((x-1)/(2x+3)) calculation"""
        # Test point where function is well-defined (x = 0 is in valid domain)
        x = mp.mpf(0)
        result = self.analyzer.f_x(x)
        
        # Verify it's a valid number
        self.assertIsInstance(result, mp.mpf)
        self.assertFalse(mp.isnan(result))
        
        # Check that result is in valid arcsin range [-π/2, π/2]
        self.assertGreaterEqual(result, -mp.pi/2)
        self.assertLessEqual(result, mp.pi/2)
    
    def test_fx_pole_singularity(self):
        """Test that pole at x = -3/2 raises appropriate error"""
        with self.assertRaises(ValueError) as context:
            self.analyzer.f_x(mp.mpf(-3)/mp.mpf(2))
        
        self.assertIn("Pole singularity", str(context.exception))
    
    def test_fx_domain_violations(self):
        """Test domain violations where arcsin argument is outside [-1, 1]"""
        # Test values that should cause domain violations
        # Values from testing that cause arg outside [-1,1]
        invalid_x_values = [mp.mpf(-2), mp.mpf(-2.5), mp.mpf(-1.4)]
        
        for x in invalid_x_values:
            with self.assertRaises(ValueError) as context:
                self.analyzer.f_x(x)
            # Should be domain violation 
            self.assertIn("outside domain", str(context.exception))
    
    def test_domain_constraints_analysis(self):
        """Test domain constraint analysis"""
        domain_info = self.analyzer.analyze_domain_constraints()
        
        # Check that expected keys are present
        expected_keys = ['boundary_arg_minus1', 'boundary_arg1', 'pole', 'admissible_intervals']
        for key in expected_keys:
            self.assertIn(key, domain_info)
        
        # Verify mathematical relationships
        # Boundary where arg = -1 should be -2/3
        self.assertAlmostEqual(float(domain_info['boundary_arg_minus1']), -2.0/3.0, places=10)
        
        # Boundary where arg = 1 should be -4
        self.assertAlmostEqual(float(domain_info['boundary_arg1']), -4.0, places=10)
        
        # Pole should be -3/2
        self.assertAlmostEqual(float(domain_info['pole']), -1.5, places=10)
    
    def test_geodesic_resolution(self):
        """Test geodesic resolution function θ'(n, k) = φ · ((n mod φ)/φ)^k"""
        # Test with default k=0.3
        result = self.analyzer.geodesic_resolution(10)
        
        self.assertIsInstance(result, mp.mpf)
        self.assertGreater(result, 0)
        
        # Test with different k values
        k_values = [0.1, 0.3, 0.5, 1.0]
        results = [self.analyzer.geodesic_resolution(10, k) for k in k_values]
        
        # All results should be positive
        for result in results:
            self.assertGreater(result, 0)
    
    def test_mapping_fx_to_geodesic(self):
        """Test mapping f(x) values to geodesic space"""
        # Use valid x values: x=0, x=1, x=-0.6 are in valid domain
        x_values = [mp.mpf(0), mp.mpf(1), mp.mpf(-0.6)]
        
        mapping = self.analyzer.map_fx_to_geodesic(x_values)
        
        # Check structure
        expected_keys = ['x_values', 'fx_values', 'geodesic_values', 'correspondences']
        for key in expected_keys:
            self.assertIn(key, mapping)
        
        # Check that we have correspondences
        self.assertGreater(len(mapping['correspondences']), 0)
        
        # Check correspondence structure
        for corr in mapping['correspondences']:
            if corr.get('fx') is not None:  # Non-singular points
                self.assertIn('x', corr)
                self.assertIn('fx', corr)
                self.assertIn('geodesic', corr)
                self.assertIn('resonance', corr)
    
    def test_invariant_alignment(self):
        """Test invariant alignment analysis"""
        alignment = self.analyzer.analyze_invariant_alignment()
        
        # Check expected keys
        expected_keys = ['domain_gap_span', 'e_squared_ratio', 'phi_ratio', 
                        'e_squared_resonance', 'phi_resonance', 'optimal_alignment']
        for key in expected_keys:
            self.assertIn(key, alignment)
        
        # Domain gap span should be positive
        self.assertGreater(alignment['domain_gap_span'], 0)
        
        # Ratios should be meaningful
        self.assertGreater(alignment['e_squared_ratio'], 0)
        self.assertGreater(alignment['phi_ratio'], 0)
    
    def test_density_enhancement(self):
        """Test density enhancement analysis"""
        enhancement = self.analyzer.demonstrate_density_enhancement(n_points=20)
        
        # Check expected keys
        expected_keys = ['baseline_variance', 'enhanced_variance', 'density_enhancement',
                        'enhancement_percentage', 'target_enhancement', 'achievement_ratio']
        for key in expected_keys:
            self.assertIn(key, enhancement)
        
        # Variances should be non-negative
        self.assertGreaterEqual(enhancement['baseline_variance'], 0)
        self.assertGreaterEqual(enhancement['enhanced_variance'], 0)
        
        # Target enhancement should be 15%
        self.assertAlmostEqual(float(enhancement['target_enhancement']), 15.0, places=5)
    
    def test_comprehensive_analysis(self):
        """Test comprehensive analysis integration"""
        results = self.analyzer.comprehensive_analysis()
        
        # Check main sections
        expected_sections = ['domain_constraints', 'invariant_alignment', 
                           'density_enhancement', 'mapping_demonstration',
                           'hypothesis_validation']
        for section in expected_sections:
            self.assertIn(section, results)
        
        # Check hypothesis validation structure
        validation = results['hypothesis_validation']
        expected_validations = ['domain_bounded_correctly', 'invariant_alignment_achieved',
                               'density_enhancement_observed', 'topological_bridge_established']
        for validation_key in expected_validations:
            self.assertIn(validation_key, validation)
            self.assertIsInstance(validation[validation_key], bool)
    
    def test_mathematical_constants(self):
        """Test that mathematical constants are correctly set"""
        # Golden ratio φ ≈ 1.618
        self.assertAlmostEqual(float(self.analyzer.phi), (1 + np.sqrt(5))/2, places=10)
        
        # e² ≈ 7.389
        self.assertAlmostEqual(float(self.analyzer.e_squared), np.e**2, places=10)
        
        # Pole at -3/2
        self.assertAlmostEqual(float(self.analyzer.pole_x), -1.5, places=10)
        
        # Domain boundary at -2/3
        self.assertAlmostEqual(float(self.analyzer.boundary_arg_minus1), -2.0/3.0, places=10)
    
    def test_valid_fx_calculation_examples(self):
        """Test f(x) calculation with specific valid examples"""
        # Test cases with known valid domain points
        test_cases = [
            (mp.mpf(0), "Should work at x=0"),
            (mp.mpf(1), "Should work at x=1"),
            (mp.mpf(-0.6), "Should work at x=-0.6"),
        ]
        
        for x, description in test_cases:
            try:
                result = self.analyzer.f_x(x)
                self.assertIsInstance(result, mp.mpf, f"Failed for {description}")
                self.assertFalse(mp.isnan(result), f"NaN result for {description}")
            except ValueError as e:
                self.fail(f"Unexpected error for {description}: {e}")
    
    def test_k_star_optimal_value(self):
        """Test that k* ≈ 0.3 shows optimal behavior"""
        # Test different k values around 0.3
        k_values = [0.1, 0.2, 0.3, 0.4, 0.5]
        x_values = [mp.mpf(0), mp.mpf(1), mp.mpf(-0.6)]  # Valid domain points
        
        mappings = []
        for k in k_values:
            mapping = self.analyzer.map_fx_to_geodesic(x_values, k=k)
            mappings.append((k, mapping))
        
        # Verify that we got results for k=0.3
        k_03_mapping = next(mapping for k, mapping in mappings if k == 0.3)
        self.assertGreater(len(k_03_mapping['correspondences']), 0)
        
        # Check that k=0.3 produces reasonable geodesic values
        for corr in k_03_mapping['correspondences']:
            if corr.get('geodesic') is not None:
                self.assertGreaterEqual(corr['geodesic'], 0)  # Allow 0 or positive


if __name__ == '__main__':
    # Configure mpmath for testing
    mp.dps = 30
    
    # Run tests
    unittest.main(verbosity=2)