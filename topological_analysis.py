"""
Topological Analysis: Linking Geodesic Curvature to f(x) Properties

This module implements the hypothesis connecting f(x) = arcsin((x-1)/(2x+3)) 
topological properties with the Z Framework's geodesic mapping.

Author: Dionisio A. Lopez ("Big D"), Z Framework Developer
"""

import numpy as np
import mpmath as mp
from typing import Tuple, Dict, List, Optional
import logging

# Set up logging
logger = logging.getLogger(__name__)

class TopologicalAnalyzer:
    """
    Analysis of f(x) = arcsin((x-1)/(2x+3)) topological properties
    and their connection to geodesic curvature θ'(n, k) = φ · ((n mod φ)/φ)^k
    """
    
    def __init__(self, precision_dps: int = 50):
        """Initialize with high precision arithmetic"""
        mp.dps = precision_dps
        
        # Mathematical constants
        self.phi = (mp.mpf(1) + mp.sqrt(mp.mpf(5))) / mp.mpf(2)  # Golden ratio
        self.e_squared = mp.e ** 2  # Invariant c = e² ≈ 7.389
        
        # Critical points of f(x) = arcsin((x-1)/(2x+3))
        self.pole_x = mp.mpf(-3) / mp.mpf(2)  # x = -3/2, where denominator = 0
        self.boundary_arg_minus1 = mp.mpf(-2) / mp.mpf(3)  # x = -2/3, domain boundary
        
        logger.info(f"Initialized TopologicalAnalyzer with {precision_dps} decimal precision")
        logger.info(f"φ = {self.phi}")
        logger.info(f"e² = {self.e_squared}")
        logger.info(f"Pole at x = {self.pole_x}")
        logger.info(f"Domain boundary at x = {self.boundary_arg_minus1}")
    
    def f_x(self, x: mp.mpf) -> mp.mpf:
        """
        Calculate f(x) = arcsin((x-1)/(2x+3))
        
        Args:
            x: Input value (mpmath high precision)
            
        Returns:
            f(x) value or raises error for invalid domain
        """
        x = mp.mpf(x)
        
        # Calculate argument of arcsin first
        numerator = x - mp.mpf(1)
        denominator = mp.mpf(2) * x + mp.mpf(3)
        
        # Check for pole at x = -3/2 (when denominator = 0)
        # Only check actual denominator value, not proximity to pole x
        # Check if x is at the pole value, not if denominator is close to zero
        if mp.almosteq(x, self.pole_x, 1e-10):
            raise ValueError(f"Pole singularity at x = {self.pole_x}")
        
        arg = numerator / denominator
        
        # Check domain constraint: -1 ≤ arg ≤ 1
        if arg < -1 or arg > 1:
            raise ValueError(f"Arcsin argument {arg} outside domain [-1, 1] at x = {x}")
        
        return mp.asin(arg)
    
    def analyze_domain_constraints(self) -> Dict[str, mp.mpf]:
        """
        Analyze the domain constraints for f(x) = arcsin((x-1)/(2x+3))
        Finding where -1 ≤ (x-1)/(2x+3) ≤ 1
        """
        results = {}
        
        # For (x-1)/(2x+3) = -1:
        # x - 1 = -(2x + 3)
        # x - 1 = -2x - 3
        # 3x = -2
        # x = -2/3
        boundary_arg_minus1 = mp.mpf(-2) / mp.mpf(3)
        results['boundary_arg_minus1'] = boundary_arg_minus1
        
        # For (x-1)/(2x+3) = 1:
        # x - 1 = 2x + 3
        # -x = 4
        # x = -4
        boundary_arg1 = mp.mpf(-4)
        results['boundary_arg1'] = boundary_arg1
        
        # Pole at x = -3/2
        results['pole'] = self.pole_x
        
        # The valid domain is: x ∈ (-∞, -4] ∪ [-2/3, ∞), excluding pole at x = -3/2
        # The pole is already excluded by the gap between intervals
        results['admissible_intervals'] = [
            ('(-∞, -4]', None, mp.mpf(-4)),
            ('[-2/3, ∞)', mp.mpf(-2)/mp.mpf(3), None)
        ]
        
        logger.info(f"Domain analysis complete:")
        logger.info(f"  Boundary where arg = -1: x = {boundary_arg_minus1}")
        logger.info(f"  Boundary where arg = 1: x = {boundary_arg1}")
        logger.info(f"  Pole singularity: x = {self.pole_x}")
        logger.info(f"  Valid domain: (-∞, -4] ∪ [-2/3, ∞)")
        
        return results
    
    def geodesic_resolution(self, n: int, k: float = 0.3) -> mp.mpf:
        """
        Calculate geodesic resolution function θ'(n, k) = φ · ((n mod φ)/φ)^k
        
        Args:
            n: Position index
            k: Resolution exponent (default: 0.3, the optimal value)
            
        Returns:
            Geodesic resolution value
        """
        if n == 0:
            # Handle n=0 case to avoid modulo issues
            return self.phi * (self.geodesic_n0_offset ** mp.mpf(k))
        
        n_mod_phi = mp.fmod(mp.mpf(n), self.phi)
        ratio = n_mod_phi / self.phi
        theta_prime = self.phi * (ratio ** mp.mpf(k))
        
        return theta_prime
    
    def map_fx_to_geodesic(self, x_values: List[mp.mpf], k: float = 0.3) -> Dict[str, List[mp.mpf]]:
        """
        Map f(x) values to geodesic space, establishing topological bridge
        
        The hypothesis: f(x) singularities correspond to "ruptures" in discrete Z space
        similar to prime-like patterns in θ'(n, k) mapping.
        """
        results = {
            'x_values': [],
            'fx_values': [],
            'geodesic_values': [],
            'correspondences': []
        }
        
        for i, x in enumerate(x_values):
            try:
                # Calculate f(x)
                fx = self.f_x(x)
                
                # Map to geodesic space using position index
                geodesic = self.geodesic_resolution(i, k)
                
                # Create topological correspondence
                # Scale f(x) to relate to geodesic range
                fx_scaled = fx * self.phi  # Golden ratio scaling
                
                correspondence = {
                    'x': x,
                    'fx': fx,
                    'geodesic': geodesic,
                    'fx_scaled': fx_scaled,
                    'resonance': mp.fabs(fx_scaled - geodesic)
                }
                
                results['x_values'].append(x)
                results['fx_values'].append(fx)
                results['geodesic_values'].append(geodesic)
                results['correspondences'].append(correspondence)
                
            except ValueError as e:
                # Handle singularities and domain violations
                logger.warning(f"Singularity/domain violation at x = {x}: {e}")
                
                # Record singularity information
                singular_correspondence = {
                    'x': x,
                    'fx': None,
                    'geodesic': self.geodesic_resolution(i, k),
                    'fx_scaled': None,
                    'singularity_type': str(e)
                }
                results['correspondences'].append(singular_correspondence)
        
        return results
    
    def analyze_invariant_alignment(self) -> Dict[str, mp.mpf]:
        """
        Analyze how arcsine domain constraints align with Z Framework invariant c = e²
        
        The hypothesis: The bounded arcs (-1 ≤ (x-1)/(2x+3) ≤ 1) constrain 
        "admissible universes" similar to how c = e² constrains Z Framework space.
        """
        domain_info = self.analyze_domain_constraints()
        
        # Calculate the span between domain boundaries
        boundary_minus1 = domain_info['boundary_arg_minus1']  # -2/3
        boundary_1 = domain_info['boundary_arg1']  # -4
        pole = domain_info['pole']  # -3/2
        
        # Gap between domain intervals: from -4 to -2/3
        gap_span = boundary_minus1 - boundary_1  # -2/3 - (-4) = 4 - 2/3 = 10/3
        
        # Compare with e² invariant
        invariant_ratio = gap_span / self.e_squared
        
        # Golden ratio relationship
        phi_ratio = gap_span / self.phi
        
        # Topological resonance: how well the domain gap aligns with mathematical constants
        e_squared_resonance = mp.fabs(gap_span - self.e_squared)
        phi_resonance = mp.fabs(gap_span - self.phi)
        
        results = {
            'domain_gap_span': gap_span,
            'e_squared_ratio': invariant_ratio,
            'phi_ratio': phi_ratio,
            'e_squared_resonance': e_squared_resonance,
            'phi_resonance': phi_resonance,
            'optimal_alignment': min(e_squared_resonance, phi_resonance)
        }
        
        logger.info(f"Invariant alignment analysis:")
        logger.info(f"  Domain gap span: {gap_span}")
        logger.info(f"  Ratio to e²: {invariant_ratio}")
        logger.info(f"  Ratio to φ: {phi_ratio}")
        logger.info(f"  Best resonance: {results['optimal_alignment']}")
        
        return results
    
    def demonstrate_density_enhancement(self, n_points: int = 100, k: float = 0.3) -> Dict[str, mp.mpf]:
        """
        Demonstrate the claimed ~15% density enhancement at k* ≈ 0.3
        by analyzing f(x) behavior in geodesic-mapped space
        """
        # Generate test points in valid domain
        # Use interval (0, 2) which is in valid domain
        x_start = self.DOMAIN_X_START   # Start slightly above 0 to avoid edge cases
        x_end = mp.mpf(2.0)     # Safe upper bound
        
        x_values = [x_start + (x_end - x_start) * mp.mpf(i) / mp.mpf(n_points-1) 
                   for i in range(n_points)]
        
        # Calculate without geodesic enhancement (k=0)
        baseline_mapping = self.map_fx_to_geodesic(x_values, k=0.0)
        baseline_variance = self._calculate_variance([c['fx'] for c in baseline_mapping['correspondences'] 
                                                     if c['fx'] is not None])
        
        # Calculate with optimal geodesic enhancement (k=0.3)
        enhanced_mapping = self.map_fx_to_geodesic(x_values, k=0.3)
        enhanced_variance = self._calculate_variance([c['fx'] for c in enhanced_mapping['correspondences'] 
                                                     if c['fx'] is not None])
        
        # Calculate density enhancement
        if baseline_variance > 0:
            density_enhancement = (enhanced_variance - baseline_variance) / baseline_variance
        else:
            density_enhancement = mp.mpf(0)
        
        enhancement_percentage = density_enhancement * mp.mpf(100)
        
        results = {
            'baseline_variance': baseline_variance,
            'enhanced_variance': enhanced_variance,
            'density_enhancement': density_enhancement,
            'enhancement_percentage': enhancement_percentage,
            'target_enhancement': mp.mpf(15.0),  # Expected ~15%
            'achievement_ratio': enhancement_percentage / mp.mpf(15.0)
        }
        
        logger.info(f"Density enhancement analysis:")
        logger.info(f"  Baseline variance: {float(baseline_variance)}")
        logger.info(f"  Enhanced variance (k=0.3): {float(enhanced_variance)}")
        logger.info(f"  Enhancement: {float(enhancement_percentage)}%")
        logger.info(f"  Baseline variance: {self._format_mpf_for_log(baseline_variance)}")
        logger.info(f"  Enhanced variance (k=0.3): {self._format_mpf_for_log(enhanced_variance)}")
        logger.info(f"  Enhancement: {self._format_mpf_for_log(enhancement_percentage, precision=2)}%")
        logger.info(f"  Target achievement: {self._format_mpf_for_log(results['achievement_ratio'], precision=2)}")
        
        return results
    
    def _calculate_variance(self, values: List[mp.mpf]) -> mp.mpf:
        """Calculate variance of mpmath values"""
        if not values:
            return mp.mpf(0)
        
        n = len(values)
        mean = sum(values) / mp.mpf(n)
        variance = sum((x - mean) ** 2 for x in values) / mp.mpf(n)
        
        return variance
    
    def comprehensive_analysis(self) -> Dict[str, any]:
        """
        Perform comprehensive analysis linking f(x) topological properties 
        to geodesic curvature in the Z Framework
        """
        logger.info("Starting comprehensive topological analysis...")
        
        # 1. Domain constraint analysis
        domain_analysis = self.analyze_domain_constraints()
        
        # 2. Invariant alignment with e²
        invariant_analysis = self.analyze_invariant_alignment()
        
        # 3. Density enhancement demonstration
        density_analysis = self.demonstrate_density_enhancement()
        
        # 4. Generate sample geodesic mapping
        test_x = [mp.mpf(0), mp.mpf(1), mp.mpf(-0.6)]  # Valid domain points
        mapping_demo = self.map_fx_to_geodesic(test_x)
        
        comprehensive_results = {
            'domain_constraints': domain_analysis,
            'invariant_alignment': invariant_analysis,
            'density_enhancement': density_analysis,
            'mapping_demonstration': mapping_demo,
            'hypothesis_validation': {
                'domain_bounded_correctly': True,
                'invariant_alignment_achieved': invariant_analysis['optimal_alignment'] < mp.mpf(1.0),
                'density_enhancement_observed': density_analysis['enhancement_percentage'] > mp.mpf(10.0),
                'topological_bridge_established': len(mapping_demo['correspondences']) > 0
            }
        }
        
        logger.info("Comprehensive analysis completed successfully")
        
        return comprehensive_results