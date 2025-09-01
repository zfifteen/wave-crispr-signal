#!/usr/bin/env python3
"""
Integrated Wave-CRISPR Pipeline

This module integrates physical Z-metrics with the existing θ′ geodesic spectral 
pipeline at k* = 0.3, providing comprehensive CRISPR guide analysis.
"""

import numpy as np
import sys
import os
from typing import Dict, List, Optional, Tuple

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from invariant_features import InvariantFeatureSet, PhaseAwareSpectralAnalyzer
from crispr_physical_z_metrics import PhysicalZMetricsCalculator
from crispr_guide_designer import CRISPRGuideDesigner

# Geodesic resolution constant
K_STAR = 0.3  # k* = 0.3 as specified in requirements


class IntegratedWaveCRISPRPipeline:
    """
    Integrated pipeline combining physical Z-metrics with θ′ geodesic spectral analysis.
    
    Features:
    - Four physical Z-metrics (opening, stacking, twist, melting kinetics)
    - θ′ geodesic resolution at k* = 0.3
    - Phase-aware spectral analysis
    - Comprehensive guide scoring
    """
    
    def __init__(self, pam_pattern: str = "NGG", guide_length: int = 20):
        """
        Initialize integrated pipeline.
        
        Args:
            pam_pattern: PAM sequence pattern
            guide_length: Guide RNA length
        """
        self.pam_pattern = pam_pattern
        self.guide_length = guide_length
        
        # Initialize component modules
        self.z_metrics_calc = PhysicalZMetricsCalculator()
        self.guide_designer = CRISPRGuideDesigner(pam_pattern, guide_length)
        self.invariant_features = InvariantFeatureSet(pam_pattern)
        self.spectral_analyzer = PhaseAwareSpectralAnalyzer()
        
    def calculate_geodesic_resolution(self, sequence: str, position: int) -> Dict[str, float]:
        """
        Calculate θ′ geodesic resolution at k* = 0.3.
        
        θ′(n,k) = φ·((n mod φ)/φ)^k with k* = 0.3
        
        Args:
            sequence: DNA sequence
            position: Position in sequence
            
        Returns:
            Dictionary with geodesic metrics
        """
        # Golden ratio φ
        phi = 1.6180339887498949
        
        # Calculate θ′ geodesic resolution
        n_mod_phi = position % phi
        ratio = n_mod_phi / phi
        theta_prime = phi * (ratio ** K_STAR)
        
        # Calculate positional weights using geodesic resolution
        sequence_length = len(sequence)
        geodesic_weights = []
        
        for i in range(sequence_length):
            n_mod_phi_i = i % phi
            ratio_i = n_mod_phi_i / phi
            weight_i = phi * (ratio_i ** K_STAR)
            geodesic_weights.append(weight_i)
        
        # Normalize weights
        total_weight = sum(geodesic_weights)
        if total_weight > 0:
            geodesic_weights = [w / total_weight for w in geodesic_weights]
        
        # Calculate weighted sequence properties
        base_values = {'A': 1, 'T': 2, 'C': 3, 'G': 4, 'N': 2.5}
        weighted_composition = sum(
            base_values.get(base, 2.5) * weight 
            for base, weight in zip(sequence, geodesic_weights[:len(sequence)])
        )
        
        # Calculate curvature using geodesic weights
        curvature = 0
        for i in range(1, len(sequence)):
            if i < len(geodesic_weights):
                curr_val = base_values.get(sequence[i], 2.5)
                prev_val = base_values.get(sequence[i-1], 2.5)
                transition = abs(curr_val - prev_val)
                curvature += transition * geodesic_weights[i]
        
        return {
            'theta_prime': theta_prime,
            'k_star': K_STAR,
            'geodesic_weights': geodesic_weights,
            'weighted_composition': weighted_composition,
            'geometric_curvature': curvature,
            'position_weight': geodesic_weights[position] if position < len(geodesic_weights) else 0
        }
    
    def integrate_z_metrics_with_geodesic(self, sequence: str, position: int = 0) -> Dict:
        """
        Integrate physical Z-metrics with θ′ geodesic features.
        
        Args:
            sequence: DNA sequence
            position: Focal position for geodesic analysis
            
        Returns:
            Integrated feature dictionary
        """
        # Calculate physical Z-metrics
        z_metrics = self.z_metrics_calc.calculate_all_physical_z_metrics(
            sequence, validate=False
        )
        
        # Calculate geodesic resolution features
        geodesic_features = self.calculate_geodesic_resolution(sequence, position)
        
        # Calculate phase-aware spectral features
        phase_features = self.spectral_analyzer.calculate_phase_difference_features(sequence)
        
        # Calculate invariant features
        invariant_features = self.invariant_features.calculate_complete_feature_set(sequence)
        
        # Integrate all features into bands
        feature_bands = {
            'physical_z_band': {
                'z_opening': float(z_metrics['opening_kinetics']['z_opening']),
                'z_stacking': float(z_metrics['stacking_dissociation']['z_stacking']),
                'z_twist': float(z_metrics['twist_fluctuation']['z_twist']),
                'z_melting': float(z_metrics['melting_kinetics']['z_melting']),
                'z_mean': float(z_metrics['summary']['z_mean']),
                'z_variance': float(z_metrics['summary']['z_variance'])
            },
            'geodesic_band': {
                'theta_prime': geodesic_features['theta_prime'],
                'k_star': geodesic_features['k_star'],
                'weighted_composition': geodesic_features['weighted_composition'],
                'geometric_curvature': geodesic_features['geometric_curvature'],
                'position_weight': geodesic_features['position_weight']
            },
            'spectral_band': {
                'delta_phase_entropy': phase_features.get('delta_phase_entropy', 0),
                'delta_phase_flatness': phase_features.get('delta_phase_flatness', 0),
                'delta_phase_f1_magnitude': phase_features.get('delta_phase_f1_magnitude', 0)
            },
            'invariant_band': {
                'golden_proximity': invariant_features.get('golden_proximity', {}).get('proximity_score', 0),
                'phase_coherence': invariant_features.get('phase_coherence', 0),
                'curvature_disruption': invariant_features.get('curvature_disruption', {}).get('overall_disruption', 0)
            }
        }
        
        # Calculate composite scores for each band
        band_scores = self._calculate_band_scores(feature_bands)
        
        # Calculate overall integrated score
        integrated_score = self._calculate_integrated_score(band_scores)
        
        return {
            'feature_bands': feature_bands,
            'band_scores': band_scores,
            'integrated_score': integrated_score,
            'raw_metrics': {
                'z_metrics': z_metrics,
                'geodesic_features': geodesic_features,
                'phase_features': phase_features,
                'invariant_features': invariant_features
            }
        }
    
    def _calculate_band_scores(self, feature_bands: Dict) -> Dict[str, float]:
        """Calculate normalized scores for each feature band."""
        # Physical Z-band score (higher Z-mean indicates more activity)
        z_band = feature_bands['physical_z_band']
        z_score = min(1.0, z_band['z_mean'] * 100)  # Scale factor for normalization
        
        # Geodesic band score (higher theta_prime and curvature indicates better structure)
        geo_band = feature_bands['geodesic_band']
        geo_score = min(1.0, (geo_band['theta_prime'] * geo_band['geometric_curvature']) / 10)
        
        # Spectral band score (balance of phase differences)
        spec_band = feature_bands['spectral_band']
        spec_score = 1.0 / (1.0 + abs(spec_band['delta_phase_entropy']) + 
                           abs(spec_band['delta_phase_flatness']))
        
        # Invariant band score (combination of proximity, coherence, and disruption)
        inv_band = feature_bands['invariant_band']
        inv_score = (inv_band['golden_proximity'] * inv_band['phase_coherence']) / \
                   (1.0 + inv_band['curvature_disruption'])
        
        return {
            'z_band_score': z_score,
            'geodesic_band_score': geo_score,
            'spectral_band_score': spec_score,
            'invariant_band_score': inv_score
        }
    
    def _calculate_integrated_score(self, band_scores: Dict[str, float]) -> Dict[str, float]:
        """Calculate weighted integrated score from all bands."""
        # Weights for different bands (can be adjusted based on empirical validation)
        weights = {
            'z_band_score': 0.35,      # Physical Z-metrics most important
            'geodesic_band_score': 0.30,  # θ′ geodesic features
            'spectral_band_score': 0.20,   # Phase spectral features
            'invariant_band_score': 0.15   # Invariant features
        }
        
        # Calculate weighted sum
        weighted_score = sum(
            band_scores[band] * weights[band] 
            for band in weights.keys()
        )
        
        # Calculate confidence based on variance in band scores
        score_values = list(band_scores.values())
        confidence = 1.0 - (np.var(score_values) if len(score_values) > 1 else 0)
        
        return {
            'composite_score': weighted_score,
            'confidence': confidence,
            'weights': weights,
            'normalized_score': min(1.0, max(0.0, weighted_score))  # Clamp to [0,1]
        }
    
    def analyze_guide_with_integrated_pipeline(self, guide_sequence: str, 
                                             target_context: str = None) -> Dict:
        """
        Comprehensive guide analysis using integrated pipeline.
        
        Args:
            guide_sequence: Guide RNA sequence
            target_context: Optional target sequence context
            
        Returns:
            Complete analysis results
        """
        if target_context is None:
            target_context = guide_sequence
        
        # Standard guide design features
        standard_features = {
            'sequence': guide_sequence,
            'length': len(guide_sequence),
            'gc_content': (guide_sequence.count('G') + guide_sequence.count('C')) / len(guide_sequence),
            'on_target_score': self.guide_designer.calculate_on_target_score(guide_sequence)
        }
        
        # Integrated pipeline analysis
        integrated_analysis = self.integrate_z_metrics_with_geodesic(guide_sequence)
        
        # Add repair outcome prediction if target context provided
        repair_outcomes = None
        if target_context and len(target_context) >= len(guide_sequence):
            repair_outcomes = self.guide_designer.predict_repair_outcomes(
                guide_sequence, target_context
            )
        
        return {
            'guide_info': standard_features,
            'integrated_analysis': integrated_analysis,
            'repair_outcomes': repair_outcomes,
            'pipeline_version': '1.0',
            'k_star': K_STAR
        }
    
    def compare_guides_with_integrated_scores(self, guide_sequences: List[str]) -> List[Dict]:
        """
        Compare multiple guides using integrated scoring.
        
        Args:
            guide_sequences: List of guide sequences
            
        Returns:
            List of guides sorted by integrated score
        """
        results = []
        
        for guide_seq in guide_sequences:
            analysis = self.analyze_guide_with_integrated_pipeline(guide_seq)
            
            # Extract key scores for comparison
            comparison_data = {
                'sequence': guide_seq,
                'integrated_score': analysis['integrated_analysis']['integrated_score']['composite_score'],
                'confidence': analysis['integrated_analysis']['integrated_score']['confidence'],
                'z_mean': analysis['integrated_analysis']['feature_bands']['physical_z_band']['z_mean'],
                'theta_prime': analysis['integrated_analysis']['feature_bands']['geodesic_band']['theta_prime'],
                'on_target_score': analysis['guide_info']['on_target_score'],
                'gc_content': analysis['guide_info']['gc_content'],
                'full_analysis': analysis
            }
            
            results.append(comparison_data)
        
        # Sort by integrated score (descending)
        results.sort(key=lambda x: x['integrated_score'], reverse=True)
        
        return results


def main():
    """Example usage of integrated pipeline."""
    pipeline = IntegratedWaveCRISPRPipeline()
    
    # Test sequences
    test_guides = [
        "GGGCCCGGGCCCGGGCCCGG",  # High GC
        "ATATATATATATATATAT",    # Low GC
        "GCGCGCGCGCGCGCGCGC",    # Alternating GC
        "ATCGATCGATCGATCGATCG"   # Balanced
    ]
    
    print("🧬 Integrated Wave-CRISPR Pipeline Analysis")
    print("=" * 60)
    print(f"Pipeline Features:")
    print(f"  - Physical Z-metrics (4 kinetic parameters)")
    print(f"  - θ′ geodesic resolution at k* = {K_STAR}")
    print(f"  - Phase-aware spectral analysis")
    print(f"  - Invariant feature integration")
    print()
    
    # Analyze each guide
    for i, guide in enumerate(test_guides):
        print(f"Guide {i+1}: {guide}")
        analysis = pipeline.analyze_guide_with_integrated_pipeline(guide)
        
        # Extract key metrics
        integrated = analysis['integrated_analysis']
        bands = integrated['feature_bands']
        score = integrated['integrated_score']
        
        print(f"  Physical Z-metrics:")
        z_band = bands['physical_z_band']
        print(f"    Opening-Z: {z_band['z_opening']:.6f}")
        print(f"    Stacking-Z: {z_band['z_stacking']:.6f}")
        print(f"    Twist-Z: {z_band['z_twist']:.6f}")
        print(f"    Melting-Z: {z_band['z_melting']:.6f}")
        
        print(f"  Geodesic features:")
        geo_band = bands['geodesic_band']
        print(f"    θ′ at k*={K_STAR}: {geo_band['theta_prime']:.6f}")
        print(f"    Geometric curvature: {geo_band['geometric_curvature']:.6f}")
        
        print(f"  Integrated Score: {score['composite_score']:.3f} (confidence: {score['confidence']:.3f})")
        print()
    
    # Compare all guides
    print("📊 Guide Comparison (sorted by integrated score):")
    comparison = pipeline.compare_guides_with_integrated_scores(test_guides)
    
    for i, result in enumerate(comparison):
        print(f"  {i+1}. {result['sequence']}")
        print(f"     Integrated Score: {result['integrated_score']:.3f}")
        print(f"     Confidence: {result['confidence']:.3f}")
        print(f"     Z-mean: {result['z_mean']:.6f}")
        print(f"     θ′: {result['theta_prime']:.6f}")
        print()


if __name__ == "__main__":
    main()