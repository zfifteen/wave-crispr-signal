"""
Z Framework Implementation for DNA Sequence Analysis

This module implements the Z Framework discrete domain form (Z = n(Δₙ/Δₘₐₓ))
with high-precision calculations for empirical validation of biological invariants.

Key Features:
- High-precision calculations using mpmath (dps=50)
- DNA sequence mapping (A=1, T=2, C=3, G=4)
- Convergence validation to golden ratio conjugate (φ-1 ≈ 0.607)
- Variance analysis with target σ ≈ 0.113
- Geodesic resolution function θ'(n, k) = φ·((n mod φ)/φ)^k
"""

import mpmath as mp
import numpy as np
from typing import List, Dict, Tuple, Optional
from collections import Counter
import logging

# Configure high precision
mp.dps = 50

# Mathematical constants with high precision
PHI = mp.mpf('1.618033988749894848204586834365638117720309179805762862135')  # Golden ratio
PHI_CONJUGATE = PHI - 1  # φ-1 ≈ 0.618...
E_SQUARED = mp.e ** 2  # e² ≈ 7.389

# DNA base mapping for Z Framework
DNA_MAPPING = {'A': mp.mpf('1'), 'T': mp.mpf('2'), 'C': mp.mpf('3'), 'G': mp.mpf('4')}

# Target variance for biological invariants
TARGET_VARIANCE = mp.mpf('0.113')

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ZFrameworkCalculator:
    """
    High-precision calculator for Z Framework analysis of DNA sequences.
    
    Implements the discrete domain form Z = n(Δₙ/Δₘₐₓ) with empirical
    validation of biological invariants and convergence properties.
    """
    
    def __init__(self, precision_dps: int = 50):
        """
        Initialize the Z Framework calculator.
        
        Args:
            precision_dps: Decimal precision for mpmath calculations
        """
        mp.dps = precision_dps
        self.precision = precision_dps
        self.phi = PHI
        self.phi_conjugate = PHI_CONJUGATE
        self.e_squared = E_SQUARED
        
        logger.info(f"Initialized Z Framework Calculator with {precision_dps} decimal precision")
        logger.info(f"φ = {self.phi}")
        logger.info(f"φ-1 = {self.phi_conjugate}")
    
    def map_dna_sequence(self, sequence: str) -> List[mp.mpf]:
        """
        Map DNA sequence to numerical values using the standard mapping.
        
        Args:
            sequence: DNA sequence string (A, T, C, G)
            
        Returns:
            List of high-precision numerical values
        """
        if not sequence:
            raise ValueError("DNA sequence cannot be empty")
        
        # Validate sequence
        valid_bases = set('ATCG')
        invalid_bases = set(sequence.upper()) - valid_bases
        if invalid_bases:
            raise ValueError(f"Invalid bases found: {invalid_bases}")
        
        # Map sequence to numbers
        mapped_values = [DNA_MAPPING[base.upper()] for base in sequence]
        
        logger.debug(f"Mapped sequence of length {len(sequence)} to numerical values")
        return mapped_values
    
    def calculate_geodesic_resolution(self, n: int, k: float = 0.3) -> mp.mpf:
        """
        Calculate geodesic resolution function θ'(n, k) = φ·((n mod φ)/φ)^k
        
        Args:
            n: Position index
            k: Resolution exponent (default: 0.3)
            
        Returns:
            High-precision geodesic resolution value
        """
        n_mod_phi = mp.fmod(mp.mpf(n), self.phi)
        ratio = n_mod_phi / self.phi
        theta_prime = self.phi * (ratio ** mp.mpf(k))
        
        return theta_prime
    
    def calculate_delta_n(self, sequence_values: List[mp.mpf], position: int) -> mp.mpf:
        """
        Calculate Δₙ for position n in the sequence.
        
        Args:
            sequence_values: Mapped DNA sequence values
            position: Position index (0-based)
            
        Returns:
            High-precision delta value at position n
        """
        if position >= len(sequence_values):
            raise IndexError(f"Position {position} exceeds sequence length {len(sequence_values)}")
        
        n = len(sequence_values)
        base_value = sequence_values[position]
        
        # Calculate positional weighting with geodesic resolution
        theta_prime = self.calculate_geodesic_resolution(position)
        
        # Enhanced delta calculation incorporating Z Framework principles
        # Δₙ = base_value * θ'(n,k) * (1 + position_weight)
        position_weight = mp.mpf(position) / mp.mpf(n)
        delta_n = base_value * theta_prime * (1 + position_weight)
        
        return delta_n
    
    def calculate_delta_max(self, sequence_values: List[mp.mpf]) -> mp.mpf:
        """
        Calculate Δₜₐₓ (maximum theoretical delta) for the sequence.
        
        Args:
            sequence_values: Mapped DNA sequence values
            
        Returns:
            High-precision maximum delta value
        """
        if not sequence_values:
            raise ValueError("Sequence values cannot be empty")
        
        n = len(sequence_values)
        max_base_value = mp.mpf('4')  # Maximum base value (G)
        
        # Calculate theoretical maximum using golden ratio scaling
        # Δₘₐₓ = max_base * φ * (1 + φ-1) = max_base * φ²
        delta_max = max_base_value * (self.phi ** 2)
        
        # Scale by sequence length for normalization
        delta_max *= mp.sqrt(mp.mpf(n))
        
        return delta_max
    
    def calculate_z_values(self, sequence: str) -> Dict[str, mp.mpf]:
        """
        Calculate Z values for the entire DNA sequence using Z = n(Δₙ/Δₘₐₓ).
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Dictionary containing Z values and related metrics
        """
        # Map sequence to numerical values
        sequence_values = self.map_dna_sequence(sequence)
        n = len(sequence_values)
        
        # Calculate Δₘₐₓ
        delta_max = self.calculate_delta_max(sequence_values)
        
        # Calculate Z values for each position
        z_values = []
        delta_values = []
        
        for i in range(n):
            delta_n = self.calculate_delta_n(sequence_values, i)
            z_n = mp.mpf(n) * (delta_n / delta_max)
            
            z_values.append(z_n)
            delta_values.append(delta_n)
        
        # Calculate statistical measures
        z_mean = sum(z_values) / mp.mpf(len(z_values))
        z_variance = sum((z - z_mean) ** 2 for z in z_values) / mp.mpf(len(z_values))
        z_std = mp.sqrt(z_variance)
        
        # Calculate convergence metrics
        phi_conjugate_convergence = abs(z_mean - self.phi_conjugate)
        variance_convergence = abs(z_variance - TARGET_VARIANCE)
        
        results = {
            'sequence_length': mp.mpf(n),
            'z_values': z_values,
            'delta_values': delta_values,
            'delta_max': delta_max,
            'z_mean': z_mean,
            'z_variance': z_variance,
            'z_std': z_std,
            'phi_conjugate_target': self.phi_conjugate,
            'phi_conjugate_convergence': phi_conjugate_convergence,
            'target_variance': TARGET_VARIANCE,
            'variance_convergence': variance_convergence,
            'converges_to_phi_conjugate': phi_conjugate_convergence < mp.mpf('0.01'),
            'converges_to_target_variance': variance_convergence < mp.mpf('0.01')
        }
        
        logger.info(f"Calculated Z values for sequence of length {n}")
        logger.info(f"Z mean: {z_mean}, convergence to φ-1: {phi_conjugate_convergence}")
        logger.info(f"Z variance: {z_variance}, convergence to target: {variance_convergence}")
        
        return results
    
    def calculate_density_enhancement(self, sequence: str, window_size: int = 10) -> Dict[str, mp.mpf]:
        """
        Calculate density enhancement across sliding windows.
        
        Args:
            sequence: DNA sequence string
            window_size: Size of sliding window for analysis
            
        Returns:
            Dictionary containing density enhancement metrics
        """
        sequence_values = self.map_dna_sequence(sequence)
        n = len(sequence_values)
        
        if window_size > n:
            window_size = n
        
        density_values = []
        enhancement_values = []
        
        for i in range(n - window_size + 1):
            window_values = sequence_values[i:i + window_size]
            
            # Calculate local density
            window_sum = sum(window_values)
            local_density = window_sum / mp.mpf(window_size)
            
            # Calculate enhancement using Z Framework principles
            z_results = self.calculate_z_values(''.join(['ATCG'[int(v)-1] for v in window_values]))
            enhancement = z_results['z_mean'] * local_density / self.phi
            
            density_values.append(local_density)
            enhancement_values.append(enhancement)
        
        # Calculate global statistics
        mean_density = sum(density_values) / mp.mpf(len(density_values))
        mean_enhancement = sum(enhancement_values) / mp.mpf(len(enhancement_values))
        
        results = {
            'window_size': mp.mpf(window_size),
            'num_windows': mp.mpf(len(density_values)),
            'density_values': density_values,
            'enhancement_values': enhancement_values,
            'mean_density': mean_density,
            'mean_enhancement': mean_enhancement,
            'density_variance': sum((d - mean_density) ** 2 for d in density_values) / mp.mpf(len(density_values)),
            'enhancement_variance': sum((e - mean_enhancement) ** 2 for e in enhancement_values) / mp.mpf(len(enhancement_values))
        }
        
        return results
    
    def perform_falsification_test(self, sequence: str, num_perturbations: int = 100, 
                                 perturbation_rate: float = 0.1) -> Dict[str, any]:
        """
        Perform falsification tests with random perturbations.
        
        Args:
            sequence: Original DNA sequence
            num_perturbations: Number of random perturbations to test
            perturbation_rate: Fraction of bases to randomly change
            
        Returns:
            Dictionary containing falsification test results
        """
        import random
        random.seed(42)  # For reproducibility
        
        # Calculate baseline metrics
        baseline_results = self.calculate_z_values(sequence)
        baseline_mean = baseline_results['z_mean']
        baseline_variance = baseline_results['z_variance']
        
        # Storage for perturbation results
        perturbed_means = []
        perturbed_variances = []
        convergence_preserved = []
        
        bases = ['A', 'T', 'C', 'G']
        
        for i in range(num_perturbations):
            # Create perturbed sequence
            seq_list = list(sequence)
            num_changes = max(1, int(len(sequence) * perturbation_rate))
            
            positions_to_change = random.sample(range(len(sequence)), num_changes)
            for pos in positions_to_change:
                # Change to a different random base
                current_base = seq_list[pos]
                available_bases = [b for b in bases if b != current_base]
                seq_list[pos] = random.choice(available_bases)
            
            perturbed_sequence = ''.join(seq_list)
            
            try:
                # Calculate metrics for perturbed sequence
                perturbed_results = self.calculate_z_values(perturbed_sequence)
                
                perturbed_means.append(perturbed_results['z_mean'])
                perturbed_variances.append(perturbed_results['z_variance'])
                
                # Check if convergence properties are preserved
                phi_convergence = abs(perturbed_results['z_mean'] - self.phi_conjugate) < mp.mpf('0.05')
                var_convergence = abs(perturbed_results['z_variance'] - TARGET_VARIANCE) < mp.mpf('0.05')
                convergence_preserved.append(phi_convergence and var_convergence)
                
            except Exception as e:
                logger.warning(f"Perturbation {i} failed: {e}")
                continue
        
        # Calculate statistics
        if perturbed_means:
            mean_stability = sum(perturbed_means) / mp.mpf(len(perturbed_means))
            variance_stability = sum(perturbed_variances) / mp.mpf(len(perturbed_variances))
            convergence_rate = sum(convergence_preserved) / len(convergence_preserved)
            
            # Deviation from baseline
            mean_deviation = abs(mean_stability - baseline_mean)
            variance_deviation = abs(variance_stability - baseline_variance)
        else:
            mean_stability = variance_stability = convergence_rate = mp.mpf('0')
            mean_deviation = variance_deviation = mp.mpf('0')
        
        results = {
            'baseline_mean': baseline_mean,
            'baseline_variance': baseline_variance,
            'num_perturbations': num_perturbations,
            'perturbation_rate': perturbation_rate,
            'perturbed_means': perturbed_means,
            'perturbed_variances': perturbed_variances,
            'mean_stability': mean_stability,
            'variance_stability': variance_stability,
            'convergence_rate': convergence_rate,
            'mean_deviation': mean_deviation,
            'variance_deviation': variance_deviation,
            'falsification_passed': convergence_rate > 0.7  # 70% threshold for passing
        }
        
        logger.info(f"Falsification test completed: {len(perturbed_means)}/{num_perturbations} successful")
        logger.info(f"Convergence rate: {convergence_rate:.3f}")
        
        return results


def format_mpmath_for_json(value):
    """Convert mpmath values to float for JSON serialization."""
    if isinstance(value, mp.mpf):
        return float(value)
    elif isinstance(value, list):
        return [format_mpmath_for_json(item) for item in value]
    elif isinstance(value, dict):
        return {key: format_mpmath_for_json(val) for key, val in value.items()}
    else:
        return value

def format_mpmath_for_display(value, precision=6):
    """Format mpmath values for display with specified precision."""
    if isinstance(value, mp.mpf):
        return f"{float(value):.{precision}f}"
    else:
        return str(value)