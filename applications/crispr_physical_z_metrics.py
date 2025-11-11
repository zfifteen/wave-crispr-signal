"""
Wave-CRISPR Discrete Biological Z-Metrics Implementation

This module implements four sequence-derivable discrete biological Z-metrics for Wave-CRISPR:
1. Base-pair opening kinetics
2. Base stacking dissociation  
3. Helical twist fluctuation rate
4. Denaturation/melting kinetics near Tm

Uses the discrete biological Z form: Z = A * (B / c) with c = e² ≈ 7.389
where A is sequence-dependent mean and B is a rate of change between adjacent bases.
"""

import mpmath as mp
import numpy as np
import re
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import logging

# Configure high precision
mp.dps = 50

# Discrete biological constants
E_SQUARED = mp.e ** 2  # c = e² ≈ 7.389
BOLTZMANN_KB = mp.mpf('1.380649e-23')  # J/K
ROOM_TEMP = mp.mpf('298.15')  # K
RT = BOLTZMANN_KB * ROOM_TEMP

# DNA physical parameters (empirical values from literature)
BASE_PAIR_ENERGIES = {
    'AT': mp.mpf('-2.0'),  # kcal/mol
    'TA': mp.mpf('-2.0'),
    'GC': mp.mpf('-3.5'),
    'CG': mp.mpf('-3.5'),
    'AA': mp.mpf('-1.2'), 'TT': mp.mpf('-1.2'),
    'GG': mp.mpf('-1.8'), 'CC': mp.mpf('-1.8'),
    'AG': mp.mpf('-1.3'), 'GA': mp.mpf('-1.3'),
    'CT': mp.mpf('-1.4'), 'TC': mp.mpf('-1.4'),
    'AC': mp.mpf('-1.1'), 'CA': mp.mpf('-1.1'),
    'GT': mp.mpf('-1.1'), 'TG': mp.mpf('-1.1'),
}

STACKING_ENERGIES = {
    'AA': mp.mpf('-1.0'), 'TT': mp.mpf('-1.0'),
    'AT': mp.mpf('-0.88'), 'TA': mp.mpf('-0.58'),
    'GG': mp.mpf('-1.3'), 'CC': mp.mpf('-1.3'),
    'GC': mp.mpf('-2.17'), 'CG': mp.mpf('-1.44'),
    'AG': mp.mpf('-1.28'), 'GA': mp.mpf('-1.30'),
    'CT': mp.mpf('-1.28'), 'TC': mp.mpf('-1.30'),
    'AC': mp.mpf('-1.45'), 'CA': mp.mpf('-1.45'),
    'GT': mp.mpf('-1.28'), 'TG': mp.mpf('-1.28'),
}

# Melting temperature parameters
TM_SALT_CORRECTION = mp.mpf('16.6')  # log10([Na+]) correction
DEFAULT_SALT_CONC = mp.mpf('0.05')  # 50 mM Na+

logger = logging.getLogger(__name__)


class DNAValidationError(ValueError):
    """Raised when DNA sequence validation fails."""
    pass


def validate_human_dna_sequence(sequence: str, header: str = "") -> None:
    """
    Validate DNA sequence for human origin and nucleotide composition.
    
    Args:
        sequence: DNA sequence string
        header: FASTA header for source validation
        
    Raises:
        DNAValidationError: If validation fails
    """
    # Clean sequence
    sequence = sequence.upper().strip()
    
    # Check for valid nucleotide alphabet (A/C/G/T/N only)
    valid_bases = set('ACGTN')
    if not set(sequence).issubset(valid_bases):
        invalid_chars = set(sequence) - valid_bases
        raise DNAValidationError(
            f"Non-nucleotide or ambiguous input: {invalid_chars}. "
            f"Only A/C/G/T/N allowed."
        )
    
    # Check for human genome reference or curated CRISPR benchmarks
    human_indicators = [
        'grch38', 'hg38', 'human', 'homo sapiens',
        'pcsk9', 'crispr', 'benchmark'
    ]
    
    header_lower = header.lower()
    is_human_ref = any(indicator in header_lower for indicator in human_indicators)
    
    if not is_human_ref and header:
        logger.warning(
            f"Sequence source not clearly identified as human (GRCh38/hg38) "
            f"or curated CRISPR benchmark. Header: {header[:100]}"
        )
    
    # Basic quality checks
    if len(sequence) < 10:
        raise DNAValidationError("Sequence too short (minimum 10 bp)")
    
    if len(sequence) > 10000:
        raise DNAValidationError("Sequence too long (maximum 10,000 bp)")
    
    # Check for excessive N content
    n_content = sequence.count('N') / len(sequence)
    if n_content > 0.1:
        raise DNAValidationError(f"Excessive ambiguous bases: {n_content:.2%} N content")


class PhysicalZMetricsCalculator:
    """
    Calculator for sequence-derivable physical Z-metrics in Wave-CRISPR.
    
    Implements Z = A * (B / c) where:
    - A: Frame-dependent parameter (sequence context)
    - B: Rate of change (kinetic parameter) 
    - c: Invariant constant (e² ≈ 7.389)
    """
    
    def __init__(self, precision_dps: int = 50):
        """Initialize calculator with high precision."""
        mp.dps = precision_dps
        self.e_squared = E_SQUARED
        
    def _safe_divide(self, numerator: mp.mpf, denominator: mp.mpf, 
                     default: mp.mpf = mp.mpf('0')) -> mp.mpf:
        """Safe division with zero guard."""
        if abs(denominator) < mp.mpf('1e-15'):
            logger.warning("Division by zero avoided in Z-metric calculation")
            return default
        return numerator / denominator
        
    def _clamp_rate(self, rate: mp.mpf) -> mp.mpf:
        """Clamp rates to B ≥ 0 as per guardrails."""
        if rate < 0:
            logger.warning(f"Negative rate {rate} clamped to 0")
            return mp.mpf('0')
        return rate
    
    def calculate_base_pair_opening_kinetics(self, sequence: str) -> Dict[str, mp.mpf]:
        """
        Calculate base-pair opening kinetics Z-metric.
        
        Z_opening = A_context * (B_opening_rate / e²)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with opening kinetics metrics
        """
        if len(sequence) < 2:
            return {'z_opening': mp.mpf('0'), 'opening_rate': mp.mpf('0')}
        
        # A_context: Sequence context factor (normalized GC content + position weight)
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        position_weights = [mp.mpf(i + 1) / len(sequence) for i in range(len(sequence))]
        a_context = mp.mpf(gc_content) * sum(position_weights) / len(position_weights)
        
        # B_opening_rate: Base-pair opening rate from thermal fluctuations
        opening_rates = []
        for i in range(len(sequence) - 1):
            bp = sequence[i:i+2]
            if bp in BASE_PAIR_ENERGIES:
                # Opening rate ~ exp(-ΔG / RT)
                delta_g = BASE_PAIR_ENERGIES[bp]  # Already in kcal/mol
                # Convert to J/mol: 1 kcal/mol = 4184 J/mol
                delta_g_j = delta_g * mp.mpf('4184')
                
                # Opening rate (simplified Arrhenius)
                opening_rate = mp.exp(-abs(delta_g_j) / (RT * mp.mpf('6.022e23')))
                opening_rates.append(opening_rate)
        
        if not opening_rates:
            b_opening_rate = mp.mpf('0')
        else:
            b_opening_rate = sum(opening_rates) / len(opening_rates)
        
        # Apply rate clamping
        b_opening_rate = self._clamp_rate(b_opening_rate)
        
        # Calculate Z_opening = A * (B / c)
        z_opening = a_context * self._safe_divide(b_opening_rate, self.e_squared)
        
        return {
            'z_opening': z_opening,
            'opening_rate': b_opening_rate,
            'context_factor': a_context,
            'gc_content': mp.mpf(gc_content)
        }
    
    def calculate_base_stacking_dissociation(self, sequence: str) -> Dict[str, mp.mpf]:
        """
        Calculate base stacking dissociation Z-metric.
        
        Z_stacking = A_stack * (B_dissociation_rate / e²)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with stacking dissociation metrics
        """
        if len(sequence) < 2:
            return {'z_stacking': mp.mpf('0'), 'dissociation_rate': mp.mpf('0')}
        
        # A_stack: Stacking context factor
        stack_count = 0
        total_stacking_energy = mp.mpf('0')
        
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            if dinuc in STACKING_ENERGIES:
                stack_count += 1
                total_stacking_energy += abs(STACKING_ENERGIES[dinuc])
        
        a_stack = total_stacking_energy / mp.mpf(max(1, stack_count))
        
        # B_dissociation_rate: Stacking dissociation rate
        dissociation_rates = []
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            if dinuc in STACKING_ENERGIES:
                # Dissociation rate ~ exp(-ΔG_stack / RT)
                delta_g_stack = STACKING_ENERGIES[dinuc]
                delta_g_j = delta_g_stack * mp.mpf('4184')
                
                dissoc_rate = mp.exp(-abs(delta_g_j) / (RT * mp.mpf('6.022e23')))
                dissociation_rates.append(dissoc_rate)
        
        if not dissociation_rates:
            b_dissociation_rate = mp.mpf('0')
        else:
            b_dissociation_rate = sum(dissociation_rates) / len(dissociation_rates)
        
        # Apply rate clamping
        b_dissociation_rate = self._clamp_rate(b_dissociation_rate)
        
        # Calculate Z_stacking = A * (B / c)
        z_stacking = a_stack * self._safe_divide(b_dissociation_rate, self.e_squared)
        
        return {
            'z_stacking': z_stacking,
            'dissociation_rate': b_dissociation_rate,
            'stack_factor': a_stack,
            'stacking_count': mp.mpf(stack_count)
        }
    
    def calculate_helical_twist_fluctuation(self, sequence: str) -> Dict[str, mp.mpf]:
        """
        Calculate helical twist fluctuation rate Z-metric.
        
        Z_twist = A_geometry * (B_twist_rate / e²)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with twist fluctuation metrics
        """
        if len(sequence) < 3:
            return {'z_twist': mp.mpf('0'), 'twist_rate': mp.mpf('0')}
        
        # A_geometry: Geometric context factor based on sequence periodicity
        # Calculate autocorrelation to detect helical periodicity (~10.5 bp)
        base_values = {'A': 1, 'T': 2, 'C': 3, 'G': 4, 'N': 2.5}
        numeric_seq = [base_values.get(base, 2.5) for base in sequence]
        
        # Simple autocorrelation at period 10-11
        periods = [10, 11]
        max_correlation = mp.mpf('0')
        
        for period in periods:
            if len(numeric_seq) > period:
                corr_sum = mp.mpf('0')
                for i in range(len(numeric_seq) - period):
                    corr_sum += abs(numeric_seq[i] - numeric_seq[i + period])
                # Inverse correlation (lower difference = higher correlation)
                correlation = mp.mpf('1') / (mp.mpf('1') + corr_sum / (len(numeric_seq) - period))
                max_correlation = max(max_correlation, correlation)
        
        a_geometry = max_correlation
        
        # B_twist_rate: Helical twist fluctuation rate
        # Based on local flexibility and base-pair step parameters
        twist_rates = []
        for i in range(len(sequence) - 2):
            triplet = sequence[i:i+3]
            
            # Twist rate influenced by dinucleotide flexibility
            dinuc1 = triplet[:2]
            dinuc2 = triplet[1:]
            
            # Flexibility scale (empirical, higher for AT-rich)
            flexibility_score = (
                (dinuc1.count('A') + dinuc1.count('T')) / 2 +
                (dinuc2.count('A') + dinuc2.count('T')) / 2
            ) / 2
            
            # Twist fluctuation rate ~ flexibility * thermal energy
            twist_rate = mp.mpf(flexibility_score) * mp.exp(-mp.mpf('1') / (RT * mp.mpf('1e20')))
            twist_rates.append(twist_rate)
        
        if not twist_rates:
            b_twist_rate = mp.mpf('0')
        else:
            b_twist_rate = sum(twist_rates) / len(twist_rates)
        
        # Apply rate clamping
        b_twist_rate = self._clamp_rate(b_twist_rate)
        
        # Calculate Z_twist = A * (B / c)
        z_twist = a_geometry * self._safe_divide(b_twist_rate, self.e_squared)
        
        return {
            'z_twist': z_twist,
            'twist_rate': b_twist_rate,
            'geometry_factor': a_geometry,
            'helical_correlation': max_correlation
        }
    
    def calculate_denaturation_melting_kinetics(self, sequence: str, 
                                              salt_conc: mp.mpf = DEFAULT_SALT_CONC) -> Dict[str, mp.mpf]:
        """
        Calculate denaturation/melting kinetics near Tm Z-metric.
        
        Z_melting = A_stability * (B_melting_rate / e²)
        
        Args:
            sequence: DNA sequence
            salt_conc: Salt concentration in M (default: 0.05 M)
            
        Returns:
            Dictionary with melting kinetics metrics
        """
        if len(sequence) < 2:
            return {'z_melting': mp.mpf('0'), 'melting_rate': mp.mpf('0')}
        
        # Calculate melting temperature (Tm)
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        
        # Simple Tm estimation (Wallace rule + salt correction)
        tm_base = mp.mpf('81.5') + TM_SALT_CORRECTION * mp.log10(salt_conc) + mp.mpf('41') * mp.mpf(gc_content) - mp.mpf('675') / mp.mpf(len(sequence))
        
        # A_stability: Thermodynamic stability factor
        # Based on total duplex stability relative to melting point
        total_binding_energy = mp.mpf('0')
        for i in range(len(sequence) - 1):
            bp = sequence[i:i+2]
            if bp in BASE_PAIR_ENERGIES:
                total_binding_energy += abs(BASE_PAIR_ENERGIES[bp])
        
        # Stability factor normalized by sequence length
        a_stability = total_binding_energy / mp.mpf(len(sequence))
        
        # B_melting_rate: Melting kinetics rate near Tm
        # Rate depends on proximity to melting temperature and cooperativity
        temp_diff = abs(ROOM_TEMP - tm_base)  # Distance from Tm
        
        # Melting rate ~ exp(-activation_energy / RT) 
        # Higher near Tm, lower when far from Tm
        activation_energy = total_binding_energy * mp.mpf('0.1')  # Simplified
        melting_rate_base = mp.exp(-activation_energy / mp.mpf('10'))
        
        # Cooperativity factor (melting is cooperative process)
        cooperativity = mp.mpf('1') + mp.mpf(len(sequence)) / mp.mpf('50')
        
        b_melting_rate = melting_rate_base * cooperativity / (mp.mpf('1') + temp_diff / mp.mpf('10'))
        
        # Apply rate clamping
        b_melting_rate = self._clamp_rate(b_melting_rate)
        
        # Calculate Z_melting = A * (B / c)
        z_melting = a_stability * self._safe_divide(b_melting_rate, self.e_squared)
        
        return {
            'z_melting': z_melting,
            'melting_rate': b_melting_rate,
            'stability_factor': a_stability,
            'estimated_tm': tm_base,
            'cooperativity': cooperativity
        }
    
    def calculate_all_physical_z_metrics(self, sequence: str, 
                                       header: str = "",
                                       validate: bool = True) -> Dict[str, Dict]:
        """
        Calculate all four physical Z-metrics for a sequence.
        
        Args:
            sequence: DNA sequence
            header: FASTA header for validation
            validate: Whether to perform validation (default: True)
            
        Returns:
            Dictionary with all Z-metrics
            
        Raises:
            DNAValidationError: If validation fails
        """
        # Input validation with guardrails
        if validate:
            validate_human_dna_sequence(sequence, header)
        
        # Clean sequence
        sequence = sequence.upper().replace('N', 'A')  # Replace N with A for calculations
        
        # Calculate all four Z-metrics
        results = {
            'opening_kinetics': self.calculate_base_pair_opening_kinetics(sequence),
            'stacking_dissociation': self.calculate_base_stacking_dissociation(sequence),
            'twist_fluctuation': self.calculate_helical_twist_fluctuation(sequence),
            'melting_kinetics': self.calculate_denaturation_melting_kinetics(sequence),
        }
        
        # Add sequence metadata
        results['sequence_info'] = {
            'length': mp.mpf(len(sequence)),
            'gc_content': mp.mpf((sequence.count('G') + sequence.count('C')) / len(sequence)),
            'header': header,
            'validated': validate
        }
        
        # Calculate summary metrics
        z_values = [
            results['opening_kinetics']['z_opening'],
            results['stacking_dissociation']['z_stacking'], 
            results['twist_fluctuation']['z_twist'],
            results['melting_kinetics']['z_melting']
        ]
        
        results['summary'] = {
            'z_mean': sum(z_values) / len(z_values),
            'z_variance': sum((z - sum(z_values)/len(z_values))**2 for z in z_values) / len(z_values),
            'z_max': max(z_values),
            'z_min': min(z_values)
        }
        
        return results


def read_fasta_with_validation(filepath: str) -> Dict[str, str]:
    """
    Read FASTA file with human DNA validation.
    
    Args:
        filepath: Path to FASTA file
        
    Returns:
        Dictionary mapping headers to sequences
        
    Raises:
        DNAValidationError: If validation fails
    """
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Process previous sequence
                if current_header is not None:
                    seq = ''.join(current_seq)
                    validate_human_dna_sequence(seq, current_header)
                    sequences[current_header] = seq
                
                # Start new sequence
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())
    
    # Process last sequence
    if current_header is not None:
        seq = ''.join(current_seq)
        validate_human_dna_sequence(seq, current_header)
        sequences[current_header] = seq
    
    return sequences


if __name__ == "__main__":
    # Example usage
    calculator = PhysicalZMetricsCalculator()
    
    # Test sequence (PCSK9 example)
    test_sequence = "ATGCGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"
    test_header = "PCSK9 test sequence GRCh38"
    
    try:
        results = calculator.calculate_all_physical_z_metrics(test_sequence, test_header)
        
        print("🧬 Physical Z-Metrics Results")
        print("=" * 50)
        print(f"Sequence: {test_sequence}")
        print(f"Length: {float(results['sequence_info']['length'])}")
        print(f"GC Content: {float(results['sequence_info']['gc_content']):.3f}")
        print()
        
        print("Z-Metrics:")
        print(f"  Opening Kinetics Z: {float(results['opening_kinetics']['z_opening']):.6f}")
        print(f"  Stacking Dissociation Z: {float(results['stacking_dissociation']['z_stacking']):.6f}")
        print(f"  Twist Fluctuation Z: {float(results['twist_fluctuation']['z_twist']):.6f}")
        print(f"  Melting Kinetics Z: {float(results['melting_kinetics']['z_melting']):.6f}")
        print()
        
        print("Summary:")
        print(f"  Z Mean: {float(results['summary']['z_mean']):.6f}")
        print(f"  Z Variance: {float(results['summary']['z_variance']):.6f}")
        
    except DNAValidationError as e:
        print(f"Validation Error: {e}")