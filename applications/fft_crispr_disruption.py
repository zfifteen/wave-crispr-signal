#!/usr/bin/env python3
"""
FFT-Based CRISPR Disruption Metrics with Golden-Ratio Phase Analysis

This module implements frequency-domain analysis of CRISPR-Cas9 edited DNA sequences
using FFT with golden-ratio-derived phase weighting θ′(n,k) = φ·((n mod φ_period)/φ_period)^k
to detect off-target periodicities and insertion/deletion disruptions.

Scientific Gates:
- Human DNA only (A/C/G/T for DNA sequences)
- Geometric resolution: θ′(n,k) with k ≈ 0.3 (default φ_period = 21 for 21-nt guides)
- Z invariants: Z = A(B/e²) with documented parameters
"""

import numpy as np
from typing import Dict, List
from scipy.fft import fft, fftfreq
from scipy.stats import entropy
import logging

# Configure logging
logger = logging.getLogger(__name__)

# Mathematical constants
PHI = 1.6180339887498949  # Golden ratio
E_SQUARED = 7.3890560989306495  # e²

# Default geometric period for CRISPR guides (21-nt)
DEFAULT_PHI_PERIOD = 21.0

# Default resolution exponent
DEFAULT_K = 0.3


class FFTCRISPRDisruptionAnalyzer:
    """
    Analyze CRISPR editing disruptions using FFT with golden-ratio phase weighting.
    
    This analyzer detects off-target periodicities and disruptions in DNA sequences
    by applying θ′(n,k) = φ·((n mod φ_period)/φ_period)^k phase weighting to FFT spectrum bins.
    """
    
    def __init__(self, phi_period: float = DEFAULT_PHI_PERIOD, k: float = DEFAULT_K):
        """
        Initialize FFT-based CRISPR disruption analyzer.
        
        Args:
            phi_period: Geometric period φ for resolution (default: 21 for 21-nt guides)
            k: Resolution exponent (default: 0.3)
        """
        self.phi_period = phi_period
        self.k = k
        self.phi = PHI
        self.e_squared = E_SQUARED
        
        logger.info(f"Initialized FFT CRISPR Disruption Analyzer")
        logger.info(f"  φ-period: {phi_period}")
        logger.info(f"  k: {k}")
    
    def validate_dna_sequence(self, sequence: str) -> str:
        """
        Validate DNA sequence contains only A/C/G/T/N.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Uppercase validated sequence
            
        Raises:
            ValueError: If sequence contains invalid bases
        """
        if not sequence:
            raise ValueError("DNA sequence cannot be empty")
        
        sequence = sequence.upper()
        
        # Check for U (RNA) in DNA context first - more specific error
        if 'U' in sequence:
            raise ValueError(
                "Found 'U' in sequence. For DNA sequences, use 'T'. "
                "For RNA sequences, use a separate RNA validation method."
            )
        
        valid_bases = set('ACGTN')
        invalid_bases = set(sequence) - valid_bases
        
        if invalid_bases:
            raise ValueError(
                f"Invalid DNA bases found: {invalid_bases}. "
                f"Only A/C/G/T/N allowed for DNA sequences."
            )
        
        return sequence
    
    def encode_dna_complex(self, sequence: str) -> np.ndarray:
        """
        Encode DNA sequence as complex waveform.
        
        Standard mapping:
        - A → 1 + 0j (positive real)
        - T → -1 + 0j (negative real)
        - C → 0 + 1j (positive imaginary)
        - G → 0 - 1j (negative imaginary)
        - N → 0 + 0j (unknown)
        
        Args:
            sequence: Validated DNA sequence
            
        Returns:
            Complex numpy array
        """
        mapping = {
            'A': 1.0 + 0.0j,
            'T': -1.0 + 0.0j,
            'C': 0.0 + 1.0j,
            'G': 0.0 - 1.0j,
            'N': 0.0 + 0.0j
        }
        
        return np.array([mapping[base] for base in sequence], dtype=np.complex128)
    
    def calculate_theta_prime(self, n: int) -> float:
        """
        Calculate geometric resolution θ′(n,k) = φ·((n mod φ_period)/φ_period)^k.
        
        Args:
            n: Position index (frequency bin)
            
        Returns:
            Geometric resolution weight
        """
        if n == 0:
            # Handle DC component specially - use small offset
            return self.phi * (0.01 ** self.k)
        
        n_mod_phi = n % self.phi_period
        ratio = n_mod_phi / self.phi_period
        
        # Prevent zero weight - add small offset if needed
        if ratio < 0.01:
            ratio = 0.01
        
        theta_prime = self.phi * (ratio ** self.k)
        
        return theta_prime
    
    def apply_golden_phase_weights(self, fft_spectrum: np.ndarray) -> np.ndarray:
        """
        Apply golden-ratio phase weights θ′(n,k) to FFT spectrum bins.
        
        Args:
            fft_spectrum: Raw FFT spectrum (complex or magnitude)
            
        Returns:
            Phase-weighted spectrum
        """
        n_bins = len(fft_spectrum)
        weights = np.array([self.calculate_theta_prime(n) for n in range(n_bins)])
        
        # Apply weights to magnitude spectrum
        if np.iscomplexobj(fft_spectrum):
            weighted_spectrum = np.abs(fft_spectrum) * weights
        else:
            weighted_spectrum = fft_spectrum * weights
        
        return weighted_spectrum
    
    def detect_off_target_periodicities(
        self, 
        sequence: str,
        threshold_percentile: float = 80.0
    ) -> Dict[str, any]:
        """
        Detect off-target periodicities using FFT with golden-ratio phase weighting.
        
        Args:
            sequence: DNA sequence to analyze
            threshold_percentile: Percentile threshold for peak detection
            
        Returns:
            Dictionary with periodicity analysis results
        """
        # Validate and encode sequence
        sequence = self.validate_dna_sequence(sequence)
        complex_wave = self.encode_dna_complex(sequence)
        
        # Compute FFT
        fft_result = fft(complex_wave)
        fft_magnitude = np.abs(fft_result)
        
        # Apply golden-ratio phase weights
        weighted_spectrum = self.apply_golden_phase_weights(fft_magnitude)
        
        # Calculate frequencies
        frequencies = fftfreq(len(sequence))
        
        # Detect peaks in weighted spectrum (use lower threshold to catch more peaks)
        if len(weighted_spectrum) > 0:
            threshold = np.percentile(weighted_spectrum, threshold_percentile)
        else:
            threshold = 0
        peak_indices = np.where(weighted_spectrum > threshold)[0]
        
        # Extract significant periodicities (only positive frequencies)
        n_half = len(sequence) // 2 if len(sequence) > 1 else 1
        significant_peaks = []
        
        for idx in peak_indices:
            if 0 < idx < n_half:  # Skip DC and negative frequencies
                freq = frequencies[idx]
                period = 1.0 / freq if freq != 0 else np.inf
                magnitude = fft_magnitude[idx]
                weighted_mag = weighted_spectrum[idx]
                theta_weight = self.calculate_theta_prime(idx)
                
                significant_peaks.append({
                    'frequency': freq,
                    'period': period,
                    'magnitude': magnitude,
                    'weighted_magnitude': weighted_mag,
                    'theta_prime_weight': theta_weight,
                    'bin_index': idx
                })
        
        # Sort by weighted magnitude
        significant_peaks = sorted(
            significant_peaks, 
            key=lambda x: x['weighted_magnitude'], 
            reverse=True
        )
        
        return {
            'sequence_length': len(sequence),
            'n_significant_peaks': len(significant_peaks),
            'significant_peaks': significant_peaks[:10],  # Top 10
            'fft_spectrum': fft_magnitude,
            'weighted_spectrum': weighted_spectrum,
            'frequencies': frequencies
        }
    
    def calculate_disruption_score(
        self,
        reference_seq: str,
        edited_seq: str
    ) -> Dict[str, float]:
        """
        Calculate disruption score comparing reference and edited sequences.
        
        Measures spectral disruption using:
        - ΔEntropy: Change in spectral entropy
        - Δf₁: Change in dominant frequency magnitude
        - Sidelobe count change
        - Phase coherence disruption
        
        Args:
            reference_seq: Reference (unedited) DNA sequence
            edited_seq: Edited DNA sequence (post-CRISPR)
            
        Returns:
            Dictionary with disruption metrics
        """
        # Validate sequences
        reference_seq = self.validate_dna_sequence(reference_seq)
        edited_seq = self.validate_dna_sequence(edited_seq)
        
        # Compute spectra for both sequences
        ref_analysis = self.detect_off_target_periodicities(reference_seq)
        edit_analysis = self.detect_off_target_periodicities(edited_seq)
        
        # Calculate spectral entropy for both
        ref_spectrum = ref_analysis['weighted_spectrum']
        edit_spectrum = edit_analysis['weighted_spectrum']
        
        ref_entropy = self._calculate_spectral_entropy(ref_spectrum)
        edit_entropy = self._calculate_spectral_entropy(edit_spectrum)
        delta_entropy = edit_entropy - ref_entropy
        
        # Calculate dominant frequency change
        ref_f1 = np.max(ref_spectrum[1:len(ref_spectrum)//2])  # Skip DC
        edit_f1 = np.max(edit_spectrum[1:len(edit_spectrum)//2])
        delta_f1 = edit_f1 - ref_f1
        
        # Sidelobe count change
        ref_sidelobes = len(ref_analysis['significant_peaks'])
        edit_sidelobes = len(edit_analysis['significant_peaks'])
        delta_sidelobes = edit_sidelobes - ref_sidelobes
        
        # Phase coherence disruption (alignment of periodicities)
        phase_disruption = self._calculate_phase_disruption(
            ref_analysis['significant_peaks'],
            edit_analysis['significant_peaks']
        )
        
        # Composite disruption score
        # Higher score = more disruption = worse for off-target
        disruption_score = (
            abs(delta_entropy) * 0.3 +
            abs(delta_f1) / max(ref_f1, 1e-10) * 0.3 +
            abs(delta_sidelobes) / max(ref_sidelobes, 1) * 0.2 +
            phase_disruption * 0.2
        )
        
        return {
            'disruption_score': disruption_score,
            'delta_entropy': delta_entropy,
            'delta_f1': delta_f1,
            'delta_f1_relative': delta_f1 / max(ref_f1, 1e-10),
            'delta_sidelobes': delta_sidelobes,
            'phase_disruption': phase_disruption,
            'reference_entropy': ref_entropy,
            'edited_entropy': edit_entropy,
            'reference_peaks': ref_sidelobes,
            'edited_peaks': edit_sidelobes
        }
    
    def _calculate_spectral_entropy(self, spectrum: np.ndarray) -> float:
        """Calculate normalized spectral entropy."""
        # Normalize spectrum to probability distribution
        spectrum_positive = spectrum[spectrum > 0]
        if len(spectrum_positive) == 0:
            return 0.0
        
        probabilities = spectrum_positive / np.sum(spectrum_positive)
        return entropy(probabilities, base=2)
    
    def _calculate_phase_disruption(
        self, 
        ref_peaks: List[Dict], 
        edit_peaks: List[Dict]
    ) -> float:
        """
        Calculate phase disruption between reference and edited peak sets.
        
        Measures how much the periodicity structure changed.
        """
        if len(ref_peaks) == 0 or len(edit_peaks) == 0:
            return 1.0  # Maximum disruption if no peaks in either
        
        # Extract frequencies from top peaks
        ref_freqs = np.array([p['frequency'] for p in ref_peaks[:5]])
        edit_freqs = np.array([p['frequency'] for p in edit_peaks[:5]])
        
        # Calculate minimum distance between frequency sets (normalized)
        if len(ref_freqs) == 0 or len(edit_freqs) == 0:
            return 1.0
        
        # Mean minimum distance
        disruption = 0.0
        for ref_f in ref_freqs:
            min_dist = np.min(np.abs(edit_freqs - ref_f))
            disruption += min_dist
        
        disruption /= len(ref_freqs)
        
        # Normalize to [0, 1]
        return min(disruption * 10, 1.0)
    
    def analyze_indel_disruption(
        self,
        sequence: str,
        indel_position: int,
        indel_length: int,
        indel_type: str = 'deletion'
    ) -> Dict[str, any]:
        """
        Analyze disruption caused by insertion or deletion.
        
        Args:
            sequence: Original DNA sequence
            indel_position: Position of indel
            indel_length: Length of indel
            indel_type: 'insertion' or 'deletion'
            
        Returns:
            Dictionary with indel disruption analysis
        """
        sequence = self.validate_dna_sequence(sequence)
        
        if indel_type == 'deletion':
            # Simulate deletion
            edited_seq = (
                sequence[:indel_position] + 
                sequence[indel_position + indel_length:]
            )
        elif indel_type == 'insertion':
            # Simulate insertion with random sequence
            insert_seq = 'N' * indel_length
            edited_seq = (
                sequence[:indel_position] + 
                insert_seq + 
                sequence[indel_position:]
            )
        else:
            raise ValueError(f"Unknown indel_type: {indel_type}")
        
        # Calculate disruption
        disruption = self.calculate_disruption_score(sequence, edited_seq)
        
        # Add indel-specific context
        disruption['indel_type'] = indel_type
        disruption['indel_position'] = indel_position
        disruption['indel_length'] = indel_length
        disruption['original_length'] = len(sequence)
        disruption['edited_length'] = len(edited_seq)
        
        return disruption
    
    def calculate_codon_aligned_features(
        self,
        sequence: str,
        frame: int = 0
    ) -> Dict[str, any]:
        """
        Calculate φ-structured codon-aligned disruption features.
        
        Analyzes sequence in codon triplets (3-bp windows) with golden-ratio
        weighting to detect codon-level periodicities.
        
        Args:
            sequence: DNA sequence
            frame: Reading frame (0, 1, or 2)
            
        Returns:
            Dictionary with codon-aligned features
        """
        sequence = self.validate_dna_sequence(sequence)
        
        # Adjust for reading frame
        if frame > 0:
            sequence = sequence[frame:]
        
        # Pad sequence to multiple of 3
        remainder = len(sequence) % 3
        if remainder != 0:
            sequence = sequence + 'N' * (3 - remainder)
        
        n_codons = len(sequence) // 3
        codon_values = []
        
        # Convert each codon to numeric value
        for i in range(n_codons):
            codon = sequence[i*3:(i+1)*3]
            codon_wave = self.encode_dna_complex(codon)
            # Use mean magnitude as codon value
            codon_values.append(np.mean(np.abs(codon_wave)))
        
        codon_array = np.array(codon_values)
        
        # Apply FFT to codon-level signal
        codon_fft = fft(codon_array)
        codon_magnitude = np.abs(codon_fft)
        
        # Apply golden-ratio phase weights
        # Use φ_period = 7 for codon analysis (7 codons = 21 bp)
        phi_codon = 7.0  # Geometric period in codon space
        weighted_codon_spectrum = np.array([
            codon_magnitude[n] * self.phi * ((n % phi_codon) / phi_codon) ** self.k
            for n in range(len(codon_magnitude))
        ])
        
        # Calculate codon-level entropy
        codon_entropy = self._calculate_spectral_entropy(weighted_codon_spectrum)
        
        return {
            'n_codons': n_codons,
            'frame': frame,
            'codon_values': codon_values,
            'codon_spectrum': codon_magnitude,
            'weighted_codon_spectrum': weighted_codon_spectrum,
            'codon_entropy': codon_entropy,
            'dominant_codon_period': self._find_dominant_period(
                weighted_codon_spectrum
            )
        }
    
    def _find_dominant_period(self, spectrum: np.ndarray) -> float:
        """Find dominant period from spectrum."""
        if len(spectrum) < 2:
            return 0.0
        
        # Skip DC component
        peak_idx = np.argmax(spectrum[1:len(spectrum)//2]) + 1
        
        if peak_idx == 0:
            return 0.0
        
        period = len(spectrum) / peak_idx
        return period


def calculate_grna_off_target_score(
    grna_sequence: str,
    phi_period: float = 21.0,
    k: float = 0.3
) -> Dict[str, float]:
    """
    Calculate off-target risk score for a gRNA using FFT-based golden-ratio analysis.
    
    Higher score = lower off-target risk (better guide).
    
    Args:
        grna_sequence: Guide RNA sequence (20-24 nt)
        phi_period: Geometric period (default: 21 for standard guides)
        k: Resolution exponent (default: 0.3)
        
    Returns:
        Dictionary with scoring metrics
    """
    analyzer = FFTCRISPRDisruptionAnalyzer(phi_period=phi_period, k=k)
    
    # Analyze periodicities
    periodicity_analysis = analyzer.detect_off_target_periodicities(grna_sequence)
    
    # Get peak magnitude if available
    significant_peaks = periodicity_analysis.get('significant_peaks', [])
    if len(significant_peaks) > 0:
        peak_magnitude = significant_peaks[0].get('weighted_magnitude', 1.0)
    else:
        # No significant peaks - use mean of spectrum
        weighted_spectrum = periodicity_analysis['weighted_spectrum']
        peak_magnitude = np.mean(weighted_spectrum) if len(weighted_spectrum) > 0 else 1.0
    
    # Lower entropy = more structured = higher specificity
    entropy_score = 1.0 / (1.0 + peak_magnitude)
    
    # Fewer peaks = less off-target potential
    peak_score = 1.0 / (1.0 + len(significant_peaks))
    
    # Composite score
    off_target_score = (entropy_score * 0.6 + peak_score * 0.4)
    
    return {
        'off_target_score': off_target_score,
        'entropy_component': entropy_score,
        'peak_component': peak_score,
        'n_significant_peaks': len(significant_peaks),
        'recommendation': 'good' if off_target_score > 0.6 else 'review' if off_target_score > 0.4 else 'poor'
    }
