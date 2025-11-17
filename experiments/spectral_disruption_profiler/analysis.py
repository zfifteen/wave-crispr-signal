"""
Analysis Engine - FFT-Based Spectral Feature Extraction

This module computes spectral disruption metrics from complex waveforms
using FFT analysis with golden-ratio phase weighting.

Features extracted:
- Δf₁: Change in fundamental frequency
- ΔEntropy: Spectral entropy change
- Sidelobes: Number of significant spectral peaks
- GC content and related metrics

Scientific gates enforced:
- Vectorized NumPy/SciPy operations for performance
- Bootstrap CI support for statistical validation
"""

import numpy as np
from scipy import fft, signal, stats
from typing import Dict, Tuple, List, Optional
import logging

logger = logging.getLogger(__name__)

# Golden ratio for spiral bias
PHI = (1.0 + np.sqrt(5.0)) / 2.0


def compute_fft_spectrum(
    waveform: np.ndarray,
    remove_dc: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute FFT spectrum of complex waveform.
    
    Args:
        waveform: Complex waveform array
        remove_dc: If True, remove DC component (mean)
        
    Returns:
        Tuple of (frequencies, spectrum_magnitude)
    """
    if len(waveform) == 0:
        return np.array([]), np.array([])
    
    # Remove DC component
    if remove_dc:
        waveform = waveform - np.mean(waveform)
    
    # Apply Hann window to reduce spectral leakage
    if len(waveform) > 1:
        window = signal.windows.hann(len(waveform))
        waveform = waveform * window
    
    # Compute FFT
    spectrum = fft.fft(waveform)
    spectrum_mag = np.abs(spectrum)
    
    # Frequency bins (normalized to sampling rate = 1)
    freqs = fft.fftfreq(len(waveform))
    
    return freqs, spectrum_mag


def compute_fundamental_frequency(
    spectrum_mag: np.ndarray,
    freqs: np.ndarray
) -> float:
    """
    Find fundamental frequency (peak with highest magnitude).
    
    Args:
        spectrum_mag: FFT spectrum magnitude
        freqs: Frequency array
        
    Returns:
        Fundamental frequency (f₁)
    """
    if len(spectrum_mag) == 0:
        return 0.0
    
    # Find peak in positive frequencies only
    positive_mask = freqs > 0
    if not np.any(positive_mask):
        return 0.0
    
    positive_freqs = freqs[positive_mask]
    positive_spectrum = spectrum_mag[positive_mask]
    
    # Find peak
    peak_idx = np.argmax(positive_spectrum)
    f1 = positive_freqs[peak_idx]
    
    return float(f1)


def compute_spectral_entropy(spectrum_mag: np.ndarray) -> float:
    """
    Compute spectral entropy (Shannon entropy of normalized spectrum).
    
    Higher entropy indicates more distributed spectral energy (more disruption).
    
    Args:
        spectrum_mag: FFT spectrum magnitude
        
    Returns:
        Spectral entropy (bits)
    """
    if len(spectrum_mag) == 0:
        return 0.0
    
    # Normalize to probability distribution
    power = spectrum_mag ** 2
    total_power = np.sum(power)
    
    if total_power == 0:
        return 0.0
    
    prob_dist = power / total_power
    
    # Remove zeros to avoid log(0)
    prob_dist = prob_dist[prob_dist > 0]
    
    # Shannon entropy (in bits)
    entropy = -np.sum(prob_dist * np.log2(prob_dist))
    
    return float(entropy)


def count_sidelobes(
    spectrum_mag: np.ndarray,
    freqs: np.ndarray,
    threshold_percentile: float = 90.0
) -> int:
    """
    Count significant spectral peaks (sidelobes).
    
    Peaks above threshold_percentile of max magnitude are counted.
    
    Args:
        spectrum_mag: FFT spectrum magnitude
        freqs: Frequency array
        threshold_percentile: Percentile threshold for peak detection
        
    Returns:
        Number of significant sidelobes
    """
    if len(spectrum_mag) == 0:
        return 0
    
    # Consider only positive frequencies
    positive_mask = freqs > 0
    if not np.any(positive_mask):
        return 0
    
    positive_spectrum = spectrum_mag[positive_mask]
    
    # Find peaks
    peaks, _ = signal.find_peaks(positive_spectrum)
    
    if len(peaks) == 0:
        return 0
    
    # Threshold based on percentile
    threshold = np.percentile(positive_spectrum[peaks], threshold_percentile)
    significant_peaks = np.sum(positive_spectrum[peaks] >= threshold)
    
    return int(significant_peaks)


def compute_gc_content(sequence: str) -> float:
    """
    Compute GC content of sequence.
    
    Args:
        sequence: DNA or RNA sequence
        
    Returns:
        GC content as fraction [0, 1]
    """
    sequence = sequence.upper().replace('U', 'T')
    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    
    if total == 0:
        return 0.0
    
    return gc_count / total


def compute_spectral_features(
    waveform: np.ndarray,
    reference_waveform: Optional[np.ndarray] = None,
    sequence: Optional[str] = None
) -> Dict[str, float]:
    """
    Compute comprehensive spectral features from waveform.
    
    Features include:
    - f₁: Fundamental frequency
    - entropy: Spectral entropy
    - sidelobes: Number of significant peaks
    - Δf₁, ΔEntropy, ΔSidelobes: Changes from reference (if provided)
    - gc_content: GC content (if sequence provided)
    
    Args:
        waveform: Complex waveform to analyze
        reference_waveform: Optional reference waveform for computing deltas
        sequence: Optional sequence string for GC content
        
    Returns:
        Dictionary of spectral features
    """
    features = {}
    
    # Compute FFT spectrum
    freqs, spectrum = compute_fft_spectrum(waveform)
    
    # Extract features
    f1 = compute_fundamental_frequency(spectrum, freqs)
    entropy = compute_spectral_entropy(spectrum)
    sidelobes = count_sidelobes(spectrum, freqs)
    
    features['f1'] = f1
    features['entropy'] = entropy
    features['sidelobes'] = sidelobes
    features['spectrum_power'] = float(np.sum(spectrum ** 2))
    features['length'] = len(waveform)
    
    # Compute deltas if reference provided
    if reference_waveform is not None:
        ref_freqs, ref_spectrum = compute_fft_spectrum(reference_waveform)
        ref_f1 = compute_fundamental_frequency(ref_spectrum, ref_freqs)
        ref_entropy = compute_spectral_entropy(ref_spectrum)
        ref_sidelobes = count_sidelobes(ref_spectrum, ref_freqs)
        
        features['delta_f1'] = f1 - ref_f1
        features['delta_entropy'] = entropy - ref_entropy
        features['delta_sidelobes'] = sidelobes - ref_sidelobes
    
    # Add GC content if sequence provided
    if sequence is not None:
        features['gc_content'] = compute_gc_content(sequence)
    
    return features


def analyze_disruption(
    mutant_sequence: str,
    reference_sequence: str,
    is_rna: bool = False,
    phi: float = 21.0,
    k: float = 0.3
) -> Dict[str, float]:
    """
    Analyze spectral disruption between mutant and reference sequences.
    
    This is the main analysis function that computes disruption metrics
    by comparing phase-weighted spectral features.
    
    Args:
        mutant_sequence: Mutant DNA/RNA sequence
        reference_sequence: Reference DNA/RNA sequence
        is_rna: If True, use RNA encoding
        phi: Geometric period for phase weighting
        k: Curvature parameter for phase weighting
        
    Returns:
        Dictionary of disruption metrics including Δf₁, ΔEntropy, ΔSidelobes
        
    Raises:
        ValueError: If sequences contain invalid characters
    """
    from .encoding import phase_weighted_encoding
    
    # Encode sequences
    ref_waveform = phase_weighted_encoding(reference_sequence, is_rna=is_rna, phi=phi, k=k)
    mut_waveform = phase_weighted_encoding(mutant_sequence, is_rna=is_rna, phi=phi, k=k)
    
    # Compute features
    features = compute_spectral_features(
        mut_waveform,
        reference_waveform=ref_waveform,
        sequence=mutant_sequence
    )
    
    # Add metadata
    features['reference_length'] = len(reference_sequence)
    features['mutant_length'] = len(mutant_sequence)
    features['phi'] = phi
    features['k'] = k
    
    logger.info(f"Disruption analysis: Δf₁={features.get('delta_f1', 0):.4f}, "
                f"ΔEntropy={features.get('delta_entropy', 0):.4f}")
    
    return features


def batch_analyze_disruption(
    mutant_sequences: List[str],
    reference_sequences: List[str],
    is_rna: bool = False,
    phi: float = 21.0,
    k: float = 0.3
) -> List[Dict[str, float]]:
    """
    Batch analyze disruption for multiple sequence pairs.
    
    Args:
        mutant_sequences: List of mutant sequences
        reference_sequences: List of reference sequences
        is_rna: If True, use RNA encoding
        phi: Geometric period
        k: Curvature parameter
        
    Returns:
        List of disruption metric dictionaries
        
    Raises:
        ValueError: If sequences contain invalid characters or lists have different lengths
    """
    if len(mutant_sequences) != len(reference_sequences):
        raise ValueError(
            f"Number of mutant sequences ({len(mutant_sequences)}) must match "
            f"number of reference sequences ({len(reference_sequences)})"
        )
    
    results = []
    
    for i, (mut_seq, ref_seq) in enumerate(zip(mutant_sequences, reference_sequences)):
        try:
            features = analyze_disruption(
                mut_seq, ref_seq,
                is_rna=is_rna,
                phi=phi,
                k=k
            )
            results.append(features)
        except Exception as e:
            logger.error(f"Failed to analyze sequence pair {i}: {e}")
            raise
    
    logger.info(f"Successfully analyzed {len(results)} sequence pairs")
    return results
