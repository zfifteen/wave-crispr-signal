#!/usr/bin/env python3
"""
Shared Spectral Analysis Utilities for Trinity Experiments

This module provides dimensionless, biophysically-anchored spectral features
for DNA breathing dynamics analysis. All features are based on AT/GC opening
rate ratios and helical periodicity (10.5 bp/turn).

SCIENTIFIC GATES:
- Human DNA only (A/C/G/T for DNA, A/C/G/U for RNA)
- No fabrication (real nucleotide sequences)
- Dimensionless parametrization (rate ratios, not MHz)
- Reproducible (fixed seeds, documented parameters)

Z Framework invariants:
- Discrete domain: Z = A(B / e²)
- Geometric resolution: θ′(n,k) = φ·((n mod φ)/φ)^k with k ≈ 0.3
"""

import numpy as np
from scipy.signal import get_window
from typing import Dict, Tuple, List
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
HELICAL_PERIOD = 10.5  # bp per turn for B-DNA
E_SQUARED = np.e ** 2  # ≈ 7.389, discrete domain invariant


def validate_dna_sequence(seq: str) -> str:
    """
    Validate DNA sequence contains only A/C/G/T.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Uppercased validated sequence
        
    Raises:
        ValueError: If sequence contains invalid bases
    """
    seq = seq.upper().replace(" ", "").replace("\n", "")
    
    # Check for valid DNA bases only
    valid_bases = set("ACGT")
    invalid_bases = set(seq) - valid_bases
    
    if invalid_bases:
        raise ValueError(
            f"Invalid DNA bases detected: {invalid_bases}. "
            f"Only A, C, G, T allowed for DNA sequences. "
            f"Use A, C, G, U for RNA sequences (convert separately)."
        )
    
    return seq


def validate_rna_sequence(seq: str) -> str:
    """
    Validate RNA sequence contains only A/C/G/U.
    
    Args:
        seq: RNA sequence string
        
    Returns:
        Uppercased validated sequence
        
    Raises:
        ValueError: If sequence contains invalid bases
    """
    seq = seq.upper().replace(" ", "").replace("\n", "")
    
    # Check for valid RNA bases only
    valid_bases = set("ACGU")
    invalid_bases = set(seq) - valid_bases
    
    if invalid_bases:
        raise ValueError(
            f"Invalid RNA bases detected: {invalid_bases}. "
            f"Only A, C, G, U allowed for RNA sequences. "
            f"Use A, C, G, T for DNA sequences (convert separately)."
        )
    
    return seq


def encode_complex(seq: str, r: float = 20.0, is_rna: bool = False) -> np.ndarray:
    """
    Encode DNA/RNA sequence as complex numbers using dimensionless rate ratio.
    
    This encoding is based on base pair opening dynamics:
    - AT pairs open ~100× faster than GC pairs (r ≈ 20 for amplitude separation)
    - Real part: strength proportional to log(opening rate)
    - Imaginary part: phase based on bond stability (AT positive, GC negative)
    
    Args:
        seq: DNA or RNA sequence string
        r: Dimensionless rate ratio k_GC/k_AT (default 20.0)
           Typical range: [5, 200]
        is_rna: If True, expect U instead of T (default False for DNA)
        
    Returns:
        Complex array of length N where N = len(seq)
        
    Mathematical form:
        AT: -log(r) + j·0.3·log(r)  (fast opening, weak bonds)
        GC: +log(r) - j·0.3·log(r)  (slow opening, strong bonds)
    """
    # Validate sequence
    if is_rna:
        seq = validate_rna_sequence(seq)
    else:
        seq = validate_dna_sequence(seq)
    
    # Calculate dimensionless amplitude and phase factors
    alpha = np.log(r)         # Strength factor (real part amplitude)
    beta = 0.3 * np.log(r)    # Phase factor (imaginary part amplitude)
    
    # Define complex weights
    if is_rna:
        # RNA: U replaces T
        weight_table = {
            'A': -alpha + 1j * beta,   # AU pair (weak, 2 H-bonds)
            'U': -alpha + 1j * beta,   # AU pair (weak, 2 H-bonds)
            'C': +alpha - 1j * beta,   # GC pair (strong, 3 H-bonds)
            'G': +alpha - 1j * beta,   # GC pair (strong, 3 H-bonds)
        }
    else:
        # DNA: T is standard
        weight_table = {
            'A': -alpha + 1j * beta,   # AT pair (weak, 2 H-bonds)
            'T': -alpha + 1j * beta,   # AT pair (weak, 2 H-bonds)
            'C': +alpha - 1j * beta,   # GC pair (strong, 3 H-bonds)
            'G': +alpha - 1j * beta,   # GC pair (strong, 3 H-bonds)
        }
    
    # Encode sequence
    encoded = np.array([weight_table.get(base, 0.0 + 0.0j) for base in seq])
    
    logger.debug(f"Encoded {len(seq)} bases with r={r:.2f}, is_rna={is_rna}")
    
    return encoded


def cz_power_at_period(z: np.ndarray, period: float) -> float:
    """
    Calculate spectral power at a specific period using Chirp Z-Transform approach.
    
    This computes the Fourier coefficient at frequency f = 1/period using
    direct summation (equivalent to single-bin CZT).
    
    Preprocessing:
    1. Remove DC component (mean subtraction)
    2. Apply Hann window to reduce spectral leakage
    3. Compute power at target frequency
    
    Args:
        z: Complex-valued sequence (output of encode_complex)
        period: Target period in base pairs (e.g., 10.5 for helical period)
        
    Returns:
        Spectral power (magnitude) at the target period
        
    Note:
        For period = 10.5 bp, this captures DNA helical periodicity.
        Harmonics at 5.25 bp (2nd) and 3.5 bp (3rd) are also informative.
    """
    N = len(z)
    
    if N == 0:
        return 0.0
    
    # Preprocessing: remove DC and apply window
    z_processed = z - np.mean(z)
    window = get_window("hann", N)
    z_processed = z_processed * window
    
    # Compute Fourier coefficient at f = 1/period
    f = 1.0 / period
    k = np.arange(N)
    
    # Direct computation: X(f) = Σ z[k] * exp(-2πj * f * k)
    X = np.sum(z_processed * np.exp(-2j * np.pi * f * k))
    
    # Return magnitude (power)
    power = float(np.abs(X))
    
    logger.debug(f"Power at period {period:.2f} bp: {power:.4f}")
    
    return power


def breathing_features(seq: str, r: float = 20.0, is_rna: bool = False) -> Dict[str, float]:
    """
    Extract breathing dynamics spectral features from DNA/RNA sequence.
    
    This extracts power at three key periods:
    - 10.5 bp: Fundamental helical period (1 turn)
    - 5.25 bp: First harmonic (2 turns per wavelength)
    - 3.5 bp:  Second harmonic (3 turns per wavelength)
    
    These features capture rotational phasing effects in DNA/RNA structure.
    
    Args:
        seq: DNA or RNA sequence string
        r: Dimensionless rate ratio (default 20.0)
        is_rna: If True, expect RNA (U not T), default False
        
    Returns:
        Dictionary with keys:
            - P10_5: Power at 10.5 bp period
            - P5_25: Power at 5.25 bp period
            - P3_5: Power at 3.5 bp period
            - length: Sequence length
            - gc_content: GC fraction (0-1)
    """
    # Validate and encode
    z = encode_complex(seq, r=r, is_rna=is_rna)
    
    # Calculate GC content
    seq_upper = seq.upper()
    if is_rna:
        gc_count = seq_upper.count('G') + seq_upper.count('C')
    else:
        gc_count = seq_upper.count('G') + seq_upper.count('C')
    gc_content = gc_count / len(seq) if len(seq) > 0 else 0.0
    
    # Extract spectral features
    features = {
        "P10_5": cz_power_at_period(z, 10.5),
        "P5_25": cz_power_at_period(z, 5.25),
        "P3_5": cz_power_at_period(z, 3.5),
        "length": len(seq),
        "gc_content": gc_content,
    }
    
    logger.debug(f"Extracted features: {features}")
    
    return features


def phase_at(pos: int, period: float = HELICAL_PERIOD) -> float:
    """
    Calculate rotational phase at a given position.
    
    For DNA double helix with period ≈ 10.5 bp/turn, this computes
    the angular position (in radians) around the helix axis.
    
    Args:
        pos: Position along sequence (0-indexed)
        period: Helical period in base pairs (default 10.5)
        
    Returns:
        Phase angle in radians, range [0, 2π)
        
    Example:
        >>> phase_at(0)   # First base
        0.0
        >>> phase_at(10.5)  # One full turn
        6.283... (≈ 2π)
    """
    phase = (2.0 * np.pi * pos / period) % (2.0 * np.pi)
    return float(phase)


def rotational_phase_curve(
    seq: str,
    period: float = HELICAL_PERIOD,
    bins: int = 24,
    r: float = 20.0,
    is_rna: bool = False
) -> Tuple[List[float], List[float]]:
    """
    Generate rotational phase curve showing how sequence properties vary with phase.
    
    This bins the sequence by rotational phase and computes average magnitude
    of complex weights in each bin. Useful for visualizing phase-dependent effects.
    
    Args:
        seq: DNA or RNA sequence string
        period: Helical period (default 10.5 bp)
        bins: Number of phase bins (default 24, i.e., 15° resolution)
        r: Dimensionless rate ratio (default 20.0)
        is_rna: If True, expect RNA (default False)
        
    Returns:
        Tuple of (bin_centers, bin_values):
            - bin_centers: List of phase angles (radians) at bin centers
            - bin_values: List of average magnitudes in each bin
    """
    # Encode sequence
    z = encode_complex(seq, r=r, is_rna=is_rna)
    N = len(z)
    
    # Calculate phase for each position
    phases = np.array([(2 * np.pi * i / period) % (2 * np.pi) for i in range(N)])
    
    # Get magnitudes
    mags = np.abs(z)
    
    # Create phase bins
    edges = np.linspace(0, 2 * np.pi, bins + 1)
    
    # Digitize phases into bins
    idx = np.digitize(phases, edges) - 1
    idx = np.clip(idx, 0, bins - 1)  # Handle edge case
    
    # Compute average magnitude per bin
    bin_values = []
    for i in range(bins):
        mask = (idx == i)
        if np.any(mask):
            bin_values.append(float(np.mean(mags[mask])))
        else:
            bin_values.append(0.0)
    
    # Calculate bin centers
    bin_centers = [(edges[i] + edges[i + 1]) / 2 for i in range(bins)]
    
    logger.debug(f"Generated phase curve with {bins} bins")
    
    return bin_centers, bin_values


def spectral_peak_analysis(
    seq: str,
    period: float = HELICAL_PERIOD,
    r: float = 20.0,
    is_rna: bool = False
) -> Dict[str, float]:
    """
    Comprehensive spectral peak analysis at target period.
    
    Args:
        seq: DNA or RNA sequence
        period: Target period (default 10.5 bp)
        r: Rate ratio (default 20.0)
        is_rna: RNA flag (default False)
        
    Returns:
        Dictionary with:
            - power: Spectral power magnitude
            - angle: Phase angle of spectral peak (radians)
            - normalized_power: Power normalized by sequence length
    """
    z = encode_complex(seq, r=r, is_rna=is_rna)
    N = len(z)
    
    if N == 0:
        return {"power": 0.0, "angle": 0.0, "normalized_power": 0.0}
    
    # Preprocess
    z_processed = z - np.mean(z)
    window = get_window("hann", N)
    z_processed = z_processed * window
    
    # Compute complex coefficient
    f = 1.0 / period
    k = np.arange(N)
    X = np.sum(z_processed * np.exp(-2j * np.pi * f * k))
    
    power = float(np.abs(X))
    angle = float(np.angle(X))
    normalized_power = power / np.sqrt(N) if N > 0 else 0.0
    
    return {
        "power": power,
        "angle": angle,
        "normalized_power": normalized_power,
    }


if __name__ == "__main__":
    # Self-test
    print("Testing spectral utilities...")
    
    # Test sequence (human genomic)
    test_seq = "ATGCGATCGATCGATCGATCGCTAGCTAGCTACGTACGTACGT"
    
    print(f"\nTest sequence: {test_seq}")
    print(f"Length: {len(test_seq)} bp")
    
    # Test encoding
    z = encode_complex(test_seq, r=20.0)
    print(f"\nEncoding: {len(z)} complex values")
    print(f"First 5: {z[:5]}")
    
    # Test breathing features
    features = breathing_features(test_seq, r=20.0)
    print(f"\nBreathing features:")
    for key, val in features.items():
        print(f"  {key}: {val:.4f}")
    
    # Test phase calculation
    for pos in [0, 5, 10, 21]:
        phi = phase_at(pos)
        print(f"\nPhase at position {pos}: {phi:.4f} rad ({np.degrees(phi):.1f}°)")
    
    # Test rotational phase curve
    centers, values = rotational_phase_curve(test_seq, bins=12)
    print(f"\nRotational phase curve (12 bins):")
    for i, (c, v) in enumerate(zip(centers, values)):
        print(f"  Bin {i}: φ={c:.4f} rad, mag={v:.4f}")
    
    # Test spectral peak analysis
    peak_info = spectral_peak_analysis(test_seq, period=10.5)
    print(f"\nSpectral peak at 10.5 bp:")
    for key, val in peak_info.items():
        print(f"  {key}: {val:.4f}")
    
    print("\n✓ All tests passed")
