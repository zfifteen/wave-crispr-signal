#!/usr/bin/env python3
"""
DNA Breathing Dynamics Encoding for CRISPR Activity Prediction

This module implements biophysically-grounded DNA breathing dynamics encoding
based on experimentally measured base-pair opening rates and thermodynamics.

Key Features:
- Biophysical breathing weights from opening lifetimes (AT: 1-5 ms, GC: 10-50 ms)
- CZT (Chirp Z-Transform) for precise fractional-period evaluation at 10.5 bp
- Goertzel algorithm for efficient single-frequency analysis
- Rotational phasing (helical period ~10.5 bp)
- Temperature and Mg²⁺ dependence
- DC removal and windowing

References:
- Opening lifetimes: PMC5393899 (Sequence dependency of canonical base pair opening)
- Helical period: ~10.5 bp for B-DNA
- R-loop formation: PNAS 1402597111
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Any
import logging
from scipy import signal
from scipy.fft import fft
import warnings

warnings.filterwarnings('ignore')

# Physical constants
HELICAL_PERIOD_BP = 10.5  # Base pairs per helical turn for B-DNA
AVOGADRO = 6.02214076e23
GAS_CONSTANT = 8.314462618  # J/(mol·K)

# Experimentally measured base-pair opening lifetimes (milliseconds)
# Source: PMC5393899 - Sequence dependency of canonical base pair opening
BP_OPENING_LIFETIMES_MS = {
    'AT': 1.0,   # AT pairs: fast opening (weaker, 2 H-bonds)
    'TA': 1.0,   # TA pairs: fast opening
    'GC': 50.0,  # GC pairs: slow opening (stronger, 3 H-bonds)
    'CG': 50.0,  # CG pairs: slow opening
}

# Nearest-neighbor thermodynamic parameters (ΔG° kcal/mol at 37°C)
# Source: SantaLucia thermodynamics
NN_DG_37C = {
    'AA': -1.00, 'AT': -0.88, 'AG': -1.28, 'AC': -1.44,
    'TA': -0.58, 'TT': -1.00, 'TG': -1.44, 'TC': -1.30,
    'GA': -1.30, 'GT': -1.44, 'GG': -1.84, 'GC': -2.17,
    'CA': -1.45, 'CT': -1.28, 'CG': -2.24, 'CC': -1.84,
}


class BreathingDynamicsEncoder:
    """
    Encode DNA sequences using biophysically-grounded breathing dynamics.
    
    This encoder maps DNA sequences to complex-valued representations based on:
    1. Base-pair opening rates (breathing frequencies)
    2. Thermodynamic stability (nearest-neighbor ΔG°)
    3. Rotational phasing (10.5 bp helical period)
    
    The encoding is dimensionless and normalized for numerical stability.
    """
    
    def __init__(self, 
                 temperature_c: float = 37.0,
                 mg_concentration_mm: float = 2.0,
                 helical_period: float = HELICAL_PERIOD_BP,
                 seed: int = 42):
        """
        Initialize breathing dynamics encoder.
        
        Args:
            temperature_c: Temperature in Celsius (default 37°C for physiological)
            mg_concentration_mm: Mg²⁺ concentration in mM (default 2 mM)
            helical_period: Helical period in base pairs (default 10.5 bp for B-DNA)
            seed: Random seed for reproducibility
        """
        self.temperature_c = temperature_c
        self.temperature_k = temperature_c + 273.15
        self.mg_concentration_mm = mg_concentration_mm
        self.helical_period = helical_period
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        
        # Create complex encoding based on breathing dynamics
        self.base_weights = self._create_breathing_weights()
        
        self.logger.info(f"BreathingDynamicsEncoder initialized: T={temperature_c}°C, Mg²⁺={mg_concentration_mm}mM")
    
    def _create_breathing_weights(self) -> Dict[str, complex]:
        """
        Create complex weights from breathing dynamics and thermodynamics.
        
        The weights encode both kinetic (opening rates) and thermodynamic (stability)
        information in a dimensionless, normalized form.
        
        Returns:
            Dictionary mapping bases to complex weights
        """
        weights = {}
        
        # Normalize opening lifetimes to dimensionless ratios
        # Use AT as reference (fastest opening)
        at_lifetime = BP_OPENING_LIFETIMES_MS['AT']
        
        for base in 'ATCGN':
            if base == 'N':
                # Unknown base: average of A and C
                real_part = 0.5 * (self._get_real_component('A') + self._get_real_component('C'))
                imag_part = 0.5 * (self._get_imag_component('A') + self._get_imag_component('C'))
            else:
                real_part = self._get_real_component(base)
                imag_part = self._get_imag_component(base)
            
            weights[base] = complex(real_part, imag_part)
        
        return weights
    
    def _get_real_component(self, base: str) -> float:
        """
        Calculate real component from breathing kinetics.
        
        The real component encodes the opening rate ratio (AT as reference).
        Faster opening (AT) → larger positive values
        Slower opening (GC) → smaller positive values
        
        Args:
            base: DNA base (A, T, C, G)
            
        Returns:
            Normalized real component
        """
        if base in 'AT':
            # Fast opening bases (AT pairs)
            # Use dimensionless ratio: log(reference_time / opening_time)
            pair = 'AT'
            lifetime = BP_OPENING_LIFETIMES_MS[pair]
            # Normalize to range suitable for spectral analysis
            real_val = np.log10(BP_OPENING_LIFETIMES_MS['GC'] / lifetime) / 2.0
        elif base in 'CG':
            # Slow opening bases (GC pairs)
            pair = 'GC'
            lifetime = BP_OPENING_LIFETIMES_MS[pair]
            real_val = np.log10(BP_OPENING_LIFETIMES_MS['GC'] / lifetime) / 2.0
        else:
            real_val = 0.5  # Neutral value for unknown
        
        return real_val
    
    def _get_imag_component(self, base: str) -> float:
        """
        Calculate imaginary component from thermodynamic stability.
        
        The imaginary component encodes base-pair stability and stacking energy.
        Uses nearest-neighbor thermodynamic parameters.
        
        Args:
            base: DNA base (A, T, C, G)
            
        Returns:
            Normalized imaginary component
        """
        # Temperature correction factor for thermodynamics
        # ΔG varies with temperature: ΔG(T) ≈ ΔG(37°C) × (T_K / 310.15)
        temp_factor = self.temperature_k / 310.15
        
        # Mg²⁺ concentration affects stability (higher Mg²⁺ → more stable)
        mg_factor = 1.0 + 0.1 * np.log10(max(self.mg_concentration_mm, 0.1) / 2.0)
        
        if base == 'A' or base == 'T':
            # AT-rich: less stable, more positive ΔG (less favorable)
            # Use average ΔG for AT-containing pairs
            avg_dg = (NN_DG_37C['AA'] + NN_DG_37C['AT'] + NN_DG_37C['TA'] + NN_DG_37C['TT']) / 4.0
            imag_val = -avg_dg * temp_factor * mg_factor / 2.0  # Normalize
        elif base == 'C' or base == 'G':
            # GC-rich: more stable, more negative ΔG (more favorable)
            # Use average ΔG for GC-containing pairs
            avg_dg = (NN_DG_37C['GG'] + NN_DG_37C['GC'] + NN_DG_37C['CG'] + NN_DG_37C['CC']) / 4.0
            imag_val = -avg_dg * temp_factor * mg_factor / 2.0  # Normalize
        else:
            imag_val = 0.5  # Neutral value for unknown
        
        return imag_val
    
    def encode_sequence(self, sequence: str, apply_helical_phase: bool = True) -> np.ndarray:
        """
        Encode DNA sequence with breathing dynamics.
        
        Args:
            sequence: DNA sequence string (A/C/G/T/N)
            apply_helical_phase: Whether to apply helical periodicity modulation
            
        Returns:
            Complex-valued array with breathing dynamics encoding
        """
        # Validate sequence
        sequence = sequence.upper().strip()
        valid_bases = set('ATCGN')
        
        if not all(base in valid_bases for base in sequence):
            invalid = set(sequence) - valid_bases
            raise ValueError(f"Invalid bases in sequence: {invalid}. Only A/C/G/T/N allowed.")
        
        encoded = []
        
        for i, base in enumerate(sequence):
            # Get base weight from breathing dynamics
            base_weight = self.base_weights.get(base, self.base_weights['N'])
            
            if apply_helical_phase:
                # Apply rotational phasing from helical structure
                # DNA wraps every ~10.5 bp, creating rotational phase modulation
                helical_phase = 2 * np.pi * i / self.helical_period
                
                # Also add weak positional phase for context
                positional_phase = 2 * np.pi * i / len(sequence) * 0.1
                
                # Combined phase
                total_phase = helical_phase + positional_phase
                
                # Apply phase modulation
                encoded_base = base_weight * np.exp(1j * total_phase)
            else:
                encoded_base = base_weight
            
            encoded.append(encoded_base)
        
        return np.array(encoded, dtype=complex)
    
    def get_encoding_info(self) -> Dict[str, Any]:
        """Get encoding information for documentation."""
        return {
            'type': 'biophysical_breathing_dynamics',
            'temperature_c': self.temperature_c,
            'mg_concentration_mm': self.mg_concentration_mm,
            'helical_period_bp': self.helical_period,
            'opening_lifetimes_ms': BP_OPENING_LIFETIMES_MS,
            'weights': {base: f"{val.real:.4f}{val.imag:+.4f}j" 
                       for base, val in self.base_weights.items()},
            'reference': 'PMC5393899'
        }


class ChirpZTransform:
    """
    Chirp Z-Transform (CZT) for precise fractional-period frequency evaluation.
    
    CZT allows evaluation of the Z-transform at arbitrary points on the z-plane,
    making it ideal for analyzing DNA sequences at specific biological frequencies
    like the helical period (1/10.5 bp⁻¹) that may not align with FFT bins.
    """
    
    def __init__(self, sequence_length: int):
        """
        Initialize CZT for a given sequence length.
        
        Args:
            sequence_length: Length of input sequence
        """
        self.N = sequence_length
        self.logger = logging.getLogger(__name__)
    
    def compute_czt(self, 
                   signal: np.ndarray, 
                   frequency_hz: float,
                   sampling_rate_hz: float = 1.0,
                   num_points: int = 1) -> np.ndarray:
        """
        Compute Chirp Z-Transform at specific frequency.
        
        Args:
            signal: Input signal (complex or real)
            frequency_hz: Target frequency to evaluate
            sampling_rate_hz: Sampling rate (default 1.0 for normalized frequency)
            num_points: Number of frequency points to evaluate
            
        Returns:
            CZT evaluated at specified frequency
        """
        # Convert frequency to angular frequency
        w = 2 * np.pi * frequency_hz / sampling_rate_hz
        
        # Define spiral contour in z-plane
        # For a specific frequency, we want to evaluate on the unit circle
        A = np.exp(1j * w * 0)  # Starting point (can be adjusted)
        W = np.exp(-1j * w)     # Spacing between points
        
        # CZT computation
        n = np.arange(self.N)
        k = np.arange(num_points)
        
        # Compute chirp factors
        chirp_n = W ** (n**2 / 2.0)
        chirp_k = W ** (-(k**2) / 2.0)
        
        # Apply chirp to signal
        signal_chirped = signal * chirp_n * (A ** -n)
        
        # Convolve with chirp
        # Use FFT for efficient convolution
        L = len(signal_chirped) + num_points - 1
        L_fft = 2 ** int(np.ceil(np.log2(L)))
        
        signal_fft = np.fft.fft(signal_chirped, L_fft)
        chirp_conv_fft = np.fft.fft(W ** (-(np.arange(L_fft)**2) / 2.0), L_fft)
        
        # Inverse FFT
        result = np.fft.ifft(signal_fft * chirp_conv_fft)
        result = result[:num_points] * chirp_k
        
        return result
    
    def evaluate_at_fractional_period(self,
                                     signal: np.ndarray,
                                     period_bp: float = HELICAL_PERIOD_BP,
                                     harmonics: int = 3) -> Dict[str, float]:
        """
        Evaluate CZT at fractional period (e.g., 10.5 bp) and harmonics.
        
        Args:
            signal: Complex-encoded DNA sequence
            period_bp: Target period in base pairs (default 10.5 bp)
            harmonics: Number of harmonics to include (default 3)
            
        Returns:
            Dictionary with power at fundamental and harmonics
        """
        # Fundamental frequency: 1/period (cycles per base pair)
        fundamental_freq = 1.0 / period_bp
        
        results = {}
        
        for h in range(1, harmonics + 1):
            # Harmonic frequency
            freq = h * fundamental_freq
            
            # Compute CZT at this frequency
            czt_value = self.compute_czt(signal, freq, sampling_rate_hz=1.0, num_points=1)[0]
            
            # Extract power
            power = np.abs(czt_value) ** 2
            phase = np.angle(czt_value)
            
            results[f'czt_period_{period_bp:.1f}_h{h}_power'] = power
            results[f'czt_period_{period_bp:.1f}_h{h}_phase'] = phase
        
        # Total power across all harmonics
        total_power = sum(results[k] for k in results if '_power' in k)
        results[f'czt_period_{period_bp:.1f}_total_power'] = total_power
        
        return results


class GoertzelAlgorithm:
    """
    Goertzel algorithm for efficient single-frequency DFT evaluation.
    
    More efficient than FFT when only analyzing a small number of frequencies,
    making it ideal for focused analysis of helical periodicity in DNA.
    """
    
    def __init__(self):
        """Initialize Goertzel algorithm."""
        self.logger = logging.getLogger(__name__)
    
    def compute_goertzel(self,
                        signal: np.ndarray,
                        frequency_hz: float,
                        sampling_rate_hz: float = 1.0) -> complex:
        """
        Compute Goertzel algorithm at specific frequency.
        
        Args:
            signal: Input signal (complex or real)
            frequency_hz: Target frequency to evaluate
            sampling_rate_hz: Sampling rate (default 1.0 for normalized)
            
        Returns:
            Complex value at specified frequency
        """
        N = len(signal)
        
        # Normalized frequency
        k = int(0.5 + (N * frequency_hz / sampling_rate_hz))
        w = (2.0 * np.pi * k) / N
        
        # Goertzel coefficients
        cosine = np.cos(w)
        coeff = 2.0 * cosine
        
        # Iterate through signal
        s0 = 0.0
        s1 = 0.0
        s2 = 0.0
        
        for sample in signal:
            s0 = sample + coeff * s1 - s2
            s2 = s1
            s1 = s0
        
        # Compute DFT value
        real_part = s1 - s2 * cosine
        imag_part = s2 * np.sin(w)
        
        return complex(real_part, imag_part)
    
    def evaluate_at_period(self,
                          signal: np.ndarray,
                          period_bp: float = HELICAL_PERIOD_BP,
                          harmonics: int = 3) -> Dict[str, float]:
        """
        Evaluate Goertzel at specific period and harmonics.
        
        Args:
            signal: Complex-encoded DNA sequence
            period_bp: Target period in base pairs
            harmonics: Number of harmonics to include
            
        Returns:
            Dictionary with power at fundamental and harmonics
        """
        fundamental_freq = 1.0 / period_bp
        
        results = {}
        
        for h in range(1, harmonics + 1):
            freq = h * fundamental_freq
            
            # Compute Goertzel
            goertzel_value = self.compute_goertzel(signal, freq, sampling_rate_hz=1.0)
            
            # Extract power and phase
            power = np.abs(goertzel_value) ** 2
            phase = np.angle(goertzel_value)
            
            results[f'goertzel_period_{period_bp:.1f}_h{h}_power'] = power
            results[f'goertzel_period_{period_bp:.1f}_h{h}_phase'] = phase
        
        # Total power
        total_power = sum(results[k] for k in results if '_power' in k)
        results[f'goertzel_period_{period_bp:.1f}_total_power'] = total_power
        
        return results


class BreathingSpectralAnalyzer:
    """
    Complete spectral analyzer for DNA breathing dynamics.
    
    Combines breathing encoding with CZT/Goertzel for precise frequency analysis.
    Includes windowing, DC removal, and phase-aware feature extraction.
    """
    
    def __init__(self,
                 temperature_c: float = 37.0,
                 mg_concentration_mm: float = 2.0,
                 helical_period: float = HELICAL_PERIOD_BP,
                 use_czt: bool = True,
                 seed: int = 42):
        """
        Initialize breathing spectral analyzer.
        
        Args:
            temperature_c: Temperature in Celsius
            mg_concentration_mm: Mg²⁺ concentration in mM
            helical_period: Helical period in base pairs
            use_czt: Use CZT (True) or Goertzel (False)
            seed: Random seed
        """
        self.encoder = BreathingDynamicsEncoder(
            temperature_c=temperature_c,
            mg_concentration_mm=mg_concentration_mm,
            helical_period=helical_period,
            seed=seed
        )
        self.use_czt = use_czt
        self.helical_period = helical_period
        self.logger = logging.getLogger(__name__)
    
    def remove_dc(self, signal: np.ndarray) -> np.ndarray:
        """
        Remove DC component from signal.
        
        Args:
            signal: Complex signal
            
        Returns:
            Signal with DC removed
        """
        return signal - np.mean(signal)
    
    def apply_window(self, input_signal: np.ndarray, window_type: str = 'hamming') -> np.ndarray:
        """
        Apply window function to signal.
        
        Args:
            input_signal: Complex signal
            window_type: Window type ('hamming', 'hann', 'blackman')
            
        Returns:
            Windowed signal
        """
        N = len(input_signal)
        
        if window_type == 'hamming':
            window = np.hamming(N)
        elif window_type == 'hann':
            window = np.hanning(N)
        elif window_type == 'blackman':
            window = np.blackman(N)
        else:
            window = np.ones(N)
        
        # Apply window to both real and imaginary parts
        return input_signal * window
    
    def extract_breathing_features(self,
                                   sequence: str,
                                   window_type: str = 'hamming',
                                   harmonics: int = 3) -> Dict[str, float]:
        """
        Extract complete breathing dynamics spectral features.
        
        Args:
            sequence: DNA sequence
            window_type: Window function type
            harmonics: Number of harmonics to analyze
            
        Returns:
            Dictionary of spectral features
        """
        # Encode with breathing dynamics
        encoded = self.encoder.encode_sequence(sequence, apply_helical_phase=True)
        
        # Remove DC component
        encoded_dc_removed = self.remove_dc(encoded)
        
        # Apply window
        encoded_windowed = self.apply_window(encoded_dc_removed, window_type)
        
        features = {}
        
        # Fractional-period analysis using CZT or Goertzel
        if self.use_czt:
            czt = ChirpZTransform(len(encoded_windowed))
            period_features = czt.evaluate_at_fractional_period(
                encoded_windowed,
                period_bp=self.helical_period,
                harmonics=harmonics
            )
        else:
            goertzel = GoertzelAlgorithm()
            period_features = goertzel.evaluate_at_period(
                encoded_windowed,
                period_bp=self.helical_period,
                harmonics=harmonics
            )
        
        features.update(period_features)
        
        # Standard FFT for comparison
        spectrum = np.abs(fft(encoded_windowed))
        features['fft_peak_power'] = np.max(spectrum)
        features['fft_mean_power'] = np.mean(spectrum)
        features['fft_total_power'] = np.sum(spectrum)
        
        # Breathing-specific metrics
        features['breathing_gc_content'] = (sequence.count('G') + sequence.count('C')) / len(sequence)
        features['breathing_at_content'] = (sequence.count('A') + sequence.count('T')) / len(sequence)
        
        # Phase coherence (measure of phase alignment across sequence)
        phases = np.angle(encoded_windowed)
        phase_coherence = np.abs(np.mean(np.exp(1j * phases)))
        features['breathing_phase_coherence'] = phase_coherence
        
        # Breathing amplitude variance (indicates regions of variable stability)
        amplitudes = np.abs(encoded_windowed)
        features['breathing_amplitude_var'] = np.var(amplitudes)
        features['breathing_amplitude_mean'] = np.mean(amplitudes)
        
        return features


if __name__ == "__main__":
    # Test breathing dynamics encoding
    import sys
    
    print("="*70)
    print("DNA BREATHING DYNAMICS ENCODING TEST")
    print("="*70)
    
    # Create encoder
    encoder = BreathingDynamicsEncoder(temperature_c=37.0, mg_concentration_mm=2.0)
    
    # Test sequences
    test_sequences = {
        'AT-rich': 'AAATTTAAATTTAAATTTAT',
        'GC-rich': 'GGGCCCGGGCCCGGGCCCGC',
        'Mixed': 'ATCGATCGATCGATCGATCG',
    }
    
    print("\nEncoding Information:")
    info = encoder.get_encoding_info()
    for key, value in info.items():
        if key != 'weights':
            print(f"  {key}: {value}")
    
    print("\nBase Weights:")
    for base, weight in info['weights'].items():
        print(f"  {base}: {weight}")
    
    print("\n" + "="*70)
    print("SEQUENCE ENCODING TESTS")
    print("="*70)
    
    for name, seq in test_sequences.items():
        print(f"\n{name} sequence: {seq}")
        
        # Encode
        encoded = encoder.encode_sequence(seq)
        
        print(f"  Length: {len(encoded)}")
        print(f"  Mean magnitude: {np.mean(np.abs(encoded)):.4f}")
        print(f"  Mean phase: {np.mean(np.angle(encoded)):.4f} rad")
    
    print("\n" + "="*70)
    print("FRACTIONAL PERIOD ANALYSIS")
    print("="*70)
    
    # Test CZT and Goertzel
    analyzer_czt = BreathingSpectralAnalyzer(use_czt=True)
    analyzer_goertzel = BreathingSpectralAnalyzer(use_czt=False)
    
    test_seq = 'ATCGATCGATCGATCGATCGATCG'  # 24 bp
    
    print(f"\nTest sequence: {test_seq}")
    print(f"Length: {len(test_seq)} bp")
    
    print("\nCZT Analysis:")
    features_czt = analyzer_czt.extract_breathing_features(test_seq, harmonics=3)
    for key in sorted(features_czt.keys()):
        if 'czt' in key:
            print(f"  {key}: {features_czt[key]:.6f}")
    
    print("\nGoertzel Analysis:")
    features_goertzel = analyzer_goertzel.extract_breathing_features(test_seq, harmonics=3)
    for key in sorted(features_goertzel.keys()):
        if 'goertzel' in key:
            print(f"  {key}: {features_goertzel[key]:.6f}")
    
    print("\n" + "="*70)
    print("✓ Breathing dynamics encoding test complete")
    print("="*70)
