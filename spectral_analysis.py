"""
Spectral Analysis Module for WAVE-CRISPR Framework

This module implements spectral analysis of DNA sequences using the Z Framework
discrete domain form with biological and arbitrary encoding comparisons.
"""

import numpy as np
import mpmath as mp
from typing import Dict, List, Tuple, Optional
from scipy.stats import pearsonr
import logging

# Configure high precision
mp.dps = 50

# Mathematical constants
PHI = mp.mpf("1.618033988749894848204586834365638117720309179805762862135")
E_SQUARED = mp.e ** 2

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SpectralDNAAnalyzer:
    """
    Spectral analysis of DNA sequences using Z Framework principles.
    
    Implements biological vs arbitrary encoding comparisons for empirical
    validation of spectral features in CRISPR efficiency prediction.
    """
    
    def __init__(self, k_star: float = 0.3, seed: int = 42):
        """
        Initialize spectral analyzer.
        
        Args:
            k_star: Geodesic resolution parameter (k* ≈ 0.3)
            seed: Random seed for reproducibility
        """
        self.k_star = k_star
        np.random.seed(seed)
        
        # Biological encoding based on polarizability
        self.biological_encoding = {
            'A': 1.65 + 0j,     # Adenine polarizability
            'T': 1.53 + 0j,     # Thymine polarizability  
            'C': 1.42 + 0j,     # Cytosine polarizability
            'G': 1.78 + 0j      # Guanine polarizability
        }
        
        # Generate arbitrary encoding (random complex weights)
        self.arbitrary_encoding = self._generate_arbitrary_encoding()
        
    def _generate_arbitrary_encoding(self) -> Dict[str, complex]:
        """Generate arbitrary complex encoding for comparison."""
        encoding = {}
        for base in ['A', 'T', 'C', 'G']:
            # Random amplitude [0.5, 2.0] and phase [0, 2π]
            amplitude = np.random.uniform(0.5, 2.0)
            phase = np.random.uniform(0, 2 * np.pi)
            encoding[base] = amplitude * np.exp(1j * phase)
        return encoding
    
    def encode_sequence(self, sequence: str, encoding_type: str = 'biological') -> np.ndarray:
        """
        Encode DNA sequence using specified encoding.
        
        Args:
            sequence: DNA sequence string
            encoding_type: 'biological' or 'arbitrary'
            
        Returns:
            Complex array of encoded sequence
        """
        if encoding_type == 'biological':
            encoding = self.biological_encoding
        elif encoding_type == 'arbitrary':
            encoding = self.arbitrary_encoding
        else:
            raise ValueError("encoding_type must be 'biological' or 'arbitrary'")
        
        encoded = []
        for base in sequence.upper():
            if base in encoding:
                encoded.append(encoding[base])
            else:
                # Handle N or other bases with average encoding
                avg_val = np.mean(list(encoding.values()))
                encoded.append(avg_val)
        
        return np.array(encoded)
    
    def compute_z_framework_values(self, encoded_sequence: np.ndarray) -> Dict[str, float]:
        """
        Compute Z Framework discrete domain values Z = n(Δₙ/Δₘₐₓ).
        
        Args:
            encoded_sequence: Complex encoded sequence
            
        Returns:
            Dictionary with Z Framework metrics
        """
        n = len(encoded_sequence)
        
        if n == 0:
            return {"z_value": 0.0, "delta_n": 0.0, "delta_max": 0.0, "kappa": 0.0}
        
        # Compute differences
        if n > 1:
            deltas = np.abs(np.diff(encoded_sequence))
            delta_n = np.mean(deltas)
            delta_max = np.max(deltas) if len(deltas) > 0 else 1.0
        else:
            delta_n = 0.0
            delta_max = 1.0
        
        # Avoid division by zero
        if delta_max == 0:
            delta_max = 1e-50
        
        # Z Framework calculation: Z = n(Δₙ/Δₘₐₓ)
        z_value = float(n * (delta_n / delta_max))
        
        # κ(n) = d(n) · ln(n+1) / e²
        d_n = self._divisor_count(n)
        kappa = float(d_n * mp.log(n + 1) / E_SQUARED)
        
        return {
            "z_value": z_value,
            "delta_n": float(delta_n),
            "delta_max": float(delta_max),
            "kappa": kappa
        }
    
    def _divisor_count(self, n: int) -> int:
        """Count divisors of n (simplified implementation)."""
        if n <= 0:
            return 1
        
        count = 0
        for i in range(1, int(np.sqrt(n)) + 1):
            if n % i == 0:
                count += 1
                if i != n // i:
                    count += 1
        return count
    
    def compute_spectral_entropy(self, encoded_sequence: np.ndarray) -> float:
        """
        Compute spectral entropy of encoded sequence.
        
        Args:
            encoded_sequence: Complex encoded sequence
            
        Returns:
            Spectral entropy value
        """
        if len(encoded_sequence) == 0:
            return 0.0
        
        # Compute power spectrum
        fft_vals = np.fft.fft(encoded_sequence)
        power_spectrum = np.abs(fft_vals) ** 2
        
        # Normalize to probability distribution
        ps = power_spectrum / np.sum(power_spectrum)
        ps = ps[ps > 0]  # Remove zeros
        
        if len(ps) <= 1:
            return 0.0
        
        # Compute entropy
        return float(-np.sum(ps * np.log2(ps)))
    
    def compute_geodesic_resolution(self, n: int) -> float:
        """
        Compute geodesic resolution θ'(n,k) = φ·((n mod φ)/φ)^k.
        
        Args:
            n: Sequence position/length
            
        Returns:
            Geodesic resolution value
        """
        phi_float = float(PHI)
        n_mod_phi = n % phi_float
        theta_prime = phi_float * ((n_mod_phi / phi_float) ** self.k_star)
        return float(theta_prime)
    
    def analyze_sequence(self, sequence: str, encoding_type: str = 'biological') -> Dict[str, any]:
        """
        Complete spectral analysis of DNA sequence.
        
        Args:
            sequence: DNA sequence string
            encoding_type: 'biological' or 'arbitrary'
            
        Returns:
            Comprehensive analysis results
        """
        encoded_seq = self.encode_sequence(sequence, encoding_type)
        z_metrics = self.compute_z_framework_values(encoded_seq)
        spectral_entropy = self.compute_spectral_entropy(encoded_seq)
        geodesic_res = self.compute_geodesic_resolution(len(sequence))
        
        return {
            "sequence_length": len(sequence),
            "encoding_type": encoding_type,
            "z_framework": z_metrics,
            "spectral_entropy": spectral_entropy,
            "geodesic_resolution": geodesic_res,
            "gc_content": (sequence.count('G') + sequence.count('C')) / len(sequence) if len(sequence) > 0 else 0
        }
    
    def compare_encodings(self, sequence: str) -> Dict[str, any]:
        """
        Compare biological vs arbitrary encodings for same sequence.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Comparison results
        """
        bio_analysis = self.analyze_sequence(sequence, 'biological')
        arb_analysis = self.analyze_sequence(sequence, 'arbitrary')
        
        # Compute differences
        z_diff = bio_analysis['z_framework']['z_value'] - arb_analysis['z_framework']['z_value']
        entropy_diff = bio_analysis['spectral_entropy'] - arb_analysis['spectral_entropy']
        
        return {
            "biological": bio_analysis,
            "arbitrary": arb_analysis,
            "differences": {
                "z_value_diff": z_diff,
                "spectral_entropy_diff": entropy_diff
            }
        }


def load_zeta_zeros(filepath: str = "data/zeta.txt") -> List[float]:
    """
    Load pre-computed zeta zeros for geodesic integration.
    
    Args:
        filepath: Path to zeta zeros file
        
    Returns:
        List of zeta zero values
    """
    try:
        with open(filepath, 'r') as f:
            zeros = [float(line.strip()) for line in f if line.strip()]
        logger.info(f"Loaded {len(zeros)} zeta zeros from {filepath}")
        return zeros
    except FileNotFoundError:
        logger.warning(f"Zeta zeros file not found: {filepath}")
        return []


def compute_zeta_correlations(spectral_features: List[float], 
                            zeta_zeros: List[float]) -> Tuple[float, float]:
    """
    Compute correlations between spectral features and zeta spacings.
    
    Args:
        spectral_features: List of spectral feature values
        zeta_zeros: List of zeta zero values
        
    Returns:
        Tuple of (correlation coefficient, p-value)
    """
    if len(spectral_features) == 0 or len(zeta_zeros) == 0:
        return 0.0, 1.0
    
    # Use minimum length for comparison
    min_len = min(len(spectral_features), len(zeta_zeros))
    features = spectral_features[:min_len]
    zeros = zeta_zeros[:min_len]
    
    if len(features) < 2:
        return 0.0, 1.0
    
    try:
        correlation, p_value = pearsonr(features, zeros)
        return float(correlation), float(p_value)
    except Exception as e:
        logger.error(f"Error computing correlation: {e}")
        return 0.0, 1.0


def demo_spectral_analysis():
    """Demonstration of spectral analysis functionality."""
    print("WAVE-CRISPR Spectral Analysis Demo")
    print("=" * 40)
    
    analyzer = SpectralDNAAnalyzer()
    
    # Test sequences
    test_sequences = [
        "GGGGGGGGGGGGGGGGGG",  # High GC
        "ATCGATCGATCGATCGAT",  # Alternating
        "AAAAAAAAAAAAAAAAAAA",  # Homopolymer
    ]
    
    for seq in test_sequences:
        print(f"\nSequence: {seq}")
        comparison = analyzer.compare_encodings(seq)
        
        bio = comparison['biological']
        arb = comparison['arbitrary']
        
        print(f"Biological encoding:")
        print(f"  Z-value: {bio['z_framework']['z_value']:.4f}")
        print(f"  Spectral entropy: {bio['spectral_entropy']:.4f}")
        
        print(f"Arbitrary encoding:")
        print(f"  Z-value: {arb['z_framework']['z_value']:.4f}")
        print(f"  Spectral entropy: {arb['spectral_entropy']:.4f}")


if __name__ == "__main__":
    demo_spectral_analysis()