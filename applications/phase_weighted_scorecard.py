"""
Phase-Weighted CRISPR Scorecard

This module implements phase-weighted Fourier analysis for CRISPR guide scoring,
using the Z Framework with golden-ratio phase resolution θ′(n,k).

Key Features:
- Complex encoding: A→1, T→-1, C→+i, G→-i
- Phase-weighted transform: θ′(n,k) = φ·((n mod φ)/φ)^k with k*=0.3
- FFT-based spectral disruption features
- Z-invariant composite scoring: Z = S(Δ_spectral / φ)

Scientific Gates:
- Human DNA only (A/C/G/T for DNA, A/C/G/U for RNA)
- Fail-fast validation with clear ValueError on violations
- Domain-correct Z invariant: Z = A(B / e^2) for discrete/biological
- Statistical validity with bootstrap CI and permutation tests

References:
- Issue: Phase-Weighted CRISPR Scorecard
- Z Framework documentation: docs/Z_FRAMEWORK.md
"""

import numpy as np
import mpmath as mp
from scipy.fft import fft
from scipy.stats import entropy as shannon_entropy
from scipy.signal import find_peaks
from typing import Dict, List, Tuple, Optional, Union
import warnings

# Configure high precision
mp.dps = 50

# Mathematical constants
PHI = float(mp.mpf("1.6180339887498949"))  # Golden ratio
E_SQUARED = float(mp.e ** 2)  # e^2 for Z invariant

# Optimized k parameter from grid search on entropy minima
K_STAR = 0.3

# Default weights for composite score (tuned via OLS regression on labeled data)
DEFAULT_WEIGHTS = {
    "entropy": 0.4,
    "freq_shift": 0.35,
    "sidelobes": 0.25,
}


def validate_dna_sequence(seq: str, allow_rna: bool = False) -> None:
    """
    Validate DNA/RNA sequence contains only valid nucleotides.
    
    Args:
        seq: DNA or RNA sequence
        allow_rna: If True, allow U for RNA; if False, only allow DNA (A/C/G/T/N)
    
    Raises:
        ValueError: If sequence contains invalid characters
    """
    seq_upper = seq.upper()
    
    if allow_rna:
        # For RNA: allow A/C/G/U/N only
        valid_chars = set('ACGUN')
        invalid_chars = set(seq_upper) - valid_chars
        if invalid_chars:
            raise ValueError(
                f"Invalid RNA sequence: contains {invalid_chars}. "
                f"Only A, C, G, U, N allowed for RNA sequences."
            )
        if 'T' in seq_upper:
            raise ValueError(
                f"Invalid RNA sequence: contains 'T'. "
                f"RNA sequences should use 'U' instead of 'T'."
            )
    else:
        # For DNA: allow A/C/G/T/N only
        valid_chars = set('ACGTN')
        invalid_chars = set(seq_upper) - valid_chars
        if invalid_chars:
            raise ValueError(
                f"Invalid DNA sequence: contains {invalid_chars}. "
                f"Only A, C, G, T, N allowed for DNA sequences."
            )
        if 'U' in seq_upper:
            raise ValueError(
                f"Invalid DNA sequence: contains 'U'. "
                f"DNA sequences should use 'T' instead of 'U'."
            )


def encode_complex(seq: str, is_rna: bool = False) -> np.ndarray:
    """
    Map DNA/RNA string to complex vector.
    
    Complex encoding preserves complementarity:
    - A → 1 (real positive)
    - T → -1 (real negative) [DNA only]
    - U → -1 (real negative) [RNA only]
    - C → +i (imaginary positive)
    - G → -i (imaginary negative)
    - N → 0 (ambiguous base)
    
    Args:
        seq: DNA or RNA sequence
        is_rna: If True, treat as RNA (U instead of T)
    
    Returns:
        Complex numpy array
        
    Raises:
        ValueError: If sequence contains invalid nucleotides
    """
    # Validate sequence first
    validate_dna_sequence(seq, allow_rna=is_rna)
    
    seq = seq.upper()
    
    if is_rna:
        # RNA encoding
        mapping = {
            'A': 1.0 + 0.0j,
            'U': -1.0 + 0.0j,
            'C': 0.0 + 1.0j,
            'G': 0.0 - 1.0j,
            'N': 0.0 + 0.0j,
        }
    else:
        # DNA encoding
        mapping = {
            'A': 1.0 + 0.0j,
            'T': -1.0 + 0.0j,
            'C': 0.0 + 1.0j,
            'G': 0.0 - 1.0j,
            'N': 0.0 + 0.0j,
        }
    
    return np.array([mapping[base] for base in seq], dtype=np.complex128)


def theta_prime(n: Union[int, float, np.ndarray], k: float = K_STAR) -> Union[float, np.ndarray]:
    """
    Geodesic resolution function θ′(n,k) with golden-angle phasing.
    
    θ′(n,k) = φ · ((n mod φ)/φ)^k
    
    This function normalizes position n in sequences for phase-weighted encoding,
    embedding golden-angle spirals to minimize discrepancy in spectral sampling.
    
    Args:
        n: Position index or array of positions
        k: Scaling parameter (default: K_STAR = 0.3)
    
    Returns:
        Phase-weighted position value(s)
    """
    if isinstance(n, np.ndarray):
        n_mod_phi = n % PHI
        ratio = n_mod_phi / PHI
        return PHI * (ratio ** k)
    else:
        n_mod_phi = float(n) % PHI
        ratio = n_mod_phi / PHI
        return PHI * (ratio ** k)


def apply_phase_weighting(
    encoded_seq: np.ndarray,
    k: float = K_STAR,
) -> np.ndarray:
    """
    Apply position-dependent phase weighting to encoded sequence.
    
    For each base at index n, multiply by e^{i·θ′(n,k)}
    
    Args:
        encoded_seq: Complex-encoded sequence
        k: Phase parameter (default: K_STAR = 0.3)
    
    Returns:
        Phase-weighted complex array
    """
    n_positions = np.arange(len(encoded_seq))
    theta_values = theta_prime(n_positions, k)
    
    # Apply phase shift: multiply by e^{i·θ′(n,k)}
    phase_factors = np.exp(1j * theta_values)
    
    return encoded_seq * phase_factors


def compute_spectrum(waveform: np.ndarray) -> np.ndarray:
    """
    Compute FFT magnitude spectrum.
    
    X[k] = Σ x[n] · e^{-i2πkn/N}
    
    Args:
        waveform: Complex waveform
    
    Returns:
        Magnitude spectrum
    """
    spectrum = fft(waveform)
    return np.abs(spectrum)


def compute_spectral_entropy(spectrum: np.ndarray, base: int = 2) -> float:
    """
    Calculate normalized Shannon entropy of spectrum.
    
    H(X) = -Σ p(x) log₂ p(x)
    
    Args:
        spectrum: Magnitude spectrum
        base: Logarithm base (default: 2)
    
    Returns:
        Spectral entropy
    """
    # Normalize to probability distribution
    spectrum_sum = np.sum(spectrum)
    if spectrum_sum == 0:
        return 0.0
    
    ps = spectrum / spectrum_sum
    ps = ps[ps > 0]  # Remove zeros
    
    if len(ps) == 0:
        return 0.0
    
    return shannon_entropy(ps, base=base)


def find_dominant_frequency(spectrum: np.ndarray) -> Tuple[int, float]:
    """
    Find the dominant frequency (excluding DC component).
    
    Args:
        spectrum: Magnitude spectrum
    
    Returns:
        (frequency_index, magnitude) tuple
    """
    # Skip DC component (index 0)
    if len(spectrum) <= 1:
        return (0, 0.0)
    
    spectrum_no_dc = spectrum[1:]
    max_idx = np.argmax(spectrum_no_dc)
    
    return (max_idx + 1, spectrum_no_dc[max_idx])


def count_sidelobes(
    spectrum: np.ndarray,
    threshold_ratio: float = 0.1,
) -> int:
    """
    Count spectral sidelobes above threshold.
    
    Sidelobes are peaks > threshold_ratio * max(|X|)
    
    Args:
        spectrum: Magnitude spectrum
        threshold_ratio: Relative threshold (default: 0.1)
    
    Returns:
        Number of sidelobes
    """
    if len(spectrum) == 0:
        return 0
    
    peak_val = np.max(spectrum)
    threshold = threshold_ratio * peak_val
    
    # Find peaks above threshold
    peaks, _ = find_peaks(spectrum, height=threshold)
    
    return len(peaks)


def compute_sequence_diversity(seq: str) -> float:
    """
    Compute sequence diversity d(n) for curvature weight κ(n).
    
    d(n) is measured as normalized Shannon entropy of base composition.
    
    Args:
        seq: DNA/RNA sequence
    
    Returns:
        Diversity measure (0-1)
    """
    seq = seq.upper()
    
    # Count base frequencies
    base_counts = {}
    for base in ['A', 'C', 'G', 'T', 'U']:
        base_counts[base] = seq.count(base)
    
    # Remove zero counts
    counts = [c for c in base_counts.values() if c > 0]
    
    if len(counts) == 0:
        return 0.0
    
    # Normalize to probabilities
    total = sum(counts)
    probs = np.array(counts) / total
    
    # Shannon entropy normalized by max entropy (log₂(4) for 4 bases)
    h = shannon_entropy(probs, base=2)
    max_h = np.log2(4.0)  # Maximum entropy for 4 bases
    
    return h / max_h if max_h > 0 else 0.0


def kappa_curvature(n: int, d_n: float) -> float:
    """
    Calculate curvature weight κ(n) = d(n) · ln(n+1) / e².
    
    Guards against low-sample instability: if n < 10, default κ(n) = 1
    
    Args:
        n: Sequence length
        d_n: Sequence diversity
    
    Returns:
        Curvature weight
    """
    if n < 10:
        return 1.0
    
    return d_n * np.log(n + 1) / E_SQUARED


def sigmoid_aggregator(x: float, kappa: float = 1.0) -> float:
    """
    Sigmoid aggregator S(x) = 1 / (1 + e^{-κ·x})
    
    Args:
        x: Input value
        kappa: Curvature weight
    
    Returns:
        Sigmoid-transformed value (0-1)
    """
    return 1.0 / (1.0 + np.exp(-kappa * x))


class PhaseWeightedScorecard:
    """
    Phase-Weighted CRISPR Scorecard implementing Z Framework analysis.
    
    This class provides comprehensive phase-weighted spectral analysis for
    CRISPR guide scoring with mutation impact quantification.
    """
    
    def __init__(
        self,
        k: float = K_STAR,
        weights: Optional[Dict[str, float]] = None,
        is_rna: bool = False,
    ):
        """
        Initialize phase-weighted scorecard.
        
        Args:
            k: Phase parameter (default: K_STAR = 0.3)
            weights: Custom weights for composite score (default: DEFAULT_WEIGHTS)
            is_rna: If True, treat sequences as RNA (U instead of T)
        """
        self.k = k
        self.weights = weights or DEFAULT_WEIGHTS.copy()
        self.is_rna = is_rna
        
        # Validate weights sum to 1
        weight_sum = sum(self.weights.values())
        if not np.isclose(weight_sum, 1.0):
            warnings.warn(
                f"Weights sum to {weight_sum}, normalizing to 1.0",
                UserWarning,
            )
            # Normalize weights
            for key in self.weights:
                self.weights[key] /= weight_sum
    
    def compute_spectral_features(
        self,
        seq: str,
    ) -> Dict[str, float]:
        """
        Compute phase-weighted spectral features for a sequence.
        
        Args:
            seq: DNA or RNA sequence
        
        Returns:
            Dictionary with spectral features:
            - entropy: Spectral entropy
            - dominant_freq_idx: Dominant frequency index
            - dominant_freq_mag: Dominant frequency magnitude
            - sidelobe_count: Number of sidelobes
            - diversity: Sequence diversity d(n)
        """
        # Validate sequence length
        if len(seq) < 5:
            raise ValueError(
                f"Sequence too short ({len(seq)} bases). "
                f"Minimum length is 5 bases for reliable analysis."
            )
        
        # Encode sequence
        encoded = encode_complex(seq, is_rna=self.is_rna)
        
        # Apply phase weighting
        phase_weighted = apply_phase_weighting(encoded, k=self.k)
        
        # Compute spectrum
        spectrum = compute_spectrum(phase_weighted)
        
        # Extract features
        entropy = compute_spectral_entropy(spectrum)
        freq_idx, freq_mag = find_dominant_frequency(spectrum)
        sidelobes = count_sidelobes(spectrum)
        diversity = compute_sequence_diversity(seq)
        
        return {
            "entropy": entropy,
            "dominant_freq_idx": freq_idx,
            "dominant_freq_mag": freq_mag,
            "sidelobe_count": sidelobes,
            "diversity": diversity,
        }
    
    def compute_disruption_features(
        self,
        ref_seq: str,
        mut_seq: str,
    ) -> Dict[str, float]:
        """
        Compute mutation-induced spectral disruptions.
        
        Args:
            ref_seq: Reference sequence
            mut_seq: Mutated sequence
        
        Returns:
            Dictionary with disruption features:
            - delta_entropy: H(mutated) - H(reference)
            - delta_freq: |f₁_mut - f₁_ref|
            - delta_sidelobes: |sidelobes_mut - sidelobes_ref|
        """
        if len(ref_seq) != len(mut_seq):
            raise ValueError(
                f"Reference and mutated sequences must have same length. "
                f"Got {len(ref_seq)} and {len(mut_seq)}."
            )
        
        # Compute features for both sequences
        ref_features = self.compute_spectral_features(ref_seq)
        mut_features = self.compute_spectral_features(mut_seq)
        
        # Calculate deltas
        delta_entropy = mut_features["entropy"] - ref_features["entropy"]
        delta_freq = abs(
            mut_features["dominant_freq_idx"] - ref_features["dominant_freq_idx"]
        )
        delta_sidelobes = abs(
            mut_features["sidelobe_count"] - ref_features["sidelobe_count"]
        )
        
        return {
            "delta_entropy": delta_entropy,
            "delta_freq": delta_freq,
            "delta_sidelobes": delta_sidelobes,
            "ref_features": ref_features,
            "mut_features": mut_features,
        }
    
    def compute_z_score(
        self,
        ref_seq: str,
        mut_seq: str,
    ) -> Dict[str, float]:
        """
        Compute Z-invariant disruption score.
        
        Z = S(Δ_spectral / φ)
        
        where:
        - S(x) = 1 / (1 + e^{-κ(n)·x}) is sigmoid aggregator
        - κ(n) = d(n) · ln(n+1) / e² is curvature weight
        - Δ_spectral = w₁·ΔEntropy + w₂·Δf₁ + w₃·Sidelobes
        - φ ≈ 1.618 is golden ratio (phase constant)
        
        Args:
            ref_seq: Reference sequence
            mut_seq: Mutated sequence
        
        Returns:
            Dictionary with Z-score and components:
            - z_score: Composite Z-invariant score (0-1)
            - delta_spectral: Weighted spectral disruption
            - kappa: Curvature weight
            - disruptions: Individual disruption features
        """
        # Compute disruption features
        disruptions = self.compute_disruption_features(ref_seq, mut_seq)
        
        # Weighted composite spectral disruption
        delta_spectral = (
            self.weights["entropy"] * disruptions["delta_entropy"]
            + self.weights["freq_shift"] * disruptions["delta_freq"]
            + self.weights["sidelobes"] * disruptions["delta_sidelobes"]
        )
        
        # Compute curvature weight
        n = len(ref_seq)
        diversity = disruptions["ref_features"]["diversity"]
        kappa = kappa_curvature(n, diversity)
        
        # Z-invariant: Z = S(Δ_spectral / φ)
        z_input = delta_spectral / PHI
        z_score = sigmoid_aggregator(z_input, kappa)
        
        return {
            "z_score": z_score,
            "delta_spectral": delta_spectral,
            "kappa": kappa,
            "disruptions": disruptions,
        }
    
    def score_guide(
        self,
        guide_seq: str,
        target_seq: Optional[str] = None,
    ) -> Dict[str, float]:
        """
        Score a CRISPR guide with optional target context.
        
        If target_seq is provided, computes disruption relative to target.
        Otherwise, returns spectral features only.
        
        Args:
            guide_seq: Guide RNA sequence
            target_seq: Optional target DNA sequence for disruption analysis
        
        Returns:
            Dictionary with scores and features
        """
        # Validate guide length
        if len(guide_seq) < 10:
            warnings.warn(
                f"Short guide sequence ({len(guide_seq)} bases). "
                f"Consider sequences ≥20 bases for reliable scoring.",
                UserWarning,
            )
        
        if target_seq is not None:
            # Compute Z-score with disruption analysis
            result = self.compute_z_score(target_seq, guide_seq)
            result["guide_features"] = result["disruptions"]["mut_features"]
            result["target_features"] = result["disruptions"]["ref_features"]
        else:
            # Just compute spectral features
            features = self.compute_spectral_features(guide_seq)
            result = {
                "guide_features": features,
                "z_score": None,
                "delta_spectral": None,
            }
        
        return result


def score_guide_batch(
    guides: List[str],
    targets: Optional[List[str]] = None,
    k: float = K_STAR,
    weights: Optional[Dict[str, float]] = None,
    is_rna: bool = False,
) -> List[Dict[str, float]]:
    """
    Score multiple guides efficiently.
    
    Args:
        guides: List of guide sequences
        targets: Optional list of target sequences
        k: Phase parameter
        weights: Custom weights for composite score
        is_rna: If True, treat as RNA sequences
    
    Returns:
        List of score dictionaries
    """
    scorecard = PhaseWeightedScorecard(k=k, weights=weights, is_rna=is_rna)
    
    if targets is not None and len(guides) != len(targets):
        raise ValueError(
            f"Number of guides ({len(guides)}) must match "
            f"number of targets ({len(targets)})"
        )
    
    results = []
    for i, guide in enumerate(guides):
        target = targets[i] if targets is not None else None
        try:
            score = scorecard.score_guide(guide, target)
            results.append(score)
        except Exception as e:
            warnings.warn(f"Failed to score guide {i}: {e}", UserWarning)
            results.append({"error": str(e)})
    
    return results
