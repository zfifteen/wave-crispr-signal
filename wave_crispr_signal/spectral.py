"""
Spectral utilities for frequency-native DNA dynamics features.
- Complex encoding with dimensionless AT/GC opening-rate ratio r = k_GC / k_AT
- Exact evaluation at fractional periods via a one-bin DFT (CZT-lite)
- Rotational phase utilities
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Tuple, List

HELICAL_PERIOD = 10.5

def validate_dna_sequence(seq: str) -> str:
    """Validate DNA sequence and return uppercase normalized text."""
    seq = seq.upper().replace(" ", "").replace("\n", "")
    invalid_bases = set(seq) - set("ACGT")
    if invalid_bases:
        raise ValueError(
            f"Invalid DNA bases detected: {invalid_bases}. "
            f"Only A, C, G, T allowed for DNA sequences. "
            f"Use A, C, G, U for RNA sequences (convert separately)."
        )
    return seq


def validate_rna_sequence(seq: str) -> str:
    """Validate RNA sequence and return uppercase normalized text."""
    seq = seq.upper().replace(" ", "").replace("\n", "")
    invalid_bases = set(seq) - set("ACGU")
    if invalid_bases:
        raise ValueError(
            f"Invalid RNA bases detected: {invalid_bases}. "
            f"Only A, C, G, U allowed for RNA sequences. "
            f"Use A, C, G, T for DNA sequences (convert separately)."
        )
    return seq


def encode_complex(seq: str, r: float = 20.0, is_rna: bool = False) -> np.ndarray:
    """
    Map DNA string -> complex vector using dimensionless opening-rate ratio.
    A,T receive negative real + positive imag; G,C inverse-signed.
    r ~ k_GC/k_AT controls contrast; default 20 (biophysically plausible order-of-magnitude).
    """
    if is_rna:
        seq = validate_rna_sequence(seq)
    else:
        seq = validate_dna_sequence(seq)

    alpha = float(np.log(max(r, 1.0000001)))  # strength (real)
    beta = 0.3 * alpha  # phase skew (imag)
    if is_rna:
        tbl = {
            "A": -alpha + 1j * beta,
            "U": -alpha + 1j * beta,
            "G": alpha - 1j * beta,
            "C": alpha - 1j * beta,
        }
    else:
        tbl = {
            "A": -alpha + 1j * beta,
            "T": -alpha + 1j * beta,
            "G": alpha - 1j * beta,
            "C": alpha - 1j * beta,
        }
    return np.array([tbl.get(ch, 0.0 + 0.0j) for ch in seq], dtype=np.complex128)


def _single_frequency_dft(x: np.ndarray, f: float) -> complex:
    """
    Evaluate the DFT exactly at normalized frequency f (cycles per sample).
    Equivalent to a 1-bin CZT. x should already be mean-centered and windowed.
    """
    n = np.arange(x.size, dtype=np.float64)
    return np.sum(x * np.exp(-2j * np.pi * f * n))


def cz_power_at_period(z: np.ndarray, period: float) -> float:
    """Return magnitude of the complex spectrum at the given period (non-integer allowed)."""
    if len(z) == 0:
        return 0.0
    # DC removal + Hann window
    z = z - np.mean(z)
    if len(z) >= 2:
        n = np.arange(len(z), dtype=np.float64)
        hann = 0.5 - 0.5 * np.cos(2 * np.pi * n / (len(z) - 1))
        z = z * hann
    f = 1.0 / float(period)
    X = _single_frequency_dft(z, f)
    return float(np.abs(X))


def breathing_features(seq: str, r: float = 20.0, is_rna: bool = False) -> Dict[str, float]:
    """Three-term spectrum around the helical period and its first two harmonics."""
    z = encode_complex(seq, r=r, is_rna=is_rna)
    seq_upper = seq.upper()
    gc_content = (
        (seq_upper.count("G") + seq_upper.count("C")) / len(seq_upper)
        if len(seq_upper) > 0
        else 0.0
    )
    return {
        "P10_5": cz_power_at_period(z, 10.5),
        "P5_25": cz_power_at_period(z, 5.25),
        "P3_5": cz_power_at_period(z, 3.5),
        "length": len(seq_upper),
        "gc_content": gc_content,
    }


def rotational_phase_at(pos: int, period: float = HELICAL_PERIOD) -> float:
    """Return circular phase in [0, 2π) for a cut-site position index."""
    return float((2.0 * np.pi * pos / period) % (2.0 * np.pi))


def phase_at(pos: int, period: float = HELICAL_PERIOD) -> float:
    """Alias kept for compatibility with existing demo imports."""
    return rotational_phase_at(pos, period=period)


def rotational_phase_curve(
    seq: str,
    period: float = HELICAL_PERIOD,
    bins: int = 24,
    r: float = 20.0,
    is_rna: bool = False,
) -> Tuple[List[float], List[float]]:
    """
    Compute a phase-binned curve of average encoded magnitude across positions.
    Returns (phase_centers, averaged_magnitudes).
    """
    z = encode_complex(seq, r=r, is_rna=is_rna)
    N = len(z)
    if N == 0:
        centers = np.linspace(0, 2 * np.pi, bins, endpoint=False)
        return centers.tolist(), [0.0] * bins
    phases = (2 * np.pi * np.arange(N) / period) % (2 * np.pi)
    mags = np.abs(z)
    edges = np.linspace(0, 2 * np.pi, bins + 1)
    idx = np.digitize(phases, edges) - 1
    tot = np.bincount(idx, weights=mags, minlength=bins).astype(np.float64)
    cnt = np.bincount(idx, minlength=bins).astype(np.float64)
    avg = tot / np.clip(cnt, 1.0, None)
    centers = (edges[:-1] + edges[1:]) * 0.5
    return centers.tolist(), avg.tolist()


def spectral_peak_analysis(
    seq: str,
    period: float = HELICAL_PERIOD,
    r: float = 20.0,
    is_rna: bool = False,
) -> Dict[str, float]:
    """Return complex peak metrics (power, angle, normalized_power) at period."""
    z = encode_complex(seq, r=r, is_rna=is_rna)
    n = len(z)
    if n == 0:
        return {"power": 0.0, "angle": 0.0, "normalized_power": 0.0}

    z = z - np.mean(z)
    if n >= 2:
        idx = np.arange(n, dtype=np.float64)
        hann = 0.5 - 0.5 * np.cos(2 * np.pi * idx / (n - 1))
        z = z * hann

    f = 1.0 / period
    idx = np.arange(n)
    x = np.sum(z * np.exp(-2j * np.pi * f * idx))
    power = float(np.abs(x))
    return {
        "power": power,
        "angle": float(np.angle(x)),
        "normalized_power": power / np.sqrt(n),
    }


def baseline_seq_features(seq: str) -> Dict[str, float]:
    """Lightweight baseline: GC%, AT%, length, dinucleotide counts (subset)."""
    seq = seq.upper().replace("U", "T")
    n = max(1, len(seq))
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    feats = {
        "len": float(n),
        "gc_frac": float(gc) / n,
        "at_frac": float(at) / n,
        "a": float(seq.count("A")) / n,
        "t": float(seq.count("T")) / n,
        "g": float(seq.count("G")) / n,
        "c": float(seq.count("C")) / n,
    }
    # Minimal di-nucs
    for d in ["AA", "AT", "TA", "TT", "GG", "GC", "CG", "CC"]:
        feats[f"di_{d}"] = float(seq.count(d)) / max(1, n - 1)
    return feats


def make_feature_row(target_seq: str, r: float = 20.0) -> Dict[str, float]:
    feats = {}
    feats.update(baseline_seq_features(target_seq))
    feats.update(breathing_features(target_seq, r=r))
    return feats
