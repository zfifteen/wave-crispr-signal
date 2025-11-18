"""
Spectral Disruption Profiler - CRISPR Guide Efficiency Analysis

This module implements spectral disruption profiling for CRISPR guide sequences
using phase-weighted FFT analysis with geometric resolution θ'(n,k) = φ·((n mod φ)/φ)^k.

Components:
- encoding: DNA sequence to weighted waveform mapping
- analysis: FFT-based feature extraction (Δf₁, ΔEntropy, sidelobes)
- scoring: Z-invariant composite scoring with bootstrap CI
- detection: Off-target detection via entropy gradients
- visualization: Dashboard and plotting utilities
- cli: Command-line interface for batch processing

Scientific Gates:
- Human DNA only (A/C/G/T for DNA, A/C/G/U for RNA)
- No fabrication (real nucleotides only)
- Fail-fast validation with clear error messages
- Bootstrap CI (≥1,000 resamples) for statistical validation
- Seed control for reproducibility
"""

__version__ = "1.0.0"

from .encoding import encode_sequence, phase_weighted_encoding
from .analysis import compute_spectral_features, analyze_disruption
from .scoring import compute_z_score, compute_composite_score
from .detection import detect_off_targets, compute_gc_resonance

__all__ = [
    "encode_sequence",
    "phase_weighted_encoding", 
    "compute_spectral_features",
    "analyze_disruption",
    "compute_z_score",
    "compute_composite_score",
    "detect_off_targets",
    "compute_gc_resonance",
]
