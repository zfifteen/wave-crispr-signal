"""
Feature extraction modules for wave-crispr-signal.

This package provides implementations of geometric resolution (theta_prime),
curvature (kappa), and φ-geometry features for CRISPR guide scoring.
"""

from .phase_weighting import theta_prime, ThetaPrimeConfig
from .curvature import kappa, KappaConfig
from .phi_geometry import (
    PHI,
    DNA_LENGTH_DIAMETER_RATIO,
    DNA_MAJOR_MINOR_SEP_RATIO,
    DNA_HELICAL_PERIOD,
    PHI_MODES,
    dna_phi_phase,
    dna_phi_phase_vectorized,
    dna_phi_curvature,
    phi_phase_score,
    phi_curvature_score,
    compute_phi_features,
    uniform_phase_score,
    random_phase_score,
    simple_gc_content,
)

__all__ = [
    "theta_prime",
    "ThetaPrimeConfig",
    "kappa",
    "KappaConfig",
    # φ-geometry constants
    "PHI",
    "DNA_LENGTH_DIAMETER_RATIO",
    "DNA_MAJOR_MINOR_SEP_RATIO",
    "DNA_HELICAL_PERIOD",
    "PHI_MODES",
    # φ-geometry functions
    "dna_phi_phase",
    "dna_phi_phase_vectorized",
    "dna_phi_curvature",
    "phi_phase_score",
    "phi_curvature_score",
    "compute_phi_features",
    # Baseline functions for comparison
    "uniform_phase_score",
    "random_phase_score",
    "simple_gc_content",
]
