"""
Feature extraction modules for wave-crispr-signal.

This package provides implementations of geometric resolution (theta_prime)
and curvature (kappa) features for CRISPR guide scoring.
"""

from .phase_weighting import theta_prime, ThetaPrimeConfig
from .curvature import kappa, KappaConfig

__all__ = [
    "theta_prime",
    "ThetaPrimeConfig",
    "kappa",
    "KappaConfig",
]
