"""
Coupling functions for wave-crispr-signal features.

This module provides functions to combine low-level features like
curvature (kappa) and phase resolution (theta_prime) into
higher-level coupled features.

This avoids circular dependencies between the feature modules.
"""

from typing import Union, Optional, Literal

import mpmath as mp
import numpy as np

from .curvature import kappa, KappaConfig
from .phase_weighting import theta_prime, ThetaPrimeConfig


def compute_coupled_features(
    n: Union[int, np.ndarray],
    k: float = 0.3,
    coupling: Literal["additive", "multiplicative"] = "multiplicative",
    kappa_config: Optional[KappaConfig] = None,
) -> Union[mp.mpf, np.ndarray]:
    """
    Compute coupled κ(n) + θ′(n,k) features with specified coupling mode.

    This function combines curvature and phase resolution using either
    additive or multiplicative coupling as specified in the bridge derivation.

    Additive coupling:
        F(n) = κ(n) + θ′(n,k)

    Multiplicative coupling:
        F(n) = κ(n) · θ′(n,k)

    Args:
        n: Position index (must be positive).
        k: Phase resolution parameter (default 0.3).
        coupling: Coupling mode ('additive' or 'multiplicative').
        kappa_config: Configuration for kappa calculation.

    Returns:
        Coupled feature value(s).

    Examples:
        >>> # Multiplicative coupling (default, consistent with derivation)
        >>> result = compute_coupled_features(21, k=0.3, coupling='multiplicative')

        >>> # Additive coupling for ablation study
        >>> result = compute_coupled_features(21, k=0.3, coupling='additive')
    """
    # Compute both features
    kappa_val = kappa(n, config=kappa_config)
    theta_config = ThetaPrimeConfig(k=k)
    theta_val = theta_prime(n, config=theta_config)

    # Apply coupling
    if coupling == "additive":
        return kappa_val + theta_val
    else:  # multiplicative
        return kappa_val * theta_val
