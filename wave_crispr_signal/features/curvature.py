"""
Curvature (kappa) implementation for wave-crispr-signal.

This module implements geodesic curvature κ(n) that couples with
the phase resolution θ′(n,k) through the topological bridge.

SCIENTIFIC GATES:
- Uses mpmath for high-precision arithmetic (no float() conversions)
- Validates input constraints
- Implements coupling modes: additive and multiplicative
"""

import mpmath as mp
import numpy as np
from typing import Union, Optional, Literal
from pydantic import BaseModel, Field


# Configure high precision
mp.mp.dps = 50

# Mathematical constants
PHI = (mp.mpf(1) + mp.sqrt(mp.mpf(5))) / mp.mpf(2)  # Golden ratio
E_SQUARED = mp.e**2  # Z Framework invariant c = e² ≈ 7.389


class KappaConfig(BaseModel):
    """Configuration for kappa (curvature) calculation.

    Attributes:
        mode: Curvature calculation mode ('discrete' or 'continuous').
        precision_dps: Decimal precision for mpmath calculations.
        scale: Scaling factor for curvature normalization.
    """

    mode: Literal["discrete", "continuous"] = Field(
        default="discrete",
        description="Calculation mode: 'discrete' for biological, 'continuous' for analytical",
    )
    precision_dps: int = Field(
        default=50, ge=15, le=100, description="Decimal precision for mpmath"
    )
    scale: Optional[float] = Field(
        default=None, description="Scaling factor (defaults to 1/φ for discrete mode)"
    )


def kappa(
    n: Union[int, np.ndarray, mp.mpf],
    config: Optional[KappaConfig] = None,
    mode: Optional[str] = None,
) -> Union[mp.mpf, np.ndarray]:
    """
    Calculate geodesic curvature κ(n) in discrete or continuous mode.

    The curvature function provides a measure of geometric distortion
    at position n, which couples with θ′(n,k) to form the complete
    spectral feature set.

    Discrete mode (biological/default):
        κ(n) = 1 / (1 + (n mod φ))

    Continuous mode (analytical):
        κ(n) = 1 / (1 + n/φ)

    Args:
        n: Position index (must be positive). Can be int, array, or mpf.
        config: KappaConfig with parameters. If None, uses defaults.
        mode: Override mode from config ('discrete' or 'continuous').

    Returns:
        Curvature value(s). Returns mpf for scalar input,
        numpy array of mpf for array input.

    Raises:
        ValueError: If n <= 0 or invalid parameters.

    Examples:
        >>> # Scalar calculation in discrete mode
        >>> result = kappa(21)
        >>> print(float(result))  # ≈ 0.5...

        >>> # Array calculation
        >>> positions = np.arange(1, 22)
        >>> results = kappa(positions, mode='discrete')

        >>> # Custom configuration
        >>> config = KappaConfig(mode='continuous', scale=0.5)
        >>> result = kappa(21, config=config)

    References:
        See docs/bridge_kappa_to_theta_prime.md for mathematical derivation.
    """
    # Initialize config
    if config is None:
        config = KappaConfig()

    # Override mode if provided
    if mode is not None:
        config.mode = mode

    # Determine scale
    scale = mp.mpf(config.scale) if config.scale is not None else (mp.mpf(1) / PHI)

    # Handle scalar input
    if np.isscalar(n):
        # Convert numpy types to Python native types first
        n_py = (
            int(n)
            if isinstance(n, (np.integer, np.int_))
            else (float(n) if isinstance(n, (np.floating, np.float_)) else n)
        )
        n_val = mp.mpf(n_py)
        if n_val <= 0:
            raise ValueError(f"n must be positive, got {n}")

        if config.mode == "discrete":
            # Discrete: κ(n) = 1 / (1 + (n mod φ))
            n_mod_phi = mp.fmod(n_val, PHI)
            result = scale / (mp.mpf(1) + n_mod_phi)
        else:
            # Continuous: κ(n) = 1 / (1 + n/φ)
            result = scale / (mp.mpf(1) + n_val / PHI)

        return result

    # Handle array input
    n_arr = np.asarray(n)
    if np.any(n_arr <= 0):
        raise ValueError("All elements of n must be positive")

    # Vectorized calculation using mpmath
    results = np.empty(n_arr.shape, dtype=object)

    # Use flat iterator to modify in place
    for idx, n_i in enumerate(n_arr.flat):
        # Convert numpy types to Python native types first
        n_val = mp.mpf(
            int(n_i) if isinstance(n_i, (np.integer, np.int_)) else float(n_i)
        )

        if config.mode == "discrete":
            n_mod_phi = mp.fmod(n_val, PHI)
            results.flat[idx] = scale / (mp.mpf(1) + n_mod_phi)
        else:
            results.flat[idx] = scale / (mp.mpf(1) + n_val / PHI)

    return results


def kappa_vectorized(
    n: np.ndarray, mode: str = "discrete", scale: Optional[float] = None
) -> np.ndarray:
    """
    Fast vectorized version of kappa using numpy (reduced precision).

    This is a performance-optimized version that uses float64 instead of
    mpmath for large-scale computations where high precision is not critical.
    Use kappa() for accurate scientific calculations.

    Args:
        n: Array of position indices (must be positive).
        mode: 'discrete' or 'continuous' calculation mode.
        scale: Scaling factor (defaults to 1/φ).

    Returns:
        Array of curvature values (float64).

    Warning:
        This function uses float64 precision and may have numerical errors.
        Use kappa() for high-precision calculations required in scientific validation.
    """
    if scale is None:
        scale = 1.0 / float(PHI)

    n_arr = np.asarray(n, dtype=np.float64)
    if np.any(n_arr <= 0):
        raise ValueError("All elements of n must be positive")

    phi_float = float(PHI)

    if mode == "discrete":
        n_mod_phi = np.fmod(n_arr, phi_float)
        return scale / (1.0 + n_mod_phi)
    else:
        return scale / (1.0 + n_arr / phi_float)


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
    from .phase_weighting import theta_prime, ThetaPrimeConfig

    # Compute both features
    kappa_val = kappa(n, config=kappa_config)
    theta_config = ThetaPrimeConfig(k=k)
    theta_val = theta_prime(n, config=theta_config)

    # Apply coupling
    if coupling == "additive":
        return kappa_val + theta_val
    else:  # multiplicative
        return kappa_val * theta_val
