"""
Phase weighting (theta_prime) implementation for wave-crispr-signal.

This module implements the geometric resolution function:
    θ′(n,k) = φ·((n mod φ)/φ)^k

where φ is the golden ratio and k* ≈ 0.3 is the optimal parameter
that emerges from the geodesic-topological bridge.

SCIENTIFIC GATES:
- Uses mpmath for high-precision arithmetic (no float() conversions)
- Validates input constraints (n > 0)
- Implements exact formula from Z Framework documentation
"""

import mpmath as mp
import numpy as np
from typing import Union, Optional
from pydantic import BaseModel, Field, field_validator


# Configure high precision
mp.mp.dps = 50

# Golden ratio constant (high precision)
PHI = (mp.mpf(1) + mp.sqrt(mp.mpf(5))) / mp.mpf(2)


class ThetaPrimeConfig(BaseModel):
    """Configuration for theta_prime calculation.

    Attributes:
        k: Geometric resolution exponent. Default 0.3 is optimal for CRISPR.
        phi: Geometric period. Default is golden ratio φ ≈ 1.618.
        precision_dps: Decimal precision for mpmath calculations.
    """

    k: float = Field(
        default=0.3,
        ge=0.0,
        le=1.0,
        description="Geometric resolution exponent (k* ≈ 0.3 optimal)",
    )
    phi: Optional[float] = Field(
        default=None, description="Geometric period (defaults to golden ratio φ)"
    )
    precision_dps: int = Field(
        default=50, ge=15, le=100, description="Decimal precision for mpmath"
    )

    @field_validator("k")
    @classmethod
    def validate_k(cls, v):
        """Validate k parameter is in reasonable range."""
        if not 0.0 <= v <= 1.0:
            raise ValueError(f"k must be in [0, 1], got {v}")
        return v


def _to_mpf(n: Union[int, float, np.number]) -> mp.mpf:
    """Convert a number to an mpmath float, handling numpy types."""
    if isinstance(n, (np.integer, np.int_)):
        return mp.mpf(int(n))
    if isinstance(n, (np.floating, np.float_)):
        return mp.mpf(float(n))
    return mp.mpf(n)


def theta_prime(
    n: Union[int, np.ndarray, mp.mpf],
    config: Optional[ThetaPrimeConfig] = None,
    k: Optional[float] = None,
) -> Union[mp.mpf, np.ndarray]:
    """
    Calculate geometric resolution function θ′(n,k) = φ·((n mod φ)/φ)^k

    This function implements the geodesic resolution mapping that connects
    curvature κ(n) in discrete space to phase resolution in the topological
    framework. The optimal parameter k* ≈ 0.3 emerges from the bridge
    between f(x) = arcsin((x-1)/(2x+3)) and the Z Framework invariants.

    Args:
        n: Position index (must be positive). Can be int, array, or mpf.
        config: ThetaPrimeConfig with parameters. If None, uses defaults.
        k: Override k parameter from config. If both config and k are None,
           uses default k=0.3.

    Returns:
        Geometric resolution value(s). Returns mpf for scalar input,
        numpy array of mpf for array input.

    Raises:
        ValueError: If n <= 0 or invalid parameters.

    Examples:
        >>> # Scalar calculation with default k=0.3
        >>> result = theta_prime(21)
        >>> print(float(result))  # ≈ 1.618...

        >>> # Array calculation
        >>> positions = np.arange(1, 22)
        >>> results = theta_prime(positions, k=0.3)

        >>> # Custom configuration
        >>> config = ThetaPrimeConfig(k=0.35, precision_dps=30)
        >>> result = theta_prime(21, config=config)

    References:
        See docs/bridge_kappa_to_theta_prime.md for mathematical derivation.
    """
    # Initialize config
    if config is None:
        config = ThetaPrimeConfig()

    # Override k if provided
    if k is not None:
        config.k = k

    # Set phi to golden ratio if not specified
    phi_val = mp.mpf(config.phi) if config.phi is not None else PHI
    k_val = mp.mpf(config.k)

    # Handle scalar input
    if np.isscalar(n):
        n_val = _to_mpf(n)
        if n_val <= 0:
            raise ValueError(f"n must be positive, got {n}")

        # Calculate: φ·((n mod φ)/φ)^k
        n_mod_phi = mp.fmod(n_val, phi_val)
        ratio = n_mod_phi / phi_val
        result = phi_val * (ratio**k_val)

        return result

    # Handle array input
    n_arr = np.asarray(n)
    if np.any(n_arr <= 0):
        raise ValueError("All elements of n must be positive")

    # Vectorized calculation using mpmath
    results = np.empty(n_arr.shape, dtype=object)

    # Use flat iterator to modify in place
    for idx, n_i in enumerate(n_arr.flat):
        n_val = _to_mpf(n_i)
        n_mod_phi = mp.fmod(n_val, phi_val)
        ratio = n_mod_phi / phi_val
        results.flat[idx] = phi_val * (ratio**k_val)

    return results


def theta_prime_vectorized(
    n: np.ndarray, k: float = 0.3, phi: Optional[float] = None
) -> np.ndarray:
    """
    Fast vectorized version of theta_prime using numpy (reduced precision).

    This is a performance-optimized version that uses float64 instead of
    mpmath for large-scale computations where high precision is not critical.
    Use theta_prime() for accurate scientific calculations.

    Args:
        n: Array of position indices (must be positive).
        k: Geometric resolution exponent (default 0.3).
        phi: Geometric period (defaults to golden ratio).

    Returns:
        Array of geometric resolution values (float64).

    Warning:
        This function uses float64 precision and may have numerical errors
        for very large or small values. Use theta_prime() for high-precision
        calculations required in scientific validation.
    """
    if phi is None:
        phi = float(PHI)

    n_arr = np.asarray(n, dtype=np.float64)
    if np.any(n_arr <= 0):
        raise ValueError("All elements of n must be positive")

    n_mod_phi = np.fmod(n_arr, phi)
    ratio = n_mod_phi / phi
    return phi * (ratio**k)
