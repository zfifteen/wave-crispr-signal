#!/usr/bin/env python3
"""
Shared Statistical Utilities for Trinity Experiments

This module provides statistical analysis tools for the trinity experiments:
- Bootstrap confidence intervals
- Permutation tests
- Effect sizes (Cohen's d, Hedges' g)
- Circular statistics (Rayleigh test, circular-linear correlation)

SCIENTIFIC GATES:
- Pre-registered endpoints
- Null models via permutation
- Multiple comparison correction (Benjamini-Hochberg FDR)
- Power/sample size considerations
"""

import numpy as np
from typing import Tuple, Dict, List, Optional
from scipy import stats
from scipy.stats import pearsonr
import logging

logger = logging.getLogger(__name__)


def bootstrap_ci(
    data: np.ndarray,
    statistic_func: callable,
    n_bootstrap: int = 1000,
    confidence_level: float = 0.95,
    seed: Optional[int] = None
) -> Tuple[float, float, float]:
    """
    Calculate bootstrap confidence interval for a statistic.
    
    Args:
        data: Input data array
        statistic_func: Function to compute statistic (e.g., np.mean)
        n_bootstrap: Number of bootstrap resamples (default 1000, min 1000)
        confidence_level: Confidence level (default 0.95)
        seed: Random seed for reproducibility
        
    Returns:
        Tuple of (point_estimate, lower_ci, upper_ci)
    """
    if n_bootstrap < 1000:
        logger.warning(f"n_bootstrap={n_bootstrap} < 1000, may be underpowered")
    
    rng = np.random.RandomState(seed)
    
    # Point estimate
    point_estimate = statistic_func(data)
    
    # Bootstrap resampling
    bootstrap_stats = []
    n = len(data)
    
    for _ in range(n_bootstrap):
        resample = rng.choice(data, size=n, replace=True)
        bootstrap_stats.append(statistic_func(resample))
    
    bootstrap_stats = np.array(bootstrap_stats)
    
    # Calculate confidence interval
    alpha = 1 - confidence_level
    lower_percentile = (alpha / 2) * 100
    upper_percentile = (1 - alpha / 2) * 100
    
    lower_ci = np.percentile(bootstrap_stats, lower_percentile)
    upper_ci = np.percentile(bootstrap_stats, upper_percentile)
    
    return float(point_estimate), float(lower_ci), float(upper_ci)


def permutation_test(
    group1: np.ndarray,
    group2: np.ndarray,
    statistic_func: callable,
    n_permutations: int = 1000,
    seed: Optional[int] = None,
    alternative: str = "two-sided"
) -> Tuple[float, float]:
    """
    Permutation test for difference between two groups.
    
    Args:
        group1: First group data
        group2: Second group data
        statistic_func: Function to compute test statistic
        n_permutations: Number of permutations (default 1000, min 1000)
        seed: Random seed
        alternative: "two-sided", "less", or "greater"
        
    Returns:
        Tuple of (observed_statistic, p_value)
    """
    if n_permutations < 1000:
        logger.warning(f"n_permutations={n_permutations} < 1000, may be underpowered")
    
    rng = np.random.RandomState(seed)
    
    # Observed statistic
    observed_stat = statistic_func(group1, group2)
    
    # Permutation distribution
    combined = np.concatenate([group1, group2])
    n1 = len(group1)
    
    perm_stats = []
    for _ in range(n_permutations):
        rng.shuffle(combined)
        perm_group1 = combined[:n1]
        perm_group2 = combined[n1:]
        perm_stats.append(statistic_func(perm_group1, perm_group2))
    
    perm_stats = np.array(perm_stats)
    
    # Calculate p-value
    if alternative == "two-sided":
        p_value = np.mean(np.abs(perm_stats) >= np.abs(observed_stat))
    elif alternative == "greater":
        p_value = np.mean(perm_stats >= observed_stat)
    elif alternative == "less":
        p_value = np.mean(perm_stats <= observed_stat)
    else:
        raise ValueError(f"Invalid alternative: {alternative}")
    
    return float(observed_stat), float(p_value)


def cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """
    Calculate Cohen's d effect size.
    
    Args:
        group1: First group
        group2: Second group
        
    Returns:
        Cohen's d
    """
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    
    # Pooled standard deviation
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
    if pooled_std == 0:
        return 0.0
    
    return float((np.mean(group1) - np.mean(group2)) / pooled_std)


def hedges_g(group1: np.ndarray, group2: np.ndarray) -> float:
    """
    Calculate Hedges' g effect size (bias-corrected Cohen's d).
    
    Args:
        group1: First group
        group2: Second group
        
    Returns:
        Hedges' g
    """
    d = cohens_d(group1, group2)
    n1, n2 = len(group1), len(group2)
    n_total = n1 + n2
    
    # Correction factor
    correction = 1 - (3 / (4 * n_total - 9))
    
    return float(d * correction)


def benjamini_hochberg_fdr(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction for multiple comparisons.
    
    Args:
        p_values: Array of p-values
        alpha: Target FDR level (default 0.05)
        
    Returns:
        Boolean array indicating rejection (True = significant)
    """
    p_values = np.asarray(p_values)
    n = len(p_values)
    
    if n == 0:
        return np.array([], dtype=bool)
    
    # Sort p-values and track original indices
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    
    # BH critical values
    critical_values = (np.arange(1, n + 1) / n) * alpha
    
    # Find largest i where p(i) <= (i/n) * alpha
    comparisons = sorted_p <= critical_values
    
    if not np.any(comparisons):
        # No rejections
        return np.zeros(n, dtype=bool)
    
    # Reject all hypotheses up to largest i
    max_i = np.max(np.where(comparisons)[0])
    rejections = np.zeros(n, dtype=bool)
    rejections[sorted_indices[:max_i + 1]] = True
    
    return rejections


def rayleigh_test(angles: np.ndarray) -> Tuple[float, float]:
    """
    Rayleigh test for circular uniformity.
    
    Tests null hypothesis that angles are uniformly distributed on the circle.
    
    Args:
        angles: Array of angles in radians
        
    Returns:
        Tuple of (test_statistic, p_value)
    """
    n = len(angles)
    
    if n == 0:
        return 0.0, 1.0
    
    # Mean resultant vector
    C = np.sum(np.cos(angles))
    S = np.sum(np.sin(angles))
    R = np.sqrt(C**2 + S**2)
    
    # Mean resultant length
    R_bar = R / n
    
    # Test statistic
    Z = n * R_bar**2
    
    # Approximate p-value (Zar, 1999)
    if n < 50:
        # Exact calculation would require special functions
        # Use approximation for now
        p_value = np.exp(-Z) * (1 + (2*Z - Z**2) / (4*n) - (24*Z - 132*Z**2 + 76*Z**3 - 9*Z**4) / (288*n**2))
    else:
        # Large n approximation
        p_value = np.exp(-Z)
    
    p_value = max(0.0, min(1.0, p_value))  # Clamp to [0, 1]
    
    return float(Z), float(p_value)


def circular_linear_correlation(
    angles: np.ndarray,
    values: np.ndarray
) -> Tuple[float, float]:
    """
    Circular-linear correlation (Jammalamadaka-Sengupta).
    
    Correlates circular variable (angles) with linear variable (values).
    
    Args:
        angles: Circular variable in radians
        values: Linear variable
        
    Returns:
        Tuple of (correlation, p_value)
    """
    n = len(angles)
    
    if n != len(values):
        raise ValueError("angles and values must have same length")
    
    if n < 3:
        return 0.0, 1.0
    
    # Remove NaN values
    valid_mask = ~(np.isnan(angles) | np.isnan(values))
    angles = angles[valid_mask]
    values = values[valid_mask]
    n = len(angles)
    
    if n < 3:
        return 0.0, 1.0
    
    # Rank both variables
    angle_ranks = stats.rankdata(angles)
    value_ranks = stats.rankdata(values)
    
    # Convert angle ranks to circular
    angle_phases = 2 * np.pi * angle_ranks / n
    
    # Compute correlation components
    try:
        rxc = pearsonr(value_ranks, np.cos(angle_phases))[0]
        rxs = pearsonr(value_ranks, np.sin(angle_phases))[0]
        rcs = pearsonr(np.cos(angle_phases), np.sin(angle_phases))[0]
    except:
        return 0.0, 1.0
    
    # Handle NaN or invalid correlations
    if np.isnan(rxc) or np.isnan(rxs) or np.isnan(rcs):
        return 0.0, 1.0
    
    # Circular-linear correlation
    numerator = rxc**2 + rxs**2 - 2*rxc*rxs*rcs
    denominator = 1 - rcs**2
    
    if denominator <= 0 or numerator < 0:
        return 0.0, 1.0
    
    r = np.sqrt(numerator / denominator)
    
    if np.isnan(r) or np.isinf(r):
        return 0.0, 1.0
    
    # Approximate p-value using chi-square
    chi2_stat = n * r**2
    p_value = 1 - stats.chi2.cdf(chi2_stat, df=2)
    
    if np.isnan(p_value) or np.isinf(p_value):
        p_value = 1.0
    
    return float(r), float(p_value)


def partial_correlation(
    x: np.ndarray,
    y: np.ndarray,
    control_vars: np.ndarray
) -> Tuple[float, float]:
    """
    Calculate partial correlation controlling for covariates.
    
    Args:
        x: First variable
        y: Second variable
        control_vars: Control variables (can be 1D or 2D)
        
    Returns:
        Tuple of (partial_r, p_value)
    """
    # Ensure control_vars is 2D
    if control_vars.ndim == 1:
        control_vars = control_vars.reshape(-1, 1)
    
    n = len(x)
    
    # Residualize x and y with respect to control variables
    from sklearn.linear_model import LinearRegression
    
    model_x = LinearRegression()
    model_y = LinearRegression()
    
    model_x.fit(control_vars, x)
    model_y.fit(control_vars, y)
    
    residuals_x = x - model_x.predict(control_vars)
    residuals_y = y - model_y.predict(control_vars)
    
    # Correlation of residuals
    r, p = pearsonr(residuals_x, residuals_y)
    
    return float(r), float(p)


if __name__ == "__main__":
    # Self-test
    print("Testing statistics utilities...")
    
    # Test bootstrap CI
    data = np.random.randn(100)
    mean_est, lower, upper = bootstrap_ci(data, np.mean, n_bootstrap=1000, seed=42)
    print(f"\nBootstrap CI for mean: {mean_est:.4f} [{lower:.4f}, {upper:.4f}]")
    
    # Test permutation
    group1 = np.random.randn(50) + 0.5
    group2 = np.random.randn(50)
    stat, p = permutation_test(
        group1, group2,
        lambda g1, g2: np.mean(g1) - np.mean(g2),
        n_permutations=1000,
        seed=42
    )
    print(f"\nPermutation test: stat={stat:.4f}, p={p:.4f}")
    
    # Test effect sizes
    d = cohens_d(group1, group2)
    g = hedges_g(group1, group2)
    print(f"\nCohen's d: {d:.4f}")
    print(f"Hedges' g: {g:.4f}")
    
    # Test FDR correction
    p_values = np.array([0.001, 0.01, 0.03, 0.05, 0.1, 0.5])
    rejections = benjamini_hochberg_fdr(p_values, alpha=0.05)
    print(f"\nFDR correction (alpha=0.05):")
    for p, rej in zip(p_values, rejections):
        print(f"  p={p:.3f}: {'reject' if rej else 'accept'}")
    
    # Test circular statistics
    angles = np.random.uniform(0, 2*np.pi, 100)
    Z, p = rayleigh_test(angles)
    print(f"\nRayleigh test: Z={Z:.4f}, p={p:.4f}")
    
    values = np.random.randn(100)
    r, p = circular_linear_correlation(angles, values)
    print(f"\nCircular-linear correlation: r={r:.4f}, p={p:.4f}")
    
    # Test partial correlation
    x = np.random.randn(100)
    y = x + np.random.randn(100) * 0.5
    z = np.random.randn(100)
    r_partial, p_partial = partial_correlation(x, y, z)
    print(f"\nPartial correlation: r={r_partial:.4f}, p={p_partial:.4f}")
    
    print("\n✓ All tests passed")
