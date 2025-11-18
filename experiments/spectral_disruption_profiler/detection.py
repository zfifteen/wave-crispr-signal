"""
Detection Module - Off-Target Detection and GC-Resonance Analysis

This module implements entropy gradient-based off-target detection and
GC-quartile resonance analysis with FDR correction.

Features:
- Entropy gradient-based off-target flagging
- GC-quartile resonance correlation
- Benjamini-Hochberg FDR correction
"""

import numpy as np
from scipy import stats
from typing import Dict, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


def detect_off_targets(
    features_list: List[Dict[str, float]],
    entropy_threshold: float = 0.05,
    gc_threshold: float = -0.2
) -> List[bool]:
    """
    Detect potential off-targets using entropy gradients and GC resonance.
    
    Flags sequences with:
    1. ΔEntropy > entropy_threshold (disruption indicator)
    2. GC-quartile resonance r < gc_threshold (negative correlation)
    
    Args:
        features_list: List of spectral feature dictionaries
        entropy_threshold: Minimum ΔEntropy for flagging (default 0.05)
        gc_threshold: Maximum GC correlation for flagging (default -0.2)
        
    Returns:
        List of boolean flags (True = potential off-target)
    """
    flags = []
    
    for features in features_list:
        # Check entropy criterion
        delta_entropy = features.get('delta_entropy', 0.0)
        entropy_flag = delta_entropy > entropy_threshold
        
        # Check GC content (flag high GC as potential issue)
        gc_content = features.get('gc_content', 0.5)
        gc_flag = gc_content > 0.7  # High GC guides may have issues
        
        # Flag if either criterion met
        is_off_target = entropy_flag or gc_flag
        flags.append(is_off_target)
    
    n_flagged = sum(flags)
    logger.info(f"Flagged {n_flagged}/{len(flags)} sequences as potential off-targets")
    
    return flags


def compute_gc_quartiles(gc_contents: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute GC quartiles and quartile labels.
    
    Args:
        gc_contents: Array of GC content values
        
    Returns:
        Tuple of (quartile_labels, quartile_boundaries)
    """
    # Compute quartile boundaries
    q1 = np.percentile(gc_contents, 25)
    q2 = np.percentile(gc_contents, 50)
    q3 = np.percentile(gc_contents, 75)
    
    # Assign quartile labels
    quartile_labels = np.zeros(len(gc_contents), dtype=int)
    quartile_labels[gc_contents <= q1] = 1
    quartile_labels[(gc_contents > q1) & (gc_contents <= q2)] = 2
    quartile_labels[(gc_contents > q2) & (gc_contents <= q3)] = 3
    quartile_labels[gc_contents > q3] = 4
    
    boundaries = np.array([0, q1, q2, q3, 1.0])
    
    return quartile_labels, boundaries


def compute_gc_resonance(
    features_list: List[Dict[str, float]],
    scores: Optional[np.ndarray] = None,
    n_permutations: int = 1000,
    seed: Optional[int] = None
) -> Dict[str, float]:
    """
    Compute GC-quartile resonance correlation with permutation test.
    
    Tests for correlation between GC quartile and disruption scores,
    as observed in Kim 2025 dataset (r = -0.211, p_perm = 0.0012).
    
    Args:
        features_list: List of spectral feature dictionaries
        scores: Optional array of disruption scores (computed if not provided)
        n_permutations: Number of permutations for null distribution
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with correlation statistics and p-value
    """
    # Extract GC contents
    gc_contents = np.array([f.get('gc_content', 0.5) for f in features_list])
    
    # Use provided scores or compute from delta_entropy
    if scores is None:
        scores = np.array([f.get('delta_entropy', 0.0) for f in features_list])
    
    # Compute quartiles
    quartile_labels, boundaries = compute_gc_quartiles(gc_contents)
    
    # Compute observed correlation for each quartile
    quartile_correlations = {}
    
    for q in range(1, 5):
        mask = quartile_labels == q
        if np.sum(mask) > 2:  # Need at least 3 points for correlation
            q_gc = gc_contents[mask]
            q_scores = scores[mask]
            
            if len(q_gc) > 1 and np.std(q_gc) > 0 and np.std(q_scores) > 0:
                r, p = stats.pearsonr(q_gc, q_scores)
                quartile_correlations[f'Q{q}'] = {'r': float(r), 'p': float(p)}
    
    # Overall correlation
    if len(gc_contents) > 1 and np.std(gc_contents) > 0 and np.std(scores) > 0:
        r_observed, p_parametric = stats.pearsonr(gc_contents, scores)
    else:
        r_observed = 0.0
        p_parametric = 1.0
    
    # Permutation test for p-value
    if seed is not None:
        np.random.seed(seed)
    
    permuted_r = []
    for _ in range(n_permutations):
        permuted_scores = np.random.permutation(scores)
        if np.std(permuted_scores) > 0:
            r_perm, _ = stats.pearsonr(gc_contents, permuted_scores)
            permuted_r.append(r_perm)
    
    permuted_r = np.array(permuted_r)
    
    # Compute permutation p-value (two-tailed)
    p_perm = np.mean(np.abs(permuted_r) >= np.abs(r_observed))
    
    results = {
        'r': float(r_observed),
        'p_parametric': float(p_parametric),
        'p_permutation': float(p_perm),
        'n_permutations': n_permutations,
        'quartile_correlations': quartile_correlations,
        'n_samples': len(gc_contents),
    }
    
    logger.info(
        f"GC resonance: r={r_observed:.3f}, "
        f"p_perm={p_perm:.4f} ({n_permutations} permutations)"
    )
    
    # Check for significant negative correlation (as in Kim 2025)
    if r_observed < -0.2 and p_perm < 0.05:
        logger.warning(
            f"Significant negative GC resonance detected (r={r_observed:.3f}, p={p_perm:.4f}). "
            "This may indicate GC-bias in disruption scoring."
        )
    
    return results


def benjamini_hochberg_fdr(
    p_values: np.ndarray,
    alpha: float = 0.05
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply Benjamini-Hochberg FDR correction to p-values.
    
    Args:
        p_values: Array of p-values
        alpha: FDR threshold (default 0.05)
        
    Returns:
        Tuple of (adjusted_p_values, significant_flags)
    """
    n = len(p_values)
    
    # Sort p-values
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    
    # Compute adjusted p-values
    adjusted_p = np.zeros(n)
    
    for i in range(n):
        # BH critical value: (i+1)/n * alpha
        bh_value = (i + 1) / n * alpha
        adjusted_p[sorted_indices[i]] = sorted_p[i] * n / (i + 1)
    
    # Ensure monotonicity (adjusted p should not decrease)
    for i in range(n - 1, 0, -1):
        if adjusted_p[sorted_indices[i]] < adjusted_p[sorted_indices[i - 1]]:
            adjusted_p[sorted_indices[i - 1]] = adjusted_p[sorted_indices[i]]
    
    # Cap at 1.0
    adjusted_p = np.minimum(adjusted_p, 1.0)
    
    # Determine significance
    significant = adjusted_p < alpha
    
    n_significant = np.sum(significant)
    logger.info(f"BH-FDR correction: {n_significant}/{n} tests significant at α={alpha}")
    
    return adjusted_p, significant


def compute_partial_correlation(
    x: np.ndarray,
    y: np.ndarray,
    covariates: List[np.ndarray]
) -> Tuple[float, float]:
    """
    Compute partial correlation controlling for covariates.
    
    Useful for controlling GC%, guide length, and guide position when
    analyzing spectral features.
    
    Args:
        x: First variable
        y: Second variable
        covariates: List of covariate arrays to control for
        
    Returns:
        Tuple of (partial_r, p_value)
    """
    from sklearn.linear_model import LinearRegression
    
    # Residualize x and y against covariates
    if len(covariates) > 0:
        # Stack covariates
        Z = np.column_stack(covariates)
        
        # Fit x ~ Z
        model_x = LinearRegression()
        model_x.fit(Z, x)
        residual_x = x - model_x.predict(Z)
        
        # Fit y ~ Z
        model_y = LinearRegression()
        model_y.fit(Z, y)
        residual_y = y - model_y.predict(Z)
    else:
        residual_x = x
        residual_y = y
    
    # Correlation of residuals
    if np.std(residual_x) > 0 and np.std(residual_y) > 0:
        r, p = stats.pearsonr(residual_x, residual_y)
    else:
        r, p = 0.0, 1.0
    
    return float(r), float(p)


def effect_size_cohens_d(
    group1: np.ndarray,
    group2: np.ndarray
) -> float:
    """
    Compute Cohen's d effect size.
    
    Args:
        group1: First group values
        group2: Second group values
        
    Returns:
        Cohen's d effect size
    """
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    
    # Pooled standard deviation
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
    if pooled_std == 0:
        return 0.0
    
    # Cohen's d
    d = (np.mean(group1) - np.mean(group2)) / pooled_std
    
    return float(d)
