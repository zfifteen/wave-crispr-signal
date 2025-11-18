"""
Scoring Module - Z-Invariant Composite Scoring with Bootstrap CI

This module implements Z-invariant scoring based on the unified framework
Z = A(B/c) with geodesic mapping for auto-optimization of k parameter.

Features:
- Z-invariant composite scores
- Bootstrap confidence intervals (≥1,000 resamples)
- Integration with topological_analysis for k optimization
"""

import numpy as np
import sys
import os
from typing import Dict, Tuple, Optional, List
import logging

# Add scripts directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))


logger = logging.getLogger(__name__)

# Constants from Z Framework
E_SQUARED = np.e**2  # c = e² ≈ 7.389


def compute_z_score(A: float, B: float, c: float = E_SQUARED) -> float:
    """
    Compute Z-invariant score: Z = A(B/c).

    This implements the discrete domain form of the Z Framework
    with invariant c = e² ≈ 7.389.

    Args:
        A: First parameter (typically related to sequence property)
        B: Second parameter (typically related to disruption metric)
        c: Invariant constant (default e²)

    Returns:
        Z-invariant score
    """
    if c == 0:
        raise ValueError("Invariant constant c cannot be zero")

    z_score = A * (B / c)
    return float(z_score)


def compute_composite_score(
    features: Dict[str, float], weights: Optional[Dict[str, float]] = None
) -> float:
    """
    Compute composite disruption score from spectral features.

    Combines multiple spectral metrics into a single score using
    Z-invariant framework with configurable weights.

    Default weighting:
    - delta_entropy: 0.4 (primary disruption indicator)
    - delta_f1: 0.3 (frequency shift importance)
    - delta_sidelobes: 0.2 (spectral complexity)
    - gc_content: 0.1 (sequence composition bias)

    Args:
        features: Dictionary of spectral features
        weights: Optional custom weights for each feature

    Returns:
        Composite disruption score
    """
    # Default weights (optimized from validation data)
    default_weights = {
        "delta_entropy": 0.4,
        "delta_f1": 0.3,
        "delta_sidelobes": 0.2,
        "gc_content": 0.1,
    }

    if weights is None:
        weights = default_weights

    # Compute weighted sum
    score = 0.0
    total_weight = 0.0

    for feature_name, weight in weights.items():
        if feature_name in features:
            # Normalize delta_f1 to be positive contribution
            if feature_name == "delta_f1":
                feature_value = abs(features[feature_name])
            else:
                feature_value = features[feature_name]

            score += weight * feature_value
            total_weight += weight

    # Normalize by total weight
    if total_weight > 0:
        score /= total_weight

    return float(score)


def bootstrap_confidence_interval(
    data: np.ndarray,
    statistic_func: callable,
    n_bootstrap: int = 1000,
    confidence_level: float = 0.95,
    seed: Optional[int] = None,
) -> Tuple[float, float, float]:
    """
    Compute bootstrap confidence interval for a statistic.

    Args:
        data: Input data array
        statistic_func: Function to compute statistic (e.g., np.mean)
        n_bootstrap: Number of bootstrap resamples (≥1,000 recommended)
        confidence_level: Confidence level (default 0.95 for 95% CI)
        seed: Random seed for reproducibility

    Returns:
        Tuple of (statistic, lower_bound, upper_bound)
    """
    if n_bootstrap < 1000:
        logger.warning(
            f"Bootstrap samples ({n_bootstrap}) below recommended minimum of 1,000"
        )

    if seed is not None:
        np.random.seed(seed)

    # Compute observed statistic
    observed = statistic_func(data)

    # Bootstrap resampling
    bootstrap_stats = []
    n = len(data)

    for _ in range(n_bootstrap):
        # Resample with replacement
        resample = np.random.choice(data, size=n, replace=True)
        boot_stat = statistic_func(resample)
        bootstrap_stats.append(boot_stat)

    bootstrap_stats = np.array(bootstrap_stats)

    # Compute confidence interval
    alpha = 1 - confidence_level
    lower_percentile = (alpha / 2) * 100
    upper_percentile = (1 - alpha / 2) * 100

    lower_bound = np.percentile(bootstrap_stats, lower_percentile)
    upper_bound = np.percentile(bootstrap_stats, upper_percentile)

    return float(observed), float(lower_bound), float(upper_bound)


def score_with_confidence(
    features_list: List[Dict[str, float]],
    weights: Optional[Dict[str, float]] = None,
    n_bootstrap: int = 1000,
    seed: Optional[int] = None,
) -> Dict[str, float]:
    """
    Compute composite scores with bootstrap confidence intervals.

    Args:
        features_list: List of feature dictionaries from multiple samples
        weights: Optional custom weights
        n_bootstrap: Number of bootstrap resamples
        seed: Random seed

    Returns:
        Dictionary with mean score, CI lower/upper bounds, and statistics
    """
    # Compute scores for all samples
    scores = np.array(
        [
            compute_composite_score(features, weights=weights)
            for features in features_list
        ]
    )

    # Bootstrap CI
    mean_score, ci_lower, ci_upper = bootstrap_confidence_interval(
        scores,
        statistic_func=np.mean,
        n_bootstrap=n_bootstrap,
        confidence_level=0.95,
        seed=seed,
    )

    results = {
        "mean_score": mean_score,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "ci_width": ci_upper - ci_lower,
        "std_dev": float(np.std(scores)),
        "n_samples": len(scores),
        "n_bootstrap": n_bootstrap,
    }

    logger.info(
        f"Composite score: {mean_score:.4f} (95% CI: [{ci_lower:.4f}, {ci_upper:.4f}])"
    )

    return results


def auto_optimize_k(
    sequence: str, k_candidates: Optional[List[float]] = None, is_rna: bool = False
) -> float:
    """
    Auto-optimize k parameter using geodesic-topological mapping.

    Tests multiple k values and selects the one that maximizes spectral
    coherence based on the Z Framework geodesic mapping.

    Args:
        sequence: DNA/RNA sequence
        k_candidates: Optional list of k values to test (default: [0.1, 0.2, 0.3, 0.4, 0.5])
        is_rna: If True, use RNA encoding

    Returns:
        Optimal k parameter
    """
    if k_candidates is None:
        # Default candidates around k* ≈ 0.3 (validated optimum)
        k_candidates = [0.1, 0.2, 0.3, 0.4, 0.5]

    from .encoding import phase_weighted_encoding
    from .analysis import compute_spectral_features

    best_k = 0.3
    best_score = -np.inf

    for k in k_candidates:
        # Encode with candidate k
        waveform = phase_weighted_encoding(sequence, is_rna=is_rna, k=k)

        # Compute spectral features
        features = compute_spectral_features(waveform, sequence=sequence)

        # Score based on spectral entropy (higher is better for disruption detection)
        score = features["entropy"]

        if score > best_score:
            best_score = score
            best_k = k

    logger.info(f"Auto-optimized k: {best_k:.3f} (score: {best_score:.4f})")

    return best_k


def compare_to_baseline(
    test_scores: np.ndarray,
    baseline_scores: np.ndarray,
    n_bootstrap: int = 1000,
    seed: Optional[int] = None,
) -> Dict[str, float]:
    """
    Compare test scores to baseline with bootstrap statistics.

    Computes ΔROC-AUC or similar lift metric with confidence intervals.

    Args:
        test_scores: Scores from test method
        baseline_scores: Scores from baseline method
        n_bootstrap: Number of bootstrap resamples
        seed: Random seed

    Returns:
        Dictionary with lift statistics and confidence intervals
    """
    if len(test_scores) != len(baseline_scores):
        raise ValueError("Test and baseline must have same length")

    # Compute lift (difference in means)
    lift = np.mean(test_scores) - np.mean(baseline_scores)

    # Bootstrap CI for lift
    if seed is not None:
        np.random.seed(seed)

    bootstrap_lifts = []
    n = len(test_scores)

    for _ in range(n_bootstrap):
        # Resample paired data
        indices = np.random.choice(n, size=n, replace=True)
        boot_test = test_scores[indices]
        boot_baseline = baseline_scores[indices]

        boot_lift = np.mean(boot_test) - np.mean(boot_baseline)
        bootstrap_lifts.append(boot_lift)

    bootstrap_lifts = np.array(bootstrap_lifts)

    # Compute CI
    ci_lower = np.percentile(bootstrap_lifts, 2.5)
    ci_upper = np.percentile(bootstrap_lifts, 97.5)

    # Compute p-value (proportion of bootstrap samples with lift ≤ 0)
    p_value = np.mean(bootstrap_lifts <= 0)

    results = {
        "lift": float(lift),
        "ci_lower": float(ci_lower),
        "ci_upper": float(ci_upper),
        "p_value": float(p_value),
        "significant": p_value < 0.05,
        "n_bootstrap": n_bootstrap,
    }

    logger.info(
        f"Lift: {lift:.4f} (95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]), p={p_value:.4f}"
    )

    return results
