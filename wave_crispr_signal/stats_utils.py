"""Shared lightweight statistical helpers for active codepaths."""

from __future__ import annotations

from typing import Callable, Iterable, Tuple
import numpy as np


def bootstrap_ci(
    values: Iterable[float],
    func: Callable[[np.ndarray], float] = np.mean,
    n_boot: int = 1000,
    ci: float = 95.0,
    seed: int | None = None,
) -> Tuple[float, float]:
    """Return bootstrap confidence interval for a statistic."""
    x = np.asarray(list(values), dtype=float)
    if x.size == 0:
        return (0.0, 0.0)

    rng = np.random.default_rng(seed)
    stats = np.empty(n_boot, dtype=float)
    n = x.size
    for i in range(n_boot):
        sample = x[rng.integers(0, n, size=n)]
        stats[i] = float(func(sample))

    alpha = (100.0 - ci) / 2.0
    lo, hi = np.percentile(stats, [alpha, 100.0 - alpha])
    return (float(lo), float(hi))
