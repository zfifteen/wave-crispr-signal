"""
Spectral Disruption Profiler Falsification Experiment

This package implements a comprehensive falsification experiment testing whether
phase-weighted FFT analysis (θ′(n,k) with k* ≈ 0.300) provides significant
enhancement to CRISPR gRNA efficiency prediction compared to unweighted baseline.

Modules:
    spectral_disruption_profiler: Main falsification experiment implementation
    smoke_test: Fast CI validation tests

Scientific Gates:
    - Human DNA only (A/C/G/T/N)
    - Z Framework invariants: Z = A(B/e²)
    - Geometric resolution: θ′(n,k) = φ·((n mod φ)/φ)^k with k ≈ 0.3
    - Statistical rigor: Bootstrap CI (≥1,000), permutation tests, FDR correction
    - Reproducibility: Fixed seed, versioned code, metadata persistence

Experiment ID: spectral_disruption_profiler_137
Version: 1.0
Author: Z Framework Falsification Team
Date: 2025-11-17
"""

__version__ = "1.0"
__experiment_id__ = "spectral_disruption_profiler_137"
__author__ = "Z Framework Falsification Team"
__date__ = "2025-11-17"
