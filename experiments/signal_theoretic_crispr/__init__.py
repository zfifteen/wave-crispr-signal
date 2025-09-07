"""
Signal-Theoretic Modeling of CRISPR-Cas9 for Cancer-Target Editing

This package implements a two-arm benchmark experiment comparing:
- Arm A: Batch-scan baseline (Spring Batch style)
- Arm B: WAVE-CRISPR spectral-geodesic features

The experiment enforces strict scientific gates for human-only FASTA validation,
reproducibility, and statistical rigor.
"""

__version__ = "1.0.0"
__author__ = "wave-crispr-signal"

from .validation import HumanFASTAValidator, DeterminismValidator
from .baseline import BaselinePipeline, BaselineModels
from .spectral import SpectralFeatureExtractor, SpectralModels
from .statistics import ExperimentalStatistics
from .main import SignalTheoreticExperiment

__all__ = [
    'HumanFASTAValidator',
    'DeterminismValidator', 
    'BaselinePipeline',
    'BaselineModels',
    'SpectralFeatureExtractor',
    'SpectralModels',
    'ExperimentalStatistics',
    'SignalTheoreticExperiment'
]