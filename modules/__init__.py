"""
Core Modules Package

This package contains the core scientific modules for the wave-crispr-signal project,
including Z Framework implementation, topological analysis, invariant features,
and molecular dynamics framework.
"""

# Import key classes and functions explicitly to avoid * imports
from .z_framework import (
    ZFrameworkCalculator,
    format_mpmath_for_json,
    format_mpmath_for_display,
    TARGET_VARIANCE,
)
from .topological_analysis import TopologicalAnalyzer
from .invariant_features import (
    InvariantFeatureSet,
    ZetaUnfoldCalculator,
    PhaseAwareSpectralAnalyzer,
    GoldenProximityCalculator,
    CurvatureDisruptionAnalyzer,
)
from .molecular_dynamics_framework import MolecularDynamicsZFramework
from .dna_storage_hypothesis import DNAStorageHypothesis
from .bio_v_arbitrary import EmpiricalValidator, ArbitraryEncoder

__version__ = "1.0.0"
__all__ = [
    "ZFrameworkCalculator",
    "format_mpmath_for_json", 
    "format_mpmath_for_display",
    "TARGET_VARIANCE",
    "TopologicalAnalyzer",
    "InvariantFeatureSet",
    "ZetaUnfoldCalculator",
    "PhaseAwareSpectralAnalyzer",
    "GoldenProximityCalculator", 
    "CurvatureDisruptionAnalyzer",
    "MolecularDynamicsZFramework",
    "DNAStorageHypothesis",
    "EmpiricalValidator",
    "ArbitraryEncoder",
]