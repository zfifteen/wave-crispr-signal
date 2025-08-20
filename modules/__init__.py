"""
Modules Package

This package contains specialized modules for the wave-crispr-signal project,
including molecular dynamics framework, DNA storage hypothesis, and biological validators.

Note: Core modules (z_framework, topological_analysis, invariant_features) have been 
moved to root level for direct import as per repository policy.
"""

# Import specialized modules that remain in the modules directory
from .molecular_dynamics_framework import MolecularDynamicsZFramework
from .dna_storage_hypothesis import DNAStorageHypothesis
from .bio_v_arbitrary import EmpiricalValidator, ArbitraryEncoder

__version__ = "1.0.0"
__all__ = [
    "MolecularDynamicsZFramework",
    "DNAStorageHypothesis", 
    "EmpiricalValidator",
    "ArbitraryEncoder",
]
