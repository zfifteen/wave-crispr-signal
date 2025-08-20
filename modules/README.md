# Core Modules

This directory contains the core reusable modules for the wave-crispr-signal project. These modules provide the fundamental algorithms, logic, and model components that are used throughout the project.

## Module Overview

- **z_framework.py** - Core Z Framework implementation for signal processing analysis
- **topological_analysis.py** - Topological analysis functions and geometric computations  
- **invariant_features.py** - Invariant feature calculations and transformations
- **molecular_dynamics_framework.py** - Molecular dynamics simulation tools and integrations
- **dna_storage_hypothesis.py** - DNA storage and encoding hypothesis implementations
- **bio_v_arbitrary.py** - Biological vs arbitrary sequence analysis tools

## Usage

These modules are designed to be imported by experiments, tools, and scripts:

```python
from modules.z_framework import ZFrameworkCore
from modules.topological_analysis import compute_topological_features
from modules.invariant_features import calculate_invariants
```

## Dependencies

All modules require the dependencies specified in the root `requirements.txt` file.

## Testing

Tests for these modules are located in the `tests/` directory with corresponding `test_*.py` files.