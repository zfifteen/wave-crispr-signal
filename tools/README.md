# Tools

This directory contains utility scripts for conversion, preprocessing, postprocessing, demonstrations, and other support functions.

## Tool Overview

- **demo_dna_storage_hypothesis.py** - Demonstration of DNA storage hypothesis capabilities
- **demo_invariant_features.py** - Interactive demonstration of invariant feature calculations
- **demo_topological_analysis.py** - Topological analysis demonstration and visualization
- **pain_management_demo.py** - Pain management analysis demonstration
- **quick_validation_demo.py** - Quick validation and testing demonstration
- **check_repo_policy.py** - Repository policy compliance checker

## Categories

### Demonstration Tools
Scripts that showcase functionality and provide interactive examples of the core modules.

### Validation Tools  
Quick validation and testing utilities for rapid verification of functionality.

### Preprocessing Tools
Data preparation and conversion utilities (to be added as needed).

### Postprocessing Tools
Result analysis and visualization utilities (to be added as needed).

## Usage

Run tools directly for demonstrations:

```bash
python tools/demo_dna_storage_hypothesis.py
python tools/quick_validation_demo.py
python tools/check_repo_policy.py
```

## Design Principles

- Tools should be focused on specific utility functions
- Each tool should be self-contained and runnable independently
- Tools import functionality from `modules/` rather than implementing core logic
- Tools should include clear usage instructions and examples