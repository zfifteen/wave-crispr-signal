# Experiments

This directory contains specific experiments and scripts that import from `modules/` to demonstrate various applications and validations of the wave-crispr-signal framework.

## Experiment Overview

- **experimental_validation.py** - Comprehensive experimental validation framework
- **corrected_validation.py** - Corrected validation experiments 
- **invariant_validation.py** - Validation of invariant features
- **enhanced_bio_analysis.py** - Enhanced biological sequence analysis
- **falsification_experiments.py** - Hypothesis falsification experiments
- **pain_management_application.py** - Pain management application experiments
- **z5d_prime_predictor_experiment.py** - Z5D predictor performance experiments
- **validate_k_parameter_falsification.py** - K parameter validation experiments

## Design Principles

All experiments in this directory follow these principles:

1. **Import from modules/** - Never hardcode logic; always import from the core modules
2. **Self-contained** - Each experiment should be runnable independently
3. **Documented** - Each experiment includes comprehensive documentation of methodology
4. **Reproducible** - Results should be reproducible with the pinned dependencies

## Usage

Run individual experiments directly:

```bash
python experiments/experimental_validation.py
python experiments/pain_management_application.py
```

## Dependencies

All experiments use the shared dependencies from `requirements.txt` and import functionality from the `modules/` package.