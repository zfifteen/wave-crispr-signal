# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Note:** For deep research on this project using Claude Sonnet 4, see the comprehensive instructions in `.claude/claude_instructions_crispr_wave_signal.md`.

## Repository Overview

This project implements signal-theoretic analysis of DNA mutations using a complex-valued spectral framework called the "Z Framework". The codebase combines bioinformatics (CRISPR guide design, DNA sequence analysis), mathematical frameworks (high-precision calculations with golden ratio convergence), and medical applications (focused ultrasound, MRI analysis, pain management).

**Key Technologies**: Python 3.12+, mpmath (high-precision arithmetic), NumPy/SciPy, BioPython, PyTorch

## Common Build and Test Commands

### Installation
```bash
pip install -r requirements.txt
```

### Testing
```bash
# Run full test suite
python run_tests.py

# Run smoke tests (CI-friendly, <5s)
make smoke

# Run specific smoke tests
make mve-smoke           # Focused ultrasound MVE test
make mri-z5d-smoke      # MRI Z5D analysis test
make fus-enhancer-smoke # FUS enhancer test
make test-z-framework-import  # Z Framework import test

# Run individual test files
python -m pytest tests/test_z_framework.py -v
python -m pytest tests/test_crispr_guide_designer.py -v
```

### Quick Validation
```bash
# Quick validation demo (~2 min)
python proof_pack/quick_validation_demo.py

# Full validation suite (bootstraps + permutation tests)
python proof_pack/run_validation.py

# GC-Quartile resonance test
python bin/bin_resonance_test.py \
  --input data/kim2025.csv \
  --output results/gc_resonance.csv \
  --n_boot 4000 --n_perm 20000 --tail two
```

### CRISPR Guide Design
```bash
# Design guides for a sequence
python applications/crispr_cli.py design "ATGCTGCGGA..." -n 5 -o guides.json

# Design guides from FASTA file
python applications/crispr_cli.py design target.fasta -o results.json

# Design guides with physical Z-metrics
python applications/crispr_cli.py design target.fasta --physical-z-metrics --format csv -o guides.csv

# Score a specific guide
python applications/crispr_cli.py score "GACGATCGATCGATCGATCG" --target-context "AATGACGATCG..."

# Batch score guides
python applications/crispr_cli.py batch-score guides.txt -o scores.csv -f csv
```

### Experiment Execution
```bash
# Run MRI Z5D analysis with DICOM data
make run-mri-z5d           # Default parameters
make run-mri-z5d-quick     # Reduced parameters (faster)

# Run focused ultrasound MVE
make run-mve               # Default parameters
make run-mve-quick         # Quick test

# Run FUS enhancer (optimized, 10^6 trials)
make run-fus-enhancer
make run-fus-enhancer-quick
```

## Architecture Overview

### Core Framework Components

**Z Framework** (`scripts/z_framework.py`):
- Implements discrete domain form: Z = n(Δₙ/Δₘₐₓ)
- High-precision calculations using mpmath (dps=50)
- DNA sequence mapping (A=1, T=2, C=3, G=4)
- Geodesic resolution function: θ'(n,k) = φ·((n mod φ)/φ)^k where k ≈ 0.3
- Golden ratio convergence (φ-1 ≈ 0.618) with target variance σ ≈ 0.113

**Key Mathematical Constants**:
- PHI (φ) = 1.618033988749894848... (golden ratio)
- PHI_CONJUGATE = φ-1 ≈ 0.618... (convergence target)
- E_SQUARED = e² ≈ 7.389 (used in discrete domain calculations)
- TARGET_VARIANCE = 0.113 (biological invariant target)

### Directory Structure

```
wave-crispr-signal/
├── applications/           # CRISPR guide design, metrics, visualization
│   ├── crispr_cli.py      # Command-line interface
│   ├── crispr_guide_designer.py
│   ├── crispr_physical_z_metrics.py  # Physical Z-metrics calculator
│   └── wave_crispr_metrics.py
├── scripts/               # Core scientific modules (importable)
│   ├── z_framework.py     # Z Framework calculator (main framework)
│   ├── topological_analysis.py
│   └── invariant_features.py
├── modules/               # Specialized analysis modules
│   ├── bio_v_arbitrary.py
│   ├── dna_storage_hypothesis.py
│   └── molecular_dynamics_framework.py
├── experiments/           # Self-contained experiments
│   ├── focused_ultrasound_mve.py
│   ├── mri_z5d_analysis.py
│   ├── pain_management_application.py
│   └── signal_theoretic_crispr/
├── proof_pack/            # Validation and verification tools
│   ├── quick_validation_demo.py
│   ├── run_validation.py
│   └── generate_synthetic_data.py
├── tests/                 # Test infrastructure
├── data/                  # Research datasets (FASTA, CSV, DICOM)
│   ├── DICOM/            # Medical imaging data
│   ├── fasta/            # DNA sequences
│   └── *.csv             # Experimental datasets
├── docs/                  # Comprehensive documentation
└── notebooks/             # Jupyter analysis notebooks
```

### Import Patterns

**Core Framework Import**:
```python
from scripts.z_framework import ZFrameworkCalculator
from scripts.topological_analysis import topological_analysis_function
from scripts.invariant_features import calculate_invariant_features
```

**Application Import**:
```python
from applications.crispr_guide_designer import CRISPRGuideDesigner
from applications.crispr_physical_z_metrics import PhysicalZMetricsCalculator
```

**Module Import**:
```python
from modules.bio_v_arbitrary import BioAnalyzer
from modules.molecular_dynamics_framework import MDSimulator
```

**Experiments import from scripts/modules**:
```python
# In experiments/, always import from core modules
from scripts.z_framework import ZFrameworkCalculator
from modules.bio_v_arbitrary import BioAnalyzer
```

### Key Design Principles

1. **Scientific Rigor**: All mathematical claims require empirical validation with bootstrap confidence intervals (≥1,000 resamples) and permutation tests
2. **High Precision**: Use mpmath with dps=50 for all Z Framework calculations
3. **Reproducibility**: Fixed seeds (42), exact version pinning, environment persistence
4. **Human DNA Only**: GRCh38/hg38 reference genome, strict nucleotide validation (A/C/G/T/N for DNA, A/C/G/U/N for RNA)
5. **Fail-Fast Validation**: Start every pipeline with nucleotide-only checks

## Scientific Gates (MANDATORY)

These gates are enforced in `.github/copilot-instructions.md` and must be followed:

### 1. Data Validation
- **Human DNA only**: GRCh38/hg38 reference
- **DNA sequences**: A/C/G/T/N only (reject U, IUPAC codes)
- **RNA sequences**: A/C/G/U/N only (reject T, IUPAC codes)
- **Fail-fast**: Raise clear `ValueError` on validation failure

### 2. Z Framework Domain Constraints
- **Discrete/biological domain (DEFAULT)**: Z = A(B / e²)
- **Geodesic resolution**: θ'(n,k) = φ·((n mod φ)/φ)^k with k ≈ 0.3
- **Document deviations**: Any changes to default k must be justified

### 3. Statistical Validity
- **Pre-register endpoints**: Pearson r with 95% bootstrap CI (≥1,000 resamples)
- **Null models**: ≥1,000× permutation/shuffle for empirical p-values
- **Leakage control**: Split-by-screen, split-by-gene (no entity in both train/eval)
- **Multiple comparisons**: Apply Benjamini-Hochberg FDR when comparing ≥3 metrics

### 4. Reproducibility Requirements
- **CLI contract**: `--seed`, `--bootstrap`, `--permutation`, `--splits`, `--domain`, `--pam`
- **Persist metadata**: seed, git commit, dataset name/version, SHA256, runtime, pip freeze
- **Precision**: mpmath with mp.dps = 50 for high-precision calculations
- **Pinned environment**: Python 3.12.*, exact version pins in `requirements.txt`

### 5. CI/CD Requirements
- **Smoke tests**: Must complete <5s (use `make smoke`)
- **CI runs**: `pytest -q` and `make smoke`
- **Experiment structure**: Each experiment needs `README.md` + `manifest.yml`

## Working with Tests

### Test Organization
- **Unit tests**: `tests/test_*.py` (test individual modules)
- **Smoke tests**: Fast validation tests (<5s total)
- **Integration tests**: Full validation suite in `proof_pack/`
- **Performance tests**: Benchmark tests in experiments

### Adding New Tests
1. Create `tests/test_<module>.py`
2. Add smoke test target to `Makefile` if needed
3. Ensure tests work with pinned dependencies
4. Include both unit and integration validation

### Test Execution Pattern
```python
# Standard test structure
import pytest
from scripts.z_framework import ZFrameworkCalculator

def test_z_framework_basic():
    calc = ZFrameworkCalculator(precision_dps=50)
    sequence = "ATCGATCG"
    results = calc.calculate_z_values(sequence)

    assert "z_mean" in results
    assert results["converges_to_phi_conjugate"] is not None
```

## File Naming Conventions

- **Python modules**: `snake_case.py`
- **Domain prefixes**: `crispr_*.py`, `wave_*.py`, `z5d_*.py`
- **Test files**: `test_*.py`
- **Major docs**: `UPPERCASE_WITH_UNDERSCORES.md`
- **Datasets**: `lowercase_with_underscores.csv`

## Dependency Management

All dependencies use **exact version pinning** in `requirements.txt`:

**Core Scientific**:
- `mpmath==1.3.0` (high-precision mathematics)
- `numpy==1.26.4` (numerical computing)
- `scipy==1.16.1` (scientific computing)
- `sympy==1.13.1` (symbolic mathematics)

**Bioinformatics**:
- `biopython==1.83` (biological sequence analysis)
- `pyfaidx==0.8.1.1` (FASTA indexing)

**Machine Learning**:
- `scikit-learn==1.5.1` (ML algorithms)
- `torch==2.8.0` (vectorized optimization)

**Visualization**:
- `matplotlib==3.10.5`
- `plotly==6.3.0`
- `seaborn==0.13.2`

**Medical Imaging**:
- `pydicom==3.0.1` (DICOM file handling)

## Special Considerations

### Working with DICOM Data
MRI/medical imaging data is in `data/DICOM/`. Use PyDICOM for reading:
```python
import pydicom
from experiments.mri_z5d_analysis import load_dicom_signals

signals = load_dicom_signals(signal_length=256, max_files_per_series=20)
```

### High-Precision Calculations
Always configure mpmath precision before calculations:
```python
import mpmath as mp
mp.dps = 50  # 50 decimal places

# Use mp.mpf for high-precision values
phi = mp.mpf("1.618033988749894848204586834365638117720309179805762862135")
```

### JSON Serialization
mpmath and numpy types need special handling for JSON:
```python
def format_mpmath_for_json(value):
    if isinstance(value, mp.mpf):
        return float(value)
    elif isinstance(value, np.integer):
        return int(value)
    elif isinstance(value, np.floating):
        return float(value)
    elif isinstance(value, np.ndarray):
        return value.tolist()
    # ... handle lists, dicts recursively
```

### Geometric Resolution (Critical Parameter)
The k parameter in θ'(n,k) is optimized across datasets:
- **Default**: k ≈ 0.3 for 21-nt CRISPR guides (φ = 21)
- **MRI/Z5D**: k ≈ 0.04449 for different domain
- **Document deviations**: Any k ≠ 0.3 must be justified

### Repository Policy Compliance
The repository follows `.github/REPOSITORY_POLICY.md`. Key points:
- Keep to established directory structure
- Every public component must have documentation
- Tests required for all core modules
- Exact version pinning in `requirements.txt`

## Common Development Patterns

### Creating a New Experiment
1. Create `experiments/<experiment_name>.py`
2. Import from `scripts/` and `modules/`, never hardcode logic
3. Add CLI arguments: `--seed`, `--bootstrap`, `--permutation`, etc.
4. Create `experiments/<experiment_name>_README.md`
5. Add smoke test to `Makefile`
6. Persist metadata to `results/<experiment_id>/`

### Adding CRISPR Analysis Features
1. Core calculations go in `applications/wave_crispr_metrics.py`
2. CLI interface updates in `applications/crispr_cli.py`
3. Physical Z-metrics in `applications/crispr_physical_z_metrics.py`
4. Tests in `tests/test_crispr_*.py`

### Working with Validation Data
Synthetic validation data lives in `proof_pack/`:
- `neural_spikes.csv` - Neural activity patterns
- `nav1.8_panel.csv` - Ion channel panel data
- `bcl11a_edits.csv` - BCL11A editing data

Generate new synthetic data:
```bash
python proof_pack/generate_synthetic_data.py
```

## Linting and Code Quality

**Configured via `pyproject.toml`**:
- Line length: 88 characters
- Target: Python 3.12
- Ruff linter with E402 allowed (module-level imports for path setup)

## Running the Full Validation Stack

```bash
# 1. Quick validation (2 min)
python proof_pack/quick_validation_demo.py

# 2. Full validation (bootstraps + permutation)
python proof_pack/run_validation.py

# 3. Smoke tests
make smoke

# 4. Full test suite
python run_tests.py

# 5. Specific experiments
make run-mri-z5d-quick
make run-mve-quick
```

## Important Notes

- **RESEARCH USE ONLY**: Not for clinical applications
- **No fabrication**: Never derive DNA from protein/RNA or fabricate nucleotides
- **Fail-fast**: Start pipelines with validation checks
- **Document everything**: Every major component needs documentation
- **Version control**: Commit includes git hash, dataset version, SHA256 in results
