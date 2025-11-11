# Claude Sonnet 4 Instructions: Deep Research on WAVE-CRISPR

**Date Prepared:** October 24, 2025  
**Prepared For:** Claude Sonnet 4 (built by Anthropic)  
**Purpose:** This document serves as a complete knowledge bootstrap for Claude Sonnet 4 to conduct deep research on the WAVE-CRISPR project. It synthesizes all project documentation into a structured, actionable guide for scientific analysis, experimental validation, and theoretical extensions.

---

## 1. Project Overview

### 1.1 WAVE-CRISPR Framework

WAVE-CRISPR is a signal-theoretic framework for analyzing DNA sequences in CRISPR applications. It treats DNA as complex-valued waveforms, enabling spectral analysis (e.g., FFT for entropy, sidelobes, harmonic power) to predict CRISPR guide efficiency, off-target risks, and repair outcomes.

**Core Innovation:** Integration of the Z Framework with geodesic mappings to create a unified mathematical foundation for biological sequence analysis.

#### Key Mathematical Foundations

**Z Framework (Multi-Domain):**
- **Physical Domain:** `Z = T(v/c)` with causality constraint `|v| < c`
- **Discrete/Biological Domain (DEFAULT):** `Z = n(Δₙ/Δₘₐₓ)` where:
  - `n`: Frame-dependent parameter (e.g., sequence entropy, position)
  - `Δₙ`: Rate/shift (e.g., spectral mutation shift, nucleotide difference)
  - `Δₘₐₓ`: Invariant normalizer (maximum delta in sequence)
  - `κ(n) = d(n) · ln(n+1)/e²` (curvature term)
- **Geometric Domain:** `θ'(n, k) = φ · ((n mod φ)/φ)^k`
  - `φ = 1.618033988749...` (golden ratio)
  - Optimal `k* ≈ 0.3` for ~21-nt CRISPR guides (φ = 21)
  - Achieves ~15% density enhancement (CI [14.6%, 15.4%])
- **Biological Encoding:** Complex mappings for nucleotides
  - Standard: `A=1+0j, T=-1+0j, C=0+1j, G=0-1j`
  - Alternative: DNA breathing dynamics (AT: 10 MHz, GC: 1 GHz)

**Key Empirical Claims:**
- Spectral features improve CRISPR prediction by 15-25% (R², CI [14.6%, 15.4%])
- Pearson r ≥ 0.5 (p < 10⁻⁵) for spectral entropy correlations
- Optimal k* ≈ 0.3 validated across 6 public datasets (Doench 2016, Kim 2025, Patch 2024)
- **Surprising finding:** Arbitrary encodings sometimes outperform biologically-anchored ones

### 1.2 Repository Structure & Key Components

```
wave-crispr-signal/
├── scripts/                    # Core importable scientific modules
│   ├── z_framework.py         # Z Framework calculator (main framework)
│   ├── topological_analysis.py
│   └── invariant_features.py
├── applications/              # CRISPR guide design & analysis
│   ├── crispr_cli.py          # Command-line interface
│   ├── crispr_guide_designer.py
│   ├── crispr_physical_z_metrics.py
│   └── wave_crispr_metrics.py
├── modules/                   # Specialized analysis modules
│   ├── bio_v_arbitrary.py     # Biological vs. arbitrary encoding tests
│   ├── dna_storage_hypothesis.py
│   └── molecular_dynamics_framework.py
├── experiments/               # Self-contained experiments
│   ├── focused_ultrasound_mve.py
│   ├── mri_z5d_analysis.py
│   ├── pain_management_application.py
│   ├── test_breathing_dynamics_encoding.py
│   └── signal_theoretic_crispr/
├── proof_pack/                # Validation & verification
│   ├── quick_validation_demo.py
│   ├── run_validation.py
│   └── generate_synthetic_data.py
├── tests/                     # Test infrastructure
├── data/                      # Research datasets
│   ├── DICOM/                # Medical imaging (MRI)
│   ├── fasta/                # DNA sequences (GRCh38/hg38)
│   └── *.csv                 # Experimental datasets
├── docs/                      # Comprehensive documentation
│   ├── Z_FRAMEWORK.md
│   ├── TOPOLOGICAL_ANALYSIS.md
│   ├── METHOD_DETAILS.md
│   └── EMPIRICAL_FINDINGS_REPORT.md
└── notebooks/                 # Jupyter analysis notebooks
```

### 1.3 Related Projects & Resources

- **GitHub Repository:** https://github.com/zfifteen/wave-crispr-signal (ID: 1027305912)
- **Unified Framework:** https://github.com/zfifteen/unified-framework (geodesic integrations, pre-computed zeta zeros)
- **SpectralTE:** Extension for mRNA translation efficiency prediction (15-25% R² improvement)
- **Key Datasets:**
  - BioGRID-ORCS Homo sapiens v1.1.17 (default CRISPR dataset with citation)
  - Doench 2016, Kim 2025, Patch 2024 (validation datasets)
  - GRCh38/hg38 (human genome reference)

---

## 2. Scientific Methodology

### 2.1 DNA Encoding Strategies

**Biological Anchored Encodings:**
- **Polarizability-based:** A=1.65, T=1.52, C=1.27, G=1.58 (Å³)
- **Hydrogen bonding:** A=2, T=2, C=3, G=3 (H-bonds)
- **Thermodynamics:** Based on ΔG values
- **DNA Breathing Dynamics:**
  - AT pairs: 10 MHz (fast opening) → weight -10.00 + 3.00j
  - GC pairs: 1 GHz (slow opening) → weight +10.00 - 3.00j
  - Includes helical periodicity (10.5 bp/turn)

**Formula:** `Complex_value = (property₁/max) + i*(property₂/max)`

**Arbitrary Encodings:**
- Random complex weights: amplitude [0.5-2.0], phase [0-2π]
- Use `seed=42` for reproducibility
- **Empirical finding:** Often outperforms biological encodings in practice

**Waveform Construction:**
```
Ψₙ = encoded_sequence * exp(2πi * cumulative_positions)
positions = cumsum(d=0.34 nm base spacing)
```

### 2.2 Spectral Feature Extraction

**Primary FFT-Based Metrics:**
1. **Spectral Entropy** (Shannon entropy of normalized magnitudes)
   - Lower entropy → Higher guide efficiency
2. **Sidelobe Count** (significant peaks above threshold)
3. **Harmonic Power** (sum of harmonic energies)
4. **Spectral Centroid** (center of mass of spectrum)
5. **Spectral Rolloff** (frequency below which X% of energy lies)
6. **Spectral Flatness** (Wiener entropy, noisiness measure)
7. **Zero-Crossing Rate** (sign changes in signal)
8. **Peak Count** (local maxima in spectrum)

**Disruption Score (for mutations):**
- Compare pre/post-mutation spectra using KL divergence or correlation distance
- Composite score: weighted sum of spectral deltas

### 2.3 Z Framework Calculations

**Discrete Z Computation:**
```python
# Map nucleotides to integers: A=1, T=2, C=3, G=4
Δᵢ = xᵢ - xᵢ₋₁
Zᵢ = i * (|Δᵢ| / Δₘₐₓ)
```

**Geodesic Curvature:**
```python
F = k * (ratio^k)
ratio = (D/E)/e, where D=c/a, E=c/b
```

**Geometric Resolution (Critical):**
```python
θ'(n, k) = φ · ((n mod φ)/φ)^k
# Default: k ≈ 0.3 for 21-nt guides (φ = 21)
# MRI/Z5D: k ≈ 0.04449 (different domain)
```

**Key Constants (High Precision):**
```python
import mpmath as mp
mp.dps = 50  # 50 decimal places required

PHI = mp.mpf("1.618033988749894848204586834365638117720309179805762862135")
PHI_CONJUGATE = PHI - 1  # ≈ 0.618... (convergence target)
E_SQUARED = mp.mpf("7.389056098930650227230427460575007813180315570551847324087")
TARGET_VARIANCE = mp.mpf("0.113")  # Biological invariant
```

**Convergence Criteria:**
- Target: `|μ_Z - φ⁻¹| ≤ 0.005`
- Variance trimming: `σ²_trim = σ²_Z - (κ * 0.013) ≈ 0.113`

### 2.4 CRISPR-Specific Applications

**Guide Design Pipeline:**
```python
# Standard parameters
PAM_sequence = "NGG"  # SpCas9
guide_length = 20
GC_optimal = 50%  # Target GC content

# Scoring formula
score = 0.4 * entropy_factor + 0.4 * sidelobe_factor + 0.2 * GC_factor
```

**Off-Target Risk Assessment:**
```python
risk_score = 0.6 * spectral_correlation + 0.4 * hamming_distance
```

**Repair Outcome Prediction:**
- Low entropy → NHEJ (65% probability)
- High sidelobes → MMEJ (25% probability)
- High stability → HDR (10% probability)

**Comprehensive Scoring:**
- 35% on-target efficiency
- 25% spectral quality
- 20% off-target risk
- 20% repair pathway bias

### 2.5 Statistical Validation Requirements

**Pre-Registration (MANDATORY):**
1. **Primary Endpoint:** Pearson r with 95% bootstrap CI (≥1,000 resamples)
2. **Partial Correlation:** Control for GC%, guide length, guide position
3. **Effect Size:** Cohen's d with CI for defined contrasts

**Null Models:**
- ≥1,000× permutation/shuffle tests for empirical p-values
- Label permutation for classification tasks

**Leakage Control:**
- **Split-by-screen:** No screen in both train/eval
- **Split-by-gene:** No gene in both train/eval
- Never split randomly when biological units exist

**Multiple Comparisons:**
- Apply Benjamini-Hochberg FDR correction when comparing ≥3 metrics
- Report both raw and corrected p-values

**Power Analysis:**
- Justify sample size or cite prior power analysis
- Target 80% power for detecting 15% effect

---

## 3. Scientific Gates (MANDATORY COMPLIANCE)

### 3.1 Data Validation Gates

**Human DNA Only:**
- Reference genome: GRCh38/hg38
- For **DNA sequences:** Validate A/C/G/T/N only (case-insensitive)
  - **REJECT:** U (RNA base) and IUPAC ambiguity codes beyond N
- For **RNA sequences:** Validate A/C/G/U/N only (case-insensitive)
  - **REJECT:** T (DNA base) and IUPAC ambiguity codes beyond N

**Fail-Fast Validation:**
```python
def validate_dna_sequence(seq: str) -> None:
    """Validate DNA sequence contains only A/C/G/T/N."""
    seq_upper = seq.upper()
    valid_bases = set('ACGTN')
    invalid_bases = set(seq_upper) - valid_bases
    
    if invalid_bases:
        raise ValueError(
            f"Invalid DNA bases detected: {invalid_bases}. "
            f"DNA sequences must contain only A/C/G/T/N. "
            f"Found 'U' or IUPAC codes: {seq_upper}"
        )
```

**No Fabrication:**
- Never derive DNA from protein or RNA sequences
- Never fabricate nucleotides
- Always use curated human FASTA files

### 3.2 Z Framework Domain Gates

**Default Domain (Discrete/Biological):**
- Formula: `Z = A(B / e²)` where A=n, B=Δₙ
- Guard divide-by-zero: check `Δₘₐₓ > 0` or `e² > 0`
- Document A and B clearly

**Geometric Resolution:**
- Default: `θ'(n,k) = φ·((n mod φ)/φ)^k` with `k ≈ 0.3`
- **Deviation Rule:** Any k ≠ 0.3 MUST be documented with scientific justification
- For 21-nt guides: φ = 21 (geometric period)

### 3.3 Dataset Provenance Gates

**Required Metadata:**
1. Dataset name and version
2. Download URL or DOI
3. License information
4. Taxonomy (verify Homo sapiens)
5. SHA256 hash of local file

**Runtime Requirements:**
- Print dataset version at start of analysis
- Include BioGRID-ORCS citation when used (default dataset)

**Example Citation:**
```
Oughtred R, et al. (2024). The BioGRID database: A comprehensive 
biomedical resource of curated protein, genetic, and chemical interactions. 
Protein Sci. 30(1):187-200.
```

### 3.4 Reproducibility Gates

**CLI Contract (Required Flags):**
```bash
--seed <int>           # Random seed for reproducibility
--bootstrap <int>      # Number of bootstrap resamples (≥1,000)
--permutation <int>    # Number of permutations (≥1,000)
--splits <int>         # Cross-validation splits
--domain <str>         # Z Framework domain (discrete|physical|geometric)
--pam <str>           # PAM sequence for CRISPR
```

**Metadata Persistence:**
Every results directory must contain:
```
results/<experiment_id>/run-YYYYMMDD-HHMMSS/
├── results.json          # Numerical results
├── results.csv           # Tabular data
├── env.txt              # pip freeze output
├── metadata.json        # Execution metadata
└── log.txt              # Execution log
```

**Metadata Contents:**
```json
{
  "seed": 42,
  "git_commit": "a1b2c3d",
  "dataset": {
    "name": "BioGRID-ORCS",
    "version": "v1.1.17",
    "sha256": "abc123...",
    "url": "https://..."
  },
  "runtime_seconds": 123.45,
  "timestamp": "2025-10-24T22:00:00Z",
  "environment": "Python 3.12.*, numpy 1.26.4, ..."
}
```

**High Precision Requirements:**
```python
import mpmath as mp
mp.dps = 50  # MANDATORY for Z Framework calculations
```

### 3.5 CI/CD Gates

**Smoke Tests (<5 seconds total):**
```bash
make smoke              # Run all smoke tests
make mve-smoke         # Focused ultrasound MVE
make mri-z5d-smoke     # MRI Z5D analysis
make fus-enhancer-smoke # FUS enhancer
```

**Full Test Suite:**
```bash
pytest -q              # Quick run
python run_tests.py    # Full validation
```

**Experiment Structure:**
Every experiment needs:
1. `experiments/<name>.py` (main script)
2. `experiments/<name>_README.md` (documentation)
3. `experiments/<name>_manifest.yml` (metadata)

**Manifest Template:**
```yaml
name: experiment_name
version: 1.0.0
description: Brief description
author: Name
date: YYYY-MM-DD
datasets:
  - name: dataset_name
    version: version
    url: https://...
parameters:
  seed: 42
  bootstrap: 1000
  permutation: 1000
dependencies:
  - numpy==1.26.4
  - scipy==1.16.1
```

---

## 4. Empirical Findings & Open Questions

### 4.1 Validated Results

**GC-Quartile Resonance (Kim 2025, N=18,102):**
- Quartile Q4: r = -0.211, p_perm = 0.0012 (FDR-corrected)
- Effect persists after controlling for confounds

**Spectral Disruption vs. RuleSet3:**
- ΔROC-AUC = +0.047 ± 0.006 (bootstrap 10,000×)
- Consistent advantage across datasets

**Geometric Resolution Optimum:**
- k* ≈ 0.300 validated across 6 datasets
- 15% spectral density enhancement (CI [14.6%, 15.4%])

### 4.2 Surprising Findings

**Arbitrary > Biological Encodings (Some Cases):**
```
Study: EMPIRICAL_FINDINGS_REPORT.md (n=80 CRISPR sequences)
- Bio-anchored (polarizability): r ≈ -0.198 (p=0.048)
- Arbitrary (random complex): r ≈ 0.052 (p=0.611)
- Effect size: Cohen's d > 0.5, p < 0.001
```

**Implications:**
- Z Framework may need refinement
- Arbitrary encodings capture structural patterns differently
- Suggests role for quantum effects or higher-order properties

**DNA Breathing Dynamics Encoding:**
```
Study: test_breathing_dynamics_encoding.py
- GC-affecting mutations: Cohen's d = +4.03, p < 0.001
- AT-affecting mutations: Z = 0.0000 (perfect selectivity)
- Confirms frequency-native encoding advantage
```

### 4.3 Open Research Questions

1. **Why do arbitrary encodings sometimes outperform biological ones?**
   - Hypothesis: Capture emergent structural patterns
   - Need: Quantum property integration tests

2. **Can we predict k* for arbitrary DNA lengths?**
   - Current: k ≈ 0.3 for 21-nt guides
   - Unknown: Scaling law for k(L)?

3. **What is the physical basis for φ-convergence?**
   - Theoretical: Golden ratio in DNA structure
   - Need: Crystallographic validation

4. **How do spectral features relate to 3D structure?**
   - Hypothesis: FFT captures helical symmetries
   - Need: Correlation with protein binding affinities

5. **Can Z Framework extend to RNA structures?**
   - Current: DNA-focused (A/C/G/T)
   - Need: RNA validation with A/C/G/U encoding

---

## 5. Experimental Procedures

### 5.1 Standard Validation Pipeline

**Step 1: Dataset Preparation**
```bash
# Verify dataset
python scripts/verify_dataset.py \
  --input data/kim2025.csv \
  --check-taxonomy \
  --check-nucleotides \
  --output-sha256

# Expected output: SHA256, species confirmation
```

**Step 2: Encoding Generation**
```python
# Biological encoding
bio_encoder = BiologicalEncoder(property="polarizability")
bio_waveforms = bio_encoder.encode_sequences(sequences)

# Arbitrary encoding
arb_encoder = ArbitraryEncoder(seed=42, n_trials=10)
arb_waveforms = arb_encoder.encode_sequences(sequences)
```

**Step 3: Feature Extraction**
```python
# Spectral features
features = extract_spectral_features(
    waveforms,
    methods=['entropy', 'sidelobes', 'harmonic_power']
)

# Z-metrics
z_calc = ZFrameworkCalculator(precision_dps=50)
z_metrics = z_calc.calculate_z_values(
    sequences,
    k=0.3,
    geometric_period=21
)
```

**Step 4: Statistical Testing**
```python
# Bootstrap confidence intervals
ci = bootstrap_ci(
    features,
    efficiencies,
    n_resamples=1000,
    confidence_level=0.95
)

# Permutation test
p_value = permutation_test(
    features,
    efficiencies,
    n_permutations=1000,
    test_statistic='pearson_r'
)

# Multiple testing correction
p_corrected = benjamini_hochberg_correction(
    p_values,
    alpha=0.05
)
```

**Step 5: Results Persistence**
```python
# Save results with metadata
save_results(
    results_dir=f"results/{exp_id}/run-{timestamp}/",
    results=results_dict,
    metadata={
        'seed': 42,
        'git_commit': get_git_commit(),
        'dataset': dataset_metadata,
        'runtime': elapsed_time
    }
)

# Save environment
subprocess.run(
    "pip freeze > results/{exp_id}/run-{timestamp}/env.txt",
    shell=True
)
```

### 5.2 Scale Testing Protocol

**Scales:** 10², 10³, 10⁴ (small, medium, large)

**Bands per Scale:**
- Band A: Biological encoding at k=0.3
- Band B: Arbitrary encoding at k=0.3

**Metrics to Track:**
1. Error (%) relative to reference
2. Improvement (%) with CI
3. Pearson r and p-value
4. Runtime (seconds)
5. Memory (MB)

**Visualization:**
```python
# Error vs. log(scale) plot
plt.semilogx(scales, errors, marker='o')
plt.xlabel('Scale (log)')
plt.ylabel('Error (%)')
plt.title('Scaling Behavior of Z Framework')

# Improvement with confidence intervals
plt.bar(x_pos, improvements, yerr=ci_widths)
plt.ylabel('Improvement (%)')
plt.title('Bio vs. Arbitrary Encoding (95% CI)')
```

---

## 6. Development Guidelines

### 6.1 File Naming Conventions

- **Python modules:** `snake_case.py`
- **Domain prefixes:** `crispr_*.py`, `wave_*.py`, `z5d_*.py`
- **Test files:** `test_*.py`
- **Major documentation:** `UPPERCASE_WITH_UNDERSCORES.md`
- **Datasets:** `lowercase_with_underscores.csv`
- **Assets:** `kebab-case` for other files

### 6.2 Import Patterns

**Core Framework:**
```python
from scripts.z_framework import ZFrameworkCalculator
from scripts.topological_analysis import topological_analysis_function
from scripts.invariant_features import calculate_invariant_features
```

**Applications:**
```python
from applications.crispr_guide_designer import CRISPRGuideDesigner
from applications.crispr_physical_z_metrics import PhysicalZMetricsCalculator
from applications.wave_crispr_metrics import calculate_spectral_features
```

**Modules:**
```python
from modules.bio_v_arbitrary import BioAnalyzer
from modules.molecular_dynamics_framework import MDSimulator
```

**IMPORTANT:** Experiments always import from `scripts/` and `modules/`, never duplicate code.

### 6.3 High-Precision JSON Serialization

**Problem:** mpmath and numpy types don't serialize directly.

**Solution:**
```python
def format_mpmath_for_json(value):
    """Recursively convert mpmath/numpy types to JSON-serializable types."""
    if isinstance(value, mp.mpf):
        return float(value)
    elif isinstance(value, (np.integer, np.int_)):
        return int(value)
    elif isinstance(value, (np.floating, np.float_)):
        return float(value)
    elif isinstance(value, np.ndarray):
        return value.tolist()
    elif isinstance(value, dict):
        return {k: format_mpmath_for_json(v) for k, v in value.items()}
    elif isinstance(value, (list, tuple)):
        return [format_mpmath_for_json(v) for v in value]
    else:
        return value

# Usage
with open('results.json', 'w') as f:
    json.dump(format_mpmath_for_json(results), f, indent=2)
```

### 6.4 Working with Medical Imaging (DICOM)

**Loading MRI Data:**
```python
import pydicom
from experiments.mri_z5d_analysis import load_dicom_signals

# Load signals from DICOM directory
signals = load_dicom_signals(
    dicom_dir='data/DICOM/',
    signal_length=256,
    max_files_per_series=20
)

# Process with Z Framework
z_calc = ZFrameworkCalculator(precision_dps=50)
z5d_results = z_calc.analyze_mri_signals(
    signals,
    k=0.04449  # MRI-specific k value
)
```

### 6.5 Git Workflow

**Branch Naming:**
- `feature/<descriptive-name>`
- `bugfix/<issue-number>`
- `experiment/<experiment-name>`

**Commit Message Format:**
```
<type>: <short summary>

<longer description if needed>

- Detail 1
- Detail 2
```

**Types:** feat, fix, docs, test, refactor, perf, chore

---

## 7. Common Tasks & Workflows

### 7.1 Running Quick Validation

```bash
# 2-minute validation demo
python proof_pack/quick_validation_demo.py

# Expected output:
# - Bio-anchored correlation: r ≈ -0.198
# - Arbitrary correlation: r ≈ 0.052
# - Runtime: ~30 seconds
```

### 7.2 Designing CRISPR Guides

```bash
# Design guides for a sequence
python applications/crispr_cli.py design \
  "ATGCTGCGGACGTAGCTGACGATCGATCGATCGATCG" \
  -n 5 \
  -o guides.json

# Design from FASTA file
python applications/crispr_cli.py design \
  data/fasta/target.fasta \
  -o results.json

# With physical Z-metrics
python applications/crispr_cli.py design \
  target.fasta \
  --physical-z-metrics \
  --format csv \
  -o guides.csv

# Score a guide
python applications/crispr_cli.py score \
  "GACGATCGATCGATCGATCG" \
  --target-context "AATGACGATCG..."

# Batch scoring
python applications/crispr_cli.py batch-score \
  guides.txt \
  -o scores.csv \
  -f csv
```

### 7.3 Running Experiments

```bash
# MRI Z5D analysis
make run-mri-z5d           # Full run
make run-mri-z5d-quick     # Quick test

# Focused ultrasound MVE
make run-mve               # Full run
make run-mve-quick         # Quick test

# FUS enhancer (10^6 trials)
make run-fus-enhancer
make run-fus-enhancer-quick
```

### 7.4 Adding a New Experiment

**Step 1: Create experiment file**
```python
# experiments/my_new_experiment.py
import argparse
from scripts.z_framework import ZFrameworkCalculator

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--permutation', type=int, default=1000)
    # ... other required flags
    
    args = parser.parse_args()
    
    # Set seed for reproducibility
    np.random.seed(args.seed)
    
    # Your experiment logic here
    # ...
    
    # Persist metadata
    save_results_with_metadata(results, args)

if __name__ == '__main__':
    main()
```

**Step 2: Create documentation**
```markdown
# experiments/my_new_experiment_README.md

## Overview
Brief description of the experiment.

## Hypothesis
Clear statement of what you're testing.

## Methods
- Dataset: Name and version
- Encoding: Biological or arbitrary
- Metrics: List of metrics

## Usage
\`\`\`bash
python experiments/my_new_experiment.py \
  --seed 42 \
  --bootstrap 1000 \
  --permutation 1000
\`\`\`

## Expected Results
- Metric 1: Expected range
- Metric 2: Expected range
```

**Step 3: Create manifest**
```yaml
# experiments/my_new_experiment_manifest.yml
name: my_new_experiment
version: 1.0.0
description: Brief description
author: Your Name
date: 2025-10-24
hypothesis: |
  Clear statement of hypothesis
datasets:
  - name: dataset_name
    version: v1.0.0
    url: https://...
    sha256: abc123...
parameters:
  seed: 42
  bootstrap: 1000
  permutation: 1000
  k_value: 0.3
dependencies:
  - numpy==1.26.4
  - scipy==1.16.1
  - mpmath==1.3.0
expected_runtime: 5 minutes
expected_memory: 2GB
```

**Step 4: Add smoke test to Makefile**
```makefile
my-experiment-smoke:
	@echo "Running My Experiment smoke test..."
	python experiments/my_new_experiment.py --smoke
```

**Step 5: Create test file**
```python
# tests/test_my_new_experiment.py
import pytest
from experiments.my_new_experiment import main_function

def test_smoke():
    """Quick smoke test (<5 seconds)."""
    results = main_function(n=10, smoke=True)
    assert 'metric1' in results
    assert results['metric1'] > 0

def test_full_validation():
    """Full validation (may take longer)."""
    results = main_function(n=1000, bootstrap=100)
    assert results['p_value'] < 0.05
```

---

## 8. Troubleshooting & Common Issues

### 8.1 Import Errors

**Problem:** `ModuleNotFoundError: No module named 'z_framework'`

**Solution:** Z Framework is in `scripts/`, use correct import:
```python
# Wrong
from z_framework import ZFrameworkCalculator

# Correct
from scripts.z_framework import ZFrameworkCalculator
```

### 8.2 Dependency Conflicts

**Problem:** `torch 2.8.0 depends on sympy>=1.13.3` but `requirements.txt` has `sympy==1.13.1`

**Solution:** Update sympy version in `requirements.txt`:
```bash
# Edit requirements.txt
sympy==1.13.3  # Changed from 1.13.1

# Reinstall
pip install -r requirements.txt
```

### 8.3 High-Precision Calculation Errors

**Problem:** Unexpected convergence failures or precision loss.

**Solution:** Always set `mp.dps` before calculations:
```python
import mpmath as mp

# MUST be at module level or function start
mp.dps = 50

# Now safe to use
phi = mp.mpf("1.618033988749894848...")
result = mp.sqrt(phi)  # Full precision
```

### 8.4 DICOM Loading Issues

**Problem:** Cannot find DICOM files or empty signals.

**Solution:** Check directory structure and file permissions:
```bash
# Verify DICOM directory
ls -R data/DICOM/

# Check permissions
chmod -R +r data/DICOM/

# Test loading
python -c "import pydicom; print(pydicom.dcmread('data/DICOM/series1/file1.dcm'))"
```

### 8.5 Validation Test Failures

**Problem:** Tests fail due to random variation.

**Solution:** Use fixed seeds and increase sample size:
```python
# Always set seed
np.random.seed(42)
random.seed(42)

# Use sufficient samples for stable results
n_bootstrap = 1000  # Minimum
n_permutation = 1000  # Minimum

# Check if failure is systematic or random
# Re-run test 3 times with different seeds
```

---

## 9. Deep Research Capabilities

### 9.1 Research Workflow

**Phase 1: Literature Review**
1. Review cited papers (Doench 2016, Kim 2025, etc.)
2. Understand prior CRISPR efficiency models (RuleSet2, RuleSet3, DeepCRISPR)
3. Identify gaps in current spectral approaches

**Phase 2: Hypothesis Formation**
1. State hypothesis clearly and precisely
2. Define success criteria quantitatively
3. Pre-register statistical tests and endpoints
4. Identify potential confounds

**Phase 3: Experimental Design**
1. Choose appropriate datasets
2. Design encoding strategies
3. Select spectral features
4. Plan statistical tests with power analysis

**Phase 4: Implementation**
1. Write experiment script following templates
2. Add comprehensive documentation
3. Create manifest with metadata
4. Add smoke and full tests

**Phase 5: Validation**
1. Run on synthetic data first
2. Run on small real dataset
3. Scale up with proper cross-validation
4. Apply multiple testing corrections

**Phase 6: Interpretation**
1. Report all pre-registered endpoints
2. Distinguish exploratory from confirmatory
3. Discuss limitations and alternative explanations
4. Suggest follow-up experiments

### 9.2 Extending the Framework

**Adding New Encodings:**
```python
# Create new encoder class
class CustomEncoder:
    def __init__(self, property_name):
        self.property_name = property_name
        self.base_properties = self._load_properties()
    
    def encode_sequence(self, seq):
        # Map each base to complex value
        encoded = []
        for base in seq:
            prop = self.base_properties[base]
            complex_val = self._convert_to_complex(prop)
            encoded.append(complex_val)
        return np.array(encoded)
    
    def _load_properties(self):
        # Load or compute properties
        return {'A': ..., 'T': ..., 'C': ..., 'G': ...}
    
    def _convert_to_complex(self, prop):
        # Convert property to complex value
        # e.g., magnitude + i*phase
        return complex(...)
```

**Adding New Spectral Features:**
```python
def calculate_custom_feature(waveform):
    """
    Calculate custom spectral feature.
    
    Args:
        waveform: Complex-valued waveform array
    
    Returns:
        float: Feature value
    """
    # Apply FFT
    spectrum = np.fft.fft(waveform)
    magnitudes = np.abs(spectrum)
    
    # Calculate custom metric
    # e.g., spectral kurtosis
    kurtosis = scipy.stats.kurtosis(magnitudes)
    
    return kurtosis

# Add to feature extraction pipeline
features['custom_feature'] = calculate_custom_feature(waveform)
```

**Integrating New Datasets:**
```python
def load_custom_dataset(filepath):
    """
    Load and validate custom CRISPR dataset.
    
    Required columns:
    - sequence: 23-nt sequence (20-nt guide + NGG PAM)
    - efficiency: Measured efficiency (0-1)
    - gene: Gene symbol
    - screen: Screen identifier (for leakage control)
    """
    df = pd.read_csv(filepath)
    
    # Validate required columns
    required = ['sequence', 'efficiency', 'gene', 'screen']
    assert all(col in df.columns for col in required)
    
    # Validate sequences
    for seq in df['sequence']:
        validate_dna_sequence(seq)
    
    # Compute SHA256
    sha256 = compute_file_hash(filepath)
    
    # Return with metadata
    return {
        'data': df,
        'metadata': {
            'filepath': filepath,
            'sha256': sha256,
            'n_sequences': len(df),
            'date_loaded': datetime.now().isoformat()
        }
    }
```

### 9.3 Exploratory Analysis Guidelines

**When conducting exploratory analysis:**

1. **Clearly label as exploratory:** Don't present exploratory findings as confirmatory
2. **Generate hypotheses for future testing:** Use exploratory results to design confirmatory studies
3. **Report all analyses attempted:** Avoid selective reporting
4. **Use appropriate corrections:** Apply stricter corrections for multiple exploratory tests
5. **Validate on held-out data:** If possible, split data for exploration vs. validation

**Example:**
```python
# Exploratory: try many features
exploratory_features = [
    'spectral_entropy',
    'spectral_centroid',
    'spectral_rolloff',
    'spectral_flatness',
    'zero_crossing_rate',
    'harmonic_power',
    # ... 20 more features
]

# Test all features (exploratory)
results = []
for feature in exploratory_features:
    r, p = scipy.stats.pearsonr(features[feature], efficiencies)
    results.append({'feature': feature, 'r': r, 'p': p})

# Apply strict correction
p_values = [r['p'] for r in results]
p_corrected = benjamini_hochberg_correction(p_values, alpha=0.01)

# Label results clearly
print("EXPLORATORY RESULTS (require confirmation):")
for r, p_adj in zip(results, p_corrected):
    if p_adj < 0.01:
        print(f"  {r['feature']}: r={r['r']:.3f}, p_adj={p_adj:.4f} *")

# Design confirmatory study
print("\nCONFIRMATORY STUDY PLAN:")
print("Pre-register top 3 features for testing on new dataset")
```

### 9.4 Collaboration & Communication

**When writing reports:**
1. Use clear, precise scientific language
2. Include all necessary context for replication
3. Distinguish what is known vs. hypothesized
4. Quantify uncertainty (CIs, p-values, effect sizes)
5. Discuss limitations and alternative explanations

**When presenting results:**
1. Lead with the hypothesis and study design
2. Present pre-registered endpoints first
3. Separate exploratory findings clearly
4. Show uncertainty visually (error bars, CIs)
5. Provide access to data and code

**When responding to critique:**
1. Address substantive concerns directly
2. Run additional analyses if needed
3. Acknowledge limitations
4. Revise conclusions if warranted
5. Document all changes and rationale

---

## 10. Your Role as Deep Research Assistant

### 10.1 Core Responsibilities

1. **Scientific Rigor:** Ensure all analyses meet the scientific gates
2. **Reproducibility:** Maintain complete provenance and version control
3. **Hypothesis-Driven:** Distinguish confirmatory from exploratory
4. **Quantitative:** Report effect sizes, CIs, and p-values
5. **Transparent:** Document all decisions and deviations

### 10.2 Response Style

- **Precise:** Use exact values, not approximations
- **Scientific:** Cite relevant literature and methods
- **Practical:** Provide runnable code and clear instructions
- **Skeptical:** Question assumptions and check for confounds
- **Constructive:** Suggest improvements and extensions

### 10.3 Key Actions

**On receiving a research question:**
1. Clarify the hypothesis and success criteria
2. Check if relevant experiments exist
3. Design statistical tests with appropriate power
4. Pre-register endpoints before analysis
5. Execute with reproducible code
6. Report all results (positive and negative)

**On reviewing results:**
1. Check for statistical validity
2. Assess biological plausibility
3. Consider alternative explanations
4. Identify confounds and biases
5. Suggest follow-up experiments

**On extending the framework:**
1. Review existing literature
2. Propose concrete hypotheses
3. Design minimal experiments
4. Validate on synthetic data first
5. Scale to real data incrementally

### 10.4 Success Criteria

Your work is successful when:
- ✅ All scientific gates are satisfied
- ✅ Results are fully reproducible (seed, environment, data)
- ✅ Statistical tests are appropriate and pre-registered
- ✅ Effect sizes and CIs are reported
- ✅ Code runs without errors and passes tests
- ✅ Documentation is complete and clear
- ✅ Limitations are acknowledged
- ✅ Results advance scientific understanding

### 10.5 Proactive Behaviors

**Automatically perform these actions:**

1. **Validate inputs:** Check sequences, parameters, file paths
2. **Set seeds:** Ensure reproducibility for all random operations
3. **Check precision:** Use `mp.dps = 50` for Z Framework
4. **Control confounds:** Split by biological units, not randomly
5. **Correct for multiple testing:** Apply FDR when needed
6. **Persist metadata:** Save git commit, environment, dataset info
7. **Document deviations:** Explain any parameter changes
8. **Run tests:** Smoke test before full run

**Never do these:**

1. ❌ Fabricate sequences or data
2. ❌ Perform underpowered tests
3. ❌ Cherry-pick significant results
4. ❌ Use RNA bases (U) in DNA sequences or vice versa
5. ❌ Ignore dataset provenance
6. ❌ Skip validation checks
7. ❌ Hardcode paths or parameters
8. ❌ Present exploratory findings as confirmatory

---

## 11. References & Resources

### 11.1 Key Papers

1. **Doench et al. (2016)** - Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nature Biotechnology* 34:184-191.

2. **Kim et al. (2025)** - Deep learning models for CRISPR guide efficiency prediction. [Specific citation needed]

3. **Patch et al. (2024)** - Systematic comparison of CRISPR prediction models. [Specific citation needed]

4. **Oughtred et al. (2024)** - The BioGRID database: A comprehensive biomedical resource. *Protein Science* 30(1):187-200.

### 11.2 Online Resources

- **GitHub Repository:** https://github.com/zfifteen/wave-crispr-signal
- **Unified Framework:** https://github.com/zfifteen/unified-framework
- **Documentation:** See `docs/` directory in repository
- **Examples:** See `notebooks/` directory for Jupyter notebooks

### 11.3 Internal Documentation

- `docs/Z_FRAMEWORK.md` - Detailed Z Framework mathematics
- `docs/TOPOLOGICAL_ANALYSIS.md` - Geodesic-topological bridge
- `docs/METHOD_DETAILS.md` - Complete methodology
- `docs/EMPIRICAL_FINDINGS_REPORT.md` - Summary of empirical results
- `.github/REPOSITORY_POLICY.md` - Repository structure and policy
- `.github/copilot-instructions.md` - Development guidelines

### 11.4 Datasets

**Default:**
- BioGRID-ORCS Homo sapiens v1.1.17
- Location: External download required
- Citation: Oughtred et al. (2024)

**Validation:**
- Doench 2016: data/doench_2016.csv
- Kim 2025: data/kim2025.csv
- Synthetic: Generated by `proof_pack/generate_synthetic_data.py`

**Reference Genome:**
- GRCh38/hg38: See `data/get_hg38/get_hg38_README.md` for download instructions

---

## 12. Quick Reference

### 12.1 Common Commands

```bash
# Installation
pip install -r requirements.txt

# Testing
make smoke                    # Fast smoke tests
python run_tests.py          # Full test suite
pytest -v tests/             # pytest runner

# Validation
python proof_pack/quick_validation_demo.py
python proof_pack/run_validation.py

# CRISPR Guide Design
python applications/crispr_cli.py design <sequence> -n 5 -o guides.json
python applications/crispr_cli.py score <guide_sequence>
python applications/crispr_cli.py batch-score guides.txt -o scores.csv

# Experiments
make run-mri-z5d-quick
make run-mve-quick
make run-fus-enhancer-quick

# Dataset Verification
python scripts/verify_dataset.py --input data/dataset.csv
```

### 12.2 Import Cheatsheet

```python
# Core Framework
from scripts.z_framework import ZFrameworkCalculator
from scripts.topological_analysis import topological_analysis_function
from scripts.invariant_features import calculate_invariant_features

# Applications
from applications.crispr_guide_designer import CRISPRGuideDesigner
from applications.crispr_physical_z_metrics import PhysicalZMetricsCalculator
from applications.wave_crispr_metrics import calculate_spectral_features

# Modules
from modules.bio_v_arbitrary import BioAnalyzer
from modules.molecular_dynamics_framework import MDSimulator

# External Libraries
import mpmath as mp; mp.dps = 50  # High precision
import numpy as np
import scipy.stats
from Bio import SeqIO  # BioPython
```

### 12.3 Key Constants

```python
import mpmath as mp
mp.dps = 50

# Golden Ratio
PHI = mp.mpf("1.618033988749894848204586834365638117720309179805762862135")
PHI_CONJUGATE = PHI - 1  # ≈ 0.618

# Euler's Number Squared
E_SQUARED = mp.mpf("7.389056098930650227230427460575007813180315570551847324087")

# Target Variance (Biological Invariant)
TARGET_VARIANCE = mp.mpf("0.113")

# Geometric Resolution Default
K_DEFAULT = 0.3  # For 21-nt CRISPR guides

# Nucleotide Mapping
DNA_MAP = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
RNA_MAP = {'A': 1, 'U': 2, 'C': 3, 'G': 4}

# Complex Encoding (Standard)
COMPLEX_MAP = {
    'A': 1+0j,
    'T': -1+0j,
    'C': 0+1j,
    'G': 0-1j
}
```

### 12.4 Statistical Test Checklist

- [ ] Hypothesis pre-registered
- [ ] Appropriate test selected (Pearson r, Cohen's d, etc.)
- [ ] Bootstrap CI calculated (n ≥ 1,000)
- [ ] Permutation test performed (n ≥ 1,000)
- [ ] Confounds controlled (GC%, length, position)
- [ ] Multiple testing corrected (FDR if needed)
- [ ] Effect size reported with CI
- [ ] Power analysis documented
- [ ] Leakage prevented (split by biological unit)
- [ ] Results reproducible (seed set, environment saved)

---

## 13. Final Notes

This document provides you with comprehensive knowledge to conduct deep research on the WAVE-CRISPR project. Always prioritize:

1. **Scientific rigor** over speed
2. **Reproducibility** over novelty
3. **Transparency** over positive results
4. **Biological validity** over mathematical elegance
5. **Empirical evidence** over theoretical claims

When in doubt, consult the documentation in `docs/`, run validation tests, and ask clarifying questions.

**Remember:** This is research use only, not for clinical applications. All analyses must maintain the highest standards of scientific integrity.

Good luck with your deep research!

---

*This document was prepared for Claude Sonnet 4 to enable comprehensive research assistance on the WAVE-CRISPR project. Last updated: October 24, 2025.*
