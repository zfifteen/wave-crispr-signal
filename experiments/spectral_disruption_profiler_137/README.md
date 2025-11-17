# Technical Design Document: Falsification Experiment for Spectral Disruption Profiler Hypothesis

**Experiment ID**: `spectral_disruption_profiler_137`  
**Version**: 1.0  
**Date Created**: 2025-11-17  
**Status**: Experimental  

---

## 0. Observations

- **DNA waveform encoding** (A→1+0i, etc.) and **FFT features** (Δf₁, ΔEntropy, sidelobes) remain foundational, but potential failures in phase weighting could arise in high-noise or unbalanced GC sequences, as seen in synthetic data with GC >70%.

- **Optimal k\* ≈ 0.300** derived from geodesic bridges (topological_analysis.py) may not generalize if curvature κ(n) = d(n) \* ln(n+1) / e² overfits to Kim 2025, per z-sandbox QMC variance reduction insights.

- **Bootstrap ΔROC-AUC +0.047 ± 0.006** may lack significance under FDR if entropy gradients inflate false positives, contradicting Doench 2016 benchmarks.

- **Unified-framework's θ ≈ 0.525** (Stadlmann 2023) and conical_flow.py enhancements suggest density corrections, but gist EXPLAIN.txt's θ=0.71 implies variability that could nullify CRISPR-specific lifts.

- **Z-invariant Z = A(B/c)** links to 5D geodesic mapping, but z-sandbox GVA on semiprimes highlights embedding failures in non-prime-like disruptions.

- No existing experiments/ folder for this specific study; new subfolder prepares for falsification runs, integrating proof_pack and run_validation.py.

- Gist kappa_signal_demo.py provides κ(n) batch processing, testable for over-reliance in auto-optimization of k.

---

## 1. Hypothesis to Falsify

### Primary Hypotheses

1. **Phase Weighting Null Hypothesis**: θ′(n,k) with k\* ≈ 0.300 does **not** enhance FFT resolution or yield >0% lift in gRNA prediction over RuleSet3, falsifiable if bootstrap CI includes zero on adversarial holdouts.

2. **Performance Null Hypothesis**: Vectorized FFT fails performance (<30s/100 sequences) under scaled noise or high-dimensional torus embedding from z-sandbox GVA.

3. **Geodesic Bridge Null Hypothesis**: Geodesic bridges (kappa_geo=0.3) do **not** mitigate k over-reliance, falsifiable if manual adjustments outperform auto-optimization in GC-biased synthetics.

4. **FDR Correction Null Hypothesis**: FDR-corrected entropy gradients do **not** reduce off-target false positives, verifiable if metrics exceed Patch 2024 thresholds in unified-framework conical flow simulations.

### Success Criteria for Falsification

- Bootstrap 95% CI for ΔROC-AUC includes or crosses zero
- Permutation test p-value > 0.05 for performance difference
- Manual k selection outperforms auto-optimization by >5%
- FDR-corrected metrics show >10% false positive rate

---

## 2. Experimentation Plan

### 2.1 Setup

**Folder Structure**:
```
experiments/spectral_disruption_profiler_137/
├── README.md                          # This document
├── manifest.yml                       # Experiment metadata
├── spectral_disruption_profiler.py    # Main implementation
├── smoke_test.py                      # Fast CI test (<5s)
└── golden_outputs/                    # Expected results for validation
    ├── smoke_test_expected.json
    └── validation_checksums.txt
```

**Dependencies**:
- NumPy/SciPy for FFT
- Bootstrap CI (≥1,000 resamples)
- z-sandbox QMC (Sobol/Owen) for sampling (if available)
- unified-framework conical_flow.py for density tests (if available)
- mpmath (dps=50) for high-precision calculations

**Datasets**:
- Kim 2025 holdouts (if available)
- Adversarial synthetics (GC 80-90%, noise-added via generate_synthetic_data.py)
- Doench 2016/Patch 2024 benchmarks (if available)

**Null Baselines**:
- RuleSet3 without weighting (baseline)
- Random k variations
- Unweighted FFT controls

### 2.2 Core Modules

#### Encoding Module
- **Purpose**: Map sequences to waveforms; apply θ′(n,k) but include unweighted controls
- **Input**: Batch CSV with noise injections
- **Output**: Comparative waveforms (weighted vs unweighted)
- **Test**: Invalid sequences trigger failures to falsify robustness

#### Analysis Module
- **Purpose**: Compute FFT features with/without golden-ratio bias
- **Features**: Batch processing with QMC variance to probe inefficiencies
- **Integration**: spectral_features.py + GeodesicMapper (if available)

#### Scoring Module
- **Purpose**: Z-invariants with CI; force k mismatches to test over-reliance
- **Implementation**: Use z_framework.py; auto-optimize vs. fixed k=0.3
- **Metrics**: ΔROC-AUC, partial correlation controlling for GC%, guide length, position

#### Detection Module
- **Purpose**: Entropy gradients with FDR; flag if ΔEntropy >0.05 persists in negatives
- **Thresholds**: Tuned to provoke false positives
- **Output**: FDR-corrected p-values

#### Visualization Module
- **Purpose**: Dashboards for null vs. alternative ROC, highlighting failures
- **Export**: JSON with p-values, bootstrap CIs, effect sizes

#### API Module
- **Purpose**: Endpoints for batch falsification runs, logging discrepancies
- **Features**: CLI interface with required flags (--seed, --bootstrap, --permutation, etc.)

### 2.3 Falsification Experiments

#### Functional Tests
- **Test**: Process 100 noisy sequences; falsify if weighted time >30s or scores match unweighted
- **Comparison**: compare_ruleset3_wave.ipynb approach
- **Criterion**: Reject if ΔROC-AUC CI crosses zero

#### Performance Tests
- **Test**: 10,000 resamples; falsify if >5s/sequence under GVA embedding
- **Scalability**: Parallel processing with adversarial GC

#### Empirical Tests
- **Test**: A/B on holdouts; falsify if lift ≤0% (p>0.05)
- **Synthetic**: run_validation.py with conical flow noise; reproduce no-bias
- **Edge Cases**: GC>80%; k adjustments show no improvement

#### Metrics Collection
- **Primary**: p-values for null hypothesis (no lift)
- **Secondary**: Effect size (Cohen's d), partial correlations
- **Output**: results/falsify_scores.csv

### 2.4 Deployment Plan

**MVP**: 
- Integrate phase_weighted_scorecard_cli.py with falsification flags
- Command-line interface for batch processing

**Cloud**: 
- For large-scale null testing (if needed)

**Coverage**: 
- 80% tests targeting failure modes
- Smoke tests for CI (<5s completion)

---

## 3. Evidence-Based Conclusions

- Proceed only if validations show potential for falsification
- Update based on empirical null rejections
- Incorporate z-sandbox QMC and unified-framework invariants to ensure rigorous testing of geometric assumptions
- Document all results with proper statistical rigor (bootstrap CIs, permutation tests, FDR correction)

---

## 4. Scientific Gates Compliance

### Human DNA Only
- Validate **A/C/G/T/N only** (case-insensitive) for DNA sequences
- Reject `U` and IUPAC ambiguity codes beyond `N`
- For RNA sequences (if needed): Validate **A/C/G/U/N only**

### Z Framework Invariants
- **Discrete/biological (DEFAULT)**: Z = A(B / e²)
- Guard divide-by-zero; document A and B parameters
- **Geometric resolution**: θ′(n,k) = φ·((n mod φ)/φ)^k with default k ≈ 0.3

### Dataset Provenance
- Record **dataset name, version, URL, license, taxonomy**, and **SHA256** of local file(s)
- **Human filter required**: Ensure **Homo sapiens**; print dataset version at runtime
- Default dataset allowed: **BioGRID-ORCS Homo sapiens v1.1.17** (with citation)

### Statistical Validity & Leakage Gates
- **Pre-register endpoints**: Pearson r with 95% bootstrap CI (≥1,000 resamples)
- **Partial r** controlling for GC%, guide length, guide position
- **Effect size** (Cohen's d with CI) for defined contrasts
- **Null model**: ≥1,000× permutation/shuffle for empirical p-values
- **Leakage control**: split-by-screen and split-by-gene; no entity appears in both train/eval
- **Multiple comparisons**: apply Benjamini–Hochberg FDR when comparing ≥3 metrics
- **Power/sample size**: include brief justification

### Reproducibility & Environment Gates
- **CLI contract** (required flags): --seed, --bootstrap, --permutation, --splits, --domain, --pam
- **Persist metadata**: seed, git commit, dataset name/version, SHA256, runtime, pip freeze
- **Precision**: use mpmath with mp.dps = 50 where high precision required
- **Pinned env**: Python 3.12.*, exact requirements.txt pins

### CI & Layout Gates
- Provide **smoke** dataset + golden outputs so CI completes <5s
- CI runs: pytest -q (validators, domain guards, stats wrappers) and make smoke
- Directory structure follows repository policy

---

## 5. Execution Commands

### Run Smoke Test (CI)
```bash
python experiments/spectral_disruption_profiler_137/smoke_test.py
```

### Run Full Falsification Experiment
```bash
python experiments/spectral_disruption_profiler_137/spectral_disruption_profiler.py \
  --seed 42 \
  --bootstrap 1000 \
  --permutation 1000 \
  --splits split-by-gene \
  --domain discrete \
  --k-parameter 0.3 \
  --output results/spectral_disruption_profiler_137/
```

### Batch Processing with Adversarial Data
```bash
python experiments/spectral_disruption_profiler_137/spectral_disruption_profiler.py \
  --input data/adversarial_high_gc.csv \
  --seed 42 \
  --bootstrap 10000 \
  --permutation 10000 \
  --k-parameter auto \
  --output results/spectral_disruption_profiler_137/adversarial/
```

---

## 6. Expected Outcomes

### If Hypotheses are Falsified (Expected)
- ΔROC-AUC bootstrap CI includes zero
- p-value > 0.05 for performance improvements
- Manual k selection outperforms auto-optimization
- FDR-corrected false positive rate >10%

### If Hypotheses are Supported (Unexpected)
- ΔROC-AUC bootstrap CI excludes zero with p < 0.05
- Significant performance improvements
- Auto-optimization performs best
- FDR-corrected false positive rate <5%

---

## 7. References

- Kim 2025: gRNA efficiency dataset
- Doench 2016: RuleSet3 baseline
- Patch 2024: Off-target threshold benchmarks
- Stadlmann 2023: θ ≈ 0.525 parameter
- BioGRID-ORCS v1.1.17: Homo sapiens CRISPR dataset

---

## 8. Version History

- **v1.0** (2025-11-17): Initial technical design document created
