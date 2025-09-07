# Signal-Theoretic Modeling of CRISPR-Cas9 for Cancer-Target Editing

## Experiment ID: signal_theoretic_crispr

### Overview
A fully reproducible, falsifiable experiment plan that evaluates **signal-theoretic (spectral–geodesic)** features from `wave-crispr-signal` for **predicting on-target efficiency** in human genomes, compared against a **batch-scan baseline**.

### Objective
Quantify whether **WAVE-CRISPR spectral–geodesic features** deliver **predictive lift** over a **batch-scan baseline** (PAM/seed scanning + simple heuristics) for:
- **On-target efficiency (regression)** on human datasets

**Note**: Off-target discrimination analysis is planned for a future release (see issue #XX).

### Falsifiable Hypotheses
- **H₀ (no lift)**: WAVE-CRISPR features do **not** improve R² (on-target) vs. baseline (**Δ ≤ 0**)
- **H₁ (lift)**: WAVE-CRISPR features improve **R² by ≥ 0.05 absolute**, with **95% bootstrap CI** excluding 0

### Scientific Gates (Must-Pass)
1. **G1 – Human source only**: Every sequence must be explicitly human (hg38 windows, human screen datasets)
2. **G2 – Alphabet**: FASTA sequences must contain only **A/C/G/T/N** (uppercase). Any other symbol → **abort**
3. **G3 – Anchoring**: Each guide must anchor to an **hg38** locus. If a 201-bp window cannot be extracted → **fail gracefully with synthetic context**
4. **G4 – Determinism**: Fixed seeds; pinned library versions; log commit hashes and parameters
5. **G5 – Metrics**: Report **R²/MAE** (on-target) with **95% bootstrap CIs**
6. **G6 – Numeric Stability**: Report numerical reproducibility separately
7. **G7 – Controls**: Include shuffled controls, PAM-broken controls, reverse-complements
8. **G8 – Ethics**: Only public, human datasets or user-owned sequences

### Execution Commands

#### Setup
```bash
cd /home/runner/work/wave-crispr-signal/wave-crispr-signal
python -m experiments.signal_theoretic_crispr.validate_setup
```

#### Run Full Experiment
```bash
python -m experiments.signal_theoretic_crispr.main \
    --seed 42 \
    --bootstrap 1000 \
    --permutation 1000 \
    --splits split-by-gene \
    --domain discrete \
    --pam NGG \
    --output results/signal_theoretic_crispr/run-$(date +%Y%m%d-%H%M%S)
```

#### Quick Validation (Smoke Test)
```bash
python -m experiments.signal_theoretic_crispr.smoke_test
```

### Expected Results
- **Baseline (Arm A)**: Simple features with linear/logistic models
- **Spectral (Arm B)**: Two encoding bands (Bio-anchored vs Arbitrary) with spectral-geodesic features
- **Success Criterion**: Δ ≥ 0.05 improvement with 95% CI excluding 0

### Data Sources
- **On-target**: Doench 2016-style human CRISPR guide datasets
- **Reference**: hg38 FASTA for 201-bp window extraction (with fallback mode)
- **Gene-level**: BioGRID-ORCS Homo sapiens v1.1.17

**Future Work**: Off-target classification with GUIDE-seq-style datasets will be implemented in a separate release.

### Theoretical Framework
- **Normalization**: Z = n(Δₙ/Δₘₐₓ) with guards for Δₘₐₓ ≠ 0
- **Geodesic resolution**: θ'(n,k) = φ·((n mod φ)/φ)^k with φ = (1+√5)/2 and k ≈ 0.3
- **Discrete domain**: No velocity terms; zero-division and domain checks only

### Timeline
- **Setup validation**: ~5 minutes
- **Full experiment**: ~30-60 minutes (depending on dataset size)
- **Results generation**: ~10 minutes

### Output Structure
```
results/signal_theoretic_crispr/run-YYYYMMDD-HHMMSS/
├── results.json          # Primary metrics and statistics
├── results.csv           # Detailed per-sample results
├── env.txt              # Environment and dependencies
├── log.txt              # Execution log
├── validation_report.txt # Gate compliance report
└── figures/             # Visualization outputs
    ├── baseline_vs_spectral.png
    ├── bootstrap_ci.png
    └── feature_importance.png
```