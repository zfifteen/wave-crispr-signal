# Model Novelty Statement

Date: 2026-03-09
Status: Active, authoritative for model-positioning claims.

## 1) Purpose and Boundary

This statement defines what is claimed as novel in this repository and what is not.

Primary boundary:
- Novelty target: the mathematical and feature-engineering model.
- Non-target: software architecture, packaging, framework choice, and deployment shape.

This means the repository is not claiming novelty for:
- Python/FastAPI/CLI composition,
- standard FFT/statistical primitives,
- generic API or tool orchestration patterns.

## 2) Active Model Surfaces (Code Ground Truth)

The active model behavior is currently implemented across these surfaces:

- `applications/phase_weighted_scorecard.py`
- `applications/genomic_disruption_api.py`
- `applications/crispr_guide_designer.py` (optional phase-weighted path)
- `wave_crispr_signal/sequence_utils.py`
- `wave_crispr_signal/features/phase_weighting.py`
- `wave_crispr_signal/spectral.py`

Runtime entrypoints exposing this model behavior:
- `applications/phase_weighted_scorecard_cli.py`
- `applications/crispr_cli.py`
- `applications/genomic_disruption_api.py`
- `web_apps/demo_mvp/app.py`

## 3) Model Components Claimed as Novel (Candidate Novelty)

### 3.1 Phase-weighted disruption framework (core scoring path)

Implemented in `applications/phase_weighted_scorecard.py` with helpers from `wave_crispr_signal`:

1. Nucleotide complex encoding for scorecard path:
- DNA: `A=1+0j`, `T=-1+0j`, `C=0+1j`, `G=0-1j`, `N=0+0j`
- RNA: `A=1+0j`, `U=-1+0j`, `C=0+1j`, `G=0-1j`, `N=0+0j`

2. Position-dependent phase modulation:
- `theta_prime(n, k) = phi * ((n mod phi)/phi)^k`
- active default `k = 0.3`
- phase weighting applied as multiplicative complex factor `exp(i * theta_prime(n,k))`

3. Mutation/disruption feature extraction:
- `delta_entropy = H(mut) - H(ref)`
- `delta_freq = |f_mut - f_ref|`
- `delta_sidelobes = |sidelobes_mut - sidelobes_ref|`

4. Composite disruption scalar:
- `delta_spectral = 0.4*delta_entropy + 0.35*delta_freq + 0.25*delta_sidelobes`

5. Diversity-coupled curvature scaling:
- `kappa = diversity * ln(n+1) / e^2` for `n >= 10`, else `1.0`

6. Final bounded score:
- `z_score = sigmoid((delta_spectral / PHI), kappa)`

Candidate novelty claim in this path:
- The specific coupling of phase-weighted spectral perturbation with diversity-scaled curvature into a bounded disruption score for CRISPR context.

### 3.2 Dimensionless AT/GC rate-ratio spectral encoding (demo/breathing path)

Implemented in `wave_crispr_signal/spectral.py` and used by `web_apps/demo_mvp/app.py`:

1. Dimensionless parameterization:
- `r = k_GC / k_AT` (default `r=20.0`)
- `alpha = log(r)`
- `beta = 0.3 * alpha`

2. Complex mapping for this path:
- DNA: `A,T -> -alpha + i*beta`, `G,C -> +alpha - i*beta`
- RNA: `A,U -> -alpha + i*beta`, `G,C -> +alpha - i*beta`

3. Fractional-period spectral probing:
- one-bin DFT/CZT-lite evaluation at non-integer periods
- primary probes at `10.5`, `5.25`, and `3.5` bp

4. Returned feature set:
- period powers (`P10_5`, `P5_25`, `P3_5`), phase angle, normalized power, GC content, and rotational phase curve summaries.

Candidate novelty claim in this path:
- A dimensionless, rate-ratio-driven complex encoding tied to helical-period harmonic probes for sequence accessibility proxies.

### 3.3 Distinct but related model families

The repo currently contains two related spectral model families:
- phase-weighted disruption scorecard family,
- rate-ratio breathing harmonic family.

Both are active and implemented, but they are not yet unified into one mathematically closed single model API.

## 4) What Is Explicitly Not Claimed as Novel

Not novel by intent:
- CLI design and argument parsing patterns.
- FastAPI service structure and route wiring.
- Use of FFT, Shannon entropy, peak counting, bootstrap mechanics by themselves.
- JSON/CSV serialization, logging, and standard Python packaging patterns.

## 5) Implementation Reality and Limits

Current implementation truth that constrains novelty statements:

1. `applications/genomic_disruption_api.py` CLI:
- Implemented commands: `score`, `design`
- Parser-only but not wired in dispatcher: `batch`, `offtarget`

2. `applications/crispr_cli.py`:
- phase-weighted scoring is available for design path via `--use-phase-weighted-scorecard`
- not the default behavior unless flag is passed

3. An authoritative in-repo comparator harness now exists (`validation/ontarget_gate`, protocol v3), but current clean-holdout results do not demonstrate superiority over the required external comparator threshold.

4. Therefore, the repository can assert:
- the model formulation is implemented,
- the model can be executed reproducibly on the active surfaces,
- the novelty is in formulation/feature coupling design.

5. The repository should not assert (yet):
- proven state-of-the-art accuracy gains,
- generalizable causal superiority across datasets,
- publication-grade empirical validation completeness.
 - superiority to modern external comparators based on current v3 clean-holdout decision outputs.

## 6) Falsifiable Novelty Claims (Testable Form)

The model-level novelty can be evaluated with falsifiable propositions:

1. Phase-weighted contribution claim:
- Holding other features fixed, introducing `theta_prime` phase weighting improves discrimination on a prespecified CRISPR task versus a non-phase-weighted spectral baseline.

2. Curvature coupling claim:
- Diversity-scaled `kappa` in the sigmoid aggregator improves calibration or ranking stability versus fixed-slope aggregation.

3. Rate-ratio encoding claim:
- Dimensionless `r`-based encoding yields more stable cross-sequence spectral signatures than fixed scalar mappings under controlled perturbations.

4. Harmonic probe claim:
- Explicit 10.5/5.25/3.5-period probes carry incremental predictive signal beyond baseline composition-only features.

Each claim is disprovable by predeclared metrics, datasets, and thresholds.

## 7) Required Evidence to Upgrade from “Novel Formulation” to “Validated Novel Model”

To move from implemented hypothesis to validated novelty, a future validation package should include:

1. Locked datasets and splits:
- fixed train/validation/test definitions with leakage controls.

2. Baseline set:
- composition-only baseline,
- non-phase spectral baseline,
- at least one strong external CRISPR scoring baseline.

3. Ablation matrix:
- remove `theta_prime`, remove `kappa`, remove harmonic probes, swap encodings.

4. Metrics and uncertainty:
- ranking/classification/regression metrics as applicable,
- confidence intervals and statistical significance tests.

5. Reproducibility contract:
- seeded runs, exact commands, pinned environment, artifact capture.

Current state update:
- v3 comparator gate is implemented with decision-grade clean holdouts and strict overlap control.
- current v3 decision outcome is `NO-GO` for baseline_c threshold criteria.

Given this state, the correct positioning is:
- novel model design and implementation,
- decision-grade comparator validation infrastructure in place,
- no validated external-comparator superiority yet.

## 8) One-line Positioning for External Use

Accurate external statement:
- "This repository implements a novel phase-weighted spectral CRISPR scoring formulation (with diversity-coupled aggregation and helical-period harmonic features); novelty is model-level, while architecture is intentionally conventional."
