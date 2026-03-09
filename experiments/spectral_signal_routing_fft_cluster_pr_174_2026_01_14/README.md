# Spectral Signal Routing Hypothesis 1 Falsification Experiment

**Experiment ID**: `spectral_signal_routing_fft_cluster_pr_174_2026_01_14`  
**Status**: Design Complete, Implementation Pending  
**PR**: #174  
**Created**: January 14, 2026

---

## Executive Summary

This experiment is designed to **falsify** the strong claim that spectral channels in Emergent Doom Engine (EDE) substrate dynamics are task-linked, routing-like, and causally necessary for high performance.

### Hypothesis (H1-strong)

> For appropriately chosen tasks and substrates in the Emergent Doom Engine, successful problem-solving runs will reliably exhibit non-trivial, stable spectral channels (distinct peaks or bands in FFT analysis) that:
> 
> 1. Are absent or strongly attenuated in matched non-solving controls
> 2. Are causally necessary for maintaining high task performance

### Falsification Strategy

The experiment applies **three independent falsification axes**:

| Axis | Test | Evidence Against H1-strong |
|------|------|---------------------------|
| **Distinctiveness** | A1 (success) vs A2 (failure) | No robust spectral differences |
| **Specificity** | B1 (structured) vs B2/B3 (scrambled/IID) | Channels indistinguishable from controls |
| **Causal Necessity** | C1 (baseline) vs C2 (knock-out) vs C3 (noise) | C2 no worse than C3 |

If the experiment fails to support H1-strong on **any** of these three axes, the hypothesis is rejected in the tested regime.

---

## Key Features

- **FFT-based spectral analysis** of EDE cell-state trajectories
- **Operational channel definition** (pre-registered, no post-hoc adjustment)
- **Phase-scrambled controls** to test channel specificity
- **Spectral knock-out interventions** to test causal necessity
- **No fabrication**: All parameters justified or fixed a priori
- **Reusable infrastructure**: Applicable to factorization, SAT, and future EDE tasks

---

## Experiment Structure

```
spectral_signal_routing_fft_cluster_pr_174_2026_01_14/
├── TECH_SPEC.md                    # Authoritative technical specification (22KB)
├── TODO_IMPLEMENTATION.md          # Implementation checklist (8KB)
├── experiment_specification.json   # Machine-readable config (8KB)
├── manifest.yml                    # Experiment metadata (11KB)
├── README.md                       # This file
└── notes/                          # Optional design notes
```

---

## Documentation Guide

### Primary Documents

1. **[TECH_SPEC.md](./TECH_SPEC.md)** - Start here
   - Complete technical specification
   - Hypothesis articulation and falsification objectives
   - Condition definitions (A/B/C families)
   - FFT analysis pipeline
   - Falsification logic

2. **[TODO_IMPLEMENTATION.md](./TODO_IMPLEMENTATION.md)** - Implementation roadmap
   - Detailed checklist (9 sections, 60+ items)
   - Task wiring, trace recording, series extraction
   - Spectral analysis harness components
   - Condition orchestration
   - Validation and testing requirements

3. **[manifest.yml](./manifest.yml)** - Experiment metadata
   - Machine-readable experiment tracking
   - Hypothesis and falsification criteria
   - Scientific gates enforcement
   - Reproducibility requirements

4. **[experiment_specification.json](./experiment_specification.json)** - Configuration
   - JSON schema for task/substrate parameters
   - To be populated during implementation
   - Currently contains placeholders marked `TO_BE_DEFINED`

---

## Experimental Design

### Condition Families

#### Family A: Performance Comparison
- **A1**: High-performance runs (meeting success criterion)
- **A2**: Low-performance runs (failing success criterion, same parameters)

**Tests**: Do spectral channels distinguish success from failure?

#### Family B: Control Comparison
- **B1**: Structured dynamics (same as A1/A2)
- **B2**: Phase-scrambled (preserve power spectrum, destroy coherence)
- **B3**: IID random (preserve marginals, no correlation)

**Tests**: Are channels specific to structured dynamics or generic artifacts?

#### Family C: Intervention Comparison
- **C1**: Baseline replay (no perturbation)
- **C2**: Spectral knock-out (attenuate identified channels)
- **C3**: Non-spectral perturbation (broad noise, same variance)

**Tests**: Are channels causally necessary for performance?

---

## Observables and Time Series

### Data Collection

From EDE substrate snapshots, extract:

1. **Task-level series** x_task(t)
   - Best score, constraint satisfaction, residual error

2. **Global aggregation series** x_global(t)
   - Mean/variance of cell-level fields

3. **Spatially localized series** x_cell_region(k, t)
   - Front band, cluster cores, random samples

### FFT Analysis Pipeline

```
Raw snapshots → Time series extraction → Windowing → Detrending (optional)
→ Tapering (Hann/Hamming) → Normalization → FFT → Channel detection
→ Feature aggregation → Statistical comparison
```

### Channel Definition (Operational)

A **channel** is a contiguous frequency band [f_low, f_high] with:

1. **Elevated magnitude**: Exceeds baseline by threshold (e.g., N σ)
2. **Reproducibility**: Detected across multiple windows/runs

Thresholds **must be fixed before labeling** conditions as success/failure.

---

## Scientific Gates

Following repository policy (`.github/REPOSITORY_POLICY.md`):

| Gate | Enforcement | Implementation |
|------|-------------|----------------|
| **No fabrication** | ✓ | All parameters justified or pre-registered |
| **Pre-registration** | ✓ | Analysis params fixed before labeling |
| **Fail-fast validation** | ✓ | Clear errors on invalid input |
| **Reproducibility** | ✓ | Seeds, SHA256, environment capture |
| **Scope control** | ✓ | Work traces to H1-strong or reusable FFT infra |

---

## Current Status

**Phase**: Design Complete  
**Implementation**: Pending

### Completed
- [x] Technical specification (TECH_SPEC.md)
- [x] Implementation checklist (TODO_IMPLEMENTATION.md)
- [x] Experiment metadata (manifest.yml)
- [x] Configuration schema (experiment_specification.json)
- [x] Repository integration (folder structure, naming)

### Next Steps

Per [TODO_IMPLEMENTATION.md](./TODO_IMPLEMENTATION.md):

1. **Task selection** (Section 1)
   - Choose EDE task (SAT solving, factorization, etc.)
   - Define binary success criterion
   
2. **Trace recorder** (Section 2)
   - Implement per-step snapshot recording
   - Document schema
   
3. **Time series extraction** (Section 3)
   - Pre-declare series types
   - Implement extraction pipeline
   
4. **Spectral analysis harness** (Section 4)
   - Windowing, FFT, channel detection
   - Feature aggregation

---

## Emergent Doom Engine Context

### Doom as Convergence

"Doom" in EDE = **inevitable convergence** toward target state, not catastrophe.

If spectral channels exist, they are **"doom pathways"** in frequency space: structurally preferred modes for information propagation.

### Alignment with EDE Philosophy

- **Domain-general substrate**: Simple local rules → emergent global solutions
- **Unreliable substrates**: Error tolerance via clustering and delayed gratification
- **FFT-compatible analysis**: Reusable across tasks without fabricating outcomes

### Scope Guardrails

- No changes to core substrate semantics
- No post-hoc redefinition of "channel" or "success"
- All new ideas trace to H1-strong falsification or reusable infrastructure

---

## Reproducibility

### Environment
- **Python**: 3.12+
- **Java**: TBD (for EDE substrate)
- **Dependencies**: `requirements.txt` (root of repo)

### Required Flags
- `--seed`: RNG seed for reproducibility
- `--task`: Task identifier
- `--output`: Results directory

### Metadata Capture
- Git commit hash
- Timestamp
- pip freeze
- Dataset checksums (SHA256)
- Configuration hash

---

## Usage (Future)

Once implemented, expected usage pattern:

```bash
# Run high-performance condition (A1)
python run_spectral_experiment.py \
  --task sat_solving \
  --condition A1 \
  --seed 42 \
  --output results/spectral_signal_routing_fft_cluster_pr_174_2026_01_14/run-20260114-120000/

# Run spectral knock-out intervention (C2)
python run_spectral_experiment.py \
  --task sat_solving \
  --condition C2 \
  --seed 42 \
  --channels results/.../identified_channels.json \
  --output results/.../run-20260114-130000/
```

---

## References

### Related Experiments
- `experiments/spectral_disruption_profiler/` - Spectral analysis for CRISPR guides
- `experiments/phi_geometry_falsification/` - φ-geometry hypothesis testing
- PR #174 - SAT clustering experiment

### Documentation
- `.github/REPOSITORY_POLICY.md` - Repository structure and standards
- `.github/copilot-instructions.md` - Z Framework and scientific gates

### Emergent Doom Engine
- EDE_PM_INSTRUCTIONS.md (referenced in spec)
- Emergent Doom Engine Space context

---

## Contact

**Owner**: Dionisio Alberto Lopez III (zfifteen)  
**Repository**: [zfifteen/wave-crispr-signal](https://github.com/zfifteen/wave-crispr-signal)  
**Issues**: [GitHub Issues](https://github.com/zfifteen/wave-crispr-signal/issues)  
**Tags**: `spectral-routing`, `fft-analysis`, `emergent-doom-engine`, `falsification`

---

## License

MIT License - See repository root `LICENSE` file

---

**Last Updated**: 2026-01-14  
**Document Version**: 1.0
