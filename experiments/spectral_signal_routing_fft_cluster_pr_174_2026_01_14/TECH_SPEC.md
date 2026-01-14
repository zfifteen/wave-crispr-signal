# Technical Design Specification: Hypothesis 1 Falsification Experiment

**Experiment ID**: `spectral_signal_routing_fft_cluster_pr_174_2026_01_14`  
**Associated PR**: #174  
**Date**: January 14, 2026 (10:36 AM EST, Pittsburgh, PA, US)  
**Owner**: Dionisio Alberto Lopez III ("zfifteen")  
**Repository**: zfifteen/wave-crispr-signal

---

## 1. Hypothesis and Falsification Target

### 1.1 Hypothesis 1 (Spectral Signal Routing via FFT-Compatible Clustering)

**Working hypothesis**: Local update rules on the Emergent Doom Engine substrate produce cell-state trajectories that, when viewed in the frequency domain, exhibit structured spectral channels that can be interpreted as routing pathways for task-relevant information. These pathways are hypothesized to be mediated by spatiotemporal clustering patterns that are "FFT-compatible" (i.e., they manifest as stable, interpretable bands or peaks under discrete Fourier analysis over space, time, or both).

This hypothesis is deliberately strong: it asserts not just that spectral structure exists, but that it is both:

1. **Task-linked**: meaningfully associated with successful problem solving rather than generic substrate behavior, and
2. **Routing-like**: consistent with the idea that information "flows" along identifiable spectral channels rather than being purely diffusive or random.

### 1.2 Strong Falsifiable Form (H1-strong)

The experiment targets the following falsifiable form:

> **H1-strong**: For appropriately chosen tasks and substrates in the Emergent Doom Engine, successful problem-solving runs will reliably exhibit non-trivial, stable spectral channels (distinct peaks or bands) that (a) are absent or strongly attenuated in matched non-solving controls, and (b) are causally necessary for maintaining high task performance.

This specification is experiment-centric rather than effect-size-centric: it defines what must be measured and how falsification pressure is applied, without assuming any particular numerical outcome.

### 1.3 Falsification Objectives

The experiment is designed to attempt to falsify H1-strong by seeking evidence against each of its core commitments:

1. **Existence/Distinctiveness**  
   Do successful runs exhibit spectral structure that is systematically different from unsuccessful or randomized controls under identical macro-parameters?

2. **Stability/Consistency**  
   Are any apparent spectral channels stable over time and reproducible across runs, rather than transient artifacts?

3. **Causal Necessity**  
   If we disrupt inferred channels via targeted spectral interventions while holding other statistics as constant as possible, does task performance degrade more than under non-channel-aligned perturbations?

If the experiment, under adequate sampling and well-controlled analysis, fails to find supporting evidence on any of these three axes, H1-strong is rejected (in the tested regime).

---

## 2. Emergent Doom Engine Context and Scope Control

### 2.1 Alignment with Emergent Doom Engine Purpose

The Emergent Doom Engine is a domain-general, morphogenesis-inspired substrate where simple local agents on unreliable substrates grind inexorably toward global solutions via emergent dynamics such as clustering, delayed gratification, and error tolerance.

Within this frame:

- **"Doom"** is interpreted as inevitable convergence toward a target state, not catastrophe.
- **Spectral channels**, if they exist, would be interpreted as emergent "doom pathways" in frequency space: structurally preferred modes along which information and influence propagate.

This experiment is explicitly scoped to:

- Test whether such "spectral doom pathways" are actually present and necessary, rather than assumed.
- Provide a reusable, FFT-compatible analysis harness that can be applied across tasks (factorization, SAT, etc.) without fabricating outcomes.

### 2.2 Scope Guardrails (Project Management Constraints)

To prevent scope creep and maintain alignment with EDE's PM discipline:

- **Top-level goal**: Evaluate H1-strong for spectral signal routing via FFT-compatible clustering on the Emergent Doom Engine substrate.
- **Current sub-goal**: Define a technical design specification and on-disk experiment skeleton to support a future implementation and analysis pass.
- **You are here**: Designing the experiment's structure and repository integration, not executing or interpreting results.

**Guardrails**:

- No fabricated metrics, p-values, or effect sizes.
- No post-hoc redefinition of "channel" or "success" to rescue the hypothesis.
- All analysis parameters (windowing, thresholds, selection of time series) must be declared before labeling conditions as "success"/"failure" during implementation.
- New ideas must be traceable back to either:
  - H1-strong falsification, or
  - Reusable infrastructure for emergent dynamical analysis (FFT harness, intervention mechanisms).

---

## 3. High-Level Experimental Design

### 3.1 Substrate and Task

**Substrate**: Emergent Doom Engine cell substrate (Java), using existing one-dimensional or two-dimensional cell arrays with local update rules. No changes to core substrate semantics are assumed at the design level.

**Task class**: Any existing EDE-integrated task with a well-defined success criterion, such as:

- Factor localization (e.g., semiprime factorization runs, as in existing factor experiments).
- SAT solving dynamics (e.g., SAT clustering experiment associated with PR #174).
- Other tasks added later, provided they define "success vs failure" at the run level.

To avoid fabricating information, this specification does not assume any particular task's performance profile; it instead defines an analysis interface that any compatible task can satisfy.

### 3.2 Observation Space

The experiment focuses on time series derived from cell-level and global metrics:

- **Per-cell scalar fields** (e.g., fitness, strategy id encoded as integer, energy).
- **Aggregated measures** (e.g., mean/variance of a field across the array, front band aggregations).
- **Task-level metrics** (e.g., best score, constraint satisfaction fraction, factor localization indicators).

These are observed at discrete time steps over the course of each run and transformed into time series suitable for FFT-based analysis.

### 3.3 Condition Families (A/B/C)

The design is organized into three condition families:

1. **Family A: High-performance vs Low-performance runs**  
   Same macro-parameters, different outcomes (success vs non-success under a fixed step budget).  
   Tests whether spectral signatures distinguish successful from unsuccessful runs.

2. **Family B: Structured vs Randomized controls**  
   Structured dynamics vs phase-scrambled and IID-random controls derived from the observed signals.  
   Tests whether any observed channels are specific to emergent structure rather than generic to the marginal distributions or analysis pipeline.

3. **Family C: Spectral disruption interventions**  
   Baseline replays vs spectral "knock-out" vs non-spectral perturbations.  
   Tests whether channels are causally necessary in the strong sense of H1-strong.

---

## 4. Condition Definitions

### 4.1 Family A: High-Performance vs Low-Performance

**Goal**: Assess existence and distinctiveness of spectral structure relative to task performance.

#### A1: High-Performance Runs

**Definition (logical)**:  
Runs that satisfy a task-specific success criterion within a predefined step budget T_max.

The success criterion must be defined strictly in the task implementation (e.g., "solution found," "all constraints satisfied," "factors localized to a given window"), with no dependence on spectral analyses.

**Data to record (per run)**:

- Full sequence of per-step snapshots as defined in Section 5.
- Task-level history (e.g., best score over time, boolean success flag).
- Run metadata: task name, configuration hash, RNG seed scheme, etc.

#### A2: Low-Performance Runs

**Definition (logical)**:  
Runs with identical macro-parameters and sampling protocols as A1, but that fail to meet the success criterion within the same T_max.

They may terminate at T_max or earlier if the underlying runner includes fail-fast logic.

**Data to record**:  
Same as A1.

**Falsification pressure**: If properly powered comparisons of A1 vs A2 show no robust differences in spectral features as defined later, H1-strong's distinctiveness requirement is undermined.

### 4.2 Family B: Structured vs Randomized Controls

**Goal**: Evaluate whether any apparent spectral channels are specific to structured EDE dynamics.

#### B1: Structured Dynamics

**Definition**:  
The same runs as A1/A2, analyzed under the same windowing and FFT pipeline.

#### B2: Phase-Scrambled Controls

**Definition (synthetic)**:  
Time series derived from B1 runs are transformed such that:

- The marginal power spectrum is preserved as closely as possible.
- Phase relationships across time and/or cells are randomized (e.g., randomizing phase in the frequency domain or permuting time indices within constraints).

The concrete implementation (e.g., real-valued phase randomization method) must be documented at implementation time.

**Purpose**:  
Preserve "how much" energy exists at each frequency but destroy coherent phase structure that could encode routing-like behavior.

#### B3: IID Random Controls

**Definition (synthetic)**:  
Time series generated by drawing each sample IID from the empirical marginal distribution of the corresponding series (e.g., same mean/variance, no temporal or spatial correlation).

**Falsification pressure**: If B1's spectral features cannot be robustly distinguished from B2/B3 under consistent criteria, then any apparent channels are likely analysis artifacts or generic properties of the marginals.

### 4.3 Family C: Spectral Disruption Interventions

**Goal**: Test causal necessity of inferred channels.

#### C1: Baseline Replay

**Definition**:  
Replays of previously recorded A1 high-performance runs (or equivalent deterministic re-executions) without any perturbation.

Used to validate determinism assumptions and as a baseline for intervention comparisons.

#### C2: Spectral Knock-Out

**Definition (conceptual)**:  
Starting from the same initial conditions and RNG seeds as an A1 high-performance configuration, introduce structured perturbations that:

- Are constructed via inverse transforms where selected spectral bands (those identified as channels under Section 7) are attenuated, removed, or phase-scrambled.
- Preserve overall amplitude distributions and variance as much as possible.

**Application modes** (to be chosen at implementation time and fixed):

- Direct modifications to time series inputs that feed into the update rules.
- Direct perturbations of cell-level state variables each step using inverse-transformed signals.

#### C3: Non-Spectral Control Perturbations

**Definition**:  
Perturbations with comparable magnitude (e.g., same variance of injected noise) that do not preferentially target the identified spectral bands.

Examples include:

- Random cell flips or swaps.
- Additive, unstructured Gaussian noise in the time domain.
- Uniform noise that is spectrally broad rather than band-limited.

**Falsification pressure**: If performance under C2 is not more degraded than under C3 (or if it remains essentially intact), then the channels are not necessary in the strong sense claimed by H1-strong.

---

## 5. Data Collection Requirements

### 5.1 Per-Step Substrate Snapshots

The spec defines what must be recorded, not how often or how many runs are needed. Power analysis and run counts must be chosen in the implementation phase using non-fabricated data and clear criteria.

For each run at time step t:

**Cell-level data** (for each cell i):

- A stable identifier (index or coordinate).
- One or more scalar state components selected for analysis (examples, not commitments):
  - Strategy or behavior ID (encoded as integer).
  - Local fitness or error proxy.
  - Any other scalar that the task implementation exposes as relevant to dynamics.

**Global step-level metadata**:

- Step index t.
- Task-level metrics at step t (e.g., best objective value, constraint satisfaction fraction).
- Run identifier and configuration hash.
- RNG seed or reproduction metadata as available.

**Format expectations**:

- Machine-readable formats such as CSV, JSON, or binary tables.
- Schema and field names documented alongside implementation; this spec only constrains semantic roles.

### 5.2 Derived Time Series for Spectral Analysis

From the raw snapshots, the implementation will derive a fixed set of time series prior to condition labeling:

1. **Task-level series** x_task(t)  
   One scalar per step capturing task progress (e.g., best fitness, residual error).

2. **Global aggregation series** x_global(t)  
   Aggregations such as mean or variance of a chosen cell-level field per step.

3. **Spatially localized series** x_cell_region(k, t)  
   Series derived from specific spatial regions or cell subsets, for example:
   - Front band (e.g., first N cells).
   - Cluster cores identified by clustering logic from previous EDE experiments.
   - Randomly sampled cells for baseline.

The precise set of series and selection rules must be:

- Declared in a configuration file or constants block.
- Fixed before analyzing "success" vs "failure" labels, to avoid cherry-picking.

---

## 6. FFT-Compatible Spectral Analysis Pipeline

### 6.1 Windowing and Preprocessing

For each time series x(t):

1. **Window segmentation**  
   Choose window length T_window and stride T_stride such that:
   - Multiple windows fit into typical run durations for the chosen task.
   - The same parameters apply across all runs and conditions.

2. **Detrending** (optional, but globally consistent)  
   Optionally subtract a simple trend (e.g., mean or linear fit) per window to remove low-frequency drift.  
   If used, the detrending method must be fixed and documented.

3. **Tapering**  
   Apply a standard window function (e.g., Hann) to each window to reduce spectral leakage.

4. **Normalization**  
   Normalize within each run/channel (e.g., z-score over time or scale to unit variance) to prevent trivial amplitude differences from dominating the analysis.

All preprocessing decisions (on/off, parameter values) must be recorded and applied identically to all condition families.

### 6.2 FFT Computation

Within each window:

- Compute the discrete Fourier transform / FFT on the preprocessed series to obtain complex coefficients X(f).
- Derive:
  - Magnitude spectrum |X(f)|.
  - Optional phase arg(X(f)) if needed for phase-based metrics.

No explicit sampling frequency or frequency resolution is fixed here; these will follow from the simulation step size and chosen T_window. Implementation must pick values compatible with typical run lengths.

### 6.3 Operational Definition of a "Channel"

This experiment defines a spectral channel operationally to avoid post-hoc re-interpretation:

A **channel** in a given windowed series is a contiguous frequency band [f_low, f_high] satisfying:

1. **Elevated magnitude**  
   The band's average magnitude exceeds a baseline (e.g., run-specific or control-derived) by a pre-specified margin (e.g., threshold in units of standard deviations or quantiles).

2. **Reproducibility**  
   The same band (within allowed frequency tolerance) is detected across multiple windows of a run or across multiple high-performance runs under identical macro-parameters.

Thresholds (e.g., min width, min strength, min number of windows) are:

- To be specified during implementation.
- To be fixed before referencing success/failure labels.

### 6.4 Run-Level Summary Features

Per run, per series type, compute:

**Channel count metrics**

- Number of channels detected.
- Distribution across frequency ranges (e.g., low, mid, high).

**Channel strength metrics**

- Mean channel magnitude or energy.
- Maximum channel magnitude.

**Stability metrics**

- Fraction of windows containing a channel within a given band.
- Run-level continuity measures (e.g., average duration of a channel's presence).

**Control metrics**

- Total spectral energy.
- Low-vs-high frequency energy ratios.
- Simple shape descriptors (e.g., spectral centroid).

These features are used in all A/B/C comparisons, with the same definitions applied across conditions.

---

## 7. Falsification Logic

### 7.1 Distinctiveness: A1 vs A2

For each series type and chosen subset of features:

Compare the distribution of channel-related features (counts, strengths, stability) between A1 (high-performance) and A2 (low-performance) runs, under identical substrate parameters and analysis settings.

**Evidence against H1-strong in this dimension**:  
No consistent, robust differences between A1 and A2 beyond pre-declared negligible effect thresholds when the implementation performs appropriate non-fabricated statistical checks.

### 7.2 Specificity: B1 vs B2/B3

For any feature that appears to differ between A1 and A2:

Compute the same features for B2 (phase-scrambled) and B3 (IID) synthetic controls derived from B1 or matched distributions.

**Evidence against H1-strong**:  
Differences observed for B1 (structured) vs A2 fail to persist when comparing B1 to B2/B3; i.e., channels are indistinguishable from those appearing in scrambled or IID series under the same detection rules.

### 7.3 Causal Necessity: C1 vs C2 vs C3

For high-performance configurations where channels have been inferred:

Compare task performance metrics across C1, C2, and C3 runs.

**Evidence against H1-strong**:

- C2 (spectral knock-out) does not degrade performance more than C3 (non-channel perturbations), or
- Performance remains essentially unchanged across C1/C2/C3, indicating that inferred channels are not necessary in the strong sense.

Across all three axes, repeated failure to support H1-strong constitutes grounds for rejecting it in the evaluated configuration.

---

## 8. Implementation Skeleton and Repository Integration

### 8.1 Experiment Folder and Naming

Create a new experiment folder under `experiments/` that includes the PR number for traceability:

```
experiments/spectral_signal_routing_fft_cluster_pr_174_2026_01_14/
├── TECH_SPEC.md              # This design specification (authoritative design doc)
├── experiment_specification.json  # Concrete configuration for task/substrate (to be implemented)
├── manifest.yml               # Experiment metadata and tracking
├── TODO_IMPLEMENTATION.md     # Implementation checklist and status tracking
└── notes/                     # Optional design notes, sketches, scratchpad
```

The folder name encodes:

- Hypothesis domain: `spectral_signal_routing_fft_cluster`
- PR linkage: `pr_174`
- Date: `2026_01_14`

This follows the existing `experiments/…` pattern in the repo (e.g., `clustering_vs_fitness_experiment_2026_01_10`), while attaching the new work to PR #174 for auditability.

### 8.2 TECH_SPEC.md Role

**TECH_SPEC.md** (this document) is the source of truth for:

- Hypothesis articulation.
- Condition definitions.
- Data collection requirements.
- Spectral analysis pipeline.
- Falsification logic.
- Scope and alignment with EDE PM constraints.

Implementation documents must reference this spec and record any deviations explicitly.

### 8.3 experiment_specification.json (To Be Implemented)

This JSON file will define:

- Task name and variant.
- Substrate parameters (array size, dimensionality, update rules).
- Run configuration (step budget, number of runs per condition, seeding strategy).
- Which cell-level fields are exported and how they map to analysis series.

No field values are specified here to avoid fabricating data; they must be defined during actual implementation in line with tasks present in the repository.

### 8.4 TODO_IMPLEMENTATION.md Checklist (Initial Stub)

The following checklist should be instantiated and refined in `TODO_IMPLEMENTATION.md`:

1. **Task Wiring**
   - [ ] Select an existing EDE task with clear success/failure definition.
   - [ ] Implement or reuse a runner capable of repeated seeded runs with per-step snapshot recording.

2. **Trace Recorder**
   - [ ] Implement a generic trace recorder that outputs:
     - Per-step cell-level fields.
     - Per-step task-level metrics.
   - [ ] Document schema for downstream analysis.

3. **Series Extraction**
   - [ ] Define which series (x_task, x_global, x_cell_region) are extracted.
   - [ ] Implement extraction without referencing outcome labels.

4. **Spectral Analysis Harness**
   - [ ] Implement windowing, detrending, tapering, and FFT.
   - [ ] Implement channel detection logic with configurable thresholds.
   - [ ] Implement per-run feature aggregation.

5. **Condition Orchestration**
   - [ ] Implement or script generation of A1/A2/B1/B2/B3/C1/C2/C3 sets.
   - [ ] Ensure reproducible configuration and seeding.

6. **Collation and Downstream Hooks**
   - [ ] Implement scripts to merge per-run features into summary tables for statistical analysis.
   - [ ] Ensure outputs are neutral (no pre-labeled "success"/"failure" columns) until analysis is ready.

7. **Documentation Hooks**
   - [ ] Add references from PR #174 description and related findings documents to this experiment folder.
   - [ ] Ensure future FINDINGS.md (if created) links back to TECH_SPEC.md.

Implementation must avoid fabricating any numeric or statistical outcomes; this checklist is structural only.

---

## 9. Current PM Snapshot (Emergent Doom Engine Space)

- **High-level goal**: Evaluate and, if appropriate, falsify strong claims about spectral signal routing via FFT-compatible clustering on EDE substrates.
- **Current sub-goal**: Establish a technically precise, repository-integrated experiment spec and folder structure keyed to PR #174.
- **You are here**: The design and skeleton folder layout are defined; actual Java implementation and analysis scripts remain to be written under this spec.

**Next concrete step**:

Create `experiments/spectral_signal_routing_fft_cluster_pr_174_2026_01_14/` in the repo on branch `feat-sat-clustering` and add:

- TECH_SPEC.md populated with this document.
- Empty `experiment_specification.json` and `TODO_IMPLEMENTATION.md` scaffolds with the checklist headers described above.

All subsequent work on this experiment should be traceable to H1-strong falsification and the Emergent Doom Engine's emergent "doom" dynamics, with scope creep flagged whenever new ideas or artifacts cannot be cleanly mapped back to this spec or to EDE_PM_INSTRUCTIONS.md in this Space.

---

**Document Version**: 1.0  
**Last Updated**: 2026-01-14  
**Status**: Design Complete, Awaiting Implementation
