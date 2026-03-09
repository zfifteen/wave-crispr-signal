# On-Target Validation Gate Protocol v3

Version: `v3`  
Status: `Active`  
Date locked: `2026-03-09`

## Purpose

Gate v3 adds a required external comparator (`baseline_c_crisprpred`) via a pluggable comparator framework.

- v1/v2 outputs are retained as historical artifacts.
- v3 is authoritative for external-comparator go/no-go claims.

## Comparator Slots

- `baseline_c`: required now (`crisprpred`)
- `baseline_d`: reserved for future DeepSpCas9 integration

## Training Manifest Sourcing (Pre-Implementation Verification)

Before decision-grade comparator execution:

1. identify CRISPRpred training sequence source,
2. produce sequence-level `training_manifest.*` with checksum,
3. update `training_manifest_status.json` to `available: true`.

Current state at lock time: training manifest unavailable.

Fail-closed policy:

- if training manifest unavailable -> decision is `INCONCLUSIVE`.
- no fallback overlap heuristic is active in v3.

## Data and Splits

Public scored datasets remain:

1. `doench2016_gecko1_lenticrispr_hg19`
2. `doench2016_hg19`
3. `hart2016_hela_lib1_avg`

Leakage-safe grouped split policy remains from v2.

Two manifests are generated:

- `decision_split_manifest_v3_clean.csv` (authoritative; overlap-sanitized)
- `exploratory_split_manifest_v3_raw.csv` (non-authoritative; original split assignment)

## Comparator Environment Isolation

Comparator execution mode is fixed:

- gate runner invokes comparator via subprocess with explicit interpreter:
  - `validation/ontarget_gate/comparators/crisprpred/.venv/bin/python`
- gate runner must not import comparator dependencies into its own process.

If comparator venv/interpreter is missing or corrupted, gate fails closed (`INCONCLUSIVE`).

## Preconditions (Hard)

All must pass before baseline_c criteria are evaluated:

1. split integrity checks pass,
2. comparator self-check passes,
3. comparator checksum verification passes,
4. overlap audit passes (exact-sequence overlap threshold = 0),
5. minimum holdout size met:
   - `primary_holdout >= 200`
   - `external_holdout >= 200`

Minimum holdout behavior:

- if either holdout has fewer than 200 scored guides -> `INCONCLUSIVE`.
- overlap audit runs before decision-grade scoring; if clean-manifest overlap fails, authoritative scoring is skipped.

## Metrics

Primary comparison metric: Spearman rank correlation.

Secondary metric: MSE (diagnostic only).

Score normalization policy:

- none required for Spearman-based gating,
- tie handling uses SciPy default average ranks,
- missing comparator predictions for decision holdouts => precondition failure (`INCONCLUSIVE`).

## Decision Logic

### Dev precondition

- if dev validation delta vs Baseline A has `CI_hi < -0.01`, decision is `NO-GO`.

### Baseline C gate criteria (fixed-delta only)

- `primary_holdout`: `delta(model - baseline_c) >= +0.01`
- `external_holdout`: `delta(model - baseline_c) >= +0.01`

Decision mapping:

- `GO`: all preconditions pass and both baseline_c thresholds pass
- `NO-GO`: dev material-negative precondition fails, or baseline_c threshold(s) fail
- `INCONCLUSIVE`: comparator/overlap/manifest/min-N preconditions fail

## Bootstrap Confidence Intervals (Diagnostic-Only Policy)

Current phase (v3):

- CIs are computed and reported for diagnostics.
- CIs do not affect gate decisions.

Potential v4 upgrade path:

- CI-aware conjunctive gating (for example `delta >= +0.01 AND CI_lo > 0`) requires protocol amendment.

## Decision vs Exploratory Lanes

- Decision-grade lane:
  - clean manifest only,
  - drives `GO` / `NO-GO`.
- Exploratory lane:
  - raw manifest only,
  - reported for context,
  - never used in decision criteria.
