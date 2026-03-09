# On-Target Validation Gate Protocol v2

Version: `v2`  
Status: `Active`  
Date locked: `2026-03-09`

## Purpose

This protocol defines the decision gate for on-target efficiency claims. It prioritizes fast risk reduction over publication-grade benchmarking.

`run_gate.py` outputs (`gate_results.json`, `gate_report.md`) are retained as historical artifacts and are non-authoritative for decisions.

## Locked Datasets

Public sources (from `maximilianh/crisporPaper/effData`):

1. `doench2016-Gecko1-LentiCrispr_hg19.scores.tab` (`doench2016_gecko1_lenticrispr_hg19`)
2. `doench2016_hg19.scores.tab` (`doench2016_hg19`)
3. `hart2016-HelaLib1Avg.scores.tab` (`hart2016_hela_lib1_avg`)

## Grouping and Split Policy (Leakage-Safe)

Grouping hierarchy is source-aware and deterministic:

1. True source metadata group when available:
- `doench2016_hg19`: group from `guide` prefix before `-` (gene-level key)
- `hart2016_hela_lib1_avg`: group from `name` prefix before `-` (gene-level key)
2. Fallback when true group is unavailable:
- SHA1 hash of `longSeq100Bp` context (`ctxsha1:<hash>`)
3. If neither is available:
- mark row as weak-group metadata

Split assignment:

- `doench2016_hg19` (strong groups only): deterministic 60/20/20 by group hash:
  - `dev_train`
  - `dev_val`
  - `primary_holdout`
- `hart2016_hela_lib1_avg` (strong groups only):
  - `external_holdout`
- Weak-group rows are excluded from decision splits and assigned to non-decision splits (`aux_train_weak` or `excluded_weak_group`).

Split integrity is fail-closed:

- no strong group overlap across `dev_train`, `dev_val`, `primary_holdout`, `external_holdout` using source-qualified keys (`source::group`) to avoid false collisions across independent datasets
- no weak-group rows in any decision split
- required decision splits must be non-empty

## Metrics and Aggregation

Primary metric: Spearman rank correlation.

Secondary metric: MSE (diagnostic only, not decision-primary).

Evaluation is per-source/per-split first. No pooled mixed-scale headline Spearman is used.

Reported deltas:

- `delta(model - baseline_a)` Spearman with 95% bootstrap CI per split
- `delta(model - baseline_b)` Spearman per split
- weighted aggregate `delta(model - baseline_a)` over decision holdouts with 95% bootstrap CI

## Baselines and Scope

Baseline A: lightweight sequence baseline (GC, homopolymer, 2-mer features + Ridge).

Baseline B: non-phase comparator (`CRISPRGuideDesigner.calculate_on_target_score`).

Model under test: `DisruptionAnalyzer.score_guide(...)[\"disruption_score\"]`.

Training scope:

- Baseline A is trained on `dev_train` plus auxiliary non-decision training rows.
- No tuning is performed on decision holdouts.

## Decision Logic

Hard precondition:

- If dev validation (`dev_val`) `delta(model - baseline_a)` is materially negative beyond CI tolerance (`CI_hi < -0.01`), decision is immediate `NO-GO`.

`GO` requires all:

1. Primary holdout delta vs Baseline A `>= +0.03` and CI excludes 0 (`CI_lo > 0`)
2. External holdout delta vs Baseline A `>= +0.01` and CI excludes 0 (`CI_lo > 0`)
3. Directional consistency (both holdout deltas > 0)
4. No underperformance vs Baseline B on either holdout
5. Aggregate holdout delta vs Baseline A > 0 and aggregate CI excludes 0

`NO-GO` if any:

1. Hard precondition fails
2. Any holdout delta vs Baseline A `<= 0`
3. Any holdout CI low bound `< -0.01`
4. Underperformance vs Baseline B on any holdout
5. Aggregate holdout CI low bound `< -0.01`

Else: `INCONCLUSIVE` (single follow-up cycle allowed; if still inconclusive, treat as no-go for further app investment).
