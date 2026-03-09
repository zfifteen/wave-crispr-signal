# Board Funding Addendum (Conditional Approval Controls)

Version: `v1`  
Status: `Active (pending board sign-off)`  
Date locked: `2026-03-09`

## Purpose

This addendum is attached to the 4-week recovery plan and defines pre-registered controls for conditional funding approval.

- Funding is authorized only when this addendum is attached to the active plan and signed.
- Scope remains technical recovery under frozen Gate v3 decision conditions.

## Success and Failure Criteria (Pre-Registered)

Primary metric: Spearman delta vs `baseline_c_crisprpred` on clean-lane holdouts.

`GO` requires all:

1. `primary_holdout delta >= +0.01`
2. `external_holdout delta >= +0.01`
3. No material subgroup harm in pre-declared decision subgroups (`PMS2`, `MSH6`) relative to baseline directionality.

`NO-GO` if any:

1. Either holdout delta is negative.
2. Subgroup harm is material.
3. Hard preconditions fail and remain unresolved inside approved cycle bounds.

`ONE FINAL SHORT CYCLE` eligibility (single exception path):

- Both holdout deltas are positive, but one or both miss `+0.01` by `<= 0.005`.
- Diagnostics isolate one intervention class with high-confidence fixability.
- Requires explicit re-approval memo before execution.

## Governance Bounds (Hard Caps)

- Duration cap: `4 weeks`
- Team cap: `2 FTE` equivalent
- Experiment cap: `4 ablations` (`A0` to `A3`)
- Compute cap: fixed budget envelope in board packet (`$CAP_USD` or `GPU_HOURS_CAP`)
- Scope expansion (new feature families, new datasets, comparator policy changes) requires re-approval.

## Week-1 Early Off-Ramp (Mandatory)

At end of Week 1, publish a checkpoint artifact with `CONTINUE` or `STOP`.

Immediate `STOP` if diagnostics show structural in-scope infeasibility, for example:

- persistent rank collapse across both holdouts,
- no directional dev lift from any approved intervention class,
- or preconditions that cannot be resolved within cycle caps.

## Short-Cycle Constraint (If Invoked)

- Maximum duration: `10 business days`
- Scope: one pre-scoped intervention class only
- Approval: board re-approval required
- Failure outcome: automatic `NO-GO` and pivot/termination recommendation

## Required Board Reporting Package

1. Baseline vs final metric table (primary/external holdouts + subgroup rows).
2. Precondition integrity table (overlap, reproducibility, fail-closed checks).
3. Budget burn vs cap table.
4. Final decision memo with criterion traceability: `GO` / `NO-GO` / `Short-cycle request`.

## Sign-Off

- Plan owner: `________________`
- Technical reviewer: `________________`
- Finance/board delegate: `________________`
- Approval date: `________________`
