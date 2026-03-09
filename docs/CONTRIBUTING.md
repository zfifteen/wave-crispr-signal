# Contributing to wave-crispr-signal

This repository uses a change-scoped workflow. Validation is introduced per change, not inherited from historical suites.

## Workflow

1. Create a branch for the change.
2. Implement scoped code updates.
3. Update docs when behavior or interfaces change.
4. Add validation for what changed.

## Validation Contract

For behavior changes, default expectation is:
- one targeted automated test, and
- one manual verification command with observed outcome.

If automated coverage is not appropriate for the change, document the manual validation clearly.

For changes affecting model scoring/evaluation claims:
- run Gate v3 (`python3 validation/ontarget_gate/scripts/run_gate_v3.py`),
- include the resulting decision and reason from `validation/ontarget_gate/outputs/gate_results_v3.json`,
- clearly separate decision-grade clean-lane conclusions from exploratory raw-lane context.

## Pull Request Checklist

- Explain what changed and why.
- Include exact validation commands run.
- Include observed outcomes.
- Call out known limitations or follow-up work.
