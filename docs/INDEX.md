# Documentation Index

This is the canonical documentation index for the active CRISPR scope.

## Core

- `WAVE_CRISPR_SIGNAL_TLDR.md` - project overview and intent.
- `PHASE_WEIGHTED_SCORECARD.md` - phase-weighted scoring method.
- `FFT_GOLDEN_RATIO_CRISPR.md` - FFT and phase-weighting details.
- `TOPOLOGICAL_ANALYSIS.md` - topology bridge and analysis notes.
- `REPRODUCIBILITY_VALIDATION.md` - reproducibility guidance.
- `INVARIANT_FEATURES.md` - invariant feature definitions.
- `CONTRIBUTING.md` - contribution workflow.
- `INVENTORY_CLASSIFICATION.md` - current top-level classification.

## Running Key Workflows

- Primary score workflow: `python applications/phase_weighted_scorecard_cli.py score --guide <GUIDE>`
- Validation proof run: `python proof_pack/run_validation.py`

## Testing and CI

- Validation is change-scoped: add tests or manual validation steps only for the behavior you change.
- CI is intentionally absent during reset and will be reintroduced later as a minimal gate.
