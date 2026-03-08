# Documentation Index

This is the canonical documentation index for the active CRISPR scope.

## Core (Current)

- `WAVE_CRISPR_SIGNAL_TLDR.md` — project overview and intent.
- `PHASE_WEIGHTED_SCORECARD.md` — phase-weighted scoring method.
- `FFT_GOLDEN_RATIO_CRISPR.md` — FFT and phase-weighting details.
- `TOPOLOGICAL_ANALYSIS.md` — topology bridge and analysis notes.
- `bridge_kappa_to_theta_prime.md` — curvature/phase mathematical bridge.
- `REPRODUCIBILITY_VALIDATION.md` — reproducibility guidance.
- `INVARIANT_FEATURES.md` — invariant feature definitions.
- `CONTRIBUTING.md` — contribution workflow.

## Structure and Policy

- `INVENTORY_CLASSIFICATION.md` — top-level subsystem classification.
- `../.github/REPOSITORY_POLICY.md` — enforced repository policy.

## Legacy (Archived)

- `../legacy/README.md` — archival scope and navigation.
- `../legacy/docs/` — non-CRISPR historical docs retained for reproducibility.

## Running Key Workflows

- Tests: `python -m pytest -q`
- Primary score workflow: `python applications/phase_weighted_scorecard_cli.py score --guide <GUIDE>`
- Validation proof run: `python proof_pack/run_validation.py`
