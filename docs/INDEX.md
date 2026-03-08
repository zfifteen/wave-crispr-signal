# Documentation Index

This is the canonical documentation index for the active CRISPR scope.

## Active Guidance

- `WAVE_CRISPR_SIGNAL_TLDR.md` - project overview and intent.
- `TOPOLOGICAL_ANALYSIS.md` - topology bridge and analysis notes.
- `REPRODUCIBILITY_VALIDATION.md` - reproducibility guidance.
- `INVARIANT_FEATURES.md` - invariant feature definitions.
- `CONTRIBUTING.md` - contribution workflow.
- `INVENTORY_CLASSIFICATION.md` - current top-level classification.
- `PRIMARY_CLI_VALIDATION.md` - manual verification commands for preserved CLIs.

## Running Key Workflows

- Phase-weighted CLI: `python applications/phase_weighted_scorecard_cli.py score --guide <GUIDE>`
- CRISPR CLI: `python applications/crispr_cli.py --help`
- Genomic disruption API CLI: `python applications/genomic_disruption_api.py --help`

## Testing and CI

- Validation is change-scoped: add tests or manual validation steps only for the behavior you change.
- CI is intentionally absent during reset and will be reintroduced later as a minimal gate.

## Historical Reference

- `ARCHIVE_INDEX.md` - archival document map (non-authoritative, non-blocking).
