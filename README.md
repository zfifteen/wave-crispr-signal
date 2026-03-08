# wave-crispr-signal

Signal-theoretic CRISPR analytics toolkit for guide scoring, disruption analysis, and reproducible validation.

## Project Charter

### Mission
Build and maintain a focused CRISPR signal-analysis toolkit that turns sequence-level signals into reproducible guide-quality and disruption metrics.

### Intended Users
- Computational biologists evaluating gRNA candidates.
- Research engineers building CRISPR scoring and validation pipelines.
- Reproducibility reviewers validating CRISPR signal claims.

### Primary Outputs
- CRISPR scoring APIs and CLIs (`applications/`, `wave_crispr_signal/`).
- Methodology and usage documentation (`docs/`).
- Demo MVP interface (`web_apps/demo_mvp/`).

### Non-Goals
- This repository is not a general cross-domain experiment lab.

## Working Model

- Keep active CRISPR work in mainline directories.
- Keep generated artifacts and machine-specific files out of versioned root.
- Introduce tests only when a specific add/remove/alter change requires them.
- Historical tests are not a blocking source of truth.
- CI is intentionally absent for now and will be reintroduced later as a minimal, change-scoped gate.

## Active Runtime Layout

- `wave_crispr_signal/`: core CRISPR feature and scoring package.
- `applications/`: primary CLIs and runtime-facing application code.
- `web_apps/demo_mvp/`: preserved web demo surface.
- `scripts/invariant_features.py`: invariant feature dependency used by active CLI flows.
- `docs/`: canonical docs for current scope.

## Archive Layout

- `archive/code/`: archived non-primary code retained for provenance.
- `docs/ARCHIVE_INDEX.md`: index of archival references.

## Contributor Decision Tree

1. Does this change advance CRISPR scoring/design/validation for active runtime surfaces?
   - Yes: implement in mainline active paths.
   - No: move to archive paths or keep out of repo.
2. Is the output generated or environment-specific?
   - Yes: keep it out of versioned root and add ignore rules if needed.
3. Does the change alter behavior?
   - Yes: add change-coupled validation (targeted automated test + manual verification command).

## Quick Start

```bash
python -m pip install -r requirements.txt

# Example: phase-weighted scorecard
python applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA
```

## Canonical Docs

- `docs/INDEX.md`: documentation index and navigation.
- `docs/CONTRIBUTING.md`: contribution workflow.
- `docs/ARCHIVE_INDEX.md`: archival reference map.

## License

MIT
