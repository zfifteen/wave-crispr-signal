# Legacy Archive

This directory contains non-CRISPR exploratory work retained for provenance and reproducibility.

## Policy

- Legacy content is **read-only by default**.
- New non-CRISPR exploration should be added here, not in mainline CRISPR paths.
- Mainline code should not depend on `legacy/` except explicit compatibility wrappers.

## Contents

- `projects/`: archived side-project directories.
- `experiments/`: non-CRISPR experiment modules and manifests.
- `docs/`: archived non-CRISPR documentation.
- `artifacts/`: generated reports and historical artifacts moved out of root.

## Compatibility

Some import paths in `experiments/` are thin wrappers to legacy modules so historical tests and scripts remain runnable.
