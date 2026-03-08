# Contributing to wave-crispr-signal

This repository uses a reset workflow: testing and validation are introduced per change, not inherited from historical suites.

## Quick Start

1. Install dependencies: `pip install -r requirements.txt`
2. Create a branch for your change.
3. Implement the change.
4. Add change-coupled validation:
   - unit/integration test when behavior is stable and testable, or
   - explicit manual validation steps for exploratory work.

## Contribution Rules

- Keep changes scoped to CRISPR functionality.
- Prefer small, reviewable commits.
- Update docs when behavior or interfaces change.
- Do not treat historical tests or historical policy docs as blockers.

## Validation Contract

For each PR, include the exact validation you ran for that change.

Examples:

```bash
# Example targeted test
python -m pytest -q path/to/new_test.py

# Example manual validation
python applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA
```

## Pull Request Checklist

- Explain what changed and why.
- Describe validation run for this change.
- Note follow-up work if validation or CI is intentionally deferred.

## CI Status

CI is intentionally absent during this reset period. A new minimal CI will be introduced later and will gate only active, change-scoped checks.
