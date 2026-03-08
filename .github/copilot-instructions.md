# Copilot Instructions (Reset Baseline)

This repository is in a reset phase.

## Working Rules

- Focus only on CRISPR-relevant changes.
- Do not enforce historical repository policy files as blockers.
- Do not assume a baseline test suite exists.
- For each add/remove/alter behavior change, add change-scoped validation:
  - targeted automated test when appropriate, or
  - explicit manual validation steps.
- Keep generated artifacts and machine-specific files out of versioned root.

## CI

CI is intentionally absent during reset. Do not add CI-specific assumptions unless explicitly requested.
