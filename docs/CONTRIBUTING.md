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

## Pull Request Checklist

- Explain what changed and why.
- Include exact validation commands run.
- Include observed outcomes.
- Call out known limitations or follow-up work.
