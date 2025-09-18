# Repository Guidelines

## Project Structure & Module Organization
- Core signal and geodesic logic lives in `z_framework.py`, `topological_analysis.py`, and the `applications/` suite for CRISPR tooling.
- Extended mathematical features are exposed via `invariant_features.py` and domain-specific modules under `modules/` and `experiments/`.
- Validation notebooks, reports, and generated outputs sit in `docs/`, `results/`, and `proof_pack/`.
- Automated checks and fixtures are in `tests/`, using matching data from `test_data/` and `proof_pack/data/`.

## Build, Test, and Development Commands
- `python -m venv .venv && source .venv/bin/activate` to create an isolated environment.
- `pip install -r requirements.txt` installs the pinned scientific stack (Python 3.12+).
- `python -m pytest -q` runs the entire regression suite.
- `python proof_pack/run_validation.py` executes the full Wave-CRISPR validation pack (long-running, uses bundled CSVs).
- `python applications/crispr_cli.py design ATCG... -n 5` demonstrates the CLI workflow end to end.

## Coding Style & Naming Conventions
- Python files use 4-space indentation, type annotations where practical, and descriptive function names.
- Follow Ruff configuration (`pyproject.toml`) with an 88-character line limit and E402 ignored only where path bootstrapping is required.
- Name modules with snake_case; classes use CapWords; CLI entry points should expose a `main()` guard.
- Prefer concise logging via the shared `logging` module; avoid print statements outside throwaway scripts.

## Testing Guidelines
- Tests are written with `pytest` and live under `tests/` using `test_<feature>.py` naming.
- When adding new functionality, create focused unit tests plus an integration scenario when the feature touches CRISPR pipelines or proof-pack assets.
- Keep fixtures lightweight; reuse sample inputs from `test_data/` rather than generating large random payloads.

## Commit & Pull Request Guidelines
- Craft commits as small, logically scoped units with imperative subject lines ("Add geodesic variance check").
- Reference relevant issue IDs in the body when available and note any validation commands executed.
- Pull requests should summarize the motivation, highlight affected modules/paths, attach screenshots or CLI transcripts for user-visible changes, and list manual or automated test results.
