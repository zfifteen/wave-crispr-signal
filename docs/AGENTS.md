# Repository Guidelines

## Project Structure
- Core CRISPR logic lives in `wave_crispr_signal/` and `applications/`.
- Preserved demo surface lives in `web_apps/demo_mvp/`.
- Non-primary code is retained under `archive/code/`.
- Canonical documentation lives in `docs/`.

## Development Commands
- `python -m venv .venv && source .venv/bin/activate`
- `pip install -r requirements.txt`
- Run change-specific validations only (tests and/or manual checks tied to the change).

## Coding Style
- Python: 4-space indentation, type annotations where practical.
- Follow Ruff config in `pyproject.toml`.
- Use snake_case modules and descriptive names.

## Validation Policy
- There is no inherited baseline test gate.
- Add tests only when your change adds/removes/alters behavior and automated coverage is appropriate.
- For exploratory work, provide explicit manual validation steps.

## PR Expectations
- Keep commits logically scoped.
- Summarize behavioral impact.
- Include exact validation evidence for the change.
