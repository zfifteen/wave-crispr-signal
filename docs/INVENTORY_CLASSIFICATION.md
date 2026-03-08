# Repository Inventory Classification

| Subsystem | Classification | Retention Rule | Notes |
|---|---|---|---|
| `wave_crispr_signal/` | `core` | Keep in mainline | Core CRISPR package APIs. |
| `applications/` | `core` | Keep in mainline | User-facing CRISPR CLIs/workflows. |
| `experiments/signal_theoretic_crispr/` | `core` | Keep in mainline | Active CRISPR experiments. |
| `proof_pack/` | `core` | Keep in mainline | Reproducibility and validation scripts. |
| `tests/` | `core` | Keep in mainline | Regression and policy checks. |
| `data/`, `test_data/` | `support` | Keep in mainline | Fixtures and datasets. |
| `bin/`, `scripts/`, `tools/`, `configs/` | `support` | Keep in mainline | Tooling and operational scripts. |
| `docs/` | `support` | Keep in mainline | Canonical docs for current scope. |
| `legacy/projects/` | `research-legacy` | Retain under `legacy/` | Archived side projects moved from root. |
| `legacy/experiments/` | `research-legacy` | Retain under `legacy/` | Non-CRISPR exploratory experiments. |
| `legacy/docs/` | `research-legacy` | Retain under `legacy/` | Non-CRISPR historical docs. |
| `legacy/artifacts/` | `research-legacy` | Retain under `legacy/` | Historical reports and generated artifacts. |
| Root generated artifacts (`*.png`, ad-hoc outputs) | `delete/regenerate` | Remove from root and ignore | Must not live in repo root. |
| IDE/system/cache artifacts | `delete/regenerate` | Remove/ignore | `.idea/`, `__pycache__/`, `.DS_Store`, etc. |

## Top-Level Purpose Statements

- `wave_crispr_signal/`: reusable CRISPR signal-analysis package.
- `applications/`: executable CRISPR tools.
- `experiments/`: active CRISPR experiments + legacy wrappers.
- `proof_pack/`: reproducibility harness.
- `docs/`: canonical project docs and policy links.
- `legacy/`: archived non-CRISPR research and side-project history.
- `tests/`: regression + repository policy checks.
