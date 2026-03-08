# Repository Inventory Classification

| Subsystem | Classification | Retention Rule | Notes |
|---|---|---|---|
| `wave_crispr_signal/` | `core` | Keep in mainline | Core CRISPR package APIs. |
| `applications/` | `core` | Keep in mainline | User-facing CRISPR CLIs/workflows. |
| `experiments/` | `core` | Keep in mainline | Active CRISPR experiments. |
| `proof_pack/` | `core` | Keep in mainline | Reproducibility and validation scripts. |
| `data/`, `test_data/` | `support` | Keep in mainline | Fixtures and datasets. |
| `bin/`, `scripts/`, `tools/`, `configs/` | `support` | Keep in mainline | Tooling and operational scripts. |
| `docs/` | `support` | Keep in mainline | Canonical docs for current scope. |
| Root generated artifacts (`*.png`, ad-hoc outputs) | `delete/regenerate` | Remove from root and ignore | Must not live in repo root. |
| IDE/system/cache artifacts | `delete/regenerate` | Remove/ignore | `.idea/`, `__pycache__/`, `.DS_Store`, etc. |

## Top-Level Purpose Statements

- `wave_crispr_signal/`: reusable CRISPR signal-analysis package.
- `applications/`: executable CRISPR tools.
- `experiments/`: active CRISPR experiments.
- `proof_pack/`: reproducibility harness.
- `docs/`: canonical project docs.
- `tools/` and `scripts/`: utility and operational helpers.

## Testing Baseline

There is no inherited baseline test suite. Tests are introduced as needed when a change adds, removes, or alters behavior.
