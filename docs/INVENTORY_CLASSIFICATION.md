# Repository Inventory Classification

| Subsystem | Classification | Retention Rule | Notes |
|---|---|---|---|
| `wave_crispr_signal/` | `core` | Keep in mainline | Core CRISPR package APIs. |
| `applications/` | `core` | Keep in mainline | Primary user-facing CRISPR CLIs/workflows. |
| `web_apps/demo_mvp/` | `core` | Keep in mainline | Preserved demo surface. |
| `scripts/invariant_features.py` | `support` | Keep in mainline | Active dependency for guide designer scoring. |
| `scripts/run_tests.py` | `support` | Keep in mainline | Change-scoped validation helper. |
| `docs/` | `support` | Keep in mainline | Canonical docs for current scope. |
| `archive/code/` | `archive` | Keep for provenance | Non-primary code moved out of active runtime paths. |
| Root generated artifacts (`*.png`, ad-hoc outputs) | `delete/regenerate` | Remove from root and ignore | Must not live in repo root. |
| IDE/system/cache artifacts | `delete/regenerate` | Remove/ignore | `.idea/`, `__pycache__/`, `.DS_Store`, etc. |

## Top-Level Purpose Statements

- `wave_crispr_signal/`: reusable CRISPR signal-analysis package.
- `applications/`: executable primary CRISPR tools.
- `web_apps/demo_mvp/`: public demo interface.
- `docs/`: canonical project docs.
- `archive/code/`: archived non-primary implementations.

## Testing Baseline

There is no inherited baseline test suite. Tests are introduced as needed when a change adds, removes, or alters behavior.
