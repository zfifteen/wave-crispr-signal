# Primary CLI Validation

Use these commands as baseline manual checks for active runtime surfaces.

## 1) Phase-Weighted Scorecard CLI

```bash
python3 applications/phase_weighted_scorecard_cli.py --help
python3 applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA
```

Expected:
- help exits 0 and lists `score`, `batch`, `analyze`.
- `score` exits 0 and prints `Phase-Weighted CRISPR Scorecard`.

## 2) CRISPR CLI

```bash
python3 applications/crispr_cli.py --help
python3 applications/crispr_cli.py design ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGGATCGATC
python3 applications/crispr_cli.py score GACGATCGATCGATCGATCG --target-context AATGACGATCGATCGATCGATCGTGG
```

Expected:
- help exits 0 and lists `design`, `score`, `batch-score`.
- `design` returns JSON with guide candidates.
- `score` returns JSON with `on_target_score`.

## 3) Genomic Disruption API CLI

```bash
python3 applications/genomic_disruption_api.py --help
python3 applications/genomic_disruption_api.py score --guide GCTGCGGAGACCTGGAGAGA
python3 applications/genomic_disruption_api.py design --target ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGGATCGATC
```

Expected:
- help exits 0 and lists `score`, `batch`, `design`, `offtarget` commands.
- `score` exits 0 with JSON output (`success: true`).
- `design` exits 0 with JSON output (`success: true`).

Current limitation (expected behavior):

```bash
python3 applications/genomic_disruption_api.py batch --guides /tmp/guides.txt
python3 applications/genomic_disruption_api.py offtarget --guide GCTGCGGAGACCTGGAGAGA
```

Expected:
- each command exits non-zero and prints `Command <name> not yet implemented in CLI`.

## 4) Demo MVP Import Surface

```bash
python3 -c "from web_apps.demo_mvp.app import app; print(sorted([r.path for r in app.router.routes]))"
```

Expected route set includes:
- `/`
- `/health`
- `/info`
- `/analyze`
