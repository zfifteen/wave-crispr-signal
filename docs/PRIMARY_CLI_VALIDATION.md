# Primary CLI Manual Validation

Use these commands for manual verification when refactoring preserved CLI surfaces.

## 1) Phase-Weighted Scorecard CLI

```bash
python applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA
```

Expected: exit 0 and output header containing `Phase-Weighted CRISPR Scorecard`.

## 2) CRISPR CLI

```bash
python applications/crispr_cli.py --help
```

Expected: exit 0 and help text describing CRISPR guide design commands.

## 3) Genomic Disruption API CLI

```bash
python applications/genomic_disruption_api.py --help
```

Expected: exit 0 and help text containing `Genomic Disruption Analyzer`.
