# Fast Go/No-Go On-Target Validation Gate

This pack runs a decision-oriented validation gate for on-target efficiency claims.

## Objective

Determine whether the current model shows real predictive lift over a simple baseline under locked data and split manifests.

The active decision protocol is `PROTOCOL_V2.md`.

## Data Sources (public)

- `doench2016-Gecko1-LentiCrispr_hg19.scores.tab` (dev pool)
- `doench2016_hg19.scores.tab` (primary holdout)
- `hart2016-HelaLib1Avg.scores.tab` (external holdout)

All are fetched from `maximilianh/crisporPaper` `effData` tabular score files.

## Locked Schema

Fields used downstream:
- `guide`
- `label`
- `source`
- `gene_or_target_group`
- `target_context` (optional)

## Run

```bash
cd /Users/velocityworks/IdeaProjects/wave-crispr-signal
python3 validation/ontarget_gate/scripts/run_gate_v2.py
```

## Outputs

Written to `validation/ontarget_gate/outputs/`:

- `locked_dataset_manifest_v2.json`
- `locked_split_manifest_v2.csv`
- `locked_schema_manifest_v2.csv`
- `gate_results_v2.json`
- `gate_report_v2.md`

## Decision Rule (v2)

Locked criteria are defined in `PROTOCOL_V2.md` and enforced by `run_gate_v2.py`.

The previous v1 outputs (`gate_results.json`, `gate_report.md`) are historical and non-authoritative.
