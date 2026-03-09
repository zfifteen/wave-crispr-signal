# Recovery Ablation v1

Generated: 2026-03-09T17:42:19.205009+00:00

Selection rule: max dev_val delta vs baseline_c.

Winner: **A2_gc_debias_plus_zblend**

## Candidate metrics (delta model - baseline_c, Spearman)
- A0_current_model: dev=0.0053, primary=-0.0165, external=-0.0166
- A1_zblend_model_plus_baseline_a: dev=0.0053, primary=-0.0165, external=-0.0166
- A2_gc_debias_plus_zblend: dev=0.0484, primary=-0.0119, external=-0.0172
- A3_stacked_ridge: dev=0.0036, primary=0.0260, external=-0.0588
- A4_pairwise_rank_logit: dev=-0.0268, primary=0.0293, external=-0.0575

## Winner holdout deltas
- primary: -0.0119
- external: -0.0172
