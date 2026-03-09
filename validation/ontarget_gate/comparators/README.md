# Comparator Adapters

Gate v3 comparator adapters live here.

- `baseline_c`: `crisprpred` (required)
- `baseline_d`: reserved for future DeepSpCas9 integration

All required comparators must pass:
1. local-env subprocess availability,
2. checksum validation,
3. self-check fixtures,
4. overlap-audit prerequisites.

Gate runner never imports comparator dependencies into its own interpreter.
