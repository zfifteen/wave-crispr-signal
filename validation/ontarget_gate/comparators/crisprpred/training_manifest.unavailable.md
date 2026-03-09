# Training Manifest Status

A sequence-level training manifest for CRISPRpred has not been obtained yet.

This repository intentionally runs fail-closed for comparator gating when the manifest is unavailable.

To unblock gate decisions, provide a checksummed `training_manifest.fasta` and update `training_manifest_status.json` to `available: true`.
