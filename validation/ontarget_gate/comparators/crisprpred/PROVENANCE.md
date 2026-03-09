# CRISPRpred Comparator Provenance

## Source
- Comparator family: CRISPRpred
- Source repository: https://github.com/khaled-rahman/CRISPRpred
- Source commit (master at acquisition): `eee7ee8a9ba61c1afcb7a56d9893f8e661f41cc0`
- Source training file: `CRISPRpred/SupplementaryFiles/Supplementary File1.csv`

## Training Manifest
- Manifest path: `training_manifest.fasta`
- Build notes: `training_manifest_build.md`
- Extraction rule: `guide20 = 30mer[4:24]` (0-based, end-exclusive)
- Records: 5310

## License Compatibility
- Repository license metadata: `NONE_DECLARED` (GitHub API)
- This comparator is for internal validation-gate use.
- If external distribution is planned, legal review is required before shipping vendored assets.

## Modifications
1. Added deterministic scaffold scorer and fixture harness for Gate v3 comparator interface.
2. Added checksum verification and fail-closed self-check behavior.
3. Added sequence-level training manifest and overlap-audit inputs.

## Integrity
- File integrity is enforced via `checksums.sha256`.
- Comparator self-check must pass before gate execution can evaluate baseline_c criteria.

## Status
- Training manifest status is tracked in `training_manifest_status.json`.
