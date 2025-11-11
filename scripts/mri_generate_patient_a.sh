#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="${1:-}"
OUTPUT_DIR="${2:-docs/mri/reports/patient_a}"

if [[ -z "$INPUT_DIR" || ! -d "$INPUT_DIR" ]]; then
  echo "ERROR: Provide the DICOM input folder (e.g., ../wave-crispr-signal/data/DICOM)." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Default anonymization ON (recommended). To disable for authorized research, you
# must export MRI_TOOL_ALLOW_PHI=AUTHORIZED_RESEARCH_USE and add --no-anonymize.
# (Enforced by CLI.)
python3 mri_enhancement_cli.py \
  -i "$INPUT_DIR" \
  -o "$OUTPUT_DIR" \
  --log-level INFO