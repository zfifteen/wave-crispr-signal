#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.lock
python self_check.py --input <(python - <<'PY'
import json
print(json.dumps({
  "checksums_path": "checksums.sha256",
  "provenance_path": "PROVENANCE.md",
  "training_manifest_path": "training_manifest.fasta",
  "training_manifest_status_path": "training_manifest_status.json"
}))
PY
) --output /tmp/crisprpred_selfcheck.json || true
cat /tmp/crisprpred_selfcheck.json
