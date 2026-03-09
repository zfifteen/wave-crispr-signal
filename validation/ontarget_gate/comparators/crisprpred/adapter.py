from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence

from ..base import ComparatorRecord, ComparatorResult


class CRISPRpredAdapter:
    name = "baseline_c_crisprpred"
    version = "v0.1-scaffold"

    def __init__(self) -> None:
        self.root = Path(__file__).resolve().parent
        self.venv_python = self.root / ".venv" / "bin" / "python"
        self.scorer = self.root / "scorer.py"
        self.self_check_script = self.root / "self_check.py"
        self.checksums = self.root / "checksums.sha256"
        self.provenance_doc = self.root / "PROVENANCE.md"
        self.training_manifest = self.root / "training_manifest.fasta"
        self.training_manifest_status = self.root / "training_manifest_status.json"

    def _run_subprocess(self, script: Path, payload: Dict[str, object]) -> ComparatorResult:
        if not self.venv_python.exists():
            return ComparatorResult(
                ok=False,
                message="comparator_venv_missing",
                details={"venv_python": str(self.venv_python)},
            )

        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as in_f:
            json.dump(payload, in_f)
            in_path = Path(in_f.name)
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as out_f:
            out_path = Path(out_f.name)

        cmd = [str(self.venv_python), str(script), "--input", str(in_path), "--output", str(out_path)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            return ComparatorResult(
                ok=False,
                message="comparator_subprocess_failed",
                details={
                    "cmd": cmd,
                    "returncode": proc.returncode,
                    "stdout": proc.stdout[-1000:],
                    "stderr": proc.stderr[-1000:],
                },
            )

        try:
            obj = json.loads(out_path.read_text(encoding="utf-8"))
        except Exception as exc:  # noqa: BLE001
            return ComparatorResult(
                ok=False,
                message="comparator_output_parse_failed",
                details={"error": str(exc), "output_path": str(out_path)},
            )

        return ComparatorResult(
            ok=bool(obj.get("ok", False)),
            message=str(obj.get("message", "unknown")),
            details=dict(obj),
        )

    def self_check(self) -> ComparatorResult:
        payload = {
            "checksums_path": str(self.checksums),
            "provenance_path": str(self.provenance_doc),
            "training_manifest_path": str(self.training_manifest),
            "training_manifest_status_path": str(self.training_manifest_status),
        }
        return self._run_subprocess(self.self_check_script, payload)

    def provenance(self) -> Dict[str, object]:
        status: Dict[str, object] = {
            "name": self.name,
            "version": self.version,
            "provenance_path": str(self.provenance_doc),
            "checksums_path": str(self.checksums),
            "training_manifest_path": str(self.training_manifest),
            "training_manifest_status_path": str(self.training_manifest_status),
            "venv_python": str(self.venv_python),
        }
        if self.training_manifest_status.exists():
            try:
                status["training_manifest_status"] = json.loads(
                    self.training_manifest_status.read_text(encoding="utf-8")
                )
            except Exception as exc:  # noqa: BLE001
                status["training_manifest_status"] = {"ok": False, "error": str(exc)}
        return status

    def predict_batch(self, records: Sequence[ComparatorRecord]) -> List[float]:
        payload = {
            "records": [
                {
                    "guide_seq": r.guide_seq,
                    "target_context": r.target_context,
                }
                for r in records
            ]
        }
        result = self._run_subprocess(self.scorer, payload)
        if not result.ok:
            raise RuntimeError(f"CRISPRpred scorer failed: {result.message} :: {result.details}")

        raw_scores = result.details.get("scores", [])
        if not isinstance(raw_scores, list) or len(raw_scores) != len(records):
            raise RuntimeError(
                "CRISPRpred scorer returned invalid score list length: "
                f"expected {len(records)}, got {len(raw_scores) if isinstance(raw_scores, list) else 'non-list'}"
            )

        scores = [float(v) for v in raw_scores]
        return scores
