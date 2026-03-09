#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(65536)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _parse_checksums(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            continue
        out[parts[1].strip()] = parts[0].strip()
    return out


def _run_fixtures(scorer_path: Path, fixtures_path: Path) -> tuple[bool, str]:
    fixtures = json.loads(fixtures_path.read_text(encoding="utf-8"))
    cases = fixtures.get("fixtures", [])
    tol = float(fixtures.get("tolerance", 1e-9))

    req = {"records": [{"guide_seq": c["seq"], "target_context": c.get("target_context", "")} for c in cases]}

    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as in_f:
        json.dump(req, in_f)
        in_path = Path(in_f.name)
    with NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as out_f:
        out_path = Path(out_f.name)

    proc = subprocess.run(["python3", str(scorer_path), "--input", str(in_path), "--output", str(out_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        return False, f"fixture_scorer_failed:{proc.stderr[-300:]}"

    obj = json.loads(out_path.read_text(encoding="utf-8"))
    if not obj.get("ok"):
        return False, f"fixture_scorer_not_ok:{obj.get('message')}"
    scores = obj.get("scores", [])
    if len(scores) != len(cases):
        return False, "fixture_score_len_mismatch"

    for i, (pred, case) in enumerate(zip(scores, cases)):
        exp = float(case["expected_score"])
        if abs(float(pred) - exp) > tol:
            return False, f"fixture_mismatch_index_{i}"
    return True, "fixtures_ok"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    req = json.loads(Path(args.input).read_text(encoding="utf-8"))
    checksums_path = Path(req["checksums_path"])
    provenance_path = Path(req["provenance_path"])
    training_manifest_path = Path(req["training_manifest_path"])
    training_manifest_status_path = Path(req["training_manifest_status_path"])

    details: dict[str, object] = {
        "checksums_path": str(checksums_path),
        "provenance_path": str(provenance_path),
        "training_manifest_path": str(training_manifest_path),
    }

    missing = [str(p) for p in [checksums_path, provenance_path, training_manifest_status_path] if not p.exists()]
    if missing:
        Path(args.output).write_text(json.dumps({"ok": False, "message": "required_files_missing", "missing": missing}), encoding="utf-8")
        return 2

    status = json.loads(training_manifest_status_path.read_text(encoding="utf-8"))
    details["training_manifest_status"] = status
    if not bool(status.get("available", False)):
        Path(args.output).write_text(json.dumps({"ok": False, "message": "training_manifest_unavailable", "details": details}), encoding="utf-8")
        return 3

    if not training_manifest_path.exists():
        Path(args.output).write_text(json.dumps({"ok": False, "message": "training_manifest_missing", "details": details}), encoding="utf-8")
        return 4

    checksums = _parse_checksums(checksums_path)
    base = checksums_path.parent
    mismatches = []
    for rel, expected in checksums.items():
        p = base / rel
        if not p.exists():
            mismatches.append({"path": rel, "error": "missing"})
            continue
        got = _sha256(p)
        if got != expected:
            mismatches.append({"path": rel, "expected": expected, "actual": got})

    details["checksum_mismatches"] = mismatches
    if mismatches:
        Path(args.output).write_text(json.dumps({"ok": False, "message": "checksum_mismatch", "details": details}), encoding="utf-8")
        return 5

    fixtures_ok, fixture_msg = _run_fixtures(base / "scorer.py", base / "test_fixtures.json")
    details["fixture_result"] = fixture_msg
    if not fixtures_ok:
        Path(args.output).write_text(json.dumps({"ok": False, "message": "fixture_validation_failed", "details": details}), encoding="utf-8")
        return 6

    Path(args.output).write_text(json.dumps({"ok": True, "message": "ok", "details": details}), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
