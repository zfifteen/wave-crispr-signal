#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def _load_holdout_sequences(split_manifest: Path) -> tuple[set[str], dict[str, int]]:
    seqs: set[str] = set()
    counts = {"primary_holdout": 0, "external_holdout": 0}
    with split_manifest.open(encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            split = row.get("split", "")
            if split not in counts:
                continue
            seq = (row.get("guide_seq", "") or "").strip().upper()
            if not seq:
                continue
            counts[split] += 1
            seqs.add(seq)
    return seqs, counts


def _load_training_sequences(training_manifest: Path) -> set[str]:
    seqs: set[str] = set()
    seq_parts = []
    with training_manifest.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_parts:
                    seqs.add("".join(seq_parts).upper())
                    seq_parts = []
            else:
                seq_parts.append(line)
        if seq_parts:
            seqs.add("".join(seq_parts).upper())
    return seqs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split-manifest", required=True)
    ap.add_argument("--training-manifest", required=True)
    ap.add_argument("--training-status", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    split_manifest = Path(args.split_manifest)
    training_manifest = Path(args.training_manifest)
    training_status = Path(args.training_status)

    out = {
        "ok": False,
        "message": "unknown",
        "overlap_count": None,
        "overlap_examples": [],
        "holdout_counts": {},
        "training_sequence_count": None,
    }

    if not split_manifest.exists():
        out["message"] = "split_manifest_missing"
        Path(args.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
        return 2

    if not training_status.exists():
        out["message"] = "training_status_missing"
        Path(args.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
        return 3

    status = json.loads(training_status.read_text(encoding="utf-8"))
    out["training_status"] = status
    if not bool(status.get("available", False)):
        out["message"] = "training_manifest_unavailable"
        Path(args.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
        return 4

    if not training_manifest.exists():
        out["message"] = "training_manifest_missing"
        Path(args.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
        return 5

    holdout_seqs, holdout_counts = _load_holdout_sequences(split_manifest)
    train_seqs = _load_training_sequences(training_manifest)

    overlap = sorted(holdout_seqs.intersection(train_seqs))
    out.update(
        {
            "holdout_counts": holdout_counts,
            "training_sequence_count": len(train_seqs),
            "overlap_count": len(overlap),
            "overlap_examples": overlap[:20],
        }
    )

    if overlap:
        out["ok"] = False
        out["message"] = "overlap_detected"
        Path(args.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
        return 6

    out["ok"] = True
    out["message"] = "no_overlap"
    Path(args.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
