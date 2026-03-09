#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def _clean_seq(seq: str) -> str:
    s = (seq or "").upper()
    return "".join(ch for ch in s if ch in "ACGT")


def _score(seq: str, context: str) -> float:
    # Lightweight deterministic scaffold scorer used only after comparator self-check passes.
    # Rank-based gating metric is Spearman; score scale is not used for thresholding.
    seq = _clean_seq(seq)
    if len(seq) != 20:
        raise ValueError("guide_seq must be 20nt after normalization")

    gc = (seq.count("G") + seq.count("C")) / len(seq)
    hom = 0
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
        else:
            hom = max(hom, cur)
            cur = 1
    hom = max(hom, cur)

    ctx = _clean_seq(context)
    ctx_bonus = 0.0
    if ctx:
        ctx_gc = (ctx.count("G") + ctx.count("C")) / max(1, len(ctx))
        ctx_bonus = 0.05 * (0.5 - abs(ctx_gc - 0.5))

    # Deterministic bounded output in [0,1]
    z = (2.0 * (0.5 - abs(gc - 0.5))) - (0.25 * (hom / 20.0)) + ctx_bonus
    return 1.0 / (1.0 + math.exp(-4.0 * z))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    inp = json.loads(Path(args.input).read_text(encoding="utf-8"))
    records = inp.get("records", [])
    if not isinstance(records, list):
        out = {"ok": False, "message": "records_not_list"}
        Path(args.output).write_text(json.dumps(out), encoding="utf-8")
        return 2

    scores = []
    try:
        for rec in records:
            seq = rec.get("guide_seq", "")
            ctx = rec.get("target_context", "")
            scores.append(_score(seq, ctx))
    except Exception as exc:  # noqa: BLE001
        Path(args.output).write_text(
            json.dumps({"ok": False, "message": "scoring_failed", "error": str(exc)}),
            encoding="utf-8",
        )
        return 3

    Path(args.output).write_text(json.dumps({"ok": True, "message": "ok", "scores": scores}), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
