#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple
from urllib.request import urlopen

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from applications.crispr_guide_designer import CRISPRGuideDesigner
from applications.genomic_disruption_api import DisruptionAnalyzer
from validation.ontarget_gate.comparators.base import ComparatorRecord
from validation.ontarget_gate.comparators.registry import build_required_comparators, required_comparator_slots


PRIMARY_GROUPED_SOURCE = "doench2016_hg19"
EXTERNAL_HOLDOUT_SOURCE = "hart2016_hela_lib1_avg"
DECISION_SPLITS = {"dev_train", "dev_val", "primary_holdout", "external_holdout"}


@dataclass(frozen=True)
class SourceSpec:
    source: str
    url: str
    role_hint: str


SOURCES: List[SourceSpec] = [
    SourceSpec(
        source="doench2016_gecko1_lenticrispr_hg19",
        url="https://raw.githubusercontent.com/maximilianh/crisporPaper/master/effData/doench2016-Gecko1-LentiCrispr_hg19.scores.tab",
        role_hint="aux_weak",
    ),
    SourceSpec(
        source="doench2016_hg19",
        url="https://raw.githubusercontent.com/maximilianh/crisporPaper/master/effData/doench2016_hg19.scores.tab",
        role_hint="primary_grouped",
    ),
    SourceSpec(
        source="hart2016_hela_lib1_avg",
        url="https://raw.githubusercontent.com/maximilianh/crisporPaper/master/effData/hart2016-HelaLib1Avg.scores.tab",
        role_hint="external_holdout",
    ),
]


def clean_seq(seq: str) -> str:
    seq = (seq or "").upper()
    return "".join(ch for ch in seq if ch in "ACGTN")


def stable_hash_bucket(value: str, mod: int = 100) -> int:
    h = int(hashlib.sha1(value.encode("utf-8")).hexdigest(), 16)
    return h % mod


def assign_primary_split(group: str) -> str:
    bucket = stable_hash_bucket(group)
    if bucket < 60:
        return "dev_train"
    if bucket < 80:
        return "dev_val"
    return "primary_holdout"


def kmer2_features(seq: str) -> np.ndarray:
    kmers = [a + b for a in "ACGT" for b in "ACGT"]
    counts = {k: 0 for k in kmers}
    for i in range(len(seq) - 1):
        k = seq[i : i + 2]
        if k in counts:
            counts[k] += 1
    denom = max(1, len(seq) - 1)
    return np.array([counts[k] / denom for k in kmers], dtype=float)


def max_homopolymer_run(seq: str) -> int:
    best = 1
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
    return best


def baseline_a_features(seq: str) -> np.ndarray:
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    hom = max_homopolymer_run(seq) / len(seq)
    return np.concatenate([np.array([gc, hom], dtype=float), kmer2_features(seq)])


def parse_gene_like_key(value: str) -> str | None:
    if not value:
        return None
    if "-" in value:
        token = value.split("-", 1)[0].strip()
        return token or None
    if "_" in value:
        token = value.split("_", 1)[0].strip()
        return token or None
    return None


def infer_group(source: str, row: dict) -> Tuple[str | None, str, str]:
    guide = row.get("guide", "") or ""
    name = row.get("name", "") or ""
    context = (row.get("longSeq100Bp", "") or "").strip().upper()

    if source == "doench2016_hg19":
        group = parse_gene_like_key(guide)
        if group:
            return group, "strong", "guide_gene_prefix"
    if source == "hart2016_hela_lib1_avg":
        group = parse_gene_like_key(name) or parse_gene_like_key(guide)
        if group:
            return group, "strong", "name_gene_prefix"

    if context:
        ctx_hash = hashlib.sha1(context.encode("utf-8")).hexdigest()
        return f"ctxsha1:{ctx_hash}", "weak", "context_hash"

    fallback = parse_gene_like_key(guide)
    if fallback:
        return fallback, "weak", "guide_prefix_weak"

    return None, "weak", "missing_group_metadata"


def download_text(url: str) -> str:
    with urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


def parse_scores_tab(text: str, source: str, role_hint: str) -> List[dict]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    reader = csv.DictReader(lines, delimiter="\t")
    out: List[dict] = []
    for row in reader:
        seq = clean_seq(row.get("seq", ""))
        if len(seq) != 20:
            continue
        try:
            label = float(row.get("modFreq", "nan"))
        except ValueError:
            continue
        if np.isnan(label):
            continue

        group_key, group_strength, group_method = infer_group(source, row)
        out.append(
            {
                "source": source,
                "role_hint": role_hint,
                "guide": row.get("guide", "") or "",
                "guide_seq": seq,
                "label": label,
                "gene_or_target_group": group_key or "",
                "group_strength": group_strength,
                "group_method": group_method,
                "target_context": row.get("longSeq100Bp", "") or "",
            }
        )
    return out


def assign_split(row: dict) -> str:
    source = row["source"]
    strength = row["group_strength"]
    group = row["gene_or_target_group"]

    if source == PRIMARY_GROUPED_SOURCE:
        if strength != "strong" or not group:
            return "excluded_weak_group"
        return assign_primary_split(group)
    if source == EXTERNAL_HOLDOUT_SOURCE:
        if strength != "strong" or not group:
            return "excluded_weak_group"
        return "external_holdout"

    if strength == "strong":
        return "aux_train_strong"
    return "aux_train_weak"


def safe_spearman(y: np.ndarray, pred: np.ndarray) -> float:
    rho = float(spearmanr(y, pred).correlation)
    if np.isnan(rho):
        return 0.0
    return rho


def bootstrap_delta_spearman(
    y: np.ndarray,
    p_model: np.ndarray,
    p_base: np.ndarray,
    n_boot: int = 1000,
    seed: int = 42,
) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(y)
    idx_all = np.arange(n)
    deltas = np.zeros(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.choice(idx_all, size=n, replace=True)
        s_m = safe_spearman(y[idx], p_model[idx])
        s_b = safe_spearman(y[idx], p_base[idx])
        deltas[i] = s_m - s_b
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return float(lo), float(hi)


def metrics(y: np.ndarray, pred: np.ndarray) -> Dict[str, float]:
    return {
        "spearman": safe_spearman(y, pred),
        "mse": float(np.mean((y - pred) ** 2)),
    }


def validate_split_integrity(rows: Sequence[dict]) -> dict:
    violations: List[str] = []
    group_to_splits: Dict[str, set] = {}
    split_counts: Dict[str, int] = {}
    weak_rows_in_decision = 0

    for r in rows:
        split = r["split"]
        split_counts[split] = split_counts.get(split, 0) + 1
        if split not in DECISION_SPLITS:
            continue
        if r["group_strength"] != "strong":
            weak_rows_in_decision += 1
            continue
        group = r["gene_or_target_group"]
        if not group:
            violations.append(f"empty_group_in_decision_split:{split}")
            continue
        source_qualified_group = f"{r['source']}::{group}"
        group_to_splits.setdefault(source_qualified_group, set()).add(split)

    if weak_rows_in_decision:
        violations.append(f"weak_rows_in_decision_splits:{weak_rows_in_decision}")

    for group, splits in group_to_splits.items():
        if len(splits) > 1:
            violations.append(f"group_overlap:{group}:{','.join(sorted(splits))}")

    for required in ["dev_train", "dev_val", "primary_holdout", "external_holdout"]:
        if split_counts.get(required, 0) == 0:
            violations.append(f"missing_required_split:{required}")

    return {
        "ok": len(violations) == 0,
        "violations": violations,
        "split_counts": split_counts,
    }


def run_overlap_audit(split_manifest_path: Path, comparator_provenance: Dict[str, object]) -> dict:
    script = Path(__file__).resolve().parent / "audit_overlap.py"
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as tmp:
        out_path = Path(tmp.name)

    training_manifest = comparator_provenance.get("training_manifest_path", "")
    training_status = comparator_provenance.get("training_manifest_status_path", "")
    cmd = [
        sys.executable,
        str(script),
        "--split-manifest",
        str(split_manifest_path),
        "--training-manifest",
        str(training_manifest),
        "--training-status",
        str(training_status),
        "--output",
        str(out_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if not out_path.exists():
        return {
            "ok": False,
            "message": "overlap_audit_output_missing",
            "returncode": proc.returncode,
            "stdout": proc.stdout[-1000:],
            "stderr": proc.stderr[-1000:],
        }

    obj = json.loads(out_path.read_text(encoding="utf-8"))
    obj["returncode"] = proc.returncode
    obj["cmd"] = cmd
    return obj


def evaluate_split(name: str, rows: Sequence[dict], pred: Dict[str, np.ndarray]) -> dict:
    y = pred["y"]
    m_a = metrics(y, pred["baseline_a"])
    m_b = metrics(y, pred["baseline_b"])
    m_c = metrics(y, pred["baseline_c"])
    m_m = metrics(y, pred["model"])
    ci_a_lo, ci_a_hi = bootstrap_delta_spearman(y, pred["model"], pred["baseline_a"], n_boot=1000, seed=42)
    ci_c_lo, ci_c_hi = bootstrap_delta_spearman(y, pred["model"], pred["baseline_c"], n_boot=1000, seed=42)
    delta_a = m_m["spearman"] - m_a["spearman"]
    delta_b = m_m["spearman"] - m_b["spearman"]
    delta_c = m_m["spearman"] - m_c["spearman"]
    return {
        "split": name,
        "source": rows[0]["source"] if rows else "",
        "n": int(len(y)),
        "metrics": {
            "baseline_a": m_a,
            "baseline_b": m_b,
            "baseline_c": m_c,
            "model": m_m,
        },
        "delta_model_minus_baseline_a_spearman": float(delta_a),
        "delta_model_minus_baseline_a_spearman_ci95": [float(ci_a_lo), float(ci_a_hi)],
        "delta_model_minus_baseline_b_spearman": float(delta_b),
        "delta_model_minus_baseline_c_spearman": float(delta_c),
        "delta_model_minus_baseline_c_spearman_ci95": [float(ci_c_lo), float(ci_c_hi)],
    }


def decide_outcome_v3(res_dev: dict, res_primary: dict, res_external: dict, preconditions: dict) -> dict:
    if not preconditions.get("comparator_self_check_ok", False):
        return {
            "decision": "INCONCLUSIVE",
            "reason": "comparator_self_check_failed",
            "criteria": {},
        }
    if not preconditions.get("overlap_audit_ok", False):
        return {
            "decision": "INCONCLUSIVE",
            "reason": "overlap_audit_failed_or_unavailable",
            "criteria": {},
        }
    if not preconditions.get("holdout_min_n_ok", False):
        return {
            "decision": "INCONCLUSIVE",
            "reason": "holdout_size_below_minimum",
            "criteria": {},
        }

    d_ci_hi = res_dev["delta_model_minus_baseline_a_spearman_ci95"][1]
    hard_fail_dev = d_ci_hi < -0.01
    p_delta_c = res_primary["delta_model_minus_baseline_c_spearman"]
    e_delta_c = res_external["delta_model_minus_baseline_c_spearman"]

    criteria = {
        "hard_precondition_dev_not_materially_negative": not hard_fail_dev,
        "primary_delta_vs_baseline_c_min_0p01": p_delta_c >= 0.01,
        "external_delta_vs_baseline_c_min_0p01": e_delta_c >= 0.01,
    }

    if hard_fail_dev:
        return {"decision": "NO-GO", "reason": "dev_material_negative", "criteria": criteria}

    if all(criteria.values()):
        return {"decision": "GO", "reason": "all_v3_criteria_passed", "criteria": criteria}

    return {"decision": "NO-GO", "reason": "baseline_c_threshold_failed", "criteria": criteria}


def write_csv(path: Path, rows: Iterable[dict], fieldnames: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    data_dir = root / "data"
    outputs_dir = root / "outputs"
    data_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    generated_at = datetime.now(timezone.utc).isoformat()

    all_rows: List[dict] = []
    source_counts: List[dict] = []
    for spec in SOURCES:
        raw_text = download_text(spec.url)
        raw_path = data_dir / f"{spec.source}.scores.tab"
        raw_path.write_text(raw_text, encoding="utf-8")

        parsed = parse_scores_tab(raw_text, spec.source, spec.role_hint)
        all_rows.extend(parsed)
        strong_rows = sum(1 for r in parsed if r["group_strength"] == "strong")
        source_counts.append(
            {
                "source": spec.source,
                "role_hint": spec.role_hint,
                "url": spec.url,
                "raw_local_path": str(raw_path),
                "row_count": len(parsed),
                "strong_group_rows": strong_rows,
                "weak_group_rows": len(parsed) - strong_rows,
            }
        )

    split_rows: List[dict] = []
    for r in all_rows:
        row = dict(r)
        row["split"] = assign_split(row)
        split_rows.append(row)

    split_manifest_path = outputs_dir / "locked_split_manifest_v3.csv"
    write_csv(
        split_manifest_path,
        split_rows,
        [
            "source",
            "guide",
            "guide_seq",
            "label",
            "gene_or_target_group",
            "group_strength",
            "group_method",
            "split",
        ],
    )

    integrity = validate_split_integrity(split_rows)

    dataset_manifest = {
        "generated_at": generated_at,
        "protocol_version": "v3",
        "schema": [
            "guide",
            "guide_seq",
            "label",
            "source",
            "gene_or_target_group",
            "group_strength",
            "group_method",
            "target_context",
            "split",
        ],
        "sources": source_counts,
        "reserved_comparator_slots": required_comparator_slots(),
    }
    (outputs_dir / "locked_dataset_manifest_v3.json").write_text(json.dumps(dataset_manifest, indent=2), encoding="utf-8")

    schema_rows = [
        {"field": "guide", "description": "Guide identifier from source"},
        {"field": "guide_seq", "description": "20nt guide sequence"},
        {"field": "label", "description": "Measured on-target activity (modFreq)"},
        {"field": "source", "description": "Dataset source identifier"},
        {"field": "gene_or_target_group", "description": "Grouping key for split integrity"},
        {"field": "group_strength", "description": "strong or weak grouping metadata confidence"},
        {"field": "group_method", "description": "Method used to derive grouping key"},
        {"field": "target_context", "description": "Optional long context sequence"},
        {"field": "split", "description": "Assigned split label"},
    ]
    write_csv(outputs_dir / "locked_schema_manifest_v3.csv", schema_rows, ["field", "description"])

    train_rows = [r for r in split_rows if r["split"] in {"dev_train", "aux_train_weak", "aux_train_strong"}]
    dev_rows = [r for r in split_rows if r["split"] == "dev_val"]
    primary_rows = [r for r in split_rows if r["split"] == "primary_holdout"]
    external_rows = [r for r in split_rows if r["split"] == "external_holdout"]

    holdout_min_n = 200
    holdout_min_n_ok = len(primary_rows) >= holdout_min_n and len(external_rows) >= holdout_min_n

    X_train = np.vstack([baseline_a_features(r["guide_seq"]) for r in train_rows])
    y_train = np.array([float(r["label"]) for r in train_rows], dtype=float)
    model_a = Ridge(alpha=1.0)
    model_a.fit(X_train, y_train)

    da = DisruptionAnalyzer(k=0.3, seed=42)
    b_model = CRISPRGuideDesigner()

    comparators = build_required_comparators()
    baseline_c = comparators["baseline_c"]
    c_self = baseline_c.self_check()
    c_prov = baseline_c.provenance()
    overlap = run_overlap_audit(split_manifest_path, c_prov)

    preconditions = {
        "split_integrity_ok": integrity["ok"],
        "comparator_self_check_ok": c_self.ok,
        "overlap_audit_ok": bool(overlap.get("ok", False)),
        "holdout_min_n_ok": holdout_min_n_ok,
        "holdout_min_n_required": holdout_min_n,
    }

    def predict_pack(rows: Sequence[dict]) -> Dict[str, np.ndarray]:
        y = np.array([float(r["label"]) for r in rows], dtype=float)
        X = np.vstack([baseline_a_features(r["guide_seq"]) for r in rows])
        p_a = model_a.predict(X)
        p_b = np.array([float(b_model.calculate_on_target_score(r["guide_seq"])) for r in rows], dtype=float)
        p_m = np.array([float(da.score_guide(r["guide_seq"])["disruption_score"]) for r in rows], dtype=float)

        if preconditions["comparator_self_check_ok"]:
            c_records = [ComparatorRecord(guide_seq=r["guide_seq"], target_context=r.get("target_context", "")) for r in rows]
            p_c = np.array(baseline_c.predict_batch(c_records), dtype=float)
        else:
            # Keep diagnostics numeric when comparator preconditions fail.
            p_c = np.zeros(len(rows), dtype=float)

        return {"y": y, "baseline_a": p_a, "baseline_b": p_b, "baseline_c": p_c, "model": p_m}

    pred_dev = predict_pack(dev_rows)
    pred_primary = predict_pack(primary_rows)
    pred_external = predict_pack(external_rows)

    res_dev = evaluate_split("dev_val", dev_rows, pred_dev)
    res_primary = evaluate_split("primary_holdout", primary_rows, pred_primary)
    res_external = evaluate_split("external_holdout", external_rows, pred_external)

    decision_pack = decide_outcome_v3(res_dev, res_primary, res_external, preconditions)

    report_json = {
        "generated_at": generated_at,
        "protocol_version": "v3",
        "decision": decision_pack["decision"],
        "decision_reason": decision_pack["reason"],
        "decision_criteria": decision_pack["criteria"],
        "preconditions": preconditions,
        "split_integrity": integrity,
        "overlap_audit": overlap,
        "baseline_c": {
            "slot": "baseline_c",
            "name": baseline_c.name,
            "version": baseline_c.version,
            "self_check": {
                "ok": c_self.ok,
                "message": c_self.message,
                "details": c_self.details,
            },
            "provenance": c_prov,
        },
        "dataset_counts": source_counts,
        "split_counts": integrity["split_counts"],
        "results": [res_dev, res_primary, res_external],
        "ci_policy": {
            "mode": "diagnostic_only",
            "note": "CI values are computed and reported, but v3 gating uses fixed delta thresholds for baseline_c only.",
        },
        "deprecated_prior_artifacts": {
            "gate_results.json": "historical_non_authoritative",
            "gate_report.md": "historical_non_authoritative",
            "gate_results_v2.json": "historical_non_authoritative_for_external_comparator_claims",
            "gate_report_v2.md": "historical_non_authoritative_for_external_comparator_claims",
        },
    }
    (outputs_dir / "gate_results_v3.json").write_text(json.dumps(report_json, indent=2), encoding="utf-8")

    md: List[str] = []
    md.append("# Fast Go/No-Go Validation Report v3 (On-Target)")
    md.append("")
    md.append(f"Generated: {generated_at}")
    md.append("")
    md.append(f"## Decision: **{decision_pack['decision']}**")
    md.append(f"- reason: {decision_pack['reason']}")
    md.append("")
    md.append("## Preconditions")
    for k, v in preconditions.items():
        md.append(f"- {k}: {v}")
    md.append("")
    md.append("## Baseline C Comparator")
    md.append(f"- slot: baseline_c")
    md.append(f"- name: {baseline_c.name}")
    md.append(f"- version: {baseline_c.version}")
    md.append(f"- self_check_ok: {c_self.ok}")
    md.append(f"- self_check_message: {c_self.message}")
    md.append(f"- provenance_path: {c_prov.get('provenance_path', '')}")
    md.append(f"- checksums_path: {c_prov.get('checksums_path', '')}")
    md.append("")
    md.append("## Overlap Audit")
    md.append(f"- status: {'PASS' if overlap.get('ok', False) else 'FAIL'}")
    md.append(f"- message: {overlap.get('message', '')}")
    md.append(f"- overlap_count: {overlap.get('overlap_count', 'n/a')}")
    md.append("")
    md.append("## Per-Split Metrics (Spearman / MSE)")
    for res in [res_dev, res_primary, res_external]:
        md.append(f"### {res['split']} ({res['source']}, n={res['n']})")
        for key in ["baseline_a", "baseline_b", "baseline_c", "model"]:
            m = res["metrics"][key]
            md.append(f"- {key}: spearman={m['spearman']:.4f}, mse={m['mse']:.4f}")
        ci_a = res["delta_model_minus_baseline_a_spearman_ci95"]
        ci_c = res["delta_model_minus_baseline_c_spearman_ci95"]
        md.append(f"- delta(model - baseline_a) spearman: {res['delta_model_minus_baseline_a_spearman']:.4f} (95% CI [{ci_a[0]:.4f}, {ci_a[1]:.4f}])")
        md.append(f"- delta(model - baseline_b) spearman: {res['delta_model_minus_baseline_b_spearman']:.4f}")
        md.append(f"- delta(model - baseline_c) spearman: {res['delta_model_minus_baseline_c_spearman']:.4f} (95% CI [{ci_c[0]:.4f}, {ci_c[1]:.4f}])")
        md.append("")

    md.append("## Rule Evaluation (CI Diagnostic-Only in v3)")
    if decision_pack["criteria"]:
        for name, passed in decision_pack["criteria"].items():
            md.append(f"- {name}: {'PASS' if passed else 'FAIL'}")
    else:
        md.append("- criteria_not_evaluated_due_to_precondition_failure: true")

    (outputs_dir / "gate_report_v3.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print(f"decision: {decision_pack['decision']}")
    print(f"reason: {decision_pack['reason']}")
    print(f"wrote {outputs_dir / 'locked_dataset_manifest_v3.json'}")
    print(f"wrote {outputs_dir / 'locked_split_manifest_v3.csv'}")
    print(f"wrote {outputs_dir / 'locked_schema_manifest_v3.csv'}")
    print(f"wrote {outputs_dir / 'gate_results_v3.json'}")
    print(f"wrote {outputs_dir / 'gate_report_v3.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
