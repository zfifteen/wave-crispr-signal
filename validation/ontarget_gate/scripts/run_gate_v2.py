#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import sys
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


def bootstrap_weighted_delta_spearman(
    payloads: Sequence[dict],
    n_boot: int = 1000,
    seed: int = 42,
) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    deltas = np.zeros(n_boot, dtype=float)
    for i in range(n_boot):
        weighted_sum = 0.0
        total_weight = 0
        for p in payloads:
            y = p["y"]
            pm = p["model"]
            pb = p["baseline_a"]
            n = len(y)
            idx = rng.choice(np.arange(n), size=n, replace=True)
            d = safe_spearman(y[idx], pm[idx]) - safe_spearman(y[idx], pb[idx])
            weighted_sum += d * n
            total_weight += n
        deltas[i] = weighted_sum / max(total_weight, 1)
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


def evaluate_split(
    name: str,
    rows: Sequence[dict],
    pred: Dict[str, np.ndarray],
) -> dict:
    y = pred["y"]
    m_a = metrics(y, pred["baseline_a"])
    m_b = metrics(y, pred["baseline_b"])
    m_m = metrics(y, pred["model"])
    ci_lo, ci_hi = bootstrap_delta_spearman(y, pred["model"], pred["baseline_a"], n_boot=1000, seed=42)
    delta_a = m_m["spearman"] - m_a["spearman"]
    delta_b = m_m["spearman"] - m_b["spearman"]
    return {
        "split": name,
        "source": rows[0]["source"] if rows else "",
        "n": int(len(y)),
        "metrics": {"baseline_a": m_a, "baseline_b": m_b, "model": m_m},
        "delta_model_minus_baseline_a_spearman": float(delta_a),
        "delta_model_minus_baseline_a_spearman_ci95": [float(ci_lo), float(ci_hi)],
        "delta_model_minus_baseline_b_spearman": float(delta_b),
    }


def decide_outcome(
    res_dev: dict,
    res_primary: dict,
    res_external: dict,
    agg_delta: float,
    agg_ci: Tuple[float, float],
) -> dict:
    d_delta = res_dev["delta_model_minus_baseline_a_spearman"]
    d_ci_lo, d_ci_hi = res_dev["delta_model_minus_baseline_a_spearman_ci95"]
    p_delta = res_primary["delta_model_minus_baseline_a_spearman"]
    p_ci_lo, p_ci_hi = res_primary["delta_model_minus_baseline_a_spearman_ci95"]
    e_delta = res_external["delta_model_minus_baseline_a_spearman"]
    e_ci_lo, e_ci_hi = res_external["delta_model_minus_baseline_a_spearman_ci95"]
    p_vs_b = res_primary["delta_model_minus_baseline_b_spearman"]
    e_vs_b = res_external["delta_model_minus_baseline_b_spearman"]
    agg_ci_lo, agg_ci_hi = agg_ci

    hard_fail_dev = d_ci_hi < -0.01

    criteria = {
        "hard_precondition_dev_not_materially_negative": not hard_fail_dev,
        "primary_delta_min": p_delta >= 0.03,
        "primary_ci_excludes_zero": p_ci_lo > 0.0,
        "external_delta_min": e_delta >= 0.01,
        "external_ci_excludes_zero": e_ci_lo > 0.0,
        "directional_consistency": (p_delta > 0.0 and e_delta > 0.0),
        "no_underperform_vs_baseline_b_primary": p_vs_b >= 0.0,
        "no_underperform_vs_baseline_b_external": e_vs_b >= 0.0,
        "aggregate_delta_positive": agg_delta > 0.0,
        "aggregate_ci_excludes_zero": agg_ci_lo > 0.0,
    }

    go = all(criteria.values())
    no_go = (
        hard_fail_dev
        or p_delta <= 0.0
        or e_delta <= 0.0
        or p_ci_lo < -0.01
        or e_ci_lo < -0.01
        or p_vs_b < 0.0
        or e_vs_b < 0.0
        or agg_ci_lo < -0.01
    )

    if go:
        decision = "GO"
    elif no_go:
        decision = "NO-GO"
    else:
        decision = "INCONCLUSIVE"

    return {
        "decision": decision,
        "criteria": criteria,
        "diagnostics": {
            "dev_delta_a": d_delta,
            "dev_ci_a": [d_ci_lo, d_ci_hi],
            "primary_delta_a": p_delta,
            "primary_ci_a": [p_ci_lo, p_ci_hi],
            "external_delta_a": e_delta,
            "external_ci_a": [e_ci_lo, e_ci_hi],
            "primary_delta_b": p_vs_b,
            "external_delta_b": e_vs_b,
            "aggregate_delta_a": agg_delta,
            "aggregate_ci_a": [agg_ci_lo, agg_ci_hi],
            "aggregate_ci_hi": agg_ci_hi,
        },
    }


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

    integrity = validate_split_integrity(split_rows)
    if not integrity["ok"]:
        out = {
            "generated_at": generated_at,
            "decision": "NO-GO",
            "fatal_error": "split_integrity_failed",
            "split_integrity": integrity,
        }
        (outputs_dir / "gate_results_v2.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
        (outputs_dir / "gate_report_v2.md").write_text(
            "# Fast Go/No-Go Validation Report v2 (On-Target)\n\n"
            "## Decision: **NO-GO**\n\n"
            "Split integrity checks failed. See `gate_results_v2.json`.\n",
            encoding="utf-8",
        )
        print("decision: NO-GO")
        print("fatal: split integrity failed")
        return 1

    dataset_manifest = {
        "generated_at": generated_at,
        "protocol_version": "v2",
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
    }
    (outputs_dir / "locked_dataset_manifest_v2.json").write_text(json.dumps(dataset_manifest, indent=2), encoding="utf-8")

    write_csv(
        outputs_dir / "locked_split_manifest_v2.csv",
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
    write_csv(outputs_dir / "locked_schema_manifest_v2.csv", schema_rows, ["field", "description"])

    train_rows = [r for r in split_rows if r["split"] in {"dev_train", "aux_train_weak", "aux_train_strong"}]
    dev_rows = [r for r in split_rows if r["split"] == "dev_val"]
    primary_rows = [r for r in split_rows if r["split"] == "primary_holdout"]
    external_rows = [r for r in split_rows if r["split"] == "external_holdout"]

    X_train = np.vstack([baseline_a_features(r["guide_seq"]) for r in train_rows])
    y_train = np.array([float(r["label"]) for r in train_rows], dtype=float)
    model_a = Ridge(alpha=1.0)
    model_a.fit(X_train, y_train)

    da = DisruptionAnalyzer(k=0.3, seed=42)
    b_model = CRISPRGuideDesigner()

    def predict_pack(rows: Sequence[dict]) -> Dict[str, np.ndarray]:
        y = np.array([float(r["label"]) for r in rows], dtype=float)
        X = np.vstack([baseline_a_features(r["guide_seq"]) for r in rows])
        p_a = model_a.predict(X)
        p_b = np.array([float(b_model.calculate_on_target_score(r["guide_seq"])) for r in rows], dtype=float)
        p_m = np.array([float(da.score_guide(r["guide_seq"])["disruption_score"]) for r in rows], dtype=float)
        return {"y": y, "baseline_a": p_a, "baseline_b": p_b, "model": p_m}

    pred_dev = predict_pack(dev_rows)
    pred_primary = predict_pack(primary_rows)
    pred_external = predict_pack(external_rows)

    res_dev = evaluate_split("dev_val", dev_rows, pred_dev)
    res_primary = evaluate_split("primary_holdout", primary_rows, pred_primary)
    res_external = evaluate_split("external_holdout", external_rows, pred_external)

    holdout_payloads = [pred_primary, pred_external]
    holdout_deltas = [
        (res_primary["delta_model_minus_baseline_a_spearman"], res_primary["n"]),
        (res_external["delta_model_minus_baseline_a_spearman"], res_external["n"]),
    ]
    total_n = sum(n for _, n in holdout_deltas)
    agg_delta = float(sum(delta * n for delta, n in holdout_deltas) / max(total_n, 1))
    agg_ci_lo, agg_ci_hi = bootstrap_weighted_delta_spearman(holdout_payloads, n_boot=1000, seed=42)

    decision_pack = decide_outcome(res_dev, res_primary, res_external, agg_delta, (agg_ci_lo, agg_ci_hi))

    report_json = {
        "generated_at": generated_at,
        "protocol_version": "v2",
        "decision": decision_pack["decision"],
        "decision_criteria": decision_pack["criteria"],
        "decision_diagnostics": decision_pack["diagnostics"],
        "split_integrity": integrity,
        "dataset_counts": source_counts,
        "split_counts": integrity["split_counts"],
        "results": [res_dev, res_primary, res_external],
        "aggregated_holdout_delta_model_minus_baseline_a_spearman": {
            "delta": agg_delta,
            "ci95": [agg_ci_lo, agg_ci_hi],
            "weights": {
                "primary_holdout_n": res_primary["n"],
                "external_holdout_n": res_external["n"],
            },
        },
        "deprecated_prior_artifacts": {
            "gate_results.json": "non-authoritative",
            "gate_report.md": "non-authoritative",
        },
    }
    (outputs_dir / "gate_results_v2.json").write_text(json.dumps(report_json, indent=2), encoding="utf-8")

    md: List[str] = []
    md.append("# Fast Go/No-Go Validation Report v2 (On-Target)")
    md.append("")
    md.append(f"Generated: {generated_at}")
    md.append("")
    md.append(f"## Decision: **{decision_pack['decision']}**")
    md.append("")
    md.append("## Split Integrity")
    md.append(f"- status: {'PASS' if integrity['ok'] else 'FAIL'}")
    for key in ["dev_train", "dev_val", "primary_holdout", "external_holdout", "aux_train_weak", "aux_train_strong", "excluded_weak_group"]:
        if key in integrity["split_counts"]:
            md.append(f"- {key}: {integrity['split_counts'][key]}")
    if integrity["violations"]:
        md.append("- violations:")
        for v in integrity["violations"]:
            md.append(f"  - {v}")
    md.append("")
    md.append("## Per-Split Metrics (Spearman / MSE)")
    for res in [res_dev, res_primary, res_external]:
        md.append(f"### {res['split']} ({res['source']}, n={res['n']})")
        for key in ["baseline_a", "baseline_b", "model"]:
            m = res["metrics"][key]
            md.append(f"- {key}: spearman={m['spearman']:.4f}, mse={m['mse']:.4f}")
        ci = res["delta_model_minus_baseline_a_spearman_ci95"]
        md.append(f"- delta(model - baseline_a) spearman: {res['delta_model_minus_baseline_a_spearman']:.4f} (95% CI [{ci[0]:.4f}, {ci[1]:.4f}])")
        md.append(f"- delta(model - baseline_b) spearman: {res['delta_model_minus_baseline_b_spearman']:.4f}")
        md.append("")
    md.append("## Aggregated Holdout Delta (No Mixed-Scale Pooled Spearman)")
    md.append(f"- weighted delta(model - baseline_a): {agg_delta:.4f} (95% CI [{agg_ci_lo:.4f}, {agg_ci_hi:.4f}])")
    md.append("")
    md.append("## Rule Evaluation")
    for name, passed in decision_pack["criteria"].items():
        md.append(f"- {name}: {'PASS' if passed else 'FAIL'}")
    md.append("")
    md.append("## One-Page Decision Summary")
    md.append(f"- outcome: {decision_pack['decision']}")
    md.append("- rationale: derived exclusively from locked v2 criteria in `PROTOCOL_V2.md`.")

    (outputs_dir / "gate_report_v2.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print(f"decision: {decision_pack['decision']}")
    print(f"wrote {outputs_dir / 'locked_dataset_manifest_v2.json'}")
    print(f"wrote {outputs_dir / 'locked_split_manifest_v2.csv'}")
    print(f"wrote {outputs_dir / 'locked_schema_manifest_v2.csv'}")
    print(f"wrote {outputs_dir / 'gate_results_v2.json'}")
    print(f"wrote {outputs_dir / 'gate_report_v2.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
