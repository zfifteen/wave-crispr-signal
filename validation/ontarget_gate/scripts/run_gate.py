#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple
from urllib.request import urlopen

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from applications.crispr_guide_designer import CRISPRGuideDesigner
from applications.genomic_disruption_api import DisruptionAnalyzer


@dataclass(frozen=True)
class SourceSpec:
    source: str
    url: str
    role: str  # dev_pool | primary_holdout | external_holdout


SOURCES: List[SourceSpec] = [
    SourceSpec(
        source="doench2016_gecko1_lenticrispr_hg19",
        url="https://raw.githubusercontent.com/maximilianh/crisporPaper/master/effData/doench2016-Gecko1-LentiCrispr_hg19.scores.tab",
        role="dev_pool",
    ),
    SourceSpec(
        source="doench2016_hg19",
        url="https://raw.githubusercontent.com/maximilianh/crisporPaper/master/effData/doench2016_hg19.scores.tab",
        role="primary_holdout",
    ),
    SourceSpec(
        source="hart2016_hela_lib1_avg",
        url="https://raw.githubusercontent.com/maximilianh/crisporPaper/master/effData/hart2016-HelaLib1Avg.scores.tab",
        role="external_holdout",
    ),
]


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
    feats = np.concatenate([
        np.array([gc, hom], dtype=float),
        kmer2_features(seq),
    ])
    return feats


def clean_seq(seq: str) -> str:
    seq = (seq or "").upper()
    seq = "".join(ch for ch in seq if ch in "ACGTN")
    return seq


def group_from_guide(guide: str) -> str:
    g = guide or ""
    if "-" in g:
        return g.split("-", 1)[0]
    if "_" in g:
        return g.split("_", 1)[0]
    return g or "unknown"


def download_text(url: str) -> str:
    with urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


def parse_scores_tab(text: str, source: str, role: str) -> List[dict]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    reader = csv.DictReader(lines, delimiter="\t")
    out = []
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

        guide = row.get("guide", "")
        out.append(
            {
                "source": source,
                "role": role,
                "guide": guide,
                "guide_seq": seq,
                "label": label,
                "gene_or_target_group": group_from_guide(guide),
                "target_context": row.get("longSeq100Bp", ""),
            }
        )
    return out


def deterministic_split(group: str) -> str:
    h = int(hashlib.sha1(group.encode("utf-8")).hexdigest(), 16)
    return "dev_train" if (h % 100) < 80 else "dev_val"


def bootstrap_delta_spearman(y: np.ndarray, p_model: np.ndarray, p_base: np.ndarray, n_boot: int = 2000, seed: int = 42) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(y)
    deltas = np.zeros(n_boot, dtype=float)
    idx_all = np.arange(n)
    for i in range(n_boot):
        idx = rng.choice(idx_all, size=n, replace=True)
        s_m = spearmanr(y[idx], p_model[idx]).correlation
        s_b = spearmanr(y[idx], p_base[idx]).correlation
        deltas[i] = float(s_m - s_b)
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return float(lo), float(hi)


def metrics(y: np.ndarray, pred: np.ndarray) -> Dict[str, float]:
    rho = float(spearmanr(y, pred).correlation)
    mse = float(np.mean((y - pred) ** 2))
    return {"spearman": rho, "mse": mse}


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    data_dir = root / "data"
    outputs_dir = root / "outputs"
    data_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    generated_at = datetime.now(timezone.utc).isoformat()

    # 1) Acquire + normalize
    all_rows: List[dict] = []
    source_counts = []
    for spec in SOURCES:
        raw_text = download_text(spec.url)
        raw_path = data_dir / f"{spec.source}.scores.tab"
        raw_path.write_text(raw_text, encoding="utf-8")

        parsed = parse_scores_tab(raw_text, spec.source, spec.role)
        all_rows.extend(parsed)
        source_counts.append({
            "source": spec.source,
            "role": spec.role,
            "url": spec.url,
            "raw_local_path": str(raw_path),
            "row_count": len(parsed),
        })

    # 2) lock schema dataset manifest
    dataset_manifest_path = outputs_dir / "locked_dataset_manifest.json"
    dataset_manifest = {
        "generated_at": generated_at,
        "schema": ["guide", "label", "source", "gene_or_target_group", "target_context"],
        "sources": source_counts,
    }
    dataset_manifest_path.write_text(json.dumps(dataset_manifest, indent=2), encoding="utf-8")

    # 3) split manifest lock
    split_rows = []
    for r in all_rows:
        split = r["role"]
        if r["role"] == "dev_pool":
            split = deterministic_split(r["gene_or_target_group"])
        split_rows.append({
            "source": r["source"],
            "guide": r["guide"],
            "guide_seq": r["guide_seq"],
            "label": r["label"],
            "gene_or_target_group": r["gene_or_target_group"],
            "split": split,
        })

    split_manifest_path = outputs_dir / "locked_split_manifest.csv"
    write_csv(split_manifest_path, split_rows, ["source", "guide", "guide_seq", "label", "gene_or_target_group", "split"])

    # 4) prepare arrays
    dev_train = [r for r in split_rows if r["split"] == "dev_train"]
    dev_val = [r for r in split_rows if r["split"] == "dev_val"]
    primary = [r for r in split_rows if r["split"] == "primary_holdout"]
    external = [r for r in split_rows if r["split"] == "external_holdout"]

    # Baseline A fit
    X_train = np.vstack([baseline_a_features(r["guide_seq"]) for r in dev_train])
    y_train = np.array([float(r["label"]) for r in dev_train], dtype=float)
    model_a = Ridge(alpha=1.0)
    model_a.fit(X_train, y_train)

    # Model + Baseline B setup
    da = DisruptionAnalyzer(k=0.3, seed=42)
    b_model = CRISPRGuideDesigner()

    def predict_pack(rows: List[dict]) -> Dict[str, np.ndarray]:
        y = np.array([float(r["label"]) for r in rows], dtype=float)

        X = np.vstack([baseline_a_features(r["guide_seq"]) for r in rows])
        p_a = model_a.predict(X)

        p_b = []
        p_m = []
        for r in rows:
            seq = r["guide_seq"]
            p_b.append(float(b_model.calculate_on_target_score(seq)))
            p_m.append(float(da.score_guide(seq)["disruption_score"]))

        return {
            "y": y,
            "baseline_a": np.array(p_a, dtype=float),
            "baseline_b": np.array(p_b, dtype=float),
            "model": np.array(p_m, dtype=float),
        }

    pred_dev = predict_pack(dev_val)
    pred_primary = predict_pack(primary)
    pred_external = predict_pack(external)

    # 5) evaluate
    def eval_split(name: str, pred: Dict[str, np.ndarray]) -> dict:
        y = pred["y"]
        m_a = metrics(y, pred["baseline_a"])
        m_b = metrics(y, pred["baseline_b"])
        m_m = metrics(y, pred["model"])

        ci_lo, ci_hi = bootstrap_delta_spearman(y, pred["model"], pred["baseline_a"], n_boot=2000, seed=42)
        delta = m_m["spearman"] - m_a["spearman"]
        delta_vs_b = m_m["spearman"] - m_b["spearman"]

        return {
            "split": name,
            "n": int(len(y)),
            "metrics": {
                "baseline_a": m_a,
                "baseline_b": m_b,
                "model": m_m,
            },
            "delta_model_minus_baseline_a_spearman": float(delta),
            "delta_model_minus_baseline_a_spearman_ci95": [ci_lo, ci_hi],
            "delta_model_minus_baseline_b_spearman": float(delta_vs_b),
        }

    res_dev = eval_split("dev_val", pred_dev)
    res_primary = eval_split("primary_holdout", pred_primary)
    res_external = eval_split("external_holdout", pred_external)

    # Decision criteria
    p_delta = res_primary["delta_model_minus_baseline_a_spearman"]
    p_ci_lo, p_ci_hi = res_primary["delta_model_minus_baseline_a_spearman_ci95"]
    e_delta = res_external["delta_model_minus_baseline_a_spearman"]
    p_vs_b = res_primary["delta_model_minus_baseline_b_spearman"]
    e_vs_b = res_external["delta_model_minus_baseline_b_spearman"]

    go = (
        (p_delta >= 0.03)
        and (p_ci_lo > 0.0)
        and (e_delta >= 0.01)
        and (p_vs_b >= 0.0)
        and (e_vs_b >= 0.0)
    )

    no_go = (
        (p_delta <= 0.0)
        or (e_delta < 0.0)
        or (p_ci_lo < -0.01)
    )

    if go:
        decision = "GO"
    elif no_go:
        decision = "NO-GO"
    else:
        decision = "INCONCLUSIVE"

    report_json = {
        "generated_at": generated_at,
        "decision": decision,
        "thresholds": {
            "primary_delta_min": 0.03,
            "primary_ci_excludes_zero": True,
            "external_delta_min": 0.01,
            "no_underperform_vs_baseline_b": True,
        },
        "dataset_counts": source_counts,
        "split_counts": {
            "dev_train": len(dev_train),
            "dev_val": len(dev_val),
            "primary_holdout": len(primary),
            "external_holdout": len(external),
        },
        "results": [res_dev, res_primary, res_external],
    }

    (outputs_dir / "gate_results.json").write_text(json.dumps(report_json, indent=2), encoding="utf-8")

    md = []
    md.append("# Fast Go/No-Go Validation Report (On-Target Efficiency)")
    md.append("")
    md.append(f"Generated: {generated_at}")
    md.append("")
    md.append(f"## Decision: **{decision}**")
    md.append("")
    md.append("## Dataset Counts")
    for s in source_counts:
        md.append(f"- `{s['source']}` ({s['role']}): {s['row_count']} rows")
    md.append("")
    md.append("## Split Counts")
    md.append(f"- dev_train: {len(dev_train)}")
    md.append(f"- dev_val: {len(dev_val)}")
    md.append(f"- primary_holdout: {len(primary)}")
    md.append(f"- external_holdout: {len(external)}")
    md.append("")
    md.append("## Metrics (Spearman / MSE)")
    for res in [res_dev, res_primary, res_external]:
        md.append(f"### {res['split']} (n={res['n']})")
        for key in ["baseline_a", "baseline_b", "model"]:
            m = res["metrics"][key]
            md.append(f"- {key}: spearman={m['spearman']:.4f}, mse={m['mse']:.4f}")
        ci = res["delta_model_minus_baseline_a_spearman_ci95"]
        md.append(
            f"- delta(model - baseline_a) spearman: {res['delta_model_minus_baseline_a_spearman']:.4f} (95% CI [{ci[0]:.4f}, {ci[1]:.4f}])"
        )
        md.append(
            f"- delta(model - baseline_b) spearman: {res['delta_model_minus_baseline_b_spearman']:.4f}"
        )
        md.append("")

    md.append("## Rule Evaluation")
    md.append(f"- Primary delta >= 0.03: {p_delta:.4f} -> {'PASS' if p_delta >= 0.03 else 'FAIL'}")
    md.append(f"- Primary CI excludes 0: [{p_ci_lo:.4f}, {p_ci_hi:.4f}] -> {'PASS' if p_ci_lo > 0 else 'FAIL'}")
    md.append(f"- External delta >= 0.01: {e_delta:.4f} -> {'PASS' if e_delta >= 0.01 else 'FAIL'}")
    md.append(f"- No underperform vs baseline_b (primary): {p_vs_b:.4f} -> {'PASS' if p_vs_b >= 0 else 'FAIL'}")
    md.append(f"- No underperform vs baseline_b (external): {e_vs_b:.4f} -> {'PASS' if e_vs_b >= 0 else 'FAIL'}")

    (outputs_dir / "gate_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    # also freeze schema manifest as csv convenience
    schema_rows = [
        {
            "field": "guide",
            "description": "20nt guide sequence (A/C/G/T only after normalization)",
        },
        {"field": "label", "description": "Measured on-target activity (modFreq from source)"},
        {"field": "source", "description": "Dataset source identifier"},
        {"field": "gene_or_target_group", "description": "Group key derived from guide identifier"},
        {"field": "target_context", "description": "Optional long context sequence from source"},
    ]
    write_csv(outputs_dir / "locked_schema_manifest.csv", schema_rows, ["field", "description"])

    print(f"decision: {decision}")
    print(f"wrote {outputs_dir / 'locked_dataset_manifest.json'}")
    print(f"wrote {outputs_dir / 'locked_split_manifest.csv'}")
    print(f"wrote {outputs_dir / 'gate_results.json'}")
    print(f"wrote {outputs_dir / 'gate_report.md'}")
    print(f"wrote {outputs_dir / 'locked_schema_manifest.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
