#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.linear_model import LogisticRegression

import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from applications.genomic_disruption_api import DisruptionAnalyzer
from validation.ontarget_gate.comparators.base import ComparatorRecord
from validation.ontarget_gate.comparators.registry import build_required_comparators


def clean_seq(seq: str) -> str:
    seq = (seq or "").upper()
    return "".join(ch for ch in seq if ch in "ACGTN")


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


def kmer2_features(seq: str) -> np.ndarray:
    kmers = [a + b for a in "ACGT" for b in "ACGT"]
    counts = {k: 0 for k in kmers}
    for i in range(len(seq) - 1):
        k = seq[i : i + 2]
        if k in counts:
            counts[k] += 1
    denom = max(1, len(seq) - 1)
    return np.array([counts[k] / denom for k in kmers], dtype=float)


def baseline_a_features(seq: str) -> np.ndarray:
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    hom = max_homopolymer_run(seq) / len(seq)
    return np.concatenate([np.array([gc, hom], dtype=float), kmer2_features(seq)])


def safe_spearman(y: np.ndarray, pred: np.ndarray) -> float:
    rho = float(spearmanr(y, pred).correlation)
    if np.isnan(rho):
        return 0.0
    return rho


def zscore(v: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    s = float(np.std(v))
    if s < eps:
        return np.zeros_like(v, dtype=float)
    return (v - float(np.mean(v))) / s


def gc_fraction(seq: str) -> float:
    return float((seq.count("G") + seq.count("C")) / len(seq))


@dataclass
class SplitPack:
    rows: List[dict]
    y: np.ndarray
    seqs: List[str]


def load_rows(manifest: Path) -> List[dict]:
    rows: List[dict] = []
    with manifest.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            seq = clean_seq(r.get("guide_seq", ""))
            if len(seq) != 20:
                continue
            rows.append(
                {
                    "source": r.get("source", ""),
                    "split": r.get("split", ""),
                    "guide": r.get("guide", ""),
                    "guide_seq": seq,
                    "label": float(r.get("label", "nan")),
                    "group": r.get("gene_or_target_group", ""),
                    "group_strength": r.get("group_strength", ""),
                }
            )
    return rows


def pack(rows: Sequence[dict], split: str) -> SplitPack:
    sub = [r for r in rows if r["split"] == split]
    return SplitPack(rows=sub, y=np.array([float(r["label"]) for r in sub], dtype=float), seqs=[r["guide_seq"] for r in sub])


def split_metrics(y: np.ndarray, pred: np.ndarray, baseline_c: np.ndarray) -> dict:
    s_model = safe_spearman(y, pred)
    s_c = safe_spearman(y, baseline_c)
    return {
        "spearman_model": float(s_model),
        "spearman_baseline_c": float(s_c),
        "delta_vs_baseline_c": float(s_model - s_c),
    }


def main() -> int:
    root = REPO_ROOT / "validation" / "ontarget_gate"
    outputs = root / "outputs"
    manifest = outputs / "locked_split_manifest_v3.csv"
    if not manifest.exists():
        manifest = outputs / "decision_split_manifest_v3_clean.csv"
    rows = load_rows(manifest)

    train_rows = [r for r in rows if r["split"] in {"dev_train", "aux_train_weak", "aux_train_strong"}]
    dev = pack(rows, "dev_val")
    primary = pack(rows, "primary_holdout")
    external = pack(rows, "external_holdout")

    X_train = np.vstack([baseline_a_features(r["guide_seq"]) for r in train_rows])
    y_train = np.array([float(r["label"]) for r in train_rows], dtype=float)
    ridge_a = Ridge(alpha=1.0)
    ridge_a.fit(X_train, y_train)

    da = DisruptionAnalyzer(k=0.3, seed=42)
    comparators = build_required_comparators()
    c = comparators["baseline_c"]
    c_records = {
        "dev": [ComparatorRecord(guide_seq=s, target_context="") for s in dev.seqs],
        "primary": [ComparatorRecord(guide_seq=s, target_context="") for s in primary.seqs],
        "external": [ComparatorRecord(guide_seq=s, target_context="") for s in external.seqs],
    }
    pred_c = {
        "dev": np.array(c.predict_batch(c_records["dev"]), dtype=float),
        "primary": np.array(c.predict_batch(c_records["primary"]), dtype=float),
        "external": np.array(c.predict_batch(c_records["external"]), dtype=float),
    }

    def preds_for(pack_: SplitPack) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        p_m = np.array([float(da.score_guide(s)["disruption_score"]) for s in pack_.seqs], dtype=float)
        X = np.vstack([baseline_a_features(s) for s in pack_.seqs])
        p_a = ridge_a.predict(X)
        gc = np.array([gc_fraction(s) for s in pack_.seqs], dtype=float)
        hom = np.array([max_homopolymer_run(s) / len(s) for s in pack_.seqs], dtype=float)
        return p_m, p_a, gc, hom

    p_dev_m, p_dev_a, p_dev_gc, p_dev_hom = preds_for(dev)
    p_primary_m, p_primary_a, p_primary_gc, p_primary_hom = preds_for(primary)
    p_external_m, p_external_a, p_external_gc, p_external_hom = preds_for(external)

    # A0: current model score
    p_a0 = {
        "dev": p_dev_m,
        "primary": p_primary_m,
        "external": p_external_m,
    }

    # A1: z-blend model + baseline_a (alpha selected on dev only)
    alphas = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
    best_a1_alpha = 0.0
    best_a1_dev = -1e9
    for a in alphas:
        pred = zscore(p_dev_m) + a * zscore(p_dev_a)
        d = safe_spearman(dev.y, pred) - safe_spearman(dev.y, pred_c["dev"])
        if d > best_a1_dev:
            best_a1_dev = d
            best_a1_alpha = a
    p_a1 = {
        "dev": zscore(p_dev_m) + best_a1_alpha * zscore(p_dev_a),
        "primary": zscore(p_primary_m) + best_a1_alpha * zscore(p_primary_a),
        "external": zscore(p_external_m) + best_a1_alpha * zscore(p_external_a),
    }

    # A2: GC debias model score, then blend with baseline_a (alpha selected on dev)
    beta_gc = float(np.polyfit(p_dev_gc, p_dev_m, deg=1)[0]) if len(p_dev_gc) > 2 else 0.0

    def debias(pred: np.ndarray, gc: np.ndarray) -> np.ndarray:
        return pred - beta_gc * (gc - float(np.mean(gc)))

    p_dev_m_db = debias(p_dev_m, p_dev_gc)
    p_primary_m_db = debias(p_primary_m, p_primary_gc)
    p_external_m_db = debias(p_external_m, p_external_gc)

    best_a2_alpha = 0.0
    best_a2_dev = -1e9
    for a in alphas:
        pred = zscore(p_dev_m_db) + a * zscore(p_dev_a)
        d = safe_spearman(dev.y, pred) - safe_spearman(dev.y, pred_c["dev"])
        if d > best_a2_dev:
            best_a2_dev = d
            best_a2_alpha = a
    p_a2 = {
        "dev": zscore(p_dev_m_db) + best_a2_alpha * zscore(p_dev_a),
        "primary": zscore(p_primary_m_db) + best_a2_alpha * zscore(p_primary_a),
        "external": zscore(p_external_m_db) + best_a2_alpha * zscore(p_external_a),
    }

    # A3: stacked ridge on model+baseline_a+gc+hom (fit on dev_train/aux only)
    p_train_m = np.array([float(da.score_guide(r["guide_seq"])["disruption_score"]) for r in train_rows], dtype=float)
    p_train_a = ridge_a.predict(np.vstack([baseline_a_features(r["guide_seq"]) for r in train_rows]))
    p_train_gc = np.array([gc_fraction(r["guide_seq"]) for r in train_rows], dtype=float)
    p_train_hom = np.array([max_homopolymer_run(r["guide_seq"]) / len(r["guide_seq"]) for r in train_rows], dtype=float)

    X_stack_train = np.vstack([p_train_m, p_train_a, p_train_gc, p_train_hom]).T
    stack = Ridge(alpha=1.0)
    stack.fit(X_stack_train, y_train)

    p_a3 = {
        "dev": stack.predict(np.vstack([p_dev_m, p_dev_a, p_dev_gc, p_dev_hom]).T),
        "primary": stack.predict(np.vstack([p_primary_m, p_primary_a, p_primary_gc, p_primary_hom]).T),
        "external": stack.predict(np.vstack([p_external_m, p_external_a, p_external_gc, p_external_hom]).T),
    }

    # A4: pairwise rank-optimized logistic model on meta features.
    # Train on sampled pairwise differences (x_i - x_j), target = 1{y_i > y_j}.
    rng = np.random.default_rng(42)
    X_train_meta = np.vstack([p_train_m, p_train_a, p_train_gc, p_train_hom]).T
    y_train_meta = y_train
    n_train = len(y_train_meta)
    max_pairs = 120000
    n_pairs = min(max_pairs, n_train * (n_train - 1) // 2)
    idx_i = rng.integers(0, n_train, size=n_pairs)
    idx_j = rng.integers(0, n_train, size=n_pairs)
    mask = idx_i != idx_j
    idx_i = idx_i[mask]
    idx_j = idx_j[mask]
    yi = y_train_meta[idx_i]
    yj = y_train_meta[idx_j]
    y_pair = (yi > yj).astype(int)
    non_tie = yi != yj
    idx_i = idx_i[non_tie]
    idx_j = idx_j[non_tie]
    y_pair = y_pair[non_tie]
    X_pair = X_train_meta[idx_i] - X_train_meta[idx_j]
    ranker = LogisticRegression(
        random_state=42,
        max_iter=1000,
        solver="lbfgs",
    )
    ranker.fit(X_pair, y_pair)
    p_a4 = {
        "dev": ranker.decision_function(np.vstack([p_dev_m, p_dev_a, p_dev_gc, p_dev_hom]).T),
        "primary": ranker.decision_function(np.vstack([p_primary_m, p_primary_a, p_primary_gc, p_primary_hom]).T),
        "external": ranker.decision_function(np.vstack([p_external_m, p_external_a, p_external_gc, p_external_hom]).T),
    }

    candidates = {
        "A0_current_model": p_a0,
        "A1_zblend_model_plus_baseline_a": p_a1,
        "A2_gc_debias_plus_zblend": p_a2,
        "A3_stacked_ridge": p_a3,
        "A4_pairwise_rank_logit": p_a4,
    }

    metrics_by_candidate: Dict[str, dict] = {}
    for name, pred in candidates.items():
        metrics_by_candidate[name] = {
            "dev_val": split_metrics(dev.y, pred["dev"], pred_c["dev"]),
            "primary_holdout": split_metrics(primary.y, pred["primary"], pred_c["primary"]),
            "external_holdout": split_metrics(external.y, pred["external"], pred_c["external"]),
        }

    winner = max(
        metrics_by_candidate.keys(),
        key=lambda n: metrics_by_candidate[n]["dev_val"]["delta_vs_baseline_c"],
    )

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "manifest_path": str(manifest),
        "n_rows": {
            "train": len(train_rows),
            "dev_val": len(dev.rows),
            "primary_holdout": len(primary.rows),
            "external_holdout": len(external.rows),
        },
        "model_config": {
            "disruption_analyzer_k": 0.3,
            "disruption_analyzer_seed": 42,
            "A1_alpha": best_a1_alpha,
            "A2_alpha": best_a2_alpha,
            "A2_beta_gc": beta_gc,
            "A3_model": "Ridge(alpha=1.0)",
            "A4_model": "LogisticRegression(pairwise differences, max_iter=1000)",
            "A4_pair_count": int(len(y_pair)),
        },
        "metrics_by_candidate": metrics_by_candidate,
        "selection_rule": "max dev_val delta_vs_baseline_c",
        "winner": winner,
        "winner_holdout_deltas": {
            "primary": metrics_by_candidate[winner]["primary_holdout"]["delta_vs_baseline_c"],
            "external": metrics_by_candidate[winner]["external_holdout"]["delta_vs_baseline_c"],
        },
    }

    out_json = outputs / "recovery_ablation_v1.json"
    out_json.write_text(json.dumps(out, indent=2), encoding="utf-8")

    md_lines = [
        "# Recovery Ablation v1",
        "",
        f"Generated: {out['generated_at']}",
        "",
        "Selection rule: max dev_val delta vs baseline_c.",
        "",
        f"Winner: **{winner}**",
        "",
        "## Candidate metrics (delta model - baseline_c, Spearman)",
    ]
    for name, m in metrics_by_candidate.items():
        md_lines.append(
            f"- {name}: dev={m['dev_val']['delta_vs_baseline_c']:.4f}, "
            f"primary={m['primary_holdout']['delta_vs_baseline_c']:.4f}, "
            f"external={m['external_holdout']['delta_vs_baseline_c']:.4f}"
        )
    md_lines.append("")
    md_lines.append("## Winner holdout deltas")
    md_lines.append(f"- primary: {out['winner_holdout_deltas']['primary']:.4f}")
    md_lines.append(f"- external: {out['winner_holdout_deltas']['external']:.4f}")

    out_md = outputs / "recovery_ablation_v1.md"
    out_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"winner: {winner}")
    print(f"wrote {out_json}")
    print(f"wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
