"""
Doench-2016 nested-model evaluation:
- Baseline features (light sequence features)
- + Breathing spectral features (10.5 / 5.25 / 3.5)
Reports ΔR^2 (regression) and, if `label` is present, ΔAUPRC.
CSV schema expected:
  guide, target, y [, label]
"""
from __future__ import annotations
import argparse, numpy as np, pandas as pd
from sklearn.model_selection import KFold
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score, average_precision_score
from scipy.stats import spearmanr

from wave_crispr_signal.spectral import baseline_seq_features, breathing_features

def build_X(df: pd.DataFrame, use_breathing: bool, r: float = 20.0) -> np.ndarray:
    rows = []
    for s in df['target']:
        feats = baseline_seq_features(s)
        if use_breathing:
            feats.update(breathing_features(s, r=r))
        rows.append(feats)
    X = pd.DataFrame(rows).astype(float)
    return X

def evaluate(df: pd.DataFrame, r: float = 20.0, folds: int = 5, alpha: float = 1.0):
    kf = KFold(n_splits=folds, shuffle=True, random_state=42)
    metrics = {"baseline": {"r2": [], "rho": [], "auprc": []},
               "+breathing": {"r2": [], "rho": [], "auprc": []}}
    y = df["y"].values.astype(float)
    y_label = df["label"].values.astype(int) if "label" in df.columns else None

    for use_breathing in [False, True]:
        X = build_X(df, use_breathing=use_breathing, r=r).values
        for tr, te in kf.split(X):
            mdl = Ridge(alpha=alpha)
            mdl.fit(X[tr], y[tr])
            yhat = mdl.predict(X[te])
            metrics["baseline" if not use_breathing else "+breathing"]["r2"].append(r2_score(y[te], yhat))
            metrics["baseline" if not use_breathing else "+breathing"]["rho"].append(spearmanr(y[te], yhat).correlation)
            if y_label is not None:
                # if labels exist, treat yhat as score
                auprc = average_precision_score(y_label[te], yhat)
                metrics["baseline" if not use_breathing else "+breathing"]["auprc"].append(auprc)
    return metrics

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="Path to Doench-2016-like CSV with columns: guide,target,y[,label]")
    ap.add_argument("--r", type=float, default=20.0, help="AT:GC opening-rate ratio (dimensionless)")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--alpha", type=float, default=1.0)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    out = evaluate(df, r=args.r, folds=args.folds, alpha=args.alpha)
    def s(x): 
        return f"{np.mean(x):.4f} ± {np.std(x):.4f}"
    print("Baseline    R2:", s(out["baseline"]["r2"]))
    print("+Breathing  R2:", s(out["+breathing"]["r2"]))
    print("ΔR2:", f"{(np.mean(out['+breathing']['r2']) - np.mean(out['baseline']['r2'])):.4f}")
    print("Baseline   rho:", s(out["baseline"]["rho"]))
    print("+Breathing rho:", s(out["+breathing"]["rho"]))
    if len(out["baseline"]["auprc"]):
        print("Baseline  AUPRC:", s(out["baseline"]["auprc"]))
        print("+Breath.  AUPRC:", s(out["+breathing"]["auprc"]))
        print("ΔAUPRC:", f"{(np.mean(out['+breathing']['auprc']) - np.mean(out['baseline']['auprc'])):.4f}")

if __name__ == "__main__":
    main()
