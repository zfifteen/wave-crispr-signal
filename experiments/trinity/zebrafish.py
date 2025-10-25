"""
Cross-species validation: Train on human (Doench), test on zebrafish (CRISPRscan).
Reports transfer performance for baseline vs. +breathing models.
"""
from __future__ import annotations
import argparse, numpy as np, pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score, spearmanr
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_csv", required=True, help="Human training data (Doench-like)")
    ap.add_argument("--test_csv", required=True, help="Zebrafish test data (CRISPRscan-like)")
    ap.add_argument("--r", type=float, default=20.0)
    ap.add_argument("--alpha", type=float, default=1.0)
    args = ap.parse_args()

    train_df = pd.read_csv(args.train_csv)
    test_df = pd.read_csv(args.test_csv)
    
    X_train = build_X(train_df, use_breathing=True, r=args.r)
    X_test = build_X(test_df, use_breathing=True, r=args.r)
    y_train = train_df['y'].values
    y_test = test_df['y'].values
    
    mdl = Ridge(alpha=args.alpha)
    mdl.fit(X_train.values, y_train)
    yhat = mdl.predict(X_test.values)
    
    r2 = r2_score(y_test, yhat)
    rho = spearmanr(y_test, yhat).correlation
    
    print(f"Cross-species R2: {r2:.4f}")
    print(f"Cross-species rho: {rho:.4f}")

if __name__ == "__main__":
    main()
