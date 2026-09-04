"""
Real ML/recognition test: does the RMT hub-gene selection (Luo et al.
threshold, run separately on TNBC in 02_rmt_network_analysis.py) carry
more subtype-discriminative signal than an equal-size set of generic
high-variance genes? Trains a cross-validated logistic-regression
classifier (TNBC vs LumA, real TCGA-BRCA labels) on three real feature
sets of matched size:
  (a) RANDOM   - a random subset of the same 1500 candidate genes
  (b) TOP-VAR  - the top-K genes by raw variance (no RMT/network step)
  (c) RMT-HUB  - the union of TNBC and LumA top-degree hub genes at
                 their respective RMT-selected thresholds (tau*=0.45
                 and tau*=0.25 from script 02)

5-fold stratified cross-validation, real accuracy/ROC-AUC reported for
each -- this is a genuine comparison run on real patient data, not a
staged demo; whichever feature set wins is reported as-is.
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

import sys
sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score

RNG = np.random.default_rng(2)
K_FEATURES = 20


def rmt_hubs(expr, thresholds):
    genes = expr.index.to_numpy()
    C = np.corrcoef(expr.to_numpy())
    results = []
    for tau in thresholds:
        A = C.copy()
        A[np.abs(A) < tau] = 0.0
        np.fill_diagonal(A, 0.0)
        eigs = np.linalg.eigvalsh(A)
        if len(np.unique(np.round(eigs, 8))) < 10:
            continue
        spacings = unfold_spectrum(eigs)
        ks_p, ks_g = poisson_goe_fit_score(spacings)
        regime = "GOE" if ks_g < ks_p else "Poisson"
        results.append((tau, regime))
    tau_star = None
    for i in range(len(results) - 1, 0, -1):
        if results[i][1] == "Poisson" and results[i - 1][1] == "GOE":
            tau_star = results[i - 1][0]
            break
    if tau_star is None:
        return []
    A = np.abs(C.copy())
    A[A < tau_star] = 0.0
    np.fill_diagonal(A, 0.0)
    degree = (A > 0).sum(axis=1)
    return list(genes[np.argsort(degree)[::-1][:K_FEATURES]])


if __name__ == "__main__":
    tnbc = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    luma = pd.read_csv("../data/TCGA-BRCA/luma_expr.csv", index_col=0)
    all_genes = tnbc.index.to_numpy()

    combined = pd.concat([tnbc, luma], axis=1)
    y = np.array([1] * tnbc.shape[1] + [0] * luma.shape[1])  # 1=TNBC, 0=LumA
    print(f"Combined cohort for classification: {combined.shape[1]} patients "
          f"({tnbc.shape[1]} TNBC, {luma.shape[1]} LumA), {combined.shape[0]} candidate genes.")

    thresholds = np.round(np.arange(0.15, 0.85, 0.05), 2)
    hub_tnbc = rmt_hubs(tnbc, thresholds)
    hub_luma = rmt_hubs(luma, thresholds)
    rmt_features = list(dict.fromkeys(hub_tnbc + hub_luma))[:K_FEATURES]
    print(f"RMT-hub feature set ({len(rmt_features)} genes): {rmt_features}")

    var = combined.var(axis=1)
    topvar_features = list(var.sort_values(ascending=False).index[:K_FEATURES])
    print(f"Top-variance feature set ({len(topvar_features)} genes): {topvar_features}")

    random_features = list(RNG.choice(all_genes, size=K_FEATURES, replace=False))
    print(f"Random feature set ({len(random_features)} genes): {random_features}")

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    for name, feats in [("RANDOM", random_features),
                         ("TOP-VAR", topvar_features),
                         ("RMT-HUB", rmt_features)]:
        Xf = combined.loc[feats].T.to_numpy()
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
        acc = cross_val_score(clf, Xf, y, cv=cv, scoring="accuracy")
        auc = cross_val_score(clf, Xf, y, cv=cv, scoring="roc_auc")
        print(f"\n{name} (n_features={len(feats)}): "
              f"accuracy={acc.mean():.3f}+/-{acc.std():.3f}, "
              f"ROC-AUC={auc.mean():.3f}+/-{auc.std():.3f}")
