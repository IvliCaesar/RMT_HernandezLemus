"""
Attacks Roadmap Line 6b directly (Sec. results-dyson only tracked the
DATA's correlation-matrix eigenvalues as sample size grows; Aarts et al.
(2024) show the same Dyson repulsion mechanism governs a neural
network's own WEIGHT-matrix eigenvalues during SGD training). This
script runs that check for real on the classifier already trained in
Sec. ml-results, rather than leaving it as a proposed extension.

A raw p x h weight matrix W is not itself real symmetric (its
eigenvalues are, generically, complex, in the Ginibre rather than GOE
universality class), so the direct analogue of the GOE-based repulsion
test used everywhere else in this paper is not W's own eigenvalues but
the real symmetric Gram matrix G = W^T W (h x h, positive
semi-definite) -- exactly the same move (raw non-symmetric data matrix
-> real symmetric correlation matrix) already used throughout this
paper. G's eigenvalues are always real and can be tracked epoch by
epoch with the identical closest-approach / reopening test used in
scripts 05 and 17.

Classifier: input = the same 20 RMT-hub genes used in Sec. ml-results
(Table tab:hubs, tau*=0.45 TNBC / tau*=0.25 LumA union), one hidden
layer of 8 units (matching K_EIGS=8 used for the data-side trajectories
in scripts 05/17, so the two are directly comparable size-for-size),
single logistic output, trained by plain SGD (no momentum) on the same
549-patient combined TNBC+LumA cohort, one epoch at a time
(warm_start=True, max_iter=1 per call) so W1 = coefs_[0] (20x8) can be
captured after every epoch.
"""
import warnings
import os
import numpy as np
import pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning

import sys
sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score

warnings.filterwarnings("ignore", category=ConvergenceWarning)

K_HIDDEN = 8
N_EPOCHS = 300


def rmt_hubs(expr, thresholds, k=20):
    """Same recipe as script 06's rmt_hubs: RMT-threshold hub genes."""
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
    return list(genes[np.argsort(degree)[::-1][:k]])


def weight_trajectory(X, y, seed, n_epochs=N_EPOCHS, n_hidden=K_HIDDEN):
    clf = MLPClassifier(hidden_layer_sizes=(n_hidden,), activation="tanh",
                         solver="sgd", learning_rate_init=0.05, alpha=0.0,
                         momentum=0.0, batch_size=64, max_iter=1,
                         warm_start=True, random_state=seed, shuffle=True)
    traj = np.full((n_epochs, n_hidden), np.nan)
    for e in range(n_epochs):
        clf.fit(X, y)
        W1 = clf.coefs_[0]           # (n_features, n_hidden), NOT symmetric
        gram = W1.T @ W1             # (n_hidden, n_hidden), real symmetric PSD
        eigs = np.sort(np.linalg.eigvalsh(gram))[::-1]
        traj[e] = eigs
    return traj, clf


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    tnbc = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    luma = pd.read_csv("../data/TCGA-BRCA/luma_expr.csv", index_col=0)
    combined = pd.concat([tnbc, luma], axis=1)
    y = np.array([1] * tnbc.shape[1] + [0] * luma.shape[1])  # 1=TNBC, 0=LumA
    print(f"Combined cohort: {combined.shape[1]} patients "
          f"({tnbc.shape[1]} TNBC, {luma.shape[1]} LumA).")

    thresholds = np.round(np.arange(0.15, 0.85, 0.05), 2)
    hub_tnbc = rmt_hubs(tnbc, thresholds)
    hub_luma = rmt_hubs(luma, thresholds)
    features = list(dict.fromkeys(hub_tnbc + hub_luma))[:20]
    print(f"Using {len(features)} RMT-hub genes as classifier input "
          f"(identical feature set to Sec. ml-results Table tab:hubs): "
          f"{features}")

    Xraw = combined.loc[features].T.to_numpy()
    X = StandardScaler().fit_transform(Xraw)

    print(f"\nTraining MLP(20 -> {K_HIDDEN} -> 1), tracking "
          f"eig(W1^T W1) for {N_EPOCHS} epochs, seed=0...")
    traj, clf = weight_trajectory(X, y, seed=0)
    np.savetxt("../results/dyson_classifier_weights.csv", traj,
               delimiter=",",
               header=",".join(f"eig{k+1}" for k in range(K_HIDDEN)),
               comments="")
    print("Final-epoch eigenvalues of W1^T W1:", np.round(traj[-1], 4))
    print(f"Final training accuracy: {clf.score(X, y):.4f}")

    gaps = np.abs(np.diff(traj, axis=1))  # (n_epochs, K_HIDDEN-1)
    min_gap = gaps.min()
    idx = np.unravel_index(gaps.argmin(), gaps.shape)
    print(f"\nSmallest adjacent-eigenvalue gap: {min_gap:.5f} at epoch "
          f"{idx[0]}, pair mu{idx[1]+1}-mu{idx[1]+2}")
    lo, hi = max(0, idx[0] - 5), min(N_EPOCHS, idx[0] + 6)
    gap_series = gaps[:, idx[1]]
    print("Gap around closest approach:")
    for e in range(lo, hi):
        marker = "  <-- closest approach" if e == idx[0] else ""
        print(f"  epoch={e:4d}  gap={gap_series[e]:.4f}{marker}")

    print("\nSaved results/dyson_classifier_weights.csv")
