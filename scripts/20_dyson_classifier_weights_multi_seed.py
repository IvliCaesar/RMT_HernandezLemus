"""
Reproducibility check for script 19's classifier-weight-matrix Dyson
trajectory, mirroring exactly what script 17 did for the data-side
trajectories of script 05: a single random seed is not evidence on its
own (script 05's own lambda2-lambda3 pair looked inconclusive under one
seed and turned out to repel in 16/20 orders once checked). This
script reruns the identical eig(W1^T W1)-trajectory scan across
N_SEEDS=20 independent MLP initializations / SGD minibatch shufflings
(same 20 RMT-hub-gene input, same 8-unit hidden layer, same 300 SGD
epochs) and classifies each adjacent eigenvalue pair's behavior with
the same closest-approach / 5-epoch-window reopening criterion used in
scripts 05 and 17.
"""
import warnings
import os
import sys
import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning

sys.path.insert(0, ".")
from importlib import import_module

warnings.filterwarnings("ignore", category=ConvergenceWarning)

N_SEEDS = 20
N_EPOCHS = 300
K_HIDDEN = 8


def gap_pattern(gaps, window=5):
    min_idx = int(np.argmin(gaps))
    min_gap = gaps[min_idx]
    lo = max(0, min_idx - window)
    hi = min(len(gaps) - 1, min_idx + window)
    reopens_after = (hi > min_idx) and (gaps[hi] > 1.5 * min_gap)
    reopens_before = (lo < min_idx) and (gaps[lo] > 1.5 * min_gap)
    at_boundary = (min_idx == len(gaps) - 1)
    if at_boundary:
        return "inconclusive (still narrowing at epoch_max)"
    elif reopens_after and reopens_before:
        return "repulsion (both sides reopen)"
    elif reopens_after or reopens_before:
        return "partial repulsion (one side reopens)"
    else:
        return "no clear repulsion"


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    m19 = import_module("19_dyson_classifier_weights")

    tnbc = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    luma = pd.read_csv("../data/TCGA-BRCA/luma_expr.csv", index_col=0)
    combined = pd.concat([tnbc, luma], axis=1)
    y = np.array([1] * tnbc.shape[1] + [0] * luma.shape[1])

    thresholds = np.round(np.arange(0.15, 0.85, 0.05), 2)
    hub_tnbc = m19.rmt_hubs(tnbc, thresholds)
    hub_luma = m19.rmt_hubs(luma, thresholds)
    features = list(dict.fromkeys(hub_tnbc + hub_luma))[:20]

    from sklearn.preprocessing import StandardScaler
    Xraw = combined.loc[features].T.to_numpy()
    X = StandardScaler().fit_transform(Xraw)

    print(f"Running {N_SEEDS} independent MLP(20->{K_HIDDEN}->1) training "
          f"runs, {N_EPOCHS} epochs each, tracking eig(W1^T W1)...")

    rows = []
    all_traj = []
    for seed in range(N_SEEDS):
        traj, clf = m19.weight_trajectory(X, y, seed=seed,
                                           n_epochs=N_EPOCHS,
                                           n_hidden=K_HIDDEN)
        all_traj.append(traj)
        gaps_all = np.abs(np.diff(traj, axis=1))
        for pair_idx in range(K_HIDDEN - 1):
            pattern = gap_pattern(gaps_all[:, pair_idx])
            rows.append(dict(seed=seed,
                              pair=f"mu{pair_idx+1}-mu{pair_idx+2}",
                              pattern=pattern,
                              final_train_acc=clf.score(X, y)))
        print(f"  seed={seed}: final train acc={clf.score(X, y):.4f}")

    df = pd.DataFrame(rows)
    df.to_csv("../results/dyson_classifier_multi_seed.csv", index=False)

    # Save one representative trajectory (seed 0) again for the figure,
    # and the full stack for a reference (not plotted, kept for the record).
    np.save("../results/dyson_classifier_all_seeds.npy", np.array(all_traj))

    print("\nReproducibility across 20 seeds, by weight-matrix eigenvalue pair:")
    summary_rows = []
    for pair in df["pair"].unique():
        sub = df[df["pair"] == pair]
        counts = sub["pattern"].value_counts()
        n_rep = counts.get("repulsion (both sides reopen)", 0)
        n_partial = counts.get("partial repulsion (one side reopens)", 0)
        n_inconclusive = counts.get("inconclusive (still narrowing at epoch_max)", 0)
        n_none = counts.get("no clear repulsion", 0)
        print(f"  {pair}: repulsion={n_rep}/{N_SEEDS}, partial={n_partial}/{N_SEEDS}, "
              f"inconclusive={n_inconclusive}/{N_SEEDS}, none={n_none}/{N_SEEDS}")
        summary_rows.append(dict(pair=pair, repulsion=n_rep, partial=n_partial,
                                  inconclusive=n_inconclusive, none=n_none))
    pd.DataFrame(summary_rows).to_csv(
        "../results/dyson_classifier_multi_seed_summary.csv", index=False)
    print("\nSaved results/dyson_classifier_multi_seed.csv, "
          "_summary.csv, and dyson_classifier_all_seeds.npy")
