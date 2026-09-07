"""
Direct test of the mechanistic prediction scripts 19-22 make but do not
themselves test: if weight-matrix Dyson repulsion is active only while
the matrix is genuinely still changing (Sec. results-dyson-classifier),
then artificially lengthening the pre-convergence transient -- without
changing the task, architecture, or feature set at all -- should
lengthen the window in which repulsion is detectable by a corresponding
amount. The original run used learning_rate_init=0.05 and plateaued by
epoch ~50 in every seed; this script uses the identical MLP(20->8->1),
identical 20 RMT-hub features, identical 549-patient cohort, and
identical 20 seeds, but a 5x smaller learning rate (0.01) and a 5x
longer run (1500 epochs instead of 300) -- the most direct, single-knob
way to slow convergence without touching anything else that could
confound the comparison.

Rather than reusing the fixed epoch-50 cutoff (which assumed the
original learning rate's convergence speed), this script detects each
seed's own empirical convergence epoch directly from its accuracy
trace (first epoch after which training accuracy never again drops
more than 0.5 points below its own final value) and splits early/late
at that per-seed point -- the fair, non-circular way to ask whether a
longer transient produces a longer repulsion window, since a fixed
global cutoff would beg the question.
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
N_EPOCHS = 1500
LR = 0.01
N_SEEDS = 20


def rmt_hubs(expr, thresholds, k=20):
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


def weight_trajectory_and_accuracy(X, y, seed):
    clf = MLPClassifier(hidden_layer_sizes=(K_HIDDEN,), activation="tanh",
                         solver="sgd", learning_rate_init=LR, alpha=0.0,
                         momentum=0.0, batch_size=64, max_iter=1,
                         warm_start=True, random_state=seed, shuffle=True)
    traj = np.full((N_EPOCHS, K_HIDDEN), np.nan)
    acc = np.full(N_EPOCHS, np.nan)
    for e in range(N_EPOCHS):
        clf.fit(X, y)
        W1 = clf.coefs_[0]
        gram = W1.T @ W1
        traj[e] = np.sort(np.linalg.eigvalsh(gram))[::-1]
        acc[e] = clf.score(X, y)
    return traj, acc


def empirical_convergence_epoch(acc, tol=0.005):
    """First epoch after which accuracy never again drops more than `tol`
    below its own final value."""
    final = acc[-1]
    below = np.where(acc < final - tol)[0]
    if len(below) == 0:
        return 0
    return int(below[-1]) + 1


def gap_pattern(gaps, window=5):
    if len(gaps) < 3:
        return "too short"
    min_idx = int(np.argmin(gaps))
    min_gap = gaps[min_idx]
    lo = max(0, min_idx - window)
    hi = min(len(gaps) - 1, min_idx + window)
    reopens_after = (hi > min_idx) and (gaps[hi] > 1.5 * min_gap)
    reopens_before = (lo < min_idx) and (gaps[lo] > 1.5 * min_gap)
    at_boundary = (min_idx == len(gaps) - 1)
    if at_boundary:
        return "inconclusive"
    elif reopens_after and reopens_before:
        return "repulsion"
    elif reopens_after or reopens_before:
        return "partial"
    else:
        return "none"


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    tnbc = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    luma = pd.read_csv("../data/TCGA-BRCA/luma_expr.csv", index_col=0)
    combined = pd.concat([tnbc, luma], axis=1)
    y = np.array([1] * tnbc.shape[1] + [0] * luma.shape[1])

    thresholds = np.round(np.arange(0.15, 0.85, 0.05), 2)
    hub_tnbc = rmt_hubs(tnbc, thresholds)
    hub_luma = rmt_hubs(luma, thresholds)
    features = list(dict.fromkeys(hub_tnbc + hub_luma))[:20]
    Xraw = combined.loc[features].T.to_numpy()
    X = StandardScaler().fit_transform(Xraw)

    print(f"MLP(20->{K_HIDDEN}->1), lr={LR} (1/5 of the original 0.05), "
          f"{N_EPOCHS} epochs (5x the original 300), {N_SEEDS} seeds.")

    all_traj = np.full((N_SEEDS, N_EPOCHS, K_HIDDEN), np.nan)
    conv_epochs = np.zeros(N_SEEDS, dtype=int)
    final_accs = np.zeros(N_SEEDS)
    for seed in range(N_SEEDS):
        traj, acc = weight_trajectory_and_accuracy(X, y, seed)
        all_traj[seed] = traj
        conv_epochs[seed] = empirical_convergence_epoch(acc)
        final_accs[seed] = acc[-1]
        print(f"  seed {seed:2d}: final acc={final_accs[seed]:.4f}, "
              f"empirical convergence epoch={conv_epochs[seed]}")

    np.save("../results/dyson_classifier_longer_transient_all_seeds.npy", all_traj)
    pd.DataFrame({"seed": range(N_SEEDS), "convergence_epoch": conv_epochs,
                  "final_accuracy": final_accs}).to_csv(
        "../results/dyson_classifier_longer_transient_convergence.csv", index=False)

    print(f"\nConvergence epoch across seeds: mean={conv_epochs.mean():.1f}, "
          f"median={np.median(conv_epochs):.1f}, min={conv_epochs.min()}, "
          f"max={conv_epochs.max()} "
          f"(original lr=0.05 run converged by ~epoch 50 in every seed; "
          f"a longer convergence epoch here is the direct test of whether "
          f"slowing the learning rate actually lengthens the transient).")

    # Early/late split at each seed's OWN empirical convergence epoch,
    # not a fixed global cutoff (which would beg the question).
    rows = []
    for seed in range(N_SEEDS):
        ce = max(conv_epochs[seed], 3)  # need >=3 points for a gap series
        traj = all_traj[seed]
        early = traj[:ce]
        late = traj[ce:]
        for pair in range(K_HIDDEN - 1):
            if len(early) >= 3:
                gaps_e = np.abs(np.diff(early, axis=1))[:, pair]
                pat_e = gap_pattern(gaps_e)
            else:
                pat_e = "too short"
            if len(late) >= 3:
                gaps_l = np.abs(np.diff(late, axis=1))[:, pair]
                pat_l = gap_pattern(gaps_l)
            else:
                pat_l = "too short"
            rows.append(dict(seed=seed, pair=f"mu{pair+1}-mu{pair+2}",
                              convergence_epoch=ce, early_pattern=pat_e,
                              late_pattern=pat_l))
    df = pd.DataFrame(rows)
    df.to_csv("../results/dyson_classifier_longer_transient_early_late.csv", index=False)

    def signal_count(col):
        return int(df[col].isin(["repulsion", "partial"]).sum())

    early_signal = signal_count("early_pattern")
    late_signal = signal_count("late_pattern")
    max_possible = N_SEEDS * (K_HIDDEN - 1)
    print(f"\nEarly (pre-own-convergence) repulsion+partial: "
          f"{early_signal}/{max_possible}")
    print(f"Late (post-own-convergence) repulsion+partial:  "
          f"{late_signal}/{max_possible}")
    ratio = (early_signal / max(late_signal, 1))
    print(f"Early:late ratio = {ratio:.2f} "
          f"(original lr=0.05/300-epoch run: 36/3 = 12.0 at a fixed "
          f"epoch-50 cutoff and a ~50-epoch transient; compare against "
          f"this run's own empirical transient length, "
          f"mean={conv_epochs.mean():.1f} epochs).")

    print("\nSaved dyson_classifier_longer_transient_all_seeds.npy, "
          "_convergence.csv, _early_late.csv")
