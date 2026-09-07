"""
Follow-up to scripts 19-20: the 20-seed reproducibility check found the
classifier's weight-matrix eigenvalue pairs far less reliably repulsive
than the data-side correlation-matrix pairs (script 17). Before
reporting that as a flat null, this script tests the obvious
alternative hypothesis directly, the same way script 12's synthetic
positive control rescued the wavelet line from a flat null: maybe
repulsion is real but confined to the early SGD transient, before
training converges (train accuracy plateaus by roughly epoch 50 in
every seed) and the weight matrix settles near a near-fixed point,
where the deterministic Dyson drift no longer dominates the SGD noise
term the way it does while the matrix is still changing rapidly.

Reruns the identical closest-approach/reopening classification
(scripts 05/17's criterion) separately on the early window
(epochs 0-49) and the late window (epochs 50-299) of the same 20
trajectories already computed and saved by script 20
(results/dyson_classifier_all_seeds.npy), rather than retraining.
"""
import numpy as np
import pandas as pd
import os


def gap_pattern(gaps, window=5):
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
    arr = np.load("../results/dyson_classifier_all_seeds.npy")  # (20, 300, 8)
    n_seeds, n_epochs, k_hidden = arr.shape
    print(f"Loaded {n_seeds} seeds x {n_epochs} epochs x {k_hidden} eigenvalues.")

    windows = {"early (epochs 0-49, pre-convergence transient)": slice(0, 50),
               "late (epochs 50-299, post-convergence)": slice(50, n_epochs)}

    rows = []
    for label, sl in windows.items():
        print(f"\n--- {label} ---")
        for pair in range(k_hidden - 1):
            counts = {"repulsion": 0, "partial": 0, "inconclusive": 0, "none": 0}
            for seed in range(n_seeds):
                traj = arr[seed, sl, :]
                gaps = np.abs(np.diff(traj, axis=1))[:, pair]
                counts[gap_pattern(gaps)] += 1
            print(f"  mu{pair+1}-mu{pair+2}: {counts}")
            rows.append(dict(window=label, pair=f"mu{pair+1}-mu{pair+2}", **counts))

    df = pd.DataFrame(rows)
    df.to_csv("../results/dyson_classifier_early_late.csv", index=False)

    early = df[df["window"].str.startswith("early")]
    late = df[df["window"].str.startswith("late")]
    early_signal = (early["repulsion"] + early["partial"]).sum()
    late_signal = (late["repulsion"] + late["partial"]).sum()
    print(f"\nTotal repulsion+partial counts, summed over all 7 pairs x 20 seeds "
          f"(max possible 140): early={early_signal}, late={late_signal}")
    print("Saved results/dyson_classifier_early_late.csv")
