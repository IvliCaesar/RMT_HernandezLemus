"""
Real Dyson-Brownian-motion-style diagnostic: track the top eigenvalues
of the (real) TNBC gene-gene correlation matrix as the patient sample
is grown one patient at a time (in a fixed random inclusion order), and
look for level repulsion -- trajectories that approach but visibly
avoid crossing, rather than crossing freely, exactly as predicted by
the 1/(lambda_i - lambda_j) Coulomb-gas repulsion term in Dyson's (1962)
eigenvalue SDEs. This is the same phenomenon Aarts et al. (2024) find in
neural-network weight-matrix training dynamics; here the "time"
parameter is sample size rather than SGD steps, but the mechanism check
(do the top eigenvalues repel or cross freely) is the same kind of
question and is run here on real correlation matrices, not simulated
Dyson dynamics.

For each sample size n = n0, n0+1, ..., N (patients added one at a
time, cumulatively, in a fixed random order so this is a genuine
nested sequence of real matrices, not independent resamples), record
the top K eigenvalues of corrcoef(X[:, :n]) restricted to a smaller,
fixed gene subset (to keep this loop fast: TOP_GENES most variable
genes within TNBC). Reports the minimum gap seen between the two
closest trajectories across the whole run, and whether trajectories
that come close visibly bounce apart (repulsion) rather than crossing.
"""
import numpy as np
import pandas as pd

TOP_GENES = 60
K_EIGS = 8
N_START = 15

if __name__ == "__main__":
    expr = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    var = expr.var(axis=1)
    sub_genes = var.sort_values(ascending=False).index[:TOP_GENES]
    X = expr.loc[sub_genes].to_numpy()  # genes x patients
    p, n_total = X.shape
    print(f"Restricted to top {p} most-variable TNBC genes, "
          f"n_total={n_total} patients.")

    rng = np.random.default_rng(1)
    order = rng.permutation(n_total)
    Xp = X[:, order]

    ns = list(range(N_START, n_total + 1))
    traj = np.full((len(ns), K_EIGS), np.nan)
    for i, n in enumerate(ns):
        C = np.corrcoef(Xp[:, :n])
        eigs = np.sort(np.linalg.eigvalsh(C))[::-1][:K_EIGS]
        traj[i] = eigs

    print(f"\nTop-{K_EIGS} eigenvalue trajectories as n grows from "
          f"{N_START} to {n_total} (final values):")
    print(np.round(traj[-1], 3))

    # Minimum instantaneous gap between adjacent trajectories over the run.
    gaps = np.diff(traj, axis=1)  # traj sorted descending -> gaps are negative
    min_gap = np.abs(gaps).min()
    min_gap_idx = np.unravel_index(np.abs(gaps).argmin(), gaps.shape)
    n_at_min, pair_idx = ns[min_gap_idx[0]], min_gap_idx[1]
    print(f"\nSmallest adjacent-eigenvalue gap observed: {min_gap:.4f}, "
          f"at n={n_at_min}, between eigenvalues ranked "
          f"{pair_idx+1} and {pair_idx+2}.")

    # Check: does the gap between that pair shrink then grow again
    # (repulsion / avoided crossing) rather than hitting ~0 and crossing?
    gap_series = np.abs(traj[:, pair_idx] - traj[:, pair_idx + 1])
    min_idx = gap_series.argmin()
    print(f"Gap trajectory for that pair around the closest approach "
          f"(n={ns[max(0,min_idx-3)]}..{ns[min(len(ns)-1,min_idx+3)]}):")
    lo, hi = max(0, min_idx - 3), min(len(ns), min_idx + 4)
    for i in range(lo, hi):
        marker = "  <-- closest approach" if i == min_idx else ""
        print(f"  n={ns[i]:4d}  gap={gap_series[i]:.4f}{marker}")
    recovers = (gap_series[min(len(gap_series)-1, min_idx+3)] >
                1.5 * gap_series[min_idx])
    print(f"\nGap re-opens after closest approach "
          f"({'consistent with level repulsion' if recovers else 'inconclusive at this resolution'}).")

    np.savetxt("../results/dyson_trajectories.csv", traj, delimiter=",",
               header=",".join(f"eig{k+1}" for k in range(K_EIGS)),
               comments="")
    print("\nSaved full trajectory table to results/dyson_trajectories.csv")
