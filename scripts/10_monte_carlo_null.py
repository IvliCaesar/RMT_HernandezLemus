"""
Monte Carlo permutation-null test for the RMT threshold (deepens
Sec. "RMT threshold selection" beyond the asymptotic Marchenko-Pastur/
Wigner-surmise comparison alone, which is what the group's own NIMEA
pipeline's bootstrap-significance step also does for its
mutual-information threshold -- see article.tex's ARACNE/NIMEA
discussion). For each cohort:

  1. Independently permute each gene's values across patients (shuffle
     within each row). This destroys every real gene-gene correlation
     while keeping each gene's own marginal distribution exactly fixed
     -- the standard permutation null for coexpression networks.
  2. Rerun the identical RMT threshold sweep on each permuted matrix.
  3. Record whether a Poisson->GOE transition appears at all under
     pure permutation, and if so at what tau.

If the real data's tau* sits outside the range ever produced by the
permutation null (or if the null never shows a transition at all), the
selected threshold reflects real correlation structure, not an
artifact of the asymptotic theoretical null used elsewhere in this
paper.
"""
import sys
import numpy as np
import pandas as pd
import os

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score

THRESHOLDS = np.round(np.arange(0.10, 0.85, 0.05), 2)
N_PERMS = 20


def find_tau_star(C, thresholds):
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
    for i in range(len(results) - 1, 0, -1):
        if results[i][1] == "Poisson" and results[i - 1][1] == "GOE":
            return results[i - 1][0]
    return None


def any_goe_at_all(C, thresholds):
    """Whether ANY threshold in the sweep shows GOE-like statistics."""
    for tau in thresholds:
        A = C.copy()
        A[np.abs(A) < tau] = 0.0
        np.fill_diagonal(A, 0.0)
        eigs = np.linalg.eigvalsh(A)
        if len(np.unique(np.round(eigs, 8))) < 10:
            continue
        spacings = unfold_spectrum(eigs)
        ks_p, ks_g = poisson_goe_fit_score(spacings)
        if ks_g < ks_p:
            return True
    return False


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    rng = np.random.default_rng(7)
    rows = []

    for name, path, real_tau_star in [
        ("TNBC", "../data/TCGA-BRCA/tnbc_expr.csv", 0.45),
        ("LumA", "../data/TCGA-BRCA/luma_expr.csv", 0.25),
    ]:
        expr = pd.read_csv(path, index_col=0)
        X = expr.to_numpy()
        p, n = X.shape
        print(f"=== {name}: p={p}, n={n}, real tau*={real_tau_star} ===")

        perm_taus = []
        perm_any_goe = 0
        for rep in range(N_PERMS):
            Xp = X.copy()
            for i in range(p):
                Xp[i] = rng.permutation(Xp[i])
            C_perm = np.corrcoef(Xp)
            tau = find_tau_star(C_perm, THRESHOLDS)
            has_goe = any_goe_at_all(C_perm, THRESHOLDS)
            perm_any_goe += int(has_goe)
            perm_taus.append(tau)
            rows.append(dict(cohort=name, rep=rep, tau_star=tau,
                              any_goe=has_goe))
            print(f"  perm {rep}: tau* = {tau}, any GOE regime found = {has_goe}")

        found = [t for t in perm_taus if t is not None]
        print(f"  Permutation null ({N_PERMS} reps): "
              f"{len(found)}/{N_PERMS} showed ANY Poisson->GOE transition; "
              f"{perm_any_goe}/{N_PERMS} showed GOE-like statistics at ANY "
              f"threshold in the sweep.")
        if found:
            print(f"  Permutation-null tau* range: "
                  f"[{min(found):.2f}, {max(found):.2f}], "
                  f"vs.\ real tau*={real_tau_star}")
        print()

    pd.DataFrame(rows).to_csv("../results/monte_carlo_null.csv", index=False)
    print("Saved results/monte_carlo_null.csv")
