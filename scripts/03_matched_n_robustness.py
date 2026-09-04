"""
Robustness check: the TNBC (n=119) and LumA (n=430) cohorts in
02_rmt_network_analysis.py have very different p/n ratios (q=12.61 vs
q=3.49), so their RMT-selected thresholds (0.45 vs 0.25) are not
directly comparable -- a larger q shifts the Marchenko-Pastur bulk
outward and changes how much correlation is needed before the spectrum
looks non-random. This script re-runs the identical method with LumA
subsampled down to n=119 (30 independent random subsamples, same seed
family) to see whether the TNBC/LumA threshold gap survives once q is
matched, or whether it was mostly a sample-size artifact.

Real data throughout (same real expression matrices, just resampled
patients), no synthetic component.
"""
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score

N_TNBC = 119
N_REPS = 30
THRESHOLDS = np.round(np.arange(0.15, 0.85, 0.05), 2)


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


if __name__ == "__main__":
    tnbc_expr = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    luma_expr = pd.read_csv("../data/TCGA-BRCA/luma_expr.csv", index_col=0)

    C_tnbc = np.corrcoef(tnbc_expr.to_numpy())
    tau_tnbc = find_tau_star(C_tnbc, THRESHOLDS)
    print(f"TNBC (n={tnbc_expr.shape[1]}, native): tau* = {tau_tnbc}")

    rng = np.random.default_rng(0)
    luma_cols = luma_expr.columns.to_numpy()
    taus = []
    for rep in range(N_REPS):
        sub_cols = rng.choice(luma_cols, size=N_TNBC, replace=False)
        C_sub = np.corrcoef(luma_expr[sub_cols].to_numpy())
        tau = find_tau_star(C_sub, THRESHOLDS)
        if tau is not None:
            taus.append(tau)
    taus = np.array(taus)
    print(f"LumA subsampled to n={N_TNBC}, {N_REPS} reps: "
          f"tau* mean={taus.mean():.3f}, std={taus.std():.3f}, "
          f"range=[{taus.min():.2f},{taus.max():.2f}]")
    print(f"\nConclusion: TNBC tau*={tau_tnbc:.2f} vs matched-n LumA "
          f"tau*={taus.mean():.2f}+/-{taus.std():.2f} "
          f"({'still separated' if abs(tau_tnbc - taus.mean()) > taus.std() else 'gap shrinks to within 1 s.d. -- mostly a sample-size artifact'}).")
