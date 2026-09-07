"""
The same Monte Carlo permutation null used for the bulk cohorts
(scripts/10_monte_carlo_null.py, Sec. montecarlo in article.tex),
applied to the single-cell result of script 23: is tau*=0.10 (the
real, sparsest-GOE threshold found on the real CID44971 single-cell
correlation matrix) actually below what pure shuffled noise can
sustain a GOE-like appearance at, or could it be a dense-matrix
universality artifact? Each of the 1500 selected genes' normalized
expression values is independently permuted across the 7986 cells, 20
times, destroying every real gene-gene correlation while keeping each
gene's own marginal distribution fixed; the identical threshold sweep
is rerun on each permuted correlation matrix.
"""
import numpy as np
import pandas as pd
from scipy.io import mmread
import sys
import os

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score

DATA_DIR = "../data/GSE176078_singlecell"
N_GENES = 1500
MIN_CELL_FRACTION = 0.03
N_PERMUTATIONS = 20
# Note: a first attempt swept tau in [0.05, 0.325] and found every single
# permutation "None" (no transition) -- not because there is no
# transition, but because the permuted null's correlations are so weak
# (SE ~ 1/sqrt(n) = 1/sqrt(7986) = 0.011 for n=7986 cells) that its
# GOE->Poisson transition sits below 0.05 and the matrix is already
# fully disconnected (degenerate, zero off-diagonal entries) by tau=0.05.
# A diagnostic single-permutation run located the real transition
# between tau=0.02 (GOE) and tau=0.04 (Poisson), so the sweep here is
# shifted and refined to actually cover it.
THRESHOLDS = np.round(np.arange(0.01, 0.09, 0.005), 3)


def load_and_select():
    mat = mmread(f"{DATA_DIR}/count_matrix_sparse.mtx").tocsr()
    genes = pd.read_csv(f"{DATA_DIR}/count_matrix_genes.tsv", header=None)[0].to_numpy()
    cell_totals = np.asarray(mat.sum(axis=0)).ravel()
    cell_totals[cell_totals == 0] = 1.0
    scale = 1e4 / cell_totals
    norm = mat.multiply(scale[np.newaxis, :]).tocsr()
    norm.data = np.log1p(norm.data)
    n_cells = mat.shape[1]
    detected_frac = np.asarray((norm > 0).sum(axis=1)).ravel() / n_cells
    keep = detected_frac >= MIN_CELL_FRACTION
    dense = np.asarray(norm[keep].todense())
    var = dense.var(axis=1)
    top_idx = np.argsort(var)[::-1][:N_GENES]
    return dense[top_idx]  # (N_GENES, n_cells)


def sweep_transition(C, thresholds):
    p = C.shape[0]
    regimes = []
    for tau in thresholds:
        A = C.copy()
        A[np.abs(A) < tau] = 0.0
        np.fill_diagonal(A, 0.0)
        eigs = np.linalg.eigvalsh(A)
        if len(np.unique(np.round(eigs, 8))) < 10:
            regimes.append(None)
            continue
        spacings = unfold_spectrum(eigs)
        ks_p, ks_g = poisson_goe_fit_score(spacings)
        regimes.append("GOE" if ks_g < ks_p else "Poisson")
    tau_star = None
    for i in range(len(regimes) - 1, 0, -1):
        if regimes[i] == "Poisson" and regimes[i - 1] == "GOE":
            tau_star = thresholds[i - 1]
            break
    return tau_star


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    print("Loading and selecting the same 1500 genes as script 23...")
    X = load_and_select()
    p, n = X.shape
    print(f"p={p}, n={n}")

    print("Real (unpermuted) threshold sweep, for reference (script 23's own "
          "wider range, 0.05-0.625, since the real matrix's transition sits "
          "well above the permuted null's -- see note on THRESHOLDS above):")
    REAL_THRESHOLDS = np.round(np.arange(0.05, 0.65, 0.025), 3)
    C_real = np.corrcoef(X)
    tau_real = sweep_transition(C_real, REAL_THRESHOLDS)
    print(f"  Real tau* = {tau_real}")

    print(f"\nRunning {N_PERMUTATIONS} independent gene-wise permutations...")
    perm_thresholds = []
    for seed in range(N_PERMUTATIONS):
        rng = np.random.default_rng(seed)
        Xp = X.copy()
        for i in range(p):
            Xp[i] = rng.permutation(Xp[i])
        Cp = np.corrcoef(Xp)
        tau_p = sweep_transition(Cp, THRESHOLDS)
        perm_thresholds.append(tau_p)
        print(f"  seed={seed}: permuted tau* = {tau_p}")

    valid = [t for t in perm_thresholds if t is not None]
    df = pd.DataFrame(dict(seed=range(N_PERMUTATIONS), tau_star=perm_thresholds))
    df.to_csv("../results/singlecell_monte_carlo_null.csv", index=False)

    if valid:
        print(f"\nPermutation tau* : mean={np.mean(valid):.4f}, "
              f"std={np.std(valid):.4f}, min={min(valid)}, max={max(valid)}, "
              f"n_valid={len(valid)}/{N_PERMUTATIONS}")
        print(f"Real tau*={tau_real} vs. permutation floor: "
              f"{'ABOVE (real structure confirmed)' if tau_real > max(valid) else 'NOT clearly above permutation floor'}")
    print("\nSaved results/singlecell_monte_carlo_null.csv")
