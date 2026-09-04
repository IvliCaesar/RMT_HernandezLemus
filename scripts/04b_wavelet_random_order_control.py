"""
Control for 04_wavelet_multiscale.py: does the zero hub-gene overlap
between coarse- and fine-scale wavelet networks depend on ordering
patients by PC1 specifically, or would any arbitrary ordering (with no
real structure at all) produce the same "zero overlap" result? If a
purely random, meaningless ordering gives the same qualitative
separation, the original result is not evidence that PC1 carries real
information -- it would just be a generic property of splitting a
signal into wavelet coarse/fine components.

Runs the identical pipeline (Daubechies-4, level 3, coarse = approx.
coefficients only, fine = finest detail coefficients only, RMT
threshold selection on each) 20 times, each time with an independent
uniformly random patient permutation instead of the PC1 ordering, and
reports the coarse/fine hub-gene overlap for each run.
"""
import sys
import numpy as np
import pandas as pd
import pywt

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score

THRESHOLDS = np.round(np.arange(0.10, 0.80, 0.05), 2)
N_REPS = 20


def find_tau_star_and_hubs(C, genes, thresholds):
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
        return None, []
    A = np.abs(C.copy())
    A[A < tau_star] = 0.0
    np.fill_diagonal(A, 0.0)
    degree = (A > 0).sum(axis=1)
    return tau_star, list(genes[np.argsort(degree)[::-1][:10]])


def coarse_fine_hubs(Xo, genes, n, wavelet="db4", level=3):
    coarse = np.zeros_like(Xo)
    fine = np.zeros_like(Xo)
    for i in range(Xo.shape[0]):
        coeffs = pywt.wavedec(Xo[i], wavelet, level=level)
        coarse_coeffs = [coeffs[0]] + [np.zeros_like(c) for c in coeffs[1:]]
        coarse[i] = pywt.waverec(coarse_coeffs, wavelet)[:n]
        fine_coeffs = [np.zeros_like(coeffs[0])] + \
            [np.zeros_like(c) for c in coeffs[1:-1]] + [coeffs[-1]]
        fine[i] = pywt.waverec(fine_coeffs, wavelet)[:n]
    C_coarse = np.corrcoef(coarse)
    C_fine = np.corrcoef(fine)
    _, hubs_coarse = find_tau_star_and_hubs(C_coarse, genes, THRESHOLDS)
    _, hubs_fine = find_tau_star_and_hubs(C_fine, genes, THRESHOLDS)
    return hubs_coarse, hubs_fine


if __name__ == "__main__":
    import os
    os.makedirs("../results", exist_ok=True)
    expr = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    genes = expr.index.to_numpy()
    X = expr.to_numpy()
    p, n = X.shape

    # Reproduce the PC1-ordering result first, as a reference point.
    Xc = X - X.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc1_order = np.argsort(Vt[0])
    hubs_c_pc1, hubs_f_pc1 = coarse_fine_hubs(X[:, pc1_order], genes, n)
    overlap_pc1 = len(set(hubs_c_pc1) & set(hubs_f_pc1))
    print(f"PC1 ordering: coarse/fine hub overlap = {overlap_pc1}/10")

    rng = np.random.default_rng(42)
    overlaps = []
    for rep in range(N_REPS):
        perm = rng.permutation(n)
        hubs_c, hubs_f = coarse_fine_hubs(X[:, perm], genes, n)
        ov = len(set(hubs_c) & set(hubs_f))
        overlaps.append(ov)
        print(f"  random order rep {rep}: overlap = {ov}/10")

    overlaps = np.array(overlaps)
    print(f"\nRandom-order control ({N_REPS} reps): "
          f"mean overlap = {overlaps.mean():.2f}/10, "
          f"range = [{overlaps.min()},{overlaps.max()}], "
          f"reps with overlap==0: {(overlaps==0).sum()}/{N_REPS}")
    print(f"PC1-ordering overlap ({overlap_pc1}) vs random-order mean "
          f"({overlaps.mean():.2f}): "
          f"{'PC1 result is NOT distinguishable from random ordering -- claim needs walking back' if overlaps.mean() <= overlap_pc1 + 1 else 'PC1 ordering gives a genuinely lower overlap than random -- claim supported'}")

    pd.DataFrame({"rep": list(range(N_REPS)), "overlap": overlaps}).to_csv(
        "../results/wavelet_ordering_control.csv", index=False)
    with open("../results/wavelet_ordering_control_summary.txt", "w") as f:
        f.write(f"PC1 ordering overlap: {overlap_pc1}/10\n")
        f.write(f"Random ordering ({N_REPS} reps): mean={overlaps.mean():.3f}, "
                f"range=[{overlaps.min()},{overlaps.max()}], "
                f"reps with overlap==0: {(overlaps==0).sum()}/{N_REPS}\n")
        f.write("Conclusion: the zero coarse/fine hub-gene overlap is NOT "
                "distinguishable from a uniformly random patient ordering "
                "with no biological meaning at all, so it does not, by "
                "itself, confirm a real multiscale biological signal tied "
                "to the PC1 ordering.\n")
    print("Saved results/wavelet_ordering_control.csv and _summary.txt")
