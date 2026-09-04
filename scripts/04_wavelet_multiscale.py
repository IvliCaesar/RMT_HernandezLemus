"""
Real multiscale (wavelet) coexpression analysis on the TNBC cohort,
following the scale-decomposition idea of Koenig et al. (BMC
Bioinformatics 2006) -- who apply a discrete wavelet transform to
ordered expression profiles and build separate coexpression networks
per scale -- adapted here since TCGA-BRCA patients have no natural time
axis: patients are ordered along PC1 of their own (real) expression
matrix (the dominant axis of variation across this real cohort, a
standard "pseudo-time" ordering also used in single-cell analysis) and
each gene's expression trace along that ordering is decomposed with a
Daubechies-4 discrete wavelet transform (pywt).

For each gene we reconstruct:
  - a COARSE signal (approximation coefficients only, i.e. low-frequency
    trend along the PC1 ordering)
  - a FINE signal (detail coefficients at the finest scale only, i.e.
    high-frequency local fluctuation)
and build a gene-gene correlation network from each, independently
running the same RMT threshold-selection method as before. Comparing
the two networks' hub genes tests whether coarse-scale (slow, cohort-
wide trend) and fine-scale (local, patient-to-patient) coexpression
structure involve different genes -- the multiscale claim from the
proposal's Linea 5, now checked against real data rather than asserted.
"""
import sys
import numpy as np
import pandas as pd
import pywt

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score

THRESHOLDS = np.round(np.arange(0.10, 0.80, 0.05), 2)


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
    hubs = list(genes[np.argsort(degree)[::-1][:10]])
    return tau_star, hubs


if __name__ == "__main__":
    expr = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
    genes = expr.index.to_numpy()
    X = expr.to_numpy()  # genes x patients
    p, n = X.shape
    print(f"TNBC cohort: p={p} genes, n={n} patients")

    # Order patients along PC1 (dominant real axis of variation).
    Xc = X - X.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc1 = Vt[0]
    order = np.argsort(pc1)
    Xo = X[:, order]
    print(f"PC1 explains {100*S[0]**2/np.sum(S**2):.1f}% of variance "
          f"(ordering axis for the wavelet decomposition).")

    coarse = np.zeros_like(Xo)
    fine = np.zeros_like(Xo)
    wavelet = "db4"
    level = pywt.dwt_max_level(n, pywt.Wavelet(wavelet).dec_len)
    level = max(1, min(level, 3))
    print(f"Wavelet: {wavelet}, decomposition level: {level}")

    for i in range(p):
        coeffs = pywt.wavedec(Xo[i], wavelet, level=level)
        # coarse: keep only approximation (cA), zero all details
        coarse_coeffs = [coeffs[0]] + [np.zeros_like(c) for c in coeffs[1:]]
        coarse[i] = pywt.waverec(coarse_coeffs, wavelet)[:n]
        # fine: keep only the finest-scale detail (cD1), zero the rest
        fine_coeffs = [np.zeros_like(coeffs[0])] + \
            [np.zeros_like(c) for c in coeffs[1:-1]] + [coeffs[-1]]
        fine[i] = pywt.waverec(fine_coeffs, wavelet)[:n]

    C_coarse = np.corrcoef(coarse)
    C_fine = np.corrcoef(fine)
    C_raw = np.corrcoef(Xo)

    tau_raw, hubs_raw = find_tau_star_and_hubs(C_raw, genes, THRESHOLDS)
    tau_coarse, hubs_coarse = find_tau_star_and_hubs(C_coarse, genes, THRESHOLDS)
    tau_fine, hubs_fine = find_tau_star_and_hubs(C_fine, genes, THRESHOLDS)

    print(f"\nRaw (unfiltered) network:    tau*={tau_raw}, hubs={hubs_raw}")
    print(f"Coarse-scale network:        tau*={tau_coarse}, hubs={hubs_coarse}")
    print(f"Fine-scale network:          tau*={tau_fine}, hubs={hubs_fine}")

    if hubs_coarse and hubs_fine:
        shared = set(hubs_coarse) & set(hubs_fine)
        print(f"\nHub genes shared between coarse and fine scale: {shared}")
        print(f"Hub genes unique to coarse scale: {set(hubs_coarse)-set(hubs_fine)}")
        print(f"Hub genes unique to fine scale: {set(hubs_fine)-set(hubs_coarse)}")
