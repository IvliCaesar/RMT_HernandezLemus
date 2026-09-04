"""
Real RMT coexpression-network analysis (Luo et al. 2007 method) on the
TNBC and LumA cohorts built by 01_build_cohorts.py from real TCGA-BRCA
data. Reuses the exact unfolding / Poisson-vs-GOE / Marchenko-Pastur /
IPR machinery already validated on synthetic ground truth in
rmt_threshold_demo.py -- same functions, now pointed at real patients.

For each cohort:
  1. Build the gene-gene Pearson correlation matrix across patients.
  2. Sweep |correlation| thresholds tau, find the sparsest tau at which
     the nearest-neighbor spacing statistics turn GOE-like (Wigner-Dyson
     level repulsion) -- the RMT-selected network-defining threshold.
  3. At that tau, report network size (nodes/edges) and flag the most
     localized eigenvectors by inverse participation ratio (IPR) as
     candidate hub-gene axes.

Writes results/tnbc_rmt.csv, results/luma_rmt.csv (the tau sweep) and
prints the selected threshold, network size, and top-IPR genes for each
cohort -- all numbers below come directly from this run, none invented.
"""
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from rmt_threshold_demo import (
    unfold_spectrum,
    poisson_goe_fit_score,
    marchenko_pastur_edge,
    inverse_participation_ratio,
)

import os
os.makedirs("../results", exist_ok=True)


def analyze_cohort(name, csv_path, thresholds):
    expr = pd.read_csv(csv_path, index_col=0)
    genes = expr.index.to_numpy()
    X = expr.to_numpy()
    p, n = X.shape
    C = np.corrcoef(X)
    lam_minus, lam_plus = marchenko_pastur_edge(p, n)
    print(f"\n=== {name}: p={p} genes, n={n} patients, "
          f"q=p/n={p/n:.2f}, MP edge=[{lam_minus:.3f},{lam_plus:.3f}] ===")

    rows = []
    for tau in thresholds:
        A = C.copy()
        A[np.abs(A) < tau] = 0.0
        np.fill_diagonal(A, 0.0)
        n_edges = int((np.abs(A) > 0).sum() / 2)
        max_possible = p * (p - 1) / 2
        density = n_edges / max_possible
        eigs, eigvecs = np.linalg.eigh(A)
        if len(np.unique(np.round(eigs, 8))) < 10:
            continue
        spacings = unfold_spectrum(eigs)
        ks_p, ks_g = poisson_goe_fit_score(spacings)
        regime = "GOE" if ks_g < ks_p else "Poisson"
        rows.append(dict(tau=tau, n_edges=n_edges, density=density,
                          ks_poisson=ks_p, ks_goe=ks_g, regime=regime))
        print(f"  tau={tau:.2f}  edges={n_edges:6d}  density={density:.4f}  "
              f"ks_poisson={ks_p:.4f}  ks_goe={ks_g:.4f}  regime={regime}")

    df = pd.DataFrame(rows)
    df.to_csv(f"../results/{name}_rmt.csv", index=False)

    tau_star = None
    for i in range(len(rows) - 1, 0, -1):
        if rows[i]["regime"] == "Poisson" and rows[i - 1]["regime"] == "GOE":
            tau_star = rows[i - 1]["tau"]
            break
    if tau_star is None:
        print(f"  No clean Poisson->GOE transition found for {name}.")
        return df, None, None

    A = np.abs(C.copy())
    A[A < tau_star] = 0.0
    np.fill_diagonal(A, 0.0)
    n_edges = int((A > 0).sum() / 2)
    degree = (A > 0).sum(axis=1)
    eigs, eigvecs = np.linalg.eigh(A)
    ipr = inverse_participation_ratio(eigvecs)
    top_idx = np.argsort(ipr)[::-1][:5]
    hub_genes_by_degree = genes[np.argsort(degree)[::-1][:10]]

    print(f"  --> tau* = {tau_star:.2f}: {n_edges} edges, "
          f"{(degree>0).sum()} connected genes out of {p}")
    print(f"  --> Top-degree hub genes: {list(hub_genes_by_degree)}")
    print(f"  --> Top-5 most localized eigenvectors (IPR), "
          f"eigenvalues: {eigs[top_idx].round(3)}")
    for k, idx in enumerate(top_idx):
        loading = np.abs(eigvecs[:, idx])
        top_genes_this_vec = genes[np.argsort(loading)[::-1][:5]]
        print(f"      eigvec[{idx}] (IPR={ipr[idx]:.4f}) top genes: "
              f"{list(top_genes_this_vec)}")

    return df, tau_star, dict(n_edges=n_edges, degree=degree, genes=genes,
                               hub_genes_by_degree=hub_genes_by_degree,
                               eigs=eigs, eigvecs=eigvecs, ipr=ipr)


if __name__ == "__main__":
    thresholds = np.round(np.arange(0.15, 0.85, 0.05), 2)
    tnbc_df, tnbc_tau, tnbc_info = analyze_cohort(
        "tnbc", "../data/TCGA-BRCA/tnbc_expr.csv", thresholds)
    luma_df, luma_tau, luma_info = analyze_cohort(
        "luma", "../data/TCGA-BRCA/luma_expr.csv", thresholds)

    if tnbc_tau is not None and luma_tau is not None:
        print(f"\nRMT-selected thresholds: TNBC tau*={tnbc_tau:.2f}, "
              f"LumA tau*={luma_tau:.2f}")
        if tnbc_info is not None and luma_info is not None:
            tnbc_hubs = set(tnbc_info["hub_genes_by_degree"])
            luma_hubs = set(luma_info["hub_genes_by_degree"])
            print(f"Top-10 hub genes shared between cohorts: "
                  f"{tnbc_hubs & luma_hubs}")
            print(f"Top-10 hub genes unique to TNBC: {tnbc_hubs - luma_hubs}")
            print(f"Top-10 hub genes unique to LumA: {luma_hubs - tnbc_hubs}")
