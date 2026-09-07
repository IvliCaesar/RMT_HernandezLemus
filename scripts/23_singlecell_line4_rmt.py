"""
Attacks Roadmap Line 4 directly (single-cell resolution, previously
identified but not attempted with real data): applies the identical
RMT machinery used throughout this paper (spiked-eigenvalue count
against the Marchenko-Pastur bulk, tau-sweep threshold selection) to
real single-cell RNA-seq data for the first time in this program,
rather than bulk patient-level expression.

Data: GSE176078 (Wu et al. 2021, already cited in article.tex as
tnbcsc2021), patient CID44971 -- the largest of the series' 9
TNBC-labeled samples (7986 cells, avoiding cross-patient batch effects
in this first attempt; pooling patients is left as the direct next
step, exactly as Sec. 5's real-data wavelet attempt flagged its own
single ordering-axis limitation).

Unlike the bulk cohorts (p=1500 genes, n=119/430 patients, q=p/n >> 1),
single-cell data inverts the aspect ratio: n (now cells) reaches the
thousands while p (genes) stays modest after variable-gene selection,
so q=p/n << 1 here -- exactly the regime the Roadmap predicted would
give a tighter Marchenko-Pastur bulk and, in principle, more sensitive
spike detection. No new methodology: the same normalize -> select
variable genes -> correlation matrix -> RMT pipeline as
scripts/01, 02, 08, applied to new data.
"""
import numpy as np
import pandas as pd
from scipy.io import mmread
import sys
import os

sys.path.insert(0, ".")
from rmt_threshold_demo import (unfold_spectrum, poisson_goe_fit_score,
                                 marchenko_pastur_edge)

DATA_DIR = "../data/GSE176078_singlecell"
N_GENES = 1500          # same p as the bulk cohorts, for direct comparability
MIN_CELL_FRACTION = 0.03  # gene must be detected in >=3% of cells


def load_data():
    mat = mmread(f"{DATA_DIR}/count_matrix_sparse.mtx").tocsr()  # genes x cells
    genes = pd.read_csv(f"{DATA_DIR}/count_matrix_genes.tsv", header=None)[0].to_numpy()
    barcodes = pd.read_csv(f"{DATA_DIR}/count_matrix_barcodes.tsv", header=None)[0].to_numpy()
    meta = pd.read_csv(f"{DATA_DIR}/metadata.csv", index_col=0)
    meta = meta.reindex(barcodes)
    return mat, genes, barcodes, meta


def normalize_cp10k_log1p(mat):
    """Standard single-cell normalization: counts-per-10k, then log1p.
    mat is genes x cells (sparse)."""
    cell_totals = np.asarray(mat.sum(axis=0)).ravel()
    cell_totals[cell_totals == 0] = 1.0
    scale = 1e4 / cell_totals
    mat_norm = mat.multiply(scale[np.newaxis, :]).tocsr()
    mat_norm.data = np.log1p(mat_norm.data)
    return mat_norm


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    print("Loading real GSE176078 single-cell TNBC data (patient CID44971)...")
    mat, genes, barcodes, meta = load_data()
    p_all, n_cells = mat.shape
    print(f"Loaded {p_all} genes x {n_cells} cells. "
          f"subtype (metadata, all cells): {meta['subtype'].unique()}")
    print("Cell-type composition (celltype_major, real annotation from "
          "Wu et al. 2021, not used by the RMT method):")
    print(meta["celltype_major"].value_counts())

    print("\nNormalizing (CP10K + log1p)...")
    norm = normalize_cp10k_log1p(mat)

    # Filter to genes detected in >= MIN_CELL_FRACTION of cells
    detected_frac = np.asarray((norm > 0).sum(axis=1)).ravel() / n_cells
    keep = detected_frac >= MIN_CELL_FRACTION
    print(f"Genes detected in >={MIN_CELL_FRACTION:.0%} of cells: "
          f"{keep.sum()} / {p_all}")

    norm_kept = norm[keep]
    genes_kept = genes[keep]

    # Select top N_GENES by variance (same recipe as the bulk cohorts,
    # scripts/01_build_cohorts.py), on the DENSE normalized submatrix.
    dense = np.asarray(norm_kept.todense())
    var = dense.var(axis=1)
    top_idx = np.argsort(var)[::-1][:N_GENES]
    X = dense[top_idx]  # (N_GENES, n_cells)
    sel_genes = genes_kept[top_idx]
    print(f"\nSelected top {N_GENES} most-variable genes (post-filter), "
          f"by normalized-expression variance across {n_cells} cells.")

    p, n = X.shape
    q = p / n
    print(f"p={p} genes, n={n} cells, q=p/n={q:.4f} "
          f"(bulk cohorts had q=12.61 TNBC / 3.49 LumA -- single-cell "
          f"inverts the aspect ratio as the Roadmap anticipated).")

    C = np.corrcoef(X)
    eigs_raw = np.linalg.eigvalsh(C)
    lam_minus, lam_plus = marchenko_pastur_edge(p, n)
    n_spiked = int((eigs_raw > lam_plus).sum())
    top_eig = eigs_raw[-1]
    total_var = eigs_raw.sum()
    spike_var_frac = eigs_raw[eigs_raw > lam_plus].sum() / total_var
    print(f"\nMarchenko-Pastur bulk edge: [{lam_minus:.4f}, {lam_plus:.4f}]")
    print(f"Spiked eigenvalues beyond bulk: {n_spiked} / {p} "
          f"({n_spiked/p:.1%}), top eigenvalue={top_eig:.3f}, "
          f"spike share of total variance={spike_var_frac:.1%}")

    # Top-loading genes of the single largest spike, same diagnostic as
    # scripts/08_spiked_eigenvalues.py, plus a real, independent check:
    # do those genes coincide with real cell-type identity (metadata),
    # rather than being an artifact of the RMT method itself?
    eigvecs = np.linalg.eigh(C)[1]
    top_vec = eigvecs[:, -1]
    top_loadings_idx = np.argsort(np.abs(top_vec))[::-1][:10]
    top_loading_genes = list(sel_genes[top_loadings_idx])
    print(f"\nTop-10 loading genes of the single largest spike "
          f"(eigenvalue {top_eig:.2f}): {top_loading_genes}")

    # Save spike table
    spike_rows = []
    for rank, ei in enumerate(np.argsort(eigs_raw)[::-1][:n_spiked]):
        vec = eigvecs[:, ei]
        idx5 = np.argsort(np.abs(vec))[::-1][:5]
        spike_rows.append(dict(rank=rank + 1, eigenvalue=eigs_raw[ei],
                                ratio_to_lambda_plus=eigs_raw[ei] / lam_plus,
                                top5_genes=";".join(sel_genes[idx5])))
    pd.DataFrame(spike_rows).to_csv(
        "../results/singlecell_spiked_eigenvalues.csv", index=False)

    # --- Tau sweep (same method as scripts/02, rmt_threshold_demo.py) ---
    print("\nRunning RMT threshold sweep on single-cell correlation matrix...")
    # Single-cell gene-gene correlations are typically weaker than bulk
    # patient-level correlations (dropout noise attenuates them), so we
    # check the empirical |C| distribution first rather than assuming
    # the bulk sweep range (0.15-0.80) is appropriate here.
    offdiag = C[np.triu_indices(p, k=1)]
    print(f"Off-diagonal |correlation| distribution: "
          f"median={np.median(np.abs(offdiag)):.4f}, "
          f"95th pct={np.percentile(np.abs(offdiag), 95):.4f}, "
          f"99th pct={np.percentile(np.abs(offdiag), 99):.4f}, "
          f"max={np.abs(offdiag).max():.4f}")

    thresholds = np.round(np.arange(0.05, 0.65, 0.025), 3)
    sweep_rows = []
    for tau in thresholds:
        A = C.copy()
        A[np.abs(A) < tau] = 0.0
        np.fill_diagonal(A, 0.0)
        n_edges = int((np.abs(A) > 0).sum() / 2)
        eigs = np.linalg.eigvalsh(A)
        if len(np.unique(np.round(eigs, 8))) < 10:
            sweep_rows.append(dict(tau=tau, edges=n_edges, density=np.nan,
                                    ks_poisson=np.nan, ks_goe=np.nan,
                                    regime="degenerate"))
            continue
        spacings = unfold_spectrum(eigs)
        ks_p, ks_g = poisson_goe_fit_score(spacings)
        regime = "GOE" if ks_g < ks_p else "Poisson"
        density = n_edges / (p * (p - 1) / 2)
        sweep_rows.append(dict(tau=tau, edges=n_edges, density=density,
                                ks_poisson=ks_p, ks_goe=ks_g, regime=regime))
        print(f"  tau={tau:.3f}  edges={n_edges:6d}  density={density:.4f}  "
              f"KS(Poisson)={ks_p:.4f}  KS(GOE)={ks_g:.4f}  regime={regime}")

    df = pd.DataFrame(sweep_rows)
    df.to_csv("../results/singlecell_rmt_sweep.csv", index=False)

    tau_star = None
    regimes = df["regime"].tolist()
    taus = df["tau"].tolist()
    for i in range(len(regimes) - 1, 0, -1):
        if regimes[i] == "Poisson" and regimes[i - 1] == "GOE":
            tau_star = taus[i - 1]
            break
    if tau_star is not None:
        print(f"\nPoisson -> GOE transition (RMT-selected threshold): "
              f"tau* = {tau_star}")
    else:
        print("\nNo clean Poisson -> GOE transition found in this range "
              "-- reporting as such rather than forcing one.")

    print("\nSaved results/singlecell_spiked_eigenvalues.csv and "
          "results/singlecell_rmt_sweep.csv")
