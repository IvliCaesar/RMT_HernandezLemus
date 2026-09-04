"""
Spiked-covariance / BBP-type diagnostic (Line 1 of the program:
spectral denoising). For each cohort's full (unthresholded) gene-gene
correlation matrix, counts how many eigenvalues exceed the
Marchenko-Pastur upper bulk edge lambda_+ = (1+sqrt(p/n))^2
\citep{marchenkopastur1967} -- candidate real, low-rank correlation
directions, as opposed to sampling noise. Also computes the inverse
participation ratio (IPR) of the corresponding eigenvectors and of the
non-spiked (bulk/noise) eigenvectors, as a check on whether "real"
components are broad multi-gene programs or dominated by a handful of
genes -- and whether IPR alone (without the Marchenko-Pastur threshold)
would even distinguish spiked from bulk eigenvectors (it does not, see
results/spiked_eigenvalues.csv and the printed comparison below).
"""
import numpy as np
import pandas as pd
import os

os.makedirs("../results", exist_ok=True)
rows = []

for name, path in [("TNBC", "../data/TCGA-BRCA/tnbc_expr.csv"),
                    ("LumA", "../data/TCGA-BRCA/luma_expr.csv")]:
    expr = pd.read_csv(path, index_col=0)
    genes = expr.index.to_numpy()
    X = expr.to_numpy()
    p, n = X.shape
    C = np.corrcoef(X)
    q = p / n
    lam_plus = (1 + np.sqrt(q)) ** 2

    eigs, eigvecs = np.linalg.eigh(C)
    order = np.argsort(eigs)[::-1]
    eigs, eigvecs = eigs[order], eigvecs[:, order]

    spiked_mask = eigs > lam_plus
    n_spiked = int(spiked_mask.sum())
    frac_var = eigs[spiked_mask].sum() / eigs.sum()

    ipr_all = np.sum(eigvecs ** 4, axis=0)
    ipr_spiked = ipr_all[spiked_mask]
    ipr_bulk = ipr_all[~spiked_mask]

    print(f"{name}: p={p} n={n} q={q:.2f} MP_+={lam_plus:.3f} "
          f"n_spiked={n_spiked} ({100*n_spiked/p:.1f}% of genes) "
          f"capturing {100*frac_var:.1f}% of total variance; "
          f"top-5 spiked eigenvalues={np.round(eigs[:5], 2)}")
    n_solid = int((eigs[spiked_mask] > 2 * lam_plus).sum())
    print(f"  spikes exceeding 2x lambda_+ (implausible as edge "
          f"fluctuation): {n_solid}/{n_spiked}")
    print(f"  IPR: spiked mean={ipr_spiked.mean():.4f}, "
          f"bulk(noise) mean={ipr_bulk.mean():.4f}, "
          f"1/p (fully delocalized)={1/p:.4f}")

    for k in range(n_spiked):
        rows.append(dict(cohort=name, rank=k + 1, eigenvalue=eigs[k],
                          lambda_plus=lam_plus, ratio=eigs[k] / lam_plus,
                          ipr=ipr_spiked[k],
                          top_genes=";".join(genes[np.argsort(
                              np.abs(eigvecs[:, k]))[::-1][:5]])))

pd.DataFrame(rows).to_csv("../results/spiked_eigenvalues.csv", index=False)
print("\nSaved results/spiked_eigenvalues.csv")
