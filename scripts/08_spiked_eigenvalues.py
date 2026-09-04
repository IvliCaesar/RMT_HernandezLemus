"""
Spiked-covariance / BBP-type diagnostic (Line 1 of the program:
spectral denoising). For each cohort's full (unthresholded) gene-gene
correlation matrix, counts how many eigenvalues exceed the
Marchenko-Pastur upper bulk edge lambda_+ = (1+sqrt(p/n))^2 -- these
are the "spiked" directions that cannot be explained by sampling noise
alone under the null of no correlation structure, and are the
components a spectral-denoising step would keep.
"""
import numpy as np
import pandas as pd

for name, path in [("TNBC", "../data/TCGA-BRCA/tnbc_expr.csv"),
                    ("LumA", "../data/TCGA-BRCA/luma_expr.csv")]:
    expr = pd.read_csv(path, index_col=0)
    X = expr.to_numpy()
    p, n = X.shape
    C = np.corrcoef(X)
    eigs = np.linalg.eigvalsh(C)[::-1]
    q = p / n
    lam_plus = (1 + np.sqrt(q)) ** 2
    n_spiked = int((eigs > lam_plus).sum())
    frac_var = eigs[eigs > lam_plus].sum() / eigs.sum()
    print(f"{name}: p={p} n={n} q={q:.2f} MP_+={lam_plus:.3f} "
          f"n_spiked={n_spiked} ({100*n_spiked/p:.1f}% of genes) "
          f"capturing {100*frac_var:.1f}% of total variance; "
          f"top-5 spiked eigenvalues={np.round(eigs[:5], 2)}")
