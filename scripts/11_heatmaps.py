"""
Real correlation-matrix heatmaps for TNBC and LumA: raw (unthresholded)
correlation matrix and the network adjacency matrix at each cohort's
own RMT-selected tau*, genes ordered by hierarchical clustering (so
block/module structure, if any, is visible) rather than left in
arbitrary variance-rank order.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

import os
os.makedirs("../figures", exist_ok=True)


def clustered_order(C):
    d = 1 - np.abs(C)
    np.fill_diagonal(d, 0.0)
    d = (d + d.T) / 2
    condensed = squareform(d, checks=False)
    Z = linkage(condensed, method="average")
    return leaves_list(Z)


fig, axes = plt.subplots(2, 2, figsize=(9, 9))

for col, (name, path, tau_star) in enumerate([
        ("TNBC", "../data/TCGA-BRCA/tnbc_expr.csv", 0.45),
        ("LumA", "../data/TCGA-BRCA/luma_expr.csv", 0.25)]):
    expr = pd.read_csv(path, index_col=0)
    X = expr.to_numpy()
    C = np.corrcoef(X)
    order = clustered_order(C)
    C_ord = C[np.ix_(order, order)]

    A = C_ord.copy()
    A[np.abs(A) < tau_star] = 0.0
    np.fill_diagonal(A, 1.0)

    ax = axes[0, col]
    im = ax.imshow(C_ord, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_title(f"{name}: raw correlation\n(hierarchically clustered order)")
    ax.set_xticks([]); ax.set_yticks([])

    ax2 = axes[1, col]
    im2 = ax2.imshow(A, cmap="RdBu_r", vmin=-1, vmax=1)
    ax2.set_title(f"{name}: thresholded at $\\tau^*={tau_star}$")
    ax2.set_xticks([]); ax2.set_yticks([])

fig.colorbar(im, ax=axes[0, :], shrink=0.7, label="Pearson correlation")
fig.colorbar(im2, ax=axes[1, :], shrink=0.7, label="Pearson correlation")
fig.savefig("../figures/fig6_correlation_heatmaps.pdf", bbox_inches="tight")
print("Saved figures/fig6_correlation_heatmaps.pdf")
