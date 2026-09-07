"""
Figures for the single-cell Line-4 attack (scripts 23-24): the tau
sweep and NNSD-at-tau* plots, in the same style as Figures 1-2 for the
bulk cohorts, plus the spiked-eigenvalue spectrum against the
Marchenko-Pastur bulk edge.
"""
import numpy as np
import pandas as pd
from scipy.io import mmread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, marchenko_pastur_edge

os.makedirs("../figures", exist_ok=True)

DATA_DIR = "../data/GSE176078_singlecell"
N_GENES = 1500
MIN_CELL_FRACTION = 0.03

# --- Recompute the real correlation matrix (same recipe as script 23) ---
mat = mmread(f"{DATA_DIR}/count_matrix_sparse.mtx").tocsr()
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
X = dense[top_idx]
p, n = X.shape
C = np.corrcoef(X)

# --- Figure 12: tau sweep (KS distance vs threshold) ---
df = pd.read_csv("../results/singlecell_rmt_sweep.csv")
fig, ax = plt.subplots(figsize=(5.5, 4))
ax.plot(df["tau"], df["ks_poisson"], "b--o", ms=3, label="KS vs Poisson")
ax.plot(df["tau"], df["ks_goe"], "r-o", ms=3, label="KS vs GOE")
ax.axvline(0.1, color="k", ls=":", lw=1.5, label="$\\tau^*=0.10$")
ax.set_xlabel("threshold $\\tau$")
ax.set_ylabel("Kolmogorov-Smirnov statistic")
ax.set_title("RMT threshold sweep, single-cell TNBC\n(patient CID44971, 7986 cells)")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("../figures/fig12_singlecell_tau_sweep.pdf")
print("Saved figures/fig12_singlecell_tau_sweep.pdf")

# --- Figure 13: NNSD at tau*=0.10 against Poisson and GOE ---
tau_star = 0.10
A = C.copy()
A[np.abs(A) < tau_star] = 0.0
np.fill_diagonal(A, 0.0)
eigs = np.linalg.eigvalsh(A)
spacings = unfold_spectrum(eigs)
spacings = spacings[spacings > 0]
spacings = spacings / spacings.mean()

s = np.linspace(0, 4, 400)
poisson_pdf = np.exp(-s)
goe_pdf = (np.pi / 2) * s * np.exp(-np.pi * s ** 2 / 4)

tail_frac = (spacings > 4).mean()
fig, ax = plt.subplots(figsize=(5.5, 4))
ax.hist(spacings[spacings <= 4], bins=40, density=False, alpha=0.5, color="grey")
# Renormalize the clipped-range histogram to a density over [0,4] only, so
# it is visually comparable to the theoretical curves on the same axis
# (the excluded tail is disclosed in the caption, not hidden silently).
counts, edges = np.histogram(spacings[spacings <= 4], bins=40, range=(0, 4))
widths = np.diff(edges)
density = counts / (counts.sum() * widths[0])
ax.clear()
ax.bar(edges[:-1], density, width=widths, align="edge", alpha=0.5, color="grey",
       label=f"Single-cell TNBC, $\\tau^*=0.10$ ({(1-tail_frac):.0%} of spacings shown)")
ax.plot(s, poisson_pdf, "b--", lw=2, label="Poisson $e^{-s}$")
ax.plot(s, goe_pdf, "r-", lw=2, label="GOE Wigner surmise")
ax.set_xlim(0, 4)
ax.set_xlabel("Unfolded nearest-neighbor spacing $s$")
ax.set_ylabel("Density")
ax.set_title("Level-spacing distribution at the RMT threshold\n(single-cell TNBC)")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig("../figures/fig13_singlecell_nnsd.pdf")
print(f"Saved figures/fig13_singlecell_nnsd.pdf (tail fraction s>4: {tail_frac:.1%})")

# --- Figure 14: raw spectrum vs Marchenko-Pastur bulk edge ---
eigs_raw = np.linalg.eigvalsh(C)
lam_minus, lam_plus = marchenko_pastur_edge(p, n)
fig, ax = plt.subplots(figsize=(5.5, 4))
ax.hist(eigs_raw, bins=80, color="steelblue", alpha=0.7)
ax.axvline(lam_plus, color="red", ls="--", lw=1.5,
           label=f"$\\lambda_+={lam_plus:.2f}$ (MP edge)")
ax.set_xlim(-1, 15)
ax.set_xlabel("eigenvalue")
ax.set_ylabel("count")
ax.set_title("Raw correlation-matrix spectrum, single-cell TNBC\n"
              "(33/1500 eigenvalues exceed $\\lambda_+$, not all shown)")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("../figures/fig14_singlecell_spiked_spectrum.pdf")
print("Saved figures/fig14_singlecell_spiked_spectrum.pdf")
