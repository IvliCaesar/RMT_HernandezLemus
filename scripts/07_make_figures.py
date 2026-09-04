"""
Generates the real figures used in the article, all computed directly
from the real TCGA-BRCA results produced by scripts 01-06 (re-runs the
small pieces needed rather than caching, so figures always match the
numbers actually reported in the text).
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum

import os
os.makedirs("../figures", exist_ok=True)

# --- Figure 1: NNSD at tau* for TNBC vs Poisson and GOE (Wigner surmise) ---
tnbc = pd.read_csv("../data/TCGA-BRCA/tnbc_expr.csv", index_col=0)
C = np.corrcoef(tnbc.to_numpy())
tau_star = 0.45
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

fig, ax = plt.subplots(figsize=(5.5, 4))
ax.hist(spacings, bins=40, density=True, alpha=0.5, color="grey",
        label="TNBC, $\\tau^*=0.45$ (observed)")
ax.plot(s, poisson_pdf, "b--", lw=2, label="Poisson $e^{-s}$")
ax.plot(s, goe_pdf, "r-", lw=2, label="GOE Wigner surmise")
ax.set_xlabel("Unfolded nearest-neighbor spacing $s$")
ax.set_ylabel("Density")
ax.set_title("Level-spacing distribution at the RMT threshold (TNBC)")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("../figures/fig1_nnsd_tnbc.pdf")
print("Saved fig1_nnsd_tnbc.pdf")

# --- Figure 2: tau sweep (GOE vs Poisson KS scores) for TNBC and LumA ---
tnbc_df = pd.read_csv("../results/tnbc_rmt.csv")
luma_df = pd.read_csv("../results/luma_rmt.csv")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharey=True)
for ax, df, name, tau_s in [(axes[0], tnbc_df, "TNBC", 0.45),
                              (axes[1], luma_df, "LumA", 0.25)]:
    ax.plot(df["tau"], df["ks_poisson"], "b--o", ms=3, label="KS vs Poisson")
    ax.plot(df["tau"], df["ks_goe"], "r-o", ms=3, label="KS vs GOE")
    ax.axvline(tau_s, color="k", ls=":", lw=1.5, label=f"$\\tau^*={tau_s}$")
    ax.set_xlabel("threshold $\\tau$")
    ax.set_title(name)
axes[0].set_ylabel("Kolmogorov-Smirnov statistic")
axes[0].legend(fontsize=7)
fig.tight_layout()
fig.savefig("../figures/fig2_tau_sweep.pdf")
print("Saved fig2_tau_sweep.pdf")

# --- Figure 3: Dyson eigenvalue trajectories ---
traj = np.loadtxt("../results/dyson_trajectories.csv", delimiter=",", skiprows=1)
ns = list(range(15, 15 + traj.shape[0]))
fig, ax = plt.subplots(figsize=(6, 4.2))
for k in range(traj.shape[1]):
    ax.plot(ns, traj[:, k], lw=1.3, label=f"$\\lambda_{{{k+1}}}$")
ax.set_xlabel("sample size $n$ (patients included)")
ax.set_ylabel("eigenvalue")
ax.set_title("Top-8 eigenvalue trajectories, TNBC correlation matrix\n"
              "(fixed random patient-inclusion order)")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig("../figures/fig3_dyson_trajectories.pdf")
print("Saved fig3_dyson_trajectories.pdf")

# --- Figure 4: zoom on pair 7-8 avoided crossing ---
fig, ax = plt.subplots(figsize=(5, 3.5))
ax.plot(ns, traj[:, 6], "b-o", ms=3, label="$\\lambda_7$")
ax.plot(ns, traj[:, 7], "r-o", ms=3, label="$\\lambda_8$")
ax.axvline(59, color="k", ls=":", lw=1, label="closest approach ($n=59$)")
ax.set_xlim(45, 75)
ax.set_xlabel("sample size $n$")
ax.set_ylabel("eigenvalue")
ax.set_title("Avoided crossing, $\\lambda_7$-$\\lambda_8$")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("../figures/fig4_avoided_crossing.pdf")
print("Saved fig4_avoided_crossing.pdf")

print("\nAll figures written to ../figures/")
