"""
Figure for the real Hi-C RMT result (script 30): three panels --
(a) the O/E-normalized bin-bin correlation heatmap for chr18, showing
the plaid/checkerboard pattern characteristic of A/B compartments;
(b) the RMT threshold sweep (KS distance to Poisson/GOE vs tau), same
style as fig2/fig12; (c) the top-eigenvector compartment track along
the chromosome, sign-colored.
"""
import numpy as np
import pandas as pd
import straw
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, ".")
from importlib import import_module
hic_mod = import_module("30_hic_line3_rmt")

HIC_FILE = "../data/GSE167150_HiC/GSM5098082_TNBC_Tissue3.hic"
CHROM = "18"
RESOLUTION = 100000

M = hic_mod.load_contact_matrix(HIC_FILE, CHROM, RESOLUTION, norm="NONE")
bin_starts_raw = np.arange(M.shape[0]) * RESOLUTION
M, bin_starts, _ = hic_mod.filter_zero_coverage_bins(M, bin_starts_raw)
OE = hic_mod.oe_normalize(M)
C = np.corrcoef(OE)
C = np.nan_to_num(C, nan=0.0)

sweep = pd.read_csv("../results/hic_chr18_rmt_sweep.csv")
compartment = pd.read_csv("../results/hic_chr18_compartment_pc1.csv")

fig = plt.figure(figsize=(12, 4.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.3])

ax0 = fig.add_subplot(gs[0])
im = ax0.imshow(C, cmap="RdBu_r", vmin=-0.5, vmax=0.5)
ax0.set_title("chr18 bin-bin correlation\n(O/E-normalized, 100kb, TNBC)")
ax0.set_xlabel("bin index")
ax0.set_ylabel("bin index")
fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)

ax1 = fig.add_subplot(gs[1])
valid = sweep.dropna(subset=["ks_poisson"])
ax1.plot(valid["tau"], valid["ks_poisson"], "o-", label="KS(Poisson)", color="#ff7f0e")
ax1.plot(valid["tau"], valid["ks_goe"], "o-", label="KS(GOE)", color="#1f77b4")
ax1.axvline(0.15, color="k", linestyle="--", linewidth=1, label=r"$\tau^*=0.15$")
ax1.set_xlabel(r"threshold $\tau$")
ax1.set_ylabel("KS distance")
ax1.set_title("RMT threshold sweep, chr18")
ax1.legend(fontsize=7)

ax2 = fig.add_subplot(gs[2])
colors = np.where(compartment["pc1_compartment_score"] > 0, "#d62728", "#1f77b4")
ax2.bar(compartment["bin_start"] / 1e6, compartment["pc1_compartment_score"],
        width=0.1, color=colors)
ax2.axhline(0, color="k", linewidth=0.5)
ax2.set_xlabel("chr18 position (Mb)")
ax2.set_ylabel("PC1 (compartment score)")
ax2.set_title("Compartment PC1 track")

fig.tight_layout()
fig.savefig("../figures/fig18_hic_compartment_rmt.pdf", bbox_inches="tight")
print("Saved figures/fig18_hic_compartment_rmt.pdf")
