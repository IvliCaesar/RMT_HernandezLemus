"""
Figures for the Hi-C compartment validation (script 32) and the
four-sample extension (script 33): (a) PC1 compartment score vs. real
GC content, TNBC_Tissue3, colored by compartment sign; (b) pairwise
compartment-sign agreement across all four real GSE167150 samples
(three TNBC tumors, one matched normal).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

val = pd.read_csv("../results/hic_chr18_compartment_validation.csv")
agree = pd.read_csv("../results/hic_multi_sample_agreement.csv")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

ax0 = axes[0]
colors = np.where(val["pc1_compartment_score"] > 0, "#d62728", "#1f77b4")
ax0.scatter(val["gc_content"], val["pc1_compartment_score"], c=colors, s=8, alpha=0.6)
ax0.axhline(0, color="k", linewidth=0.5)
r = np.corrcoef(val["gc_content"], val["pc1_compartment_score"])[0, 1]
ax0.set_xlabel("GC content (hg19 chr18, per 100kb bin)")
ax0.set_ylabel("PC1 compartment score")
ax0.set_title(f"TNBC_Tissue3: PC1 vs. real GC content\n"
              f"(r={r:.2f}, independent validation)")

names = ["TNBC_Tissue1", "TNBC_Tissue2", "TNBC_Tissue3", "Normal_Tissue"]
mat = np.full((4, 4), np.nan)
for _, row in agree.iterrows():
    i, j = names.index(row["sample_a"]), names.index(row["sample_b"])
    mat[i, j] = row["agreement"]
    mat[j, i] = row["agreement"]
np.fill_diagonal(mat, 1.0)

ax1 = axes[1]
im = ax1.imshow(mat, cmap="viridis", vmin=0.5, vmax=1.0)
ax1.set_xticks(range(4))
ax1.set_yticks(range(4))
short = ["TNBC1", "TNBC2", "TNBC3", "Normal"]
ax1.set_xticklabels(short, rotation=30, ha="right")
ax1.set_yticklabels(short)
for i in range(4):
    for j in range(4):
        ax1.text(j, i, f"{mat[i,j]:.0%}", ha="center", va="center",
                  color="white" if mat[i, j] < 0.85 else "black", fontsize=9)
ax1.set_title("Compartment-sign agreement,\nchr18, all four real samples")
fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)

fig.tight_layout()
fig.savefig("../figures/fig19_hic_validation_multisample.pdf", bbox_inches="tight")
print("Saved figures/fig19_hic_validation_multisample.pdf")
