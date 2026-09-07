"""
Figures for the three new real results added in this round:
  fig15: LumA Dyson-repulsion reproducibility (script 27), same style
         as fig8 (TNBC) for direct visual comparison.
  fig16: single-cell Cancer-Epithelial pseudo-time Dyson trajectories,
         real order vs. the 20-random-order control (script 26) --
         shows the repulsion pattern is a property of growing sample
         size, not of this specific ordering.
  fig17: classifier weight-matrix repulsion, early vs. late window, at
         5x lower learning rate / 5x more epochs (script 28) -- direct
         test of whether a longer transient produces a longer
         detectable-repulsion window.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- fig15: LumA repulsion reproducibility (same style as fig8) ---
df = pd.read_csv("../results/luma_dyson_repulsion_multi_order_summary.csv")
pairs = df["pair"]
fig, ax = plt.subplots(figsize=(7, 4))
x = np.arange(len(pairs))
width = 0.6
bottoms = np.zeros(len(pairs))
colors = {"repulsion": "#2ca02c", "partial": "#98df8a",
          "inconclusive": "#ff7f0e", "none": "#d62728"}
labels = {"repulsion": "clear repulsion (both sides reopen)",
          "partial": "partial repulsion (one side)",
          "inconclusive": "inconclusive (still narrowing)",
          "none": "no clear repulsion"}
for key in ["repulsion", "partial", "inconclusive", "none"]:
    vals = df[key].to_numpy()
    ax.bar(x, vals, width, bottom=bottoms, color=colors[key], label=labels[key])
    bottoms += vals
ax.set_xticks(x)
ax.set_xticklabels(pairs, rotation=30, ha="right")
ax.set_ylabel("count out of 20 random orders")
ax.set_title("Dyson repulsion reproducibility across 20 random\n"
              "patient-inclusion orders (LumA, top 60 genes, n=430)")
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=2)
fig.tight_layout()
fig.savefig("../figures/fig15_luma_repulsion_reproducibility.pdf", bbox_inches="tight")
print("Saved figures/fig15_luma_repulsion_reproducibility.pdf")

# --- fig16: single-cell epithelial pseudo-time, real vs random-order ---
summ = pd.read_csv("../results/singlecell_epithelial_dyson_summary.csv")
pairs2 = summ["pair"]
fig, ax = plt.subplots(figsize=(7, 4))
x2 = np.arange(len(pairs2))
w = 0.35
real_repulsive = (summ["real_pattern"] == "repulsion").astype(int) * 20
ax.bar(x2 - w/2, real_repulsive, w, label="real pseudo-time\n(20/20 if repulsion, 0 else)",
       color="#1f77b4")
ax.bar(x2 + w/2, summ["random_repulsion"], w, label="random cell order\n(count out of 20)",
       color="#ff7f0e")
ax.set_xticks(x2)
ax.set_xticklabels(pairs2, rotation=30, ha="right")
ax.set_ylabel("clear-repulsion count (out of 20 random orders)")
ax.set_title("Single-cell Cancer-Epithelial pseudo-time: real ordering\n"
              "vs. 20 random cell orderings, same criterion as Fig. repro")
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig("../figures/fig16_singlecell_epithelial_real_vs_random.pdf", bbox_inches="tight")
print("Saved figures/fig16_singlecell_epithelial_real_vs_random.pdf")

# --- fig17: classifier longer-transient early vs late ---
el = pd.read_csv("../results/dyson_classifier_longer_transient_early_late.csv")
conv = pd.read_csv("../results/dyson_classifier_longer_transient_convergence.csv")


def signal(df_, col):
    return df_[col].isin(["repulsion", "partial"]).astype(int)


early_by_pair = el.groupby("pair").apply(lambda d: signal(d, "early_pattern").sum())
late_by_pair = el.groupby("pair").apply(lambda d: signal(d, "late_pattern").sum())
pairs3 = early_by_pair.index

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
ax = axes[0]
x3 = np.arange(len(pairs3))
w = 0.35
ax.bar(x3 - w/2, early_by_pair.to_numpy(), w, label="early (pre-own-convergence)",
       color="#d62728")
ax.bar(x3 + w/2, late_by_pair.to_numpy(), w, label="late (post-own-convergence)",
       color="#1f77b4")
ax.set_xticks(x3)
ax.set_xticklabels(pairs3, rotation=30, ha="right")
ax.set_ylabel("repulsion+partial count (out of 20 seeds)")
ax.set_title("Longer-transient run (lr=0.01, 1500 epochs):\nearly vs. late, per-seed convergence split")
ax.legend(fontsize=7)

ax2 = axes[1]
ax2.hist(conv["convergence_epoch"], bins=15, color="#9467bd", edgecolor="k")
ax2.set_xlabel("empirical convergence epoch")
ax2.set_ylabel("seeds")
ax2.set_title(f"Empirical convergence epoch\n(mean={conv['convergence_epoch'].mean():.0f}, "
              f"vs. ~50 at the original lr=0.05)")
fig.tight_layout()
fig.savefig("../figures/fig17_classifier_longer_transient.pdf", bbox_inches="tight")
print("Saved figures/fig17_classifier_longer_transient.pdf")
