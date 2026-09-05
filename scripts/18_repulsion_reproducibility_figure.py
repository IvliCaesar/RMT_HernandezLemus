"""
Figure for the multi-order Dyson repulsion reproducibility check
(scripts/17_dyson_repulsion_multi_order.py,
results/dyson_repulsion_multi_order_summary.csv): stacked bar chart of
how each adjacent eigenvalue pair's gap behaves (clear repulsion,
partial, inconclusive, or none) across 20 independent random
patient-inclusion orders.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("../results/dyson_repulsion_multi_order_summary.csv")
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
              "patient-inclusion orders (TNBC, top 60 genes)")
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=2)
fig.tight_layout()
fig.savefig("../figures/fig8_repulsion_reproducibility.pdf", bbox_inches="tight")
print("Saved figures/fig8_repulsion_reproducibility.pdf")
