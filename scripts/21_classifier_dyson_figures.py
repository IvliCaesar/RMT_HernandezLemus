"""
Figures for scripts 19-20 (classifier weight-matrix Dyson-repulsion
check, Sec. results-dyson-classifier): the eig(W1^T W1) trajectory
across SGD training epochs for one representative seed, and the
20-seed reproducibility bar chart, in exactly the same style as
Figures 3 and 8 for the data-side check.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

os.makedirs("../figures", exist_ok=True)

# --- Figure 9: weight-matrix eigenvalue trajectories across epochs ---
traj = np.loadtxt("../results/dyson_classifier_weights.csv", delimiter=",",
                   skiprows=1)
epochs = np.arange(traj.shape[0])
fig, ax = plt.subplots(figsize=(6, 4.2))
for k in range(traj.shape[1]):
    ax.plot(epochs, traj[:, k], lw=1.3, label=f"$\\mu_{{{k+1}}}$")
ax.set_xlabel("SGD training epoch")
ax.set_ylabel("eigenvalue of $W_1^\\top W_1$")
ax.set_title("Classifier weight-matrix eigenvalue trajectories\n"
              "(seed 0, MLP(20$\\to$8$\\to$1), RMT-hub input genes)")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig("../figures/fig9_classifier_weight_trajectories.pdf")
print("Saved figures/fig9_classifier_weight_trajectories.pdf")

# --- Figure 10: 20-seed reproducibility, same style as fig8 ---
df = pd.read_csv("../results/dyson_classifier_multi_seed_summary.csv")
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
ax.set_ylabel("count out of 20 random seeds")
ax.set_title("Classifier weight-matrix repulsion reproducibility\n"
              "across 20 independent SGD training runs")
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=2)
fig.tight_layout()
fig.savefig("../figures/fig10_classifier_repulsion_reproducibility.pdf",
            bbox_inches="tight")
print("Saved figures/fig10_classifier_repulsion_reproducibility.pdf")

# --- Figure 11: a clean early-transient avoided crossing example ---
# seed 17, pair mu5-mu6 (columns 4,5, zero-indexed), closest approach at
# epoch 5 within the pre-convergence window (script 22).
all_seeds = np.load("../results/dyson_classifier_all_seeds.npy")  # (20,300,8)
seed, epoch_lo, epoch_hi = 17, 0, 30
mu5 = all_seeds[seed, epoch_lo:epoch_hi, 4]
mu6 = all_seeds[seed, epoch_lo:epoch_hi, 5]
epochs_zoom = np.arange(epoch_lo, epoch_hi)
fig, ax = plt.subplots(figsize=(5, 3.5))
ax.plot(epochs_zoom, mu5, "b-o", ms=3, label="$\\mu_5$")
ax.plot(epochs_zoom, mu6, "r-o", ms=3, label="$\\mu_6$")
ax.axvline(epoch_lo + 5, color="k", ls=":", lw=1,
           label="closest approach (epoch 5)")
ax.set_xlabel("SGD training epoch")
ax.set_ylabel("eigenvalue of $W_1^\\top W_1$")
ax.set_title("Classifier weight-matrix avoided crossing,\n"
              "$\\mu_5$-$\\mu_6$, seed 17, early training transient")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("../figures/fig11_classifier_avoided_crossing.pdf")
print("Saved figures/fig11_classifier_avoided_crossing.pdf")
