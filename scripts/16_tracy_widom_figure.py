"""
Figure for the Tracy-Widom edge-fluctuation scaling check
(scripts/15_entropy_and_edge_scaling.py, results/tracy_widom_scaling.csv):
log-log plot of the empirical top-eigenvalue standard deviation vs.
sample size n (aspect ratio q=p/n held fixed at 3.0), against the
fitted power law and the Tracy-Widom (n^-2/3) and Gaussian/CLT (n^-1)
reference slopes -- each reference line anchored to pass through the
first empirical point, not floating freely (an earlier draft of this
figure anchored the Tracy-Widom reference line incorrectly, placing it
far above the data despite the fit being correct; fixed here).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("../results/tracy_widom_scaling.csv")
n = df["n"].to_numpy()
std = df["top_std"].to_numpy()
b, a = np.polyfit(np.log(n), np.log(std), 1)

fig, ax = plt.subplots(figsize=(5.5, 4.2))
ax.loglog(n, std, "o", color="k", label="empirical std (this study)", ms=7)
nfit = np.array([n.min(), n.max()], dtype=float)
ax.loglog(nfit, std[0] * (nfit / n[0]) ** b, "k-", lw=1,
          label=f"fit: $n^{{{b:.3f}}}$")
ax.loglog(nfit, std[0] * (nfit / n[0]) ** (-2 / 3), "g--", lw=1.5,
          label="Tracy-Widom: $n^{-2/3}$")
ax.loglog(nfit, std[0] * (nfit / n[0]) ** (-1.0), "r:", lw=1.5,
          label="Gaussian/CLT: $n^{-1}$")
ax.set_xlabel(r"sample size $n$ ($q=p/n=3.0$ fixed)")
ax.set_ylabel("std. dev. of top eigenvalue")
ax.set_title("Edge-fluctuation scaling: synthetic i.i.d. null")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("../figures/fig7_tracy_widom_scaling.pdf")
print(f"Saved figures/fig7_tracy_widom_scaling.pdf (fit exponent b={b:.3f})")
