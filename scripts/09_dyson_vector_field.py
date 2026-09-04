"""
Real illustration of the Dyson Brownian motion drift term stated in
article.tex's Definition (Dyson Brownian motion): for two eigenvalues
lambda_i, lambda_j, the deterministic (non-noise) part of the SDE is

    d(lambda_i)/dt =  1/(lambda_i - lambda_j)
    d(lambda_j)/dt = -1/(lambda_i - lambda_j)

This is an exact vector field on the (lambda_i, lambda_j) plane, not a
decorative diagram: it is the literal drift term from the SDE already
written in the paper, restricted to a pair of eigenvalues. Plotting it
shows the repulsion mechanism directly -- arrows point away from the
diagonal lambda_i=lambda_j everywhere, diverging as the two approach
it, which is the mechanical reason two eigenvalues never cross under
this dynamics (matches the avoided crossings in Fig. 3/4, now for a
real gene-correlation matrix rather than this idealized 2-eigenvalue
toy system).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

lo, hi = -3, 3
n = 25
li = np.linspace(lo, hi, n)
lj = np.linspace(lo, hi, n)
LI, LJ = np.meshgrid(li, lj)

diff = LI - LJ
with np.errstate(divide="ignore", invalid="ignore"):
    dLI = 1.0 / diff
    dLJ = -1.0 / diff

# Mask the singular diagonal (excluded region: eigenvalues never sit
# exactly on top of each other under this dynamics).
mask = np.abs(diff) < 0.15
dLI[mask] = np.nan
dLJ[mask] = np.nan

# Normalize arrow length for display (direction matters, not magnitude,
# since the drift diverges at the diagonal); color encodes true
# magnitude on a log scale.
mag = np.sqrt(dLI ** 2 + dLJ ** 2)
dLI_n = dLI / mag
dLJ_n = dLJ / mag

fig, ax = plt.subplots(figsize=(5.5, 5))
strm = ax.quiver(LI, LJ, dLI_n, dLJ_n, np.log10(mag), cmap="viridis",
                  pivot="mid", scale=22)
ax.plot([lo, hi], [lo, hi], "r--", lw=1.5, label=r"$\lambda_i=\lambda_j$ (excluded)")
cbar = fig.colorbar(strm, ax=ax)
cbar.set_label(r"$\log_{10}|\mathrm{drift}|$")
ax.set_xlabel(r"$\lambda_i$")
ax.set_ylabel(r"$\lambda_j$")
ax.set_title("Dyson repulsion drift field\n" r"$\dot\lambda_i=1/(\lambda_i-\lambda_j)$")
ax.legend(loc="upper left", fontsize=8)
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig("../figures/fig5_dyson_vector_field.pdf")
print("Saved figures/fig5_dyson_vector_field.pdf")
