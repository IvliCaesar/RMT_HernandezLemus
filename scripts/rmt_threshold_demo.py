"""
Reference implementation of the classical RMT correlation-threshold
method (Luo et al., BMC Bioinformatics 8:299, 2007) plus the two
diagnostics used throughout the program (Marchenko-Pastur bulk edge,
inverse participation ratio), run first on SYNTHETIC data with a known
planted block structure -- so the method's correctness can be checked
against ground truth before it is ever pointed at real GEO data.

No network access, no real patient data: this script only needs numpy
and scipy, and is meant to be the first thing run (per the proposal's
"Primeros pasos sugeridos"), producing a printed threshold sweep and,
if matplotlib is available, a diagnostic figure.

Method sketch (Luo et al. 2007):
  1. Build the correlation matrix, take |C| >= tau for a threshold tau
     (this is the "adjacency-like" matrix RMT is applied to).
  2. Diagonalize, get eigenvalues.
  3. "Unfold" the spectrum: fit the cumulative eigenvalue density with
     a smooth polynomial and use it to rescale eigenvalues so the mean
     local spacing is 1 everywhere (spacing statistics are only
     meaningful after unfolding).
  4. Compute nearest-neighbor spacings s_i, compare their distribution
     against the Poisson law exp(-s) and the GOE Wigner surmise
     (pi/2) s exp(-pi s^2/4).
  5. Sweep tau; the threshold at which the fit switches from
     Poisson-like to GOE-like is the RMT-selected correlation cutoff.
"""
import numpy as np
from scipy import stats


def planted_block_correlation(n_vars=400, n_modules=5, module_size=20,
                                n_samples=120, within_rho=0.6, seed=0):
    """Synthetic data with KNOWN structure: `n_modules` disjoint blocks
    of `module_size` co-regulated variables (correlation `within_rho`,
    i.e. a real planted signal), embedded in `n_vars` otherwise
    independent variables (pure noise). Returns the sample correlation
    matrix and the true module membership (for later validation only,
    never used by the method itself)."""
    rng = np.random.default_rng(seed)
    membership = -np.ones(n_vars, dtype=int)
    latent = rng.standard_normal((n_modules, n_samples))
    X = rng.standard_normal((n_vars, n_samples))
    idx = 0
    for m in range(n_modules):
        block = slice(idx, idx + module_size)
        membership[block] = m
        # each variable in the block = shared module signal + private noise,
        # mixed to give population correlation ~= within_rho
        a = np.sqrt(within_rho)
        b = np.sqrt(1 - within_rho)
        X[block] = a * latent[m][None, :] + b * X[block]
        idx += module_size
    C = np.corrcoef(X)
    return C, membership


def unfold_spectrum(eigs):
    """Unfold via a degree-7 polynomial fit to the empirical cumulative
    distribution function of the eigenvalues (standard RMT practice),
    so nearest-neighbor spacings have mean 1 and are comparable across
    matrices/thresholds."""
    eigs = np.sort(eigs)
    n = len(eigs)
    ecdf = (np.arange(1, n + 1) - 0.5) / n
    coeffs = np.polyfit(eigs, ecdf, deg=7)
    smooth_cdf = np.polyval(coeffs, eigs)
    unfolded = smooth_cdf * n
    return np.diff(np.sort(unfolded))


def poisson_goe_fit_score(spacings):
    """Returns (ks_poisson, ks_goe): Kolmogorov-Smirnov statistic of the
    empirical spacing distribution against the Poisson law exp(-s) and
    against the GOE Wigner surmise (pi/2) s exp(-pi s^2/4). Smaller is a
    better fit; whichever is smaller "wins" at this threshold."""
    spacings = spacings[spacings > 0]
    ks_poisson = stats.kstest(spacings, "expon").statistic

    def goe_cdf(s):
        return 1.0 - np.exp(-np.pi * s ** 2 / 4.0)

    ks_goe = stats.kstest(spacings, goe_cdf).statistic
    return ks_poisson, ks_goe


def marchenko_pastur_edge(p, n):
    """Theoretical bulk edge [lambda_-, lambda_+] for the eigenvalues of
    a p x n Wishart-type sample correlation matrix under pure noise
    (aspect ratio q = p/n)."""
    q = p / n
    lam_plus = (1 + np.sqrt(q)) ** 2
    lam_minus = (1 - np.sqrt(q)) ** 2
    return lam_minus, lam_plus


def inverse_participation_ratio(eigvecs):
    """IPR_k = sum_i v_{i,k}^4 for each eigenvector k (already unit
    norm). IPR ~ 1/p for a delocalized (noise) eigenvector, IPR ~ 1 for
    a fully localized one -- the diagnostic used to flag module axes."""
    return np.sum(eigvecs ** 4, axis=0)


def rmt_threshold_sweep(C, thresholds, n_samples):
    p = C.shape[0]
    lam_minus, lam_plus = marchenko_pastur_edge(p, n_samples)
    print(f"Marchenko-Pastur bulk edge for p={p}, n={n_samples}: "
          f"[{lam_minus:.3f}, {lam_plus:.3f}]")
    print(f"{'tau':>6}{'n_edges':>10}{'ks_poisson':>13}{'ks_goe':>10}{'regime':>10}")
    results = []
    for tau in thresholds:
        A = C.copy()
        A[np.abs(A) < tau] = 0.0
        np.fill_diagonal(A, 0.0)
        n_edges = int((np.abs(A) > 0).sum() / 2)
        eigs = np.linalg.eigvalsh(A)
        if len(np.unique(np.round(eigs, 8))) < 10:
            continue
        spacings = unfold_spectrum(eigs)
        ks_p, ks_g = poisson_goe_fit_score(spacings)
        regime = "GOE" if ks_g < ks_p else "Poisson"
        print(f"{tau:6.2f}{n_edges:10d}{ks_p:13.4f}{ks_g:10.4f}{regime:>10}")
        results.append((tau, n_edges, ks_p, ks_g, regime))
    return results


if __name__ == "__main__":
    C, membership = planted_block_correlation()
    n_samples = 120

    print("Synthetic ground truth: 5 modules of 20 co-regulated "
          "variables each, embedded in 400 total (300 pure noise).\n")

    thresholds = np.arange(0.10, 0.85, 0.05)
    results = rmt_threshold_sweep(C, thresholds, n_samples)

    # Dense (low tau) matrices are generic large random matrices and are
    # GOE-like almost by default (Wigner-Dyson universality); very sparse
    # (high tau) matrices are Poisson-like almost by default (a sparse,
    # nearly-disconnected graph has no level-repulsion mechanism). Neither
    # extreme is informative. The RMT-selected threshold (Luo et al. 2007)
    # is the transition point itself: sweeping tau DOWN from sparse to
    # dense, it is the largest tau at which the statistics have just
    # turned GOE-like -- i.e. the sparsest network that already shows
    # non-trivial structure, not the fully dense, mostly-noise network at
    # the smallest tau tried.
    tau_star = None
    for i in range(len(results) - 1, 0, -1):
        if results[i][4] == "Poisson" and results[i - 1][4] == "GOE":
            tau_star = results[i - 1][0]
            break
    if tau_star is not None:
        print(f"\nPoisson -> GOE transition (sweeping tau downward from "
              f"sparse to dense) at tau = {tau_star:.2f} "
              f"(RMT-selected threshold).")
    else:
        print("\nNo clean Poisson -> GOE transition found in this range; "
              "widen or refine `thresholds`.")

    # Sanity check against ground truth: at the selected threshold, do
    # within-module edges dominate over noise edges?
    if tau_star is not None:
        A = np.abs(C.copy())
        A[A < tau_star] = 0.0
        np.fill_diagonal(A, 0.0)
        within = ((membership[:, None] == membership[None, :]) &
                   (membership[:, None] >= 0) & (A > 0))
        total_edges = (A > 0).sum()
        within_edges = within.sum()
        frac = within_edges / total_edges if total_edges else float("nan")
        print(f"At tau*={tau_star:.2f}: {frac:.1%} of surviving edges are "
              f"within a true planted module (ground-truth check only, "
              f"not used by the method itself).")
