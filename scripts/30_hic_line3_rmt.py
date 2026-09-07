"""
Attacks Roadmap Line 3 directly (3D genomic structure via Hi-C,
previously identified but genuinely unattempted with real data): builds
a bin-bin correlation matrix from a real TNBC Hi-C contact map and
applies the identical RMT threshold/spiked-eigenvalue machinery used
throughout this paper for gene-expression correlation matrices, exactly
as the Roadmap proposed -- no new methodology, only new data and the
one standard Hi-C-specific preprocessing step (observed/expected
normalization) needed before that machinery is meaningful.

Data: GSE167150 (Reyes-Gopar et al. 2025, already cited in article.tex
as reyesgopar2025), sample GSM5098082 ("TNBC_Tissue3"), the smallest of
the series' three real TNBC tumor .hic files (645MB), downloaded once
via a single targeted HTTP Range request against GSE167150_RAW.tar
(scripts/_tar_range_extract.py) rather than the full 6.2GB archive.
Chromosome 18 (one of the smaller autosomes, ~78Mb, ~780 bins at 100kb
resolution) at 100kb resolution, read with the pure-Python `straw`
reader (pip package hic-straw==0.0.6, the aidenlab/straw project before
its C++/pybind11 rewrite -- no C++ compiler was available in this
environment, and this legacy pure-Python reader works identically for
a local file).

Standard Hi-C preprocessing before any spectral analysis: raw contact
counts are dominated by genomic distance alone (nearby bins always
contact more, regardless of biology), so before correlating anything
we divide out the genome-wide expected count at each genomic distance
(observed/expected, O/E) -- the same normalization step used in
Lieberman-Aiden et al. (2009), the paper that introduced exactly this
bin-bin-correlation-of-the-O/E-matrix construction to reveal A/B
chromatin compartments. Pearson-correlating the O/E matrix's rows then
gives a bin x bin correlation matrix directly analogous to the gene x
gene correlation matrices used everywhere else in this paper, to which
the identical tau-sweep/NNSD/Monte-Carlo-null pipeline is applied
unchanged.
"""
import sys
import os
import numpy as np
import pandas as pd
import straw

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score

HIC_FILE = "../data/GSE167150_HiC/GSM5098082_TNBC_Tissue3.hic"
CHROM = "18"
RESOLUTION = 100000
THRESHOLDS = np.round(np.arange(0.05, 0.85, 0.05), 2)
N_PERMUTATIONS = 20


def load_contact_matrix(hic_file, chrom, resolution, norm="NONE"):
    x, y, counts = straw.straw(norm, hic_file, chrom, chrom, "BP", resolution)
    max_pos = max(max(x), max(y))
    n_bins = max_pos // resolution + 1
    M = np.zeros((n_bins, n_bins))
    xi = np.array(x) // resolution
    yi = np.array(y) // resolution
    c = np.array(counts, dtype=float)
    M[xi, yi] = c
    M[yi, xi] = c
    return M


def filter_zero_coverage_bins(M, bin_starts):
    """Real Hi-C matrices have bins with zero (or near-zero) total
    coverage -- unmappable/repetitive regions (centromeres, assembly
    gaps), a known technical artifact, not biological signal. Dropping
    them before O/E normalization and correlation is the direct
    genomic-bin analogue of this paper's own single-cell gene filter
    (genes detected in <3% of cells, Sec. singlecell-methods): without
    it, a handful of all-zero rows produce NaN correlations that (if
    silently zeroed) dominate the top eigenvector artificially."""
    rowsum = M.sum(axis=1)
    keep = rowsum > 0
    n_dropped = int((~keep).sum())
    return M[np.ix_(keep, keep)], bin_starts[keep], n_dropped


def oe_normalize(M):
    """Observed/expected: divide every diagonal band |i-j|=d by its own
    genome-wide mean count, removing the distance-decay trend that would
    otherwise dominate any correlation computed directly on M."""
    n = M.shape[0]
    OE = np.zeros_like(M)
    for d in range(n):
        idx = np.arange(0, n - d)
        band = M[idx, idx + d]
        mean_d = band.mean()
        if mean_d > 0:
            OE[idx, idx + d] = band / mean_d
            OE[idx + d, idx] = band / mean_d
    return OE


def find_tau_star(C, thresholds):
    results = []
    for tau in thresholds:
        A = C.copy()
        A[np.abs(A) < tau] = 0.0
        np.fill_diagonal(A, 0.0)
        eigs = np.linalg.eigvalsh(A)
        if len(np.unique(np.round(eigs, 8))) < 10:
            results.append((tau, "degenerate", np.nan, np.nan))
            continue
        spacings = unfold_spectrum(eigs)
        ks_p, ks_g = poisson_goe_fit_score(spacings)
        regime = "GOE" if ks_g < ks_p else "Poisson"
        results.append((tau, regime, ks_p, ks_g))
    tau_star = None
    for i in range(len(results) - 1, 0, -1):
        if results[i][1] == "Poisson" and results[i - 1][1] == "GOE":
            tau_star = results[i - 1][0]
            break
    return tau_star, results


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    print(f"Loading real TNBC Hi-C contact matrix: chr{CHROM} at "
          f"{RESOLUTION}bp resolution from {HIC_FILE} ...")
    M = load_contact_matrix(HIC_FILE, CHROM, RESOLUTION, norm="NONE")
    n_raw = M.shape[0]
    bin_starts = np.arange(n_raw) * RESOLUTION
    print(f"chr{CHROM}: {n_raw} bins (raw), {(M > 0).sum() // 2} nonzero "
          f"contact pairs (upper triangle), matrix density "
          f"{(M > 0).sum() / (n_raw * n_raw):.3f}.")

    M, bin_starts, n_dropped = filter_zero_coverage_bins(M, bin_starts)
    n = M.shape[0]
    print(f"Dropped {n_dropped} zero-coverage bins (unmappable/centromeric "
          f"regions, a real technical artifact -- see docstring); "
          f"{n} bins remain.")

    print("\nApplying observed/expected (O/E) normalization "
          "(Lieberman-Aiden et al. 2009) to remove genomic-distance decay...")
    OE = oe_normalize(M)

    print("Computing bin-bin Pearson correlation matrix of the O/E matrix "
          "(the standard Hi-C compartment-analysis construction)...")
    C = np.corrcoef(OE)
    n_nan = int(np.isnan(C).sum())
    if n_nan > 0:
        print(f"WARNING: {n_nan} NaN entries remain after bin filtering "
              f"(zero-variance O/E rows despite nonzero coverage) -- "
              f"zeroed rather than silently left in place.")
        C = np.nan_to_num(C, nan=0.0)

    print(f"\nRunning RMT threshold sweep (identical method to "
          f"scripts/02, 23; Sec. tau-sweep) on the {n}x{n} bin-bin "
          f"correlation matrix...")
    tau_star, results = find_tau_star(C, THRESHOLDS)
    df = pd.DataFrame(results, columns=["tau", "regime", "ks_poisson", "ks_goe"])
    df.to_csv("../results/hic_chr18_rmt_sweep.csv", index=False)
    for tau, regime, ks_p, ks_g in results:
        print(f"  tau={tau:.2f}  regime={regime:>10s}  "
              f"KS(Poisson)={ks_p if ks_p==ks_p else float('nan'):.4f}  "
              f"KS(GOE)={ks_g if ks_g==ks_g else float('nan'):.4f}")
    if tau_star is not None:
        print(f"\nPoisson -> GOE transition (RMT-selected threshold): "
              f"tau* = {tau_star}")
    else:
        print("\nNo clean Poisson -> GOE transition found in this range.")

    # --- Monte Carlo permutation null: shuffle each bin's own O/E profile
    # independently (destroys real bin-bin correlation, keeps each bin's
    # own marginal fixed) -- identical logic to scripts/10, 24. ---
    print(f"\nMonte Carlo permutation null ({N_PERMUTATIONS} replicates, "
          f"identical recipe to Sec. montecarlo)...")
    rng = np.random.default_rng(0)
    perm_taus = []
    for rep in range(N_PERMUTATIONS):
        OE_perm = np.array([rng.permutation(row) for row in OE])
        C_perm = np.corrcoef(OE_perm)
        C_perm = np.nan_to_num(C_perm, nan=0.0)
        tau_p, _ = find_tau_star(C_perm, THRESHOLDS)
        perm_taus.append(tau_p)
        print(f"  permutation {rep}: tau*={tau_p}")
    perm_taus_clean = [t for t in perm_taus if t is not None]
    pd.DataFrame({"replicate": range(N_PERMUTATIONS), "tau_star": perm_taus}).to_csv(
        "../results/hic_chr18_monte_carlo_null.csv", index=False)
    if perm_taus_clean:
        print(f"\nPermutation null tau* range: "
              f"{min(perm_taus_clean)}-{max(perm_taus_clean)} "
              f"(real tau*={tau_star}).")

    # --- Top eigenvector: is it the classic A/B compartment split? ---
    eigvals, eigvecs = np.linalg.eigh(C)
    top_vec = eigvecs[:, -1]
    n_pos = int((top_vec > 0).sum())
    n_neg = int((top_vec < 0).sum())
    print(f"\nTop eigenvalue={eigvals[-1]:.3f}. Sign split of its "
          f"eigenvector across {n} bins: {n_pos} positive, {n_neg} negative "
          f"-- the classic two-compartment (A/B) signature if this "
          f"bimodal split correlates with a real, independent genomic "
          f"feature rather than being an arbitrary sign pattern.")

    pd.DataFrame({"bin_start": bin_starts,
                  "pc1_compartment_score": top_vec}).to_csv(
        "../results/hic_chr18_compartment_pc1.csv", index=False)

    print("\nSaved hic_chr18_rmt_sweep.csv, hic_chr18_monte_carlo_null.csv, "
          "hic_chr18_compartment_pc1.csv")
