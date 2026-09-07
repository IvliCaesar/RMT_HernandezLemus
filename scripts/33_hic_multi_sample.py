"""
Extends the single-tumor, single-chromosome Line 3 result (scripts 30,
32) to the other two real TNBC tumors and the matched-normal file in
GSE167150 -- the direct next step that single-sample result flagged
for itself, since one tumor cannot distinguish real patient-to-patient
variability from one sample's own idiosyncrasy. Runs the identical
chr18 pipeline (load -> filter zero-coverage bins -> O/E normalize ->
correlate -> RMT threshold sweep -> Monte Carlo null -> GC/gene-density
validation) on all four real GSE167150 samples, then asks the concrete
comparison \citet{reyesgopar2025}'s own network measures were designed
for: does the RMT-selected threshold or the compartment structure
differ between the three TNBC tumors, and between tumor and matched
normal tissue?
"""
import sys
import os
import numpy as np
import pandas as pd
import straw
from scipy import stats

sys.path.insert(0, ".")
from rmt_threshold_demo import unfold_spectrum, poisson_goe_fit_score
from importlib import import_module
hic30 = import_module("30_hic_line3_rmt")
hic32 = import_module("32_hic_compartment_validation")

CHROM = "18"
RESOLUTION = 100000
THRESHOLDS = np.round(np.arange(0.05, 0.85, 0.05), 2)
N_PERMUTATIONS = 20

SAMPLES = {
    "TNBC_Tissue1": "../data/GSE167150_HiC/GSM5098079_TNBC_Tissue1.hic",
    "TNBC_Tissue2": "../data/GSE167150_HiC/GSM5098080_TNBC_Tissue2.hic",
    "TNBC_Tissue3": "../data/GSE167150_HiC/GSM5098082_TNBC_Tissue3.hic",
    "Normal_Tissue": "../data/GSE167150_HiC/GSM5098074_Normal_Tissue.hic",
}


def run_sample(name, hic_file):
    print(f"\n{'='*70}\n{name}: {hic_file}\n{'='*70}")
    M = hic30.load_contact_matrix(hic_file, CHROM, RESOLUTION, norm="NONE")
    n_raw = M.shape[0]
    bin_starts = np.arange(n_raw) * RESOLUTION
    M, bin_starts, n_dropped = hic30.filter_zero_coverage_bins(M, bin_starts)
    n = M.shape[0]
    print(f"chr{CHROM}: {n_raw} raw bins, {n_dropped} zero-coverage dropped, {n} remain.")

    OE = hic30.oe_normalize(M)
    C = np.corrcoef(OE)
    n_nan = int(np.isnan(C).sum())
    if n_nan > 0:
        C = np.nan_to_num(C, nan=0.0)
        print(f"  ({n_nan} residual NaNs zeroed)")

    tau_star, results = hic30.find_tau_star(C, THRESHOLDS)
    print(f"tau* = {tau_star}")

    rng = np.random.default_rng(0)
    perm_taus = []
    for rep in range(N_PERMUTATIONS):
        OE_perm = np.array([rng.permutation(row) for row in OE])
        C_perm = np.corrcoef(OE_perm)
        C_perm = np.nan_to_num(C_perm, nan=0.0)
        tau_p, _ = hic30.find_tau_star(C_perm, THRESHOLDS)
        perm_taus.append(tau_p)
    perm_clean = [t for t in perm_taus if t is not None]
    perm_summary = f"{min(perm_clean)}-{max(perm_clean)}" if perm_clean else "none clean"
    print(f"permutation null tau* range: {perm_summary} "
          f"({sum(1 for t in perm_taus if t is None)}/{N_PERMUTATIONS} no clean transition)")

    eigvals, eigvecs = np.linalg.eigh(C)
    top_vec = eigvecs[:, -1]

    # GC/gene-density validation, same recipe as script 32
    seq = hic32.load_chr18_sequence(hic32.FASTA)
    gc = hic32.gc_content_per_bin(seq, RESOLUTION)
    gene_density = hic32.gene_density_per_bin(hic32.REFGENE, "chr18", len(gc), RESOLUTION)
    bin_idx = (bin_starts // RESOLUTION).astype(int)
    gc_here = gc[bin_idx]
    gd_here = gene_density[bin_idx]
    valid = ~np.isnan(gc_here)
    r_gc, p_gc = stats.pearsonr(top_vec[valid], gc_here[valid])
    r_gd, p_gd = stats.pearsonr(top_vec[valid], gd_here[valid])
    print(f"PC1 vs GC content: r={r_gc:.3f} (p={p_gc:.1e}); "
          f"PC1 vs gene density: r={r_gd:.3f} (p={p_gd:.1e})")
    # Orient PC1 sign consistently: positive = higher GC/gene density (the A compartment)
    if r_gc < 0:
        top_vec = -top_vec
        r_gc, r_gd = -r_gc, -r_gd

    return dict(name=name, n_bins_raw=n_raw, n_dropped=n_dropped, n_bins=n,
                tau_star=tau_star, perm_tau_range=perm_summary,
                perm_n_none=sum(1 for t in perm_taus if t is None),
                r_gc=r_gc, p_gc=p_gc, r_gd=r_gd, p_gd=p_gd,
                bin_starts=bin_starts, pc1=top_vec)


if __name__ == "__main__":
    os.makedirs("../results", exist_ok=True)
    all_results = {}
    for name, path in SAMPLES.items():
        if not os.path.exists(path):
            print(f"SKIP {name}: {path} not found.")
            continue
        all_results[name] = run_sample(name, path)

    summary_rows = []
    for name, r in all_results.items():
        summary_rows.append({k: v for k, v in r.items() if k not in ("bin_starts", "pc1")})
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv("../results/hic_multi_sample_summary.csv", index=False)
    print(f"\n{'='*70}\nSummary across all samples:\n{'='*70}")
    print(summary_df.to_string(index=False))

    # Pairwise compartment-call agreement (sign match at each shared bin)
    print("\nPairwise compartment sign agreement (fraction of shared bins "
          "where both samples' PC1 signs agree, after orienting both to "
          "positive = higher-GC as above):")
    names = list(all_results.keys())
    agree_rows = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = all_results[names[i]], all_results[names[j]]
            common = np.intersect1d(a["bin_starts"], b["bin_starts"])
            ia = np.searchsorted(a["bin_starts"], common)
            ib = np.searchsorted(b["bin_starts"], common)
            sign_a = np.sign(a["pc1"][ia])
            sign_b = np.sign(b["pc1"][ib])
            agreement = (sign_a == sign_b).mean()
            print(f"  {names[i]} vs {names[j]}: {agreement:.1%} agreement "
                  f"over {len(common)} shared bins")
            agree_rows.append(dict(sample_a=names[i], sample_b=names[j],
                                    n_shared_bins=len(common), agreement=agreement))
    pd.DataFrame(agree_rows).to_csv("../results/hic_multi_sample_agreement.csv", index=False)
    print("\nSaved hic_multi_sample_summary.csv and hic_multi_sample_agreement.csv")
