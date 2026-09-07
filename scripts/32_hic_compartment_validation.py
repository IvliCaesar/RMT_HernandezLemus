"""
Independent validation of the compartment interpretation script 30
flagged as real-but-unverified: correlates the chr18 top-eigenvector
sign pattern (results/hic_chr18_compartment_pc1.csv) against two real,
independent genomic features never seen by the RMT method itself --
GC content and RefSeq gene density -- computed directly from the same
genome build the .hic file itself declares (hg19, read from the file's
own header) rather than assumed. Lieberman-Aiden et al. (2009) and the
compartment literature since predict A compartments (open, active
chromatin) should have higher GC content and higher gene density than
B compartments; this is the concrete, feature-level check that a
plaid-looking correlation matrix and a contiguous eigenvector sign
pattern, on their own, cannot provide.

Data: hg19 chr18 reference sequence (data/hg19_annotation/chr18.fa.gz,
UCSC) and hg19 RefSeq gene table (data/hg19_annotation/refGene.txt.gz,
UCSC) -- both public, downloaded directly, same build the GSE167150
.hic file's own header declares
(chrom_hg19.sizes, verified by reading the raw file header before
downloading anything).
"""
import gzip
import numpy as np
import pandas as pd
from scipy import stats

RESOLUTION = 100000
FASTA = "../data/hg19_annotation/chr18.fa.gz"
REFGENE = "../data/hg19_annotation/refGene.txt.gz"
COMPARTMENT_CSV = "../results/hic_chr18_compartment_pc1.csv"


def load_chr18_sequence(fasta_gz):
    seq_parts = []
    with gzip.open(fasta_gz, "rt") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip())
    return "".join(seq_parts).upper()


def gc_content_per_bin(seq, resolution):
    n_bins = len(seq) // resolution + 1
    gc = np.full(n_bins, np.nan)
    for i in range(n_bins):
        chunk = seq[i * resolution:(i + 1) * resolution]
        acgt = sum(chunk.count(b) for b in "ACGT")
        if acgt > 0:
            gc[i] = (chunk.count("G") + chunk.count("C")) / acgt
    return gc


def gene_density_per_bin(refgene_gz, chrom, n_bins, resolution):
    density = np.zeros(n_bins)
    with gzip.open(refgene_gz, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            # UCSC refGene.txt columns: bin, name, chrom, strand, txStart, txEnd, ...
            if fields[2] != chrom:
                continue
            tx_start, tx_end = int(fields[4]), int(fields[5])
            b0, b1 = tx_start // resolution, tx_end // resolution
            for b in range(b0, min(b1, n_bins - 1) + 1):
                density[b] += 1
    return density


if __name__ == "__main__":
    print("Loading hg19 chr18 reference sequence...")
    seq = load_chr18_sequence(FASTA)
    print(f"chr18 length: {len(seq)} bp (hg19, matches the .hic file's own "
          f"declared build, chrom_hg19.sizes, read from its raw header).")

    gc = gc_content_per_bin(seq, RESOLUTION)
    n_bins = len(gc)
    print(f"Computed GC content for {n_bins} bins at {RESOLUTION}bp resolution.")

    print("\nParsing hg19 RefSeq gene annotation (refGene.txt) for chr18...")
    gene_density = gene_density_per_bin(REFGENE, "chr18", n_bins, RESOLUTION)
    print(f"Total RefSeq transcripts overlapping chr18 bins: {gene_density.sum():.0f}")

    comp = pd.read_csv(COMPARTMENT_CSV)
    comp["bin_idx"] = (comp["bin_start"] // RESOLUTION).astype(int)

    df = pd.DataFrame({
        "bin_idx": np.arange(n_bins),
        "bin_start": np.arange(n_bins) * RESOLUTION,
        "gc_content": gc,
        "gene_density": gene_density,
    })
    merged = comp.merge(df, on=["bin_idx", "bin_start"], how="inner")
    merged = merged.dropna(subset=["gc_content"])
    print(f"\nMerged {len(merged)} bins (of {len(comp)} retained by script 30 "
          f"after zero-coverage filtering) with valid GC content.")

    pc1 = merged["pc1_compartment_score"].to_numpy()
    gc_vals = merged["gc_content"].to_numpy()
    gd_vals = merged["gene_density"].to_numpy()

    r_gc, p_gc = stats.pearsonr(pc1, gc_vals)
    r_gd, p_gd = stats.pearsonr(pc1, gd_vals)
    print(f"\nPearson correlation, PC1 compartment score vs. GC content: "
          f"r={r_gc:.3f}, p={p_gc:.2e}")
    print(f"Pearson correlation, PC1 compartment score vs. gene density: "
          f"r={r_gd:.3f}, p={p_gd:.2e}")

    pos = merged[merged["pc1_compartment_score"] > 0]
    neg = merged[merged["pc1_compartment_score"] < 0]
    print(f"\nPositive-PC1 bins (n={len(pos)}): mean GC={pos['gc_content'].mean():.4f}, "
          f"mean gene density={pos['gene_density'].mean():.3f} transcripts/bin")
    print(f"Negative-PC1 bins (n={len(neg)}): mean GC={neg['gc_content'].mean():.4f}, "
          f"mean gene density={neg['gene_density'].mean():.3f} transcripts/bin")

    t_gc, pt_gc = stats.ttest_ind(pos["gc_content"], neg["gc_content"], equal_var=False)
    t_gd, pt_gd = stats.ttest_ind(pos["gene_density"], neg["gene_density"], equal_var=False)
    print(f"\nWelch t-test (positive vs. negative PC1 bins): "
          f"GC content t={t_gc:.2f}, p={pt_gc:.2e}; "
          f"gene density t={t_gd:.2f}, p={pt_gd:.2e}")

    merged.to_csv("../results/hic_chr18_compartment_validation.csv", index=False)
    print("\nSaved results/hic_chr18_compartment_validation.csv")
