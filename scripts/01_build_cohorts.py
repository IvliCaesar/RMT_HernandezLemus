"""
Real TCGA-BRCA cohort construction from the data downloaded 2026-09-04
(UCSC Xena TCGA hub, public, no login required):
  data/TCGA-BRCA/HiSeqV2.gz           -- RNA-seq (RSEM, log2(x+1)), 1218 samples
  data/TCGA-BRCA/BRCA_clinicalMatrix  -- clinical/subtype annotation, 1247 samples

Defines two real, non-overlapping cohorts by receptor status / PAM50 call
(both fields independently recorded in the TCGA Nature 2012 clinical
supplement, not derived from each other):
  - TNBC:  ER-negative AND PR-negative AND HER2-negative
  - LumA:  PAM50Call_RNAseq == 'LumA' (luminal A, the least aggressive
           subtype -- used here purely as a contrasting, non-TNBC cohort)

Restricts to the top `N_GENES` most variable genes across ALL samples
(variance computed before cohort splitting, so gene selection cannot
leak cohort-label information into which genes are kept) -- standard
practice for RMT coexpression analysis (Luo et al. 2007 apply the same
kind of variance/expression filter before building the correlation
matrix), and keeps the eigendecomposition (N_GENES x N_GENES) tractable.

Writes data/TCGA-BRCA/tnbc_expr.csv and data/TCGA-BRCA/luma_expr.csv
(genes x samples, real values, no synthetic data anywhere in this
script) for use by the downstream RMT/wavelet/Dyson/ML scripts.
"""
import gzip
import numpy as np
import pandas as pd

N_GENES = 1500

print("Loading clinical matrix...")
clin = pd.read_csv("../data/TCGA-BRCA/BRCA_clinicalMatrix", sep="\t", low_memory=False)
clin = clin.set_index("sampleID")

is_tnbc = (
    (clin["ER_Status_nature2012"] == "Negative")
    & (clin["PR_Status_nature2012"] == "Negative")
    & (clin["HER2_Final_Status_nature2012"] == "Negative")
)
is_luma = clin["PAM50Call_RNAseq"] == "LumA"

tnbc_ids = set(clin.index[is_tnbc])
luma_ids = set(clin.index[is_luma])
overlap = tnbc_ids & luma_ids
print(f"TNBC samples (clinical): {len(tnbc_ids)}")
print(f"LumA samples (clinical): {len(luma_ids)}")
print(f"Overlap (dropped from both, should be 0 or near-0): {len(overlap)}")
tnbc_ids -= overlap
luma_ids -= overlap

print("Loading expression matrix (this reads all 1218 samples)...")
expr = pd.read_csv("../data/TCGA-BRCA/HiSeqV2.gz", sep="\t", index_col=0)
expr = expr[~expr.index.duplicated(keep="first")]  # a few duplicate gene symbols in HiSeqV2
print(f"Expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

tnbc_cols = [c for c in expr.columns if c in tnbc_ids]
luma_cols = [c for c in expr.columns if c in luma_ids]
print(f"TNBC samples with expression data: {len(tnbc_cols)}")
print(f"LumA samples with expression data: {len(luma_cols)}")

# Gene selection by variance across ALL samples (cohort-blind), before splitting.
gene_var = expr.var(axis=1)
top_genes = gene_var.sort_values(ascending=False).index[:N_GENES]
print(f"Kept top {N_GENES} most-variable genes (cohort-blind selection).")

tnbc_expr = expr.loc[top_genes, tnbc_cols]
luma_expr = expr.loc[top_genes, luma_cols]

tnbc_expr.to_csv("../data/TCGA-BRCA/tnbc_expr.csv")
luma_expr.to_csv("../data/TCGA-BRCA/luma_expr.csv")

print(f"\nSaved tnbc_expr.csv: {tnbc_expr.shape}")
print(f"Saved luma_expr.csv: {luma_expr.shape}")
print("\nSample sizes for the RMT network analysis:")
print(f"  TNBC: p={tnbc_expr.shape[0]} genes, n={tnbc_expr.shape[1]} patients")
print(f"  LumA: p={luma_expr.shape[0]} genes, n={luma_expr.shape[1]} patients")
