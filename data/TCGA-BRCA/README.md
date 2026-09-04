# TCGA-BRCA (downloaded 2026-09-04)

Raw files (not committed, see `.gitignore`; ~69MB total):

- `HiSeqV2.gz` — RNA-seq (RSEM, log2(x+1)), 20530 genes x 1218 samples.
  Source: UCSC Xena TCGA hub,
  `https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz`
- `BRCA_clinicalMatrix` — clinical/subtype annotation, 1247 samples
  (TCGA Nature 2012 supplement fields: receptor status, PAM50 call,
  survival, etc.). Source: UCSC Xena TCGA hub,
  `https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix`
- `tnbc_expr.csv`, `luma_expr.csv` — derived cohort expression matrices
  (1500 genes x 119 TNBC / 430 LumA patients), produced by
  `../../scripts/01_build_cohorts.py`.

Both URLs are public, no login required. To reproduce:

```bash
cd data/TCGA-BRCA
curl -sL -o HiSeqV2.gz "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz"
curl -sL -o BRCA_clinicalMatrix "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
cd ../../scripts
python3 01_build_cohorts.py
```

Used directly by `article.tex` — every number in that paper's Results
section comes from these two files via `scripts/01_build_cohorts.py`
through `scripts/08_spiked_eigenvalues.py`, in order.
