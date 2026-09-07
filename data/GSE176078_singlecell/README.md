# GSE176078 single-cell TNBC (downloaded 2026-09-05)

Raw files (not committed, see `.gitignore`; ~178MB total), one real
patient sample, chosen as the largest of the 9 TNBC-labeled samples in
the series (avoids cross-patient batch effects in a first single-cell
attempt; pooling multiple patients is a direct next step, see
`sec:roadmap` in `article.tex`):

- `count_matrix_sparse.mtx` — sparse gene x cell UMI count matrix,
  29733 genes x 7986 cells, MatrixMarket format.
- `count_matrix_genes.tsv`, `count_matrix_barcodes.tsv` — row/column
  labels for the matrix above.
- `metadata.csv` — per-cell annotation from Wu et al. (2021), including
  `subtype` (confirms TNBC) and `celltype_major`/`celltype_minor`
  (the original paper's own cell-type calls, used here only as an
  independent, real check on RMT spike interpretation, never as an
  input to the RMT method itself).

Source: patient CID44971, GEO accession GSM5354531, part of series
GSE176078 (\citet{tnbcsc2021} in `article.tex`; original publication
Wu et al., *Nat. Genet.* 2021, "A single-cell and spatially resolved
atlas of human breast cancers"). Public, no login required.

To reproduce:

```bash
cd data/GSE176078_singlecell
curl -sL -o CID44971.tar.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354531/suppl/GSM5354531_CID44971.tar.gz"
tar -xzf CID44971.tar.gz --strip-components=1
rm CID44971.tar.gz
cd ../../scripts
python3 23_singlecell_line4_rmt.py
```

Used directly by `article.tex` Sec. results-singlecell (Line 4,
previously "not yet attempted with real data" in the Roadmap) via
`scripts/23_singlecell_line4_rmt.py`.
