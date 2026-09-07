# hg19 reference annotation (downloaded 2026-09-07)

Raw files (not committed, see `.gitignore`; ~32MB total), used only for
independently validating the Hi-C compartment interpretation
(Sec. results-hic in `article.tex`) against real genomic features the
RMT method itself never sees:

- `chr18.fa.gz` — hg19 chromosome 18 reference sequence (UCSC,
  ~24.6MB), used to compute real GC content per 100kb bin.
- `refGene.txt.gz` — hg19 RefSeq gene table, genome-wide (UCSC,
  ~8MB), filtered to chr18 to compute real transcript density per bin.

Both public, no login required. This is the same genome build
(`chrom_hg19.sizes`) the GSE167150 `.hic` files' own raw headers
declare, verified directly before downloading anything, not assumed.

To reproduce:

```bash
cd data/hg19_annotation
curl -sL -o chr18.fa.gz "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz"
curl -sL -o refGene.txt.gz "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
cd ../../scripts
python3 32_hic_compartment_validation.py
```

Used directly by `scripts/32_hic_compartment_validation.py` and
`scripts/33_hic_multi_sample.py`.
