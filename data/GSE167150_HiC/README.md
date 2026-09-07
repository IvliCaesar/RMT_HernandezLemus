# GSE167150 real TNBC Hi-C (downloaded 2026-09-07)

Raw files (not committed, see `.gitignore`; ~3.5GB total). All four
real samples in the series are now downloaded and analyzed (extended
2026-09-07d from the original single-sample round):

- `GSM5098082_TNBC_Tissue3.hic` (645MB) — real TNBC tumor, the smallest
  of the three, analyzed first.
- `GSM5098079_TNBC_Tissue1.hic` (846MB), `GSM5098080_TNBC_Tissue2.hic`
  (755MB) — the other two real TNBC tumors.
- `GSM5098074_Normal_Tissue.hic` (2.48GB) — matched contralateral
  healthy-tissue sample.

All four: multi-resolution (10/20/40/100/500kb) Juicer `.hic` v8
format.

Source: GSE167150 (\citet{reyesgopar2025} in `article.tex`; Reyes-Gopar
et al. 2025, *Front. Cell Dev. Biol.*, "Integration of chromosome
conformation and gene expression networks reveals regulatory
mechanisms in triple negative breast cancer" — Hernández-Lemus's own
group). Public, no login required, but GEO's plain directory listing
only exposes the combined `GSE167150_RAW.tar` (6.2GB, mostly unrelated
cell-line/normal-tissue files), not per-sample downloads.

To reproduce without downloading the full 6.2GB archive:

```bash
cd scripts
python3 download_hic_tnbc_tissue3.py GSM5098082_TNBC_Tissue3 ../data/GSE167150_HiC/GSM5098082_TNBC_Tissue3.hic
python3 download_hic_tnbc_tissue3.py GSM5098079_TNBC_Tissue1 ../data/GSE167150_HiC/GSM5098079_TNBC_Tissue1.hic
python3 download_hic_tnbc_tissue3.py GSM5098080_TNBC_Tissue2 ../data/GSE167150_HiC/GSM5098080_TNBC_Tissue2.hic
python3 download_hic_tnbc_tissue3.py GSM5098074_Normal_Tissue ../data/GSE167150_HiC/GSM5098074_Normal_Tissue.hic
python3 30_hic_line3_rmt.py
python3 32_hic_compartment_validation.py
python3 33_hic_multi_sample.py
python3 31_hic_figures.py
python3 35_hic_validation_figures.py
```

`download_hic_tnbc_tissue3.py` (name kept from the original
single-sample round; it takes the target member's name and output path
as arguments and works for any member) locates the target member inside
the remote tar with a sequence of tiny HTTP Range requests against its
USTAR headers (each header's declared file size gives the exact byte
offset of the next header, so no data is downloaded while searching),
then issues one Range GET for exactly that member's data — a public
technique for extracting one file from a remote tar, not anything
GSE167150-specific.

Used directly by `article.tex` Sec. results-hic (Line 3, the one
originally proposed line left unattempted through several earlier
rounds) via `scripts/30_hic_line3_rmt.py`,
`scripts/32_hic_compartment_validation.py`,
`scripts/33_hic_multi_sample.py`, `scripts/31_hic_figures.py`, and
`scripts/35_hic_validation_figures.py`.
