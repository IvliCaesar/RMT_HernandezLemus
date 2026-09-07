# GSE167150 real TNBC Hi-C (downloaded 2026-09-07)

Raw file (not committed, see `.gitignore`; ~645MB):

- `GSM5098082_TNBC_Tissue3.hic` — real patient-tumor Hi-C contact map,
  multi-resolution (10/20/40/100/500kb), Juicer `.hic` v8 format.
  Chosen as the smallest of GSE167150's three real TNBC tumor samples
  (GSM5098079 "TNBC_Tissue1" 846MB, GSM5098080 "TNBC_Tissue2" 755MB,
  GSM5098082 "TNBC_Tissue3" 645MB).

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
python3 download_hic_tnbc_tissue3.py
python3 30_hic_line3_rmt.py
python3 31_hic_figures.py
```

`download_hic_tnbc_tissue3.py` locates the target member inside the
remote tar with a sequence of tiny HTTP Range requests against its
USTAR headers (each header's declared file size gives the exact byte
offset of the next header, so no data is downloaded while searching),
then issues one Range GET for exactly that member's ~645MB — a public
technique for extracting one file from a remote tar, not anything
GSE167150-specific.

Used directly by `article.tex` Sec. results-hic (Line 3, previously
"genuinely unattempted with real data" in the Roadmap) via
`scripts/30_hic_line3_rmt.py` and `scripts/31_hic_figures.py`.
