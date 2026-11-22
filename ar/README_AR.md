# camscast: Atmospheric Research version (2026-11-22)

**The CAMS reanalysis follows column dust in time but not surface PM10 across the Abu Dhabi
monitoring network.** Jin Yoo, Jeongyun Kim, Doo-Hee Lee, Ku Kang.

This folder holds the version submitted to *Atmospheric Research*. It starts from `ae/` (the
Atmospheric Environment version) and changes the framing and the order of the results: the
reanalysis is first evaluated against independent dust records (column versus surface), then for
the spatial pattern of its bias and the local deficit of the conurbation, the sub-grid variance,
the statistical correction and the exceedance of the WHO 24-h guideline.

## What is new relative to `ae/`

1. **Column versus surface at the AERONET site** (`code/column_surface_check.py` ->
   `results/column_surface_check.json`, SI Table S8 via `code/make_tables_revision.py`):
   on 228 days, CAMS total AOD / AERONET total AOD = 1.07 (median), CAMS dust AOD / AERONET
   coarse-mode AOD = 0.76, CAMS dust share 0.31 against AERONET coarse fraction 0.48, and
   regional PM10 per unit of column dust 0.86 of the observed value (Wilcoxon p = 0.03).
   Exploratory; added after the main results were known.
2. **Graphical abstract** from the column/surface comparison (`code/make_graphical_abstract_ar.py`).
   It reads `results/column_surface_daily.csv`, which `column_surface_check.py` writes; that
   file holds daily AERONET values and is not redistributed here.
3. **Figure 1 without national borders** (`code/make_figs_ae.py`), following the journal's map
   policy.
4. **Manuscript** (`manuscript/`): author-year references (elsarticle-harv, 81 entries), new
   discussion of why column dust and surface PM10 diverge, including the CAMS PM10 diagnostic,
   which counts 0.4 of the mass of the largest dust bin.
5. **Checks**: `code/verify_numbers_ae.py` rebuilds 147 quoted phrases from `results/` (0 missing);
   `step7/` holds the four-layer citation check (`citation_verification_ar.csv`,
   `gate7_summary_ar.json`, `relevance_ar.json`). `code/citation_layers.py` reads an optional
   Semantic Scholar key from the environment variable `S2_KEY`; no key is stored.
6. `submission_AR_20261122/` is the upload set (LaTeX source, figures, SI, highlights, graphical
   abstract, declaration, cover letter).

## Running the code

The scripts expect `../data_raw`, `../data_processed` and `../data_aux` next to `code/`. In this
repository these are links to the folders of `ae/`. The AERONET, METAR and CAMS raw files are not
redistributed: `code/fetch_aeronet.py`, `code/fetch_metar_omaa.py`,
`code/fetch_cams_eac4_pm10.py` and `code/fetch_cams_eac4_aod_wind.py` retrieve them (the last
two need an account on the Copernicus Atmosphere Data Store).

Contains modified Copernicus Atmosphere Monitoring Service information 2026.
