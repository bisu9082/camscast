# camscast: Atmospheric Environment version (2026-11-14)

**A global reanalysis reproduces when, but not where, PM10 is high across the dense monitoring
network of an arid emirate.** Jin Yoo, Jeongyun Kim, Doo-Hee Lee, Ku Kang.

## What this folder adds to the repository

1. **The archived Abu Dhabi PM10 values are 24-h running means, not hourly concentrations.**
   Evidence: `code/archive_checks.py` (hourly archive: increments anti-correlated at 24 h at
   19/19 monitors, unlike NO2/SO2/O3/CO from the same stations) and
   `code/averaging_fingerprint.py`. Every model comparison here averages CAMS over the same 24 h.
2. **The full OpenAQ archive** of the 19 locations (`data_raw/abudhabi_openaq_hourly.csv.gz`,
   98,459 hourly PM10 values, 29 Oct 2022 to 20 May 2024; `code/fetch_openaq_archive.py`).
3. **A continuous CAMS EAC4 PM10 series** for the nine occupied cells, 1 Oct 2022 to 31 May 2024
   (`data_processed/cams_eac4_pm10_cells.csv`, extracted by `code/build_cams_cells.py` from the
   ADS files fetched with `code/fetch_cams_eac4_pm10.py`). It reproduces all 20,818 model values
   of release v1.3 to within 1.5e-5 ug/m3. Contains modified Copernicus Atmosphere Monitoring
   Service information 2026 (CC BY 4.0); the raw NetCDF files are not included.
4. The correction analysis rerun on the same 32,135 values with the matched 24-h model mean and
   with the instantaneous model value (`results/audit_matched.json`,
   `results/audit_instantaneous.json`, `results/audit_matched_guard3.json`).
   Results of the earlier August-December 2023 subset are kept in `results/v13period/`.

## Run order (from code/)

    python fetch_openaq_archive.py            # optional: rebuilds data_raw/ from the archive
    python fetch_cams_eac4_pm10.py ../data_raw/cams   # needs an ADS account (~/.adsapirc)
    python build_cams_cells.py                # -> data_processed/cams_eac4_pm10_cells.csv
    python archive_checks.py                  # -> results/archive_checks.json
    python averaging_fingerprint.py           # -> results/averaging_fingerprint.json
    python build_records_full.py              # -> data_processed/full/ (3-hourly + daily 08 UTC)
    python qc_flatline_ae.py                  # -> data_processed/full/aligned_qc6/
    python station_representativeness.py      # -> results/station_representativeness.json
    python extended_observations.py           # -> results/extended_observations.json
    python exceedance_detection.py            # -> results/exceedance_detection.json
    python make_matched_records.py            # -> data_processed/full/derived/{matched,instantaneous}
    python feature_sets_loso.py ../data_processed/full/derived/instantaneous ../data_processed/full/derived/matched
    # correction analysis, from code/audit/ (about 30 min each on two cores):
    #   CAMSCAST_DATA=../../data_processed/full/derived/matched GUARD_DAYS=10 BLOCKS=quantile \
    #     CAMSCAST_OUT=../../results/audit_matched.json python audit_v7.py
    #   same with GUARD_DAYS=3  -> audit_matched_guard3.json
    #   CAMSCAST_DATA=../../data_processed/full/derived/instantaneous GUARD_DAYS=10 ... -> audit_instantaneous.json
    python make_tables_ae.py && python make_figs_ae.py
    python verify_numbers_ae.py               # 101 quoted phrases checked against results/
    cd ../manuscript && pdflatex main_ae && bibtex main_ae && pdflatex main_ae && pdflatex main_ae
    pdflatex SI_ae && pdflatex SI_ae && pdflatex cover_letter_ae

Software: Python 3.13, scikit-learn 1.9.1, NumPy 2.5.3, pandas 3.0.5, SciPy 1.18.1, h5py,
Matplotlib.
