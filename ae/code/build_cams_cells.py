"""Extract the CAMS EAC4 surface PM10 series of the occupied model cells from the ADS NetCDF files.

Input   data_raw/cams/cams_eac4_pm10_abudhabi_<year>.nc (from fetch_cams_eac4_pm10.py)
Output  data_processed/cams_eac4_pm10_cells.csv: one row per three-hourly time (UTC), one column
        per occupied cell ("<lat>_<lon>"), in ug m-3 (kg m-3 x 1e9)
        results/cams_cells_check.json: grid, period, and agreement with the model values of the
        release-v1.3 aligned records at the time stamps they share

Each monitor is assigned the grid point (round(lat/0.75)*0.75, round(lon/0.75)*0.75), the rule
of the original retrieval; the script stops if any such point is missing from the files.
"""
import glob
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd

from common_ae import GRID, RES, STATIONS, cell_of

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data_raw", "cams")
OUT = os.path.join(HERE, "..", "data_processed", "cams_eac4_pm10_cells.csv")
OLD = os.path.join(HERE, "..", "data_processed")          # release-v1.3 aligned records


def read(f, var="pm10"):
    with h5py.File(f, "r") as h:
        tn = "valid_time" if "valid_time" in h else "time"
        t = pd.to_datetime(h[tn][:], unit="s")
        lat, lon = h["latitude"][:], h["longitude"][:]
        x = h[var][:]
        if x.ndim == 4:                       # (time, level/number, lat, lon) -> squeeze
            x = x.reshape(x.shape[0], -1, x.shape[-2], x.shape[-1])[:, 0]
        units = h[var].attrs.get("units", b"")
    return t, lat, lon, x, units.decode() if isinstance(units, bytes) else str(units)


def main():
    files = sorted(glob.glob(os.path.join(RAW, "cams_eac4_pm10_abudhabi_*.nc")))
    if not files:
        sys.exit(f"no files in {RAW}")
    cells = sorted({cell_of(v[1], v[2]) for v in STATIONS.values()})
    parts, info = [], []
    for f in files:
        t, lat, lon, x, units = read(f)
        assert units.replace(" ", "") in ("kgm**-3", "kgm-3"), units
        cols = {}
        for la, lo in cells:
            i = np.where(np.isclose(lat, la, atol=1e-6))[0]
            j = np.where(np.isclose(lon, lo, atol=1e-6))[0]
            if len(i) != 1 or len(j) != 1:
                sys.exit(f"grid point {la},{lo} not in {os.path.basename(f)}")
            cols[f"{la:.2f}_{lo:.2f}"] = x[:, i[0], j[0]].astype(float) * 1e9
        parts.append(pd.DataFrame(cols, index=t))
        info.append(dict(file=os.path.basename(f), n_times=int(len(t)), first=str(t.min()),
                         last=str(t.max()), lat=[float(lat.min()), float(lat.max())],
                         lon=[float(lon.min()), float(lon.max())]))
    c = pd.concat(parts).sort_index()
    c = c[~c.index.duplicated()]
    full = pd.date_range(c.index.min(), c.index.max(), freq="3h")
    missing = full.difference(c.index)
    c.index.name = "time"
    c.to_csv(OUT, float_format="%.5f")

    # agreement with the model values used in the release-v1.3 aligned records
    diffs, n = [], 0
    for sid, (_, la, lo, _) in STATIONS.items():
        f = os.path.join(OLD, f"aligned_pm10_{sid}.csv")
        if not os.path.exists(f):
            continue
        a = pd.read_csv(f, index_col=0, parse_dates=True)
        k = f"{cell_of(la, lo)[0]:.2f}_{cell_of(la, lo)[1]:.2f}"
        m = a.cams.to_frame().join(c[k].rename("new"), how="inner")
        diffs.append((m.cams - m.new).abs().values)
        n += len(m)
    d = np.concatenate(diffs) if diffs else np.array([np.nan])
    res = dict(files=info, cells=list(c.columns), n_times=int(len(c)),
               first=str(c.index.min()), last=str(c.index.max()),
               missing_three_hourly_steps=int(len(missing)),
               v13_values_compared=int(n), v13_max_abs_diff=float(np.nanmax(d)),
               v13_median_abs_diff=float(np.nanmedian(d)),
               v13_share_within_0p01=float(np.mean(d < 0.01)))
    json.dump(res, open(os.path.join(RES, "cams_cells_check.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "files"}, indent=1))


if __name__ == "__main__":
    main()
