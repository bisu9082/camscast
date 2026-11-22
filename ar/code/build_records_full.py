"""Build the observation-model records for the whole archive.

Inputs   data_raw/abudhabi_openaq_hourly.csv.gz          (archived hourly PM10, all monitors)
         data_processed/cams_eac4_pm10_cells.csv          (continuous three-hourly CAMS, by cell)
Outputs  data_processed/full/aligned_pm10_<id>.csv  every archived value at a three-hourly time
             stamp, with the instantaneous CAMS value of the monitor's cell (cams) and
             resid = obs - cams; same layout as the release-v1.3 records
         data_processed/full/daily_08utc.csv  one row per monitor and day: the archived value at
             08 UTC (the 24-h mean ending then) and the mean of the eight CAMS values in the 24 h
             ending at 08 UTC (09, 12, ..., 06 UTC); 08 UTC is the only hour held in every month
"""
import os

import numpy as np
import pandas as pd

from common_ae import STATIONS, cell_of

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data_raw", "abudhabi_openaq_hourly.csv.gz")
CAMS = os.path.join(HERE, "..", "data_processed", "cams_eac4_pm10_cells.csv")
OUT = os.path.join(HERE, "..", "data_processed", "full")


def main():
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_csv(RAW, parse_dates=["t"])
    d = d[(d.parameter == "pm10") & (d.value >= 0)]
    d["t"] = d.t.dt.tz_convert(None)
    c = pd.read_csv(CAMS, index_col=0, parse_dates=True)
    grid = pd.date_range(c.index.min(), c.index.max(), freq="3h")
    c = c.reindex(grid)
    n3, rows = 0, []
    for sid, (_, la, lo, _) in STATIONS.items():
        k = f"{cell_of(la, lo)[0]:.2f}_{cell_of(la, lo)[1]:.2f}"
        s = d[d.location_id == int(sid)].set_index("t").value.sort_index()
        s = s[(s.index >= c.index.min()) & (s.index <= c.index.max())]
        # three-hourly record
        o3 = s[s.index.hour % 3 == 0]
        a = pd.DataFrame({"cams": c[k].reindex(o3.index), "obs": o3})
        a["resid"] = a.obs - a.cams
        a = a.dropna()
        a.to_csv(os.path.join(OUT, f"aligned_pm10_{sid}.csv"))
        n3 += len(a)
        # daily record at 08 UTC: CAMS mean over (t - 24 h, t]
        o8 = s[s.index.hour == 8]
        for t, v in o8.items():
            w = c[k][(c.index > t - pd.Timedelta(hours=24)) & (c.index <= t)]
            rows.append(dict(sid=sid, day=t.normalize(), obs=float(v),
                             cams24=float(w.mean()) if w.notna().sum() == 8 else np.nan,
                             n_cams=int(w.notna().sum())))
    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(OUT, "daily_08utc.csv"), index=False)
    print(f"three-hourly records: {n3} values; daily: {len(D)} monitor-days, "
          f"{int(D.cams24.notna().sum())} with a complete CAMS window")


if __name__ == "__main__":
    main()
