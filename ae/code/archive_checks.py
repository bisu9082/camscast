"""Checks of the analysed record against the OpenAQ archive, at the archive's hourly resolution.

1. Every value of the analysed three-hourly record is located in the archive and compared.
2. Hourly signature of the averaging applied by the reporting chain, for PM10 and for the gases
   reported by the same stations. For a trailing 24-h running mean x of hourly values h,
   x(t) - x(t-1) = [h(t) - h(t-24)] / 24, so the hourly increments of x are negatively correlated
   with themselves 24 h later; an hourly concentration with a diurnal cycle gives a positive
   correlation at that lag.
3. Which hours of the day the archive holds, month by month.

Reads  data_raw/abudhabi_openaq_hourly.csv.gz, data_processed/aligned_pm10_<id>.csv
Writes results/archive_checks.json
"""
import json
import os

import numpy as np
import pandas as pd

from common_ae import RES, STATIONS

# the release-v1.3 three-hourly records (to 31 December 2023), with CAMS at the observation
# time stamps only, exactly as analysed in the earlier versions
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data_processed")

RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data_raw",
                   "abudhabi_openaq_hourly.csv.gz")
LAGS = (1, 8, 12, 23, 24, 25, 48)


def corr_at(x, k):
    a, b = x, x.shift(k)
    ok = a.notna() & b.notna()
    return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 200 else float("nan")


def main():
    d = pd.read_csv(RAW, parse_dates=["t"])
    d["t"] = d.t.dt.tz_convert(None)
    res = {"_meta": dict(rows=int(len(d)), first=str(d.t.min()), last=str(d.t.max()),
                         rows_by_parameter={k: int(v) for k, v in
                                            d.parameter.value_counts().items()},
                         negative_values=int((d.value < 0).sum()),
                         duplicate_rows=int(d.duplicated(["location_id", "parameter", "t"]).sum()))}
    p = d[d.parameter == "pm10"]

    # 1. the analysed record against the archive
    n = same = extra = 0
    dcoord = 0.0
    for sid, (_, lat, lon, _) in STATIONS.items():
        a = pd.read_csv(os.path.join(DATA, f"aligned_pm10_{sid}.csv"), index_col=0,
                        parse_dates=True)
        g = p[p.location_id == int(sid)]
        s = g.set_index("t").value
        c = s.reindex(a.index)
        n += len(a)
        same += int((np.abs(c - a.obs) < 1e-9).sum())
        r3 = s[(s.index.hour % 3 == 0) & (s.index <= a.index.max())]
        extra += len(r3.index.difference(a.index))
        dcoord = max(dcoord, abs(g.lat.iloc[0] - lat), abs(g.lon.iloc[0] - lon))
    res["record_vs_archive"] = dict(
        n_values=n, identical_to_archive_value_at_same_time=same,
        archive_values_at_3h_stamps_not_in_record=extra,
        max_coordinate_difference_deg=float(dcoord),
        note="the analysed value is the archived hourly value at the three-hourly time stamp")

    # 2. hourly signature by parameter
    rows = []
    for (loc, par), g in d.groupby(["location_id", "parameter"]):
        s = g.set_index("t").value.sort_index()
        s = s.reindex(pd.date_range(s.index.min(), s.index.max(), freq="h"))
        dx = s.diff()
        r = dict(location_id=int(loc), parameter=par, n=int(s.notna().sum()),
                 acf_1h=corr_at(s, 1), frac_integer=float((s.dropna() % 1 == 0).mean()),
                 frac_unchanged=float((dx.dropna() == 0).mean()))
        for k in LAGS:
            r[f"increment_acf_{k}h"] = corr_at(dx, k)
        rows.append(r)
    T = pd.DataFrame(rows)
    res["hourly_signature"] = {}
    for par, g in T.groupby("parameter"):
        e = dict(n_stations=int(len(g)), n_values=int(g.n.sum()))
        for c in g.columns[3:]:
            e[c + "_median"] = float(g[c].median())
        e["increment_acf_24h_min"] = float(g.increment_acf_24h.min())
        e["increment_acf_24h_max"] = float(g.increment_acf_24h.max())
        e["n_stations_negative_at_24h"] = int((g.increment_acf_24h < 0).sum())
        res["hourly_signature"][par] = e
    res["hourly_signature_by_station_pm10"] = T[T.parameter == "pm10"].to_dict("records")

    # 3. hours of the day held by the archive, by month
    res["pm10_by_month"] = {}
    for m, g in p.groupby(p.t.dt.strftime("%Y-%m")):
        h = g.t.dt.hour
        res["pm10_by_month"][m] = dict(
            n=int(len(g)), n_at_08utc=int((h == 8).sum()),
            n_days_with_08utc=int(g[h == 8].t.dt.floor("D").nunique()),
            hours_present=int(h.nunique()))
    json.dump(res, open(os.path.join(RES, "archive_checks.json"), "w"), indent=1)

    print(res["record_vs_archive"])
    for par, e in res["hourly_signature"].items():
        print(f"{par:5s} n={e['n_values']:6d} acf1 {e['acf_1h_median']:.3f} "
              f"integer {e['frac_integer_median']:.2f} unchanged {e['frac_unchanged_median']:.3f} "
              f"incr-acf 23/24/25 h {e['increment_acf_23h_median']:+.2f}/"
              f"{e['increment_acf_24h_median']:+.2f}/{e['increment_acf_25h_median']:+.2f} "
              f"[{e['increment_acf_24h_min']:+.2f},{e['increment_acf_24h_max']:+.2f}] "
              f"neg {e['n_stations_negative_at_24h']}/{e['n_stations']}")
    for m, e in res["pm10_by_month"].items():
        print(m, e)


if __name__ == "__main__":
    main()
