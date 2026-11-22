"""Evidence that the archived Abu Dhabi PM10 values are trailing 24-h running means.

Four signatures, each computed per monitor on the regular three-hourly grid. If hourly series
from other sites are present in the data directory (aligned_<Name>.csv) they are processed as a
control; none is used in the manuscript, where the control is the gases reported by the same
stations (archive_checks.py):

  acf_3h        lag-3 h autocorrelation of the observed series
  var_ratio     Var(x_t - x_{t-3h}) / Var(x_t - x_{t-24h}); 0.125 for a 24-h running mean of
                white noise, 1 for white noise itself
  r_trailing / r_centred / r_leading
                correlation of the observed series with the 24-h mean of CAMS ending at,
                centred on, and starting at the time stamp. A trailing average in the
                observations makes the trailing window the best match.
  diurnal_amp   range of the mean within-day departure by hour of day, on complete days,
                relative to the mean concentration

Writes results/averaging_fingerprint.json
"""
import glob
import json
import os

import numpy as np
import pandas as pd

from common_ae import RES, STATIONS

# the release-v1.3 three-hourly records (to 31 December 2023), with CAMS at the observation
# time stamps only, exactly as analysed in the earlier versions
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data_processed")


def fingerprint(f):
    d = pd.read_csv(f, index_col=0, parse_dates=True)[["cams", "obs"]].astype(float)
    g = d.reindex(pd.date_range(d.index.min(), d.index.max(), freq="3h"))
    tr = g.cams.rolling(8, min_periods=7).mean()
    x = g.dropna()
    full = x.groupby(x.index.floor("D")).obs.transform("size") == 8
    y = x[full]
    dep = y.obs - y.groupby(y.index.floor("D")).obs.transform("mean")
    prof = dep.groupby(dep.index.hour).mean()
    step = g.obs - g.obs.shift(1)
    return dict(
        n=int(len(d)), mean=float(d.obs.mean()),
        acf_3h=float(g.obs.corr(g.obs.shift(1))), acf_24h=float(g.obs.corr(g.obs.shift(8))),
        acf_3h_cams=float(g.cams.corr(g.cams.shift(1))),
        var_ratio=float(step.var() / (g.obs - g.obs.shift(8)).var()),
        var_ratio_cams=float((g.cams - g.cams.shift(1)).var()
                             / (g.cams - g.cams.shift(8)).var()),
        r_instant=float(g.obs.corr(g.cams)), r_trailing=float(g.obs.corr(tr)),
        r_centred=float(g.obs.corr(tr.shift(-4))), r_leading=float(g.obs.corr(tr.shift(-7))),
        diurnal_amp_rel=float((prof.max() - prof.min()) / d.obs.mean()),
        n_complete_days=int(full.sum() // 8),
        frac_repeat=float((step == 0).sum() / step.notna().sum()),
    )


def med(rows, k):
    return float(np.median([r[k] for r in rows]))


def main():
    pm10 = {s: fingerprint(os.path.join(DATA, f"aligned_pm10_{s}.csv")) for s in STATIONS}
    ctrl = {os.path.basename(f)[8:-4]: fingerprint(f)
            for f in sorted(glob.glob(os.path.join(DATA, "aligned_[A-Z]*.csv")))}
    keys = ("acf_3h", "acf_24h", "var_ratio", "var_ratio_cams", "acf_3h_cams", "r_instant",
            "r_trailing", "r_centred", "r_leading", "diurnal_amp_rel", "frac_repeat")
    P = list(pm10.values())
    C = list(ctrl.values())
    res = dict(
        pm10=pm10, control_pm25=ctrl,
        pm10_median={k: med(P, k) for k in keys},
        control_median={k: med(C, k) for k in keys} if C else None,
        control_range={k: [float(min(r[k] for r in C)), float(max(r[k] for r in C))]
                       for k in keys} if C else None,
        pm10_range={k: [float(min(r[k] for r in P)), float(max(r[k] for r in P))]
                    for k in keys},
        n_pm10_trailing_best=int(sum(r["r_trailing"] > r["r_centred"] > r["r_leading"]
                                     for r in P)),
        n_control_trailing_best=int(sum(r["r_trailing"] > r["r_centred"] > r["r_leading"]
                                        for r in C)),
    )
    os.makedirs(RES, exist_ok=True)
    json.dump(res, open(os.path.join(RES, "averaging_fingerprint.json"), "w"), indent=1)
    for lab, m, rg in (("PM10 (19)", res["pm10_median"], res["pm10_range"]),
                       ("PM2.5 control (%d)" % len(C), res["control_median"],
                        res["control_range"])):
        if m is None:
            continue
        print(f"{lab:18s} acf3h {m['acf_3h']:.3f} [{rg['acf_3h'][0]:.2f},{rg['acf_3h'][1]:.2f}]"
              f" | var ratio {m['var_ratio']:.3f} [{rg['var_ratio'][0]:.2f},"
              f"{rg['var_ratio'][1]:.2f}] | r inst/trail/centre/lead {m['r_instant']:.2f}/"
              f"{m['r_trailing']:.2f}/{m['r_centred']:.2f}/{m['r_leading']:.2f} | diurnal amp "
              f"{100*m['diurnal_amp_rel']:.1f}% of mean | repeats {100*m['frac_repeat']:.1f}%")
    print(f"trailing > centred > leading at {res['n_pm10_trailing_best']}/19 PM10 monitors, "
          f"{res['n_control_trailing_best']}/{len(C)} control series")


if __name__ == "__main__":
    main()
