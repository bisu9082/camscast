"""Column amount versus column-to-surface ratio, in the cell of the AERONET site.

Question: is the reanalysis deficit in the regional PM10 level a deficit of column dust, or of the
surface PM10 that the model associates with a given column?  On the days with a sun-photometer value
(AERONET DEWA_ResearchCentre, V3 L1.5 SDA, mean of day d-1 and d as in aeronet_dust_check.py):
  - magnitude: CAMS dust AOD (550 nm) against AERONET coarse-mode AOD (500 nm), CAMS total AOD
    (550 nm) against AERONET total AOD (500 nm), CAMS dust share against AERONET coarse fraction 1-eta
  - surface per unit column: regional PM10 level (network median of the matched daily record) divided
    by the column dust measure, observed (obs / AERONET coarse AOD) and modelled (CAMS / CAMS dust AOD)
All by terciles of AERONET coarse-mode AOD.  The wavelength difference (500 vs 550 nm) is small for
the coarse mode and is not corrected.  Output: results/column_surface_check.json and the daily
values in results/column_surface_daily.csv
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

from common_ae import RES, STATIONS, daily, load_all
from aeronet_dust_check import read_sda
from dust_wind_analysis import read

HERE = os.path.dirname(os.path.abspath(__file__))


def cams_daily_dewa(c):
    r = c[["dewa_total", "dewa_dust"]].rolling(8, min_periods=8).mean()
    r = r[r.index.hour == 6]
    r.index = r.index.normalize()
    return r


def main():
    a = read_sda()
    full = a.reindex(pd.date_range(a.index.min() - pd.Timedelta("1D"), a.index.max()))
    aw = pd.concat([full.shift(1), full]).groupby(level=0).mean()

    rec = load_all()
    obs = pd.DataFrame({s: daily(g).obs for s, g in rec.items()})
    cams = pd.DataFrame({s: daily(g).cams24 for s, g in rec.items()})
    ok = obs.notna().sum(axis=1) >= 10
    reg = pd.DataFrame({"reg_obs": obs[ok].median(axis=1), "reg_cams": cams[ok].median(axis=1)})

    D = cams_daily_dewa(read()).join(aw, how="inner").join(reg, how="inner").dropna()
    D["ratio_dust_to_coarse"] = D.dewa_dust / D.tau_c
    D["ratio_total"] = D.dewa_total / D.tau
    D["dust_share"] = D.dewa_dust / D.dewa_total
    D["coarse_fraction"] = 1 - D.eta
    D["surf_per_col_obs"] = D.reg_obs / D.tau_c
    D["surf_per_col_cams"] = D.reg_cams / D.dewa_dust
    D["surf_per_col_ratio"] = D.surf_per_col_cams / D.surf_per_col_obs

    q1, q2 = D.tau_c.quantile([1 / 3, 2 / 3])
    t = pd.Series(np.where(D.tau_c <= q1, "low", np.where(D.tau_c <= q2, "mid", "high")), index=D.index)
    cols = ["tau_c", "dewa_dust", "ratio_dust_to_coarse", "tau", "dewa_total", "ratio_total",
            "coarse_fraction", "dust_share", "reg_obs", "reg_cams", "surf_per_col_obs",
            "surf_per_col_cams", "surf_per_col_ratio"]
    med = lambda m: {k: float(D[k][m].median()) for k in cols}
    res = dict(n_days=int(len(D)), coarse_tercile_bounds=[float(q1), float(q2)], all=med(D.index == D.index))
    res["by_coarse_tercile"] = {k: dict(n_days=int((t == k).sum()), **med(t == k)) for k in ("low", "mid", "high")}
    w = stats.wilcoxon(D.dust_share, D.coarse_fraction)
    res["wilcoxon_dust_share_vs_coarse_fraction_p"] = float(w.pvalue)
    w2 = stats.wilcoxon(D.surf_per_col_cams, D.surf_per_col_obs)
    res["wilcoxon_surf_per_col_p"] = float(w2.pvalue)
    res["share_days_dust_aod_below_coarse_aod"] = float((D.dewa_dust < D.tau_c).mean())
    res["share_days_surf_per_col_cams_below_obs"] = float((D.surf_per_col_cams < D.surf_per_col_obs).mean())
    json.dump(res, open(os.path.join(RES, "column_surface_check.json"), "w"), indent=1)
    D[["tau_c", "dewa_dust", "tau", "dewa_total", "eta", "reg_obs", "reg_cams"]].to_csv(
        os.path.join(RES, "column_surface_daily.csv"), index_label="day", float_format="%.5f")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
