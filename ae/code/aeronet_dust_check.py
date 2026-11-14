"""Independent dust indicator from AERONET (DEWA_ResearchCentre, Dubai; Version 3, Level 1.5 daily
averages, Spectral Deconvolution Algorithm) and the regional terciles of the daily PM10 record.

Input   data_raw/aeronet/DEWA_SDA15.txt (print_web_data_v3, SDA15=1, AVG=20, 2022-10-01..2024-05-31)
Output  results/aeronet_dust_check.json

For the 24-h mean ending 08 UTC on day d, the AERONET value is the mean of the daily averages of
days d-1 and d (whichever exist), because the window spans the daylight hours of both days.
Days are classed by terciles of the coarse-mode AOD at 500 nm (no external threshold), and the
metropolitan and Al Ain bias is computed by coarse-AOD tercile and, separately, by regional tercile
within the low- and high-coarse-AOD days.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

from common_ae import ALAIN, METRO, RES, SEED, STATIONS, daily, load_all

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "data_raw", "aeronet", "DEWA_SDA15.txt")


def read_sda():
    lines = open(SRC, encoding="latin-1").read().splitlines()
    h = next(i for i, l in enumerate(lines) if l.startswith("AERONET_Site,"))
    from io import StringIO
    d = pd.read_csv(StringIO("\n".join(lines[h:])), na_values=[-999, "-999", "-999.000000"])
    d["day"] = pd.to_datetime(d.iloc[:, 1], format="%d:%m:%Y")
    out = pd.DataFrame({"tau": d["Total_AOD_500nm[tau_a]"].astype(float),
                        "tau_c": d["Coarse_Mode_AOD_500nm[tau_c]"].astype(float),
                        "eta": d["FineModeFraction_500nm[eta]"].astype(float)})
    out.index = d["day"]
    return out[~out.index.duplicated()].sort_index()


def terc(x):
    q1, q2 = x.quantile([1 / 3, 2 / 3])
    return pd.Series(np.where(x <= q1, "low", np.where(x <= q2, "mid", "high")), index=x.index), [float(q1), float(q2)]


def main():
    a = read_sda()
    full = a.reindex(pd.date_range(a.index.min() - pd.Timedelta("1D"), a.index.max()))
    win = pd.concat([full.shift(1), full], axis=1, keys=["prev", "same"])
    dust = pd.DataFrame({k: win.xs(k, axis=1, level=1).mean(axis=1) for k in ("tau", "tau_c", "eta")})

    rec = load_all()
    sids = list(STATIONS)
    name = {s: STATIONS[s][0] for s in sids}
    metro = [s for s in sids if name[s] in METRO]
    alain = [s for s in sids if name[s] in ALAIN]
    obs = pd.DataFrame({s: daily(g).obs for s, g in rec.items()})
    cams = pd.DataFrame({s: daily(g).cams24 for s, g in rec.items()})
    ok = obs.notna().sum(axis=1) >= 10
    reg = obs[ok].median(axis=1)
    creg = cams[ok].median(axis=1)
    bias = (cams - obs)[ok]
    rcls, _ = terc(reg)

    D = pd.DataFrame({"reg": reg, "creg": creg}).join(dust, how="inner").dropna(subset=["tau_c"])
    res = dict(_meta=dict(site="DEWA_ResearchCentre", lat=24.766847, lon=55.369125,
                          level="Version 3 Level 1.5 (not final calibration)", product="SDA daily average",
                          n_aeronet_days=int(len(a)), n_days_matched=int(len(D)),
                          first=str(D.index.min().date()), last=str(D.index.max().date())))
    sp = lambda x, y: dict(rho=float(stats.spearmanr(x, y).statistic), p=float(stats.spearmanr(x, y).pvalue))
    res["spearman"] = dict(coarse_aod_vs_obs_regional=sp(D.tau_c, D.reg),
                           coarse_aod_vs_cams_regional=sp(D.tau_c, D.creg),
                           fine_fraction_vs_obs_regional=sp(D.eta, D.reg))
    dcls, qb = terc(D.tau_c)
    res["coarse_aod_tercile_bounds"] = qb
    res["regional_tercile_x_dust_tercile"] = {r: {k: int(((rcls.reindex(D.index) == r) & (dcls == k)).sum())
                                                 for k in ("low", "mid", "high")} for r in ("low", "mid", "high")}
    by = {}
    for k in ("low", "mid", "high"):
        d = dcls.index[dcls == k]
        bm = float(np.nanmean(bias.loc[d, metro].values))
        ba = float(np.nanmean(bias.loc[d, alain].values))
        by[k] = dict(n_days=int(len(d)), coarse_aod_median=float(D.tau_c[d].median()),
                     fine_fraction_median=float(D.eta[d].median()),
                     regional_obs_median=float(D.reg[d].median()), regional_cams_median=float(D.creg[d].median()),
                     regional_bias_median=float((D.creg - D.reg)[d].median()),
                     bias_metro=bm, bias_al_ain=ba, metro_minus_al_ain=bm - ba,
                     metro_minus_al_ain_excl_al_mafraq=float(
                         np.nanmean(bias.loc[d, [s for s in metro if name[s] != "Al Mafraq"]].values)) - ba)
    res["by_coarse_aod_tercile"] = by
    # low-dust days only: does the metropolitan deficit persist?
    low = dcls.index[dcls == "low"]
    res["low_dust_days_metro_bias_by_monitor"] = {name[s]: float(bias.loc[low, s].mean()) for s in metro}
    json.dump(res, open(os.path.join(RES, "aeronet_dust_check.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
