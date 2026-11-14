"""Diagnostics added in response to the GATE-8 review (daily record, 24-h means ending 08 UTC).

  A. season of the days in each regional tercile, and terciles formed within each season
  B. sensitivity of the metropolitan minus Al Ain bias difference ("local offset") to how it is
     aggregated, and its dependence on the regional level (daily regression, bootstrap CI)
  C. nearest grid point versus bilinear interpolation of the CAMS field to the monitors

Inputs  data_processed/full/daily_08utc.csv, data_raw/cams/cams_eac4_pm10_abudhabi_*.nc
Output  results/revision_diagnostics.json
Seasons: warm = May-September, cool = November-February; March, April and October are
transitional and enter only the all-season analysis.
"""
import glob
import json
import os

import h5py
import numpy as np
import pandas as pd
from scipy import stats

from common_ae import ALAIN, METRO, RES, SEED, STATIONS, cell_of, load_all, daily

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data_raw", "cams")
WARM, COOL = (5, 6, 7, 8, 9), (11, 12, 1, 2)
NBOOT = 10000


def season(ix):
    m = pd.DatetimeIndex(ix).month
    return np.where(np.isin(m, WARM), "warm", np.where(np.isin(m, COOL), "cool", "transition"))


def terciles(reg):
    q1, q2 = reg.quantile([1 / 3, 2 / 3])
    return pd.Series(np.where(reg <= q1, "low", np.where(reg <= q2, "mid", "high")), index=reg.index)


def offset(bias, days, metro, alain, how="pooled"):
    b = bias.loc[days]
    if how == "pooled":            # as in the manuscript: mean over all monitor-days
        return float(np.nanmean(b[metro].values) - np.nanmean(b[alain].values))
    mm = b.mean()                  # per-monitor mean bias first
    if how == "monitor_mean":
        return float(mm[metro].mean() - mm[alain].mean())
    return float(mm[metro].median() - mm[alain].median())


def bilinear_daily():
    """CAMS 24-h mean ending 08 UTC, bilinearly interpolated to each monitor."""
    parts = []
    for f in sorted(glob.glob(os.path.join(RAW, "cams_eac4_pm10_abudhabi_*.nc"))):
        with h5py.File(f, "r") as h:
            t = pd.to_datetime(h["valid_time"][:], unit="s")
            lat, lon, x = h["latitude"][:], h["longitude"][:], h["pm10"][:] * 1e9
        cols = {}
        for sid, (_, la, lo, _) in STATIONS.items():
            la_s, lat_s = np.sort(lat), np.argsort(lat)
            i = np.searchsorted(la_s, la) - 1
            j = np.searchsorted(lon, lo) - 1
            wy = (la - la_s[i]) / (la_s[i + 1] - la_s[i])
            wx = (lo - lon[j]) / (lon[j + 1] - lon[j])
            i0, i1 = lat_s[i], lat_s[i + 1]
            cols[sid] = ((1 - wy) * (1 - wx) * x[:, i0, j] + (1 - wy) * wx * x[:, i0, j + 1]
                         + wy * (1 - wx) * x[:, i1, j] + wy * wx * x[:, i1, j + 1])
        parts.append(pd.DataFrame(cols, index=t))
    c = pd.concat(parts).sort_index()
    c = c[~c.index.duplicated()]
    c = c.reindex(pd.date_range(c.index.min(), c.index.max(), freq="3h"))
    r = c.rolling(8, min_periods=8).mean()          # 24 h ending at 06 UTC ...
    # the 24 h ending 08 UTC holds the steps 09, 12, ..., 06 UTC: the mean ending at 06 UTC
    r = r[r.index.hour == 6]
    r.index = r.index.normalize()
    return r


def main():
    rng = np.random.default_rng(SEED)
    rec = load_all()
    sids = list(STATIONS)
    name = {s: STATIONS[s][0] for s in sids}
    metro = [s for s in sids if name[s] in METRO]
    alain = [s for s in sids if name[s] in ALAIN]
    obs = pd.DataFrame({s: daily(g).obs for s, g in rec.items()})
    cams = pd.DataFrame({s: daily(g).cams24 for s, g in rec.items()})
    ok = obs.notna().sum(axis=1) >= 10
    reg = obs[ok].median(axis=1)
    bias = (cams - obs)[ok]
    cls = terciles(reg)
    sea = pd.Series(season(reg.index), index=reg.index)
    res = dict(_meta=dict(seasons=dict(warm=list(WARM), cool=list(COOL)), n_days=int(ok.sum()),
                          n_boot=NBOOT, seed=SEED))

    # ---- A. season composition and within-season terciles ---------------------------------
    comp = {k: {s: int(((cls == k) & (sea == s)).sum()) for s in ("cool", "transition", "warm")}
            for k in ("low", "mid", "high")}
    within = {}
    for s in ("cool", "warm"):
        r_s = reg[sea == s]
        c_s = terciles(r_s)
        within[s] = dict(n_days=int(len(r_s)), classes={})
        for k in ("low", "mid", "high"):
            d = c_s.index[c_s == k]
            within[s]["classes"][k] = dict(
                n_days=int(len(d)), regional_median=float(r_s[d].median()),
                bias_metro=float(np.nanmean(bias.loc[d, metro].values)),
                bias_al_ain=float(np.nanmean(bias.loc[d, alain].values)),
                metro_minus_al_ain=offset(bias, d, metro, alain))
    res["A_season"] = dict(tercile_by_season=comp, within_season_terciles=within)

    # ---- B. aggregation sensitivity of the offset -------------------------------------------
    variants = {}
    for lab, how, m in (("pooled_all", "pooled", metro),
                        ("pooled_excl_al_mafraq", "pooled", [s for s in metro if name[s] != "Al Mafraq"]),
                        ("monitor_mean", "monitor_mean", metro),
                        ("monitor_median", "monitor_median", metro)):
        variants[lab] = {k: offset(bias, cls.index[cls == k], m, alain, how)
                         for k in ("low", "mid", "high")}
    loo = {}
    for s in metro:
        mm = [x for x in metro if x != s]
        loo[name[s]] = {k: offset(bias, cls.index[cls == k], mm, alain) for k in ("low", "mid", "high")}
    # daily offset versus regional level
    slopes = {}
    nomaf = [s for s in metro if name[s] != "Al Mafraq"]
    for lab, mset, sname in (("all", metro, None), ("warm", metro, "warm"), ("cool", metro, "cool"),
                             ("all_excl_al_mafraq", nomaf, None),
                             ("warm_excl_al_mafraq", nomaf, "warm")):
        dm = bias[mset].mean(axis=1, skipna=True)
        da = bias[alain].mean(axis=1, skipna=True)
        keep = (bias[mset].notna().sum(axis=1) >= 4) & (bias[alain].notna().sum(axis=1) >= 2)
        D = pd.DataFrame(dict(off=dm - da, reg=reg, sea=sea))[keep]
        sel = (D.sea == sname) if sname else (D.index == D.index)
        x, y = D.reg[sel].values, D.off[sel].values
        b = stats.linregress(x, y)
        bs = []
        for _ in range(NBOOT):
            i = rng.integers(0, len(x), len(x))
            if np.ptp(x[i]) > 0:
                bs.append(stats.linregress(x[i], y[i]).slope)
        slopes[lab] = dict(n_days=int(len(x)), slope=float(b.slope), r=float(b.rvalue),
                           ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                           offset_median=float(np.median(y)),
                           offset_on_days_below_q1=float(np.median(y[x <= reg.quantile(1 / 3)]))
                           if (x <= reg.quantile(1 / 3)).any() else None)
    res["B_offset"] = dict(variants=variants, leave_one_metro_out=loo,
                           loo_range={k: [min(v[k] for v in loo.values()), max(v[k] for v in loo.values())]
                                      for k in ("low", "mid", "high")},
                           daily_regression=dict(min_metro=4, min_al_ain=2, fits=slopes))

    # ---- C. nearest versus bilinear ---------------------------------------------------------
    bl = bilinear_daily().reindex(obs.index)
    out = {}
    for lab, C in (("nearest", cams), ("bilinear", bl[sids])):
        B = (C - obs)
        sm_o, sm_c = obs.mean(), C.where(obs.notna()).mean()
        r_each = [stats.pearsonr(*pd.concat([obs[s], C[s]], axis=1).dropna().values.T).statistic
                  for s in sids]
        b_ok = B[ok]
        out[lab] = dict(
            station_bias=dict(min=float(B.mean().min()), max=float(B.mean().max()),
                              mean=float(B.mean().mean())),
            pearson_r_median=float(np.median(r_each)),
            spearman_station_means=float(stats.spearmanr(sm_o, sm_c).statistic),
            cams_station_mean_range=[float(sm_c.min()), float(sm_c.max())],
            regional_spearman=float(stats.spearmanr(reg, C[ok].median(axis=1), nan_policy="omit").statistic),
            metro_minus_al_ain={k: offset(b_ok, cls.index[cls == k], metro, alain) for k in ("low", "mid", "high")},
            cams_between_monitor_sd_median=float(C[ok].std(axis=1, ddof=1).median()))
        out[lab]["_bias_by_station"] = {name[s]: float(B[s].mean()) for s in sids}
    bn, bb = pd.Series(out["nearest"]["_bias_by_station"]), pd.Series(out["bilinear"]["_bias_by_station"])
    out["max_abs_change_station_bias"] = float((bb - bn).abs().max())
    out["spearman_bias_rank_nearest_vs_bilinear"] = float(stats.spearmanr(bn, bb).statistic)
    res["C_interpolation"] = out

    os.makedirs(RES, exist_ok=True)
    json.dump(res, open(os.path.join(RES, "revision_diagnostics.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "C_interpolation"}, indent=1))
    for lab in ("nearest", "bilinear"):
        print(lab, json.dumps({k: v for k, v in out[lab].items() if k != "_bias_by_station"}))
    print("max abs change", out["max_abs_change_station_bias"], "rank rho",
          out["spearman_bias_rank_nearest_vs_bilinear"])


if __name__ == "__main__":
    main()
