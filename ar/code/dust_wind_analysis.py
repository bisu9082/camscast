"""Dust share of the CAMS aerosol and 10-m wind, against the local deficit of the conurbation.

Run after fetch_cams_eac4_aod_wind.py has written data_raw/cams/cams_eac4_aodwind_abudhabi_*.nc.

For each day of the daily record (24 h ending 08 UTC = the eight analysis times 09..06 UTC):
  - CAMS dust AOD, sea-salt AOD and total AOD at 550 nm, averaged over the three cells that hold
    the metropolitan monitors and in the Al Ain cell, and the dust share dust/total
  - 10-m wind speed and direction (vector mean) over the metropolitan cells
Then the metropolitan minus Al Ain bias is stratified by terciles of the CAMS dust AOD, of the
wind speed, and by wind sector (NW quadrant 270-360 deg, the shamal and sea-breeze direction,
versus the rest), with Spearman correlations.  Output: results/dust_wind_analysis.json
"""
import glob
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import stats

from common_ae import ALAIN, METRO, RES, STATIONS, daily, load_all

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data_raw", "cams")
CELLS = {"metro": [(24.00, 54.75), (24.75, 54.00), (24.75, 54.75)], "al_ain": [(24.00, 55.50)],
         "dewa": [(24.75, 55.50)]}            # cell of the AERONET site, for the consistency check
VARS = {"aod550": "total", "duaod550": "dust", "ssaod550": "seasalt", "u10": "u10", "v10": "v10"}


def read():
    files = sorted(glob.glob(os.path.join(RAW, "cams_eac4_aodwind_abudhabi_*.nc")))
    if not files:
        sys.exit("run fetch_cams_eac4_aod_wind.py first")
    parts = []
    for f in files:
        with h5py.File(f, "r") as h:
            t = pd.to_datetime(h["valid_time"][:], unit="s")
            lat, lon = h["latitude"][:], h["longitude"][:]
            missing = [v for v in VARS if v not in h]
            if missing:
                sys.exit(f"{os.path.basename(f)}: variables {missing} not found; present: {list(h.keys())}")
            cols = {}
            for cname, pts in CELLS.items():
                ij = []
                for la, lo in pts:
                    i = int(np.argmin(abs(lat - la))); j = int(np.argmin(abs(lon - lo)))
                    assert abs(lat[i] - la) < 1e-6 and abs(lon[j] - lo) < 1e-6
                    ij.append((i, j))
                for v, lab in VARS.items():
                    x = h[v][:]
                    cols[f"{cname}_{lab}"] = np.mean([x[:, i, j] for i, j in ij], axis=0).astype(float)
        parts.append(pd.DataFrame(cols, index=t))
    c = pd.concat(parts).sort_index()
    c = c[~c.index.duplicated()]
    return c.reindex(pd.date_range(c.index.min(), c.index.max(), freq="3h"))


def to_daily(c):
    r = c.rolling(8, min_periods=8).mean()            # window 09..06 UTC, ending 06 UTC
    r = r[r.index.hour == 6]
    r.index = r.index.normalize()
    out = pd.DataFrame(index=r.index)
    for k in CELLS:
        out[f"{k}_dust"] = r[f"{k}_dust"]
        out[f"{k}_seasalt"] = r[f"{k}_seasalt"]
        out[f"{k}_dust_share"] = r[f"{k}_dust"] / r[f"{k}_total"]
        out[f"{k}_seasalt_share"] = r[f"{k}_seasalt"] / r[f"{k}_total"]
    u, v = r["metro_u10"], r["metro_v10"]
    spd = np.sqrt(c["metro_u10"] ** 2 + c["metro_v10"] ** 2).rolling(8, min_periods=8).mean()
    spd = spd[spd.index.hour == 6]; spd.index = spd.index.normalize()
    out["wind_speed"] = spd
    out["wind_dir_from"] = (np.degrees(np.arctan2(-u, -v)) + 360) % 360
    return out


def terc(x):
    q1, q2 = x.quantile([1 / 3, 2 / 3])
    return pd.Series(np.where(x <= q1, "low", np.where(x <= q2, "mid", "high")), index=x.index)


def main():
    rec = load_all()
    sids = list(STATIONS)
    name = {s: STATIONS[s][0] for s in sids}
    metro = [s for s in sids if name[s] in METRO]
    alain = [s for s in sids if name[s] in ALAIN]
    obs = pd.DataFrame({s: daily(g).obs for s, g in rec.items()})
    cams = pd.DataFrame({s: daily(g).cams24 for s, g in rec.items()})
    ok = obs.notna().sum(axis=1) >= 10
    bias = (cams - obs)[ok]
    keep = (bias[metro].notna().sum(axis=1) >= 4) & (bias[alain].notna().sum(axis=1) >= 2)
    off = (bias[metro].mean(axis=1) - bias[alain].mean(axis=1))[keep]
    off_nm = (bias[[s for s in metro if name[s] != "Al Mafraq"]].mean(axis=1) - bias[alain].mean(axis=1))[keep]
    D = to_daily(read()).join(pd.DataFrame({"off": off, "off_excl_al_mafraq": off_nm,
                                            "reg": obs[ok].median(axis=1)}), how="inner").dropna()
    res = dict(_meta=dict(cells={k: [list(p) for p in v] for k, v in CELLS.items()}, n_days=int(len(D))))
    res["dust_share_median"] = {k: float(D[f"{k}_dust_share"].median()) for k in CELLS}
    res["seasalt_aod_median"] = {k: float(D[f"{k}_seasalt"].median()) for k in CELLS}
    res["seasalt_share_median"] = {k: float(D[f"{k}_seasalt_share"].median()) for k in CELLS}
    res["seasalt_share_max"] = {k: float(D[f"{k}_seasalt_share"].max()) for k in CELLS}
    res["nw_sector_days"] = int(((D.wind_dir_from >= 270) & (D.wind_dir_from < 360)).sum())
    sp = lambda a, b: dict(rho=float(stats.spearmanr(a, b).statistic), p=float(stats.spearmanr(a, b).pvalue))
    res["spearman"] = dict(off_vs_dust_aod=sp(D.off, D.metro_dust), off_vs_wind_speed=sp(D.off, D.wind_speed),
                           off_excl_vs_dust_aod=sp(D.off_excl_al_mafraq, D.metro_dust),
                           off_excl_vs_wind_speed=sp(D.off_excl_al_mafraq, D.wind_speed),
                           reg_vs_dust_aod=sp(D.reg, D.metro_dust))
    for var, lab in (("metro_dust", "by_dust_aod_tercile"), ("wind_speed", "by_wind_speed_tercile")):
        t = terc(D[var])
        res[lab] = {k: dict(n_days=int((t == k).sum()), median=float(D[var][t == k].median()),
                            off_median=float(D.off[t == k].median()),
                            off_excl_al_mafraq_median=float(D.off_excl_al_mafraq[t == k].median()))
                    for k in ("low", "mid", "high")}
    nw = (D.wind_dir_from >= 270) & (D.wind_dir_from < 360)
    res["by_sector"] = {lab: dict(n_days=int(m.sum()), off_median=float(D.off[m].median()),
                                  off_excl_al_mafraq_median=float(D.off_excl_al_mafraq[m].median()))
                        for lab, m in (("NW_270_360", nw), ("other", ~nw))}
    mw = stats.mannwhitneyu(D.off[nw], D.off[~nw])
    res["by_sector"]["mann_whitney_p"] = float(mw.pvalue)
    # consistency check against the independent sun photometer (same cell, same days)
    from aeronet_dust_check import read_sda
    a = read_sda()
    full = a.reindex(pd.date_range(a.index.min() - pd.Timedelta("1D"), a.index.max()))
    aw = pd.concat([full.shift(1), full]).groupby(level=0).mean()
    J = D[["dewa_dust", "dewa_dust_share"]].join(aw, how="inner").dropna()
    res["aeronet_consistency"] = dict(
        n_days=int(len(J)),
        cams_dust_aod_vs_aeronet_coarse_aod=sp(J.dewa_dust, J.tau_c),
        cams_dust_share_vs_aeronet_fine_fraction=sp(J.dewa_dust_share, J.eta))
    json.dump(res, open(os.path.join(RES, "dust_wind_analysis.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
