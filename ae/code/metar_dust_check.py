"""Dust indicator from the METAR reports of Abu Dhabi International Airport (OMAA), inside the
metropolitan area of the network, against the daily PM10 record.

Input   data_raw/metar/OMAA_metar_20220930_20240531.csv (fetch_metar_omaa.py)
Output  results/metar_dust_check.json

For the day d of the daily record, the reports in the 24 h ending 08 UTC on d, i.e. valid in
(d-1 08:00, d 08:00] UTC, are used. Two indicators:
  dust_report   at least one report with a dust or sand code (DU, BLDU, SA, BLSA, DS, SS, PO)
  obscured      share of reports with visibility below 10 km (6.2 statute miles reported as
                CAVOK or 9999 m) and no precipitation, fog or mist code, i.e. visibility reduced
                by lithometeors (dust, sand, haze)
Days need at least 18 reports in the window.
"""
import json
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

from common_ae import ALAIN, METRO, RES, STATIONS, daily, load_all

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "data_raw", "metar", "OMAA_metar_20220930_20240531.csv")
DUST = re.compile(r"(?:^|\s|[-+])(?:BL)?(?:DU|SA)\b|\b(?:DS|SS|PO)\b")
WET = re.compile(r"RA|DZ|SN|TS|SH|FG|BR|GR|GS|UP")


def indicators():
    m = pd.read_csv(SRC, na_values="M", parse_dates=["valid"])
    m["wx"] = m.wxcodes.fillna("")
    m["dust"] = m.wx.str.contains(DUST)
    m["wet"] = m.wx.str.contains(WET)
    m["vis_km"] = m.vsby * 1.609344
    m["obsc"] = (m.vis_km < 9.9) & ~m.wet
    m["haze"] = m.wx.str.contains(r"\bHZ\b")
    # assign each report to the day whose 24-h window (ending 08 UTC) contains it
    m["day"] = (m.valid - pd.Timedelta(hours=8) - pd.Timedelta(seconds=1)).dt.normalize() + pd.Timedelta(days=1)
    g = m.groupby("day")
    out = pd.DataFrame({"n": g.size(), "dust_report": g.dust.any(), "n_dust": g.dust.sum(),
                        "obscured": g.obsc.mean(), "haze_report": g.haze.any(),
                        "vis_min_km": g.vis_km.min(), "vis_median_km": g.vis_km.median()})
    return out[out.n >= 18], m


def main():
    I, m = indicators()
    rec = load_all()
    sids = list(STATIONS)
    name = {s: STATIONS[s][0] for s in sids}
    metro = [s for s in sids if name[s] in METRO]
    alain = [s for s in sids if name[s] in ALAIN]
    nomaf = [s for s in metro if name[s] != "Al Mafraq"]
    obs = pd.DataFrame({s: daily(g).obs for s, g in rec.items()})
    cams = pd.DataFrame({s: daily(g).cams24 for s, g in rec.items()})
    ok = obs.notna().sum(axis=1) >= 10
    reg, creg = obs[ok].median(axis=1), cams[ok].median(axis=1)
    bias = (cams - obs)[ok]
    keep = (bias[metro].notna().sum(axis=1) >= 4) & (bias[alain].notna().sum(axis=1) >= 2)
    off = (bias[metro].mean(axis=1) - bias[alain].mean(axis=1))[keep]
    off_nm = (bias[nomaf].mean(axis=1) - bias[alain].mean(axis=1))[keep]
    D = I.join(pd.DataFrame({"reg": reg, "creg": creg}), how="inner").dropna(subset=["reg"])
    D = D.join(pd.DataFrame({"off": off, "off_nm": off_nm}), how="left")
    sp = lambda a, b: dict(rho=float(stats.spearmanr(a, b, nan_policy="omit").statistic),
                           p=float(stats.spearmanr(a, b, nan_policy="omit").pvalue))
    codes = m.loc[m.dust, "wx"].value_counts().to_dict()
    res = dict(_meta=dict(station="OMAA", lat=24.433, lon=54.6511, n_reports=int(len(m)),
                          n_days_with_record=int(len(D)), min_reports=18, dust_codes_found=codes,
                          n_reports_dust=int(m.dust.sum()), n_reports_haze=int(m.haze.sum())))
    res["spearman"] = dict(obscured_vs_obs_regional=sp(D.obscured, D.reg),
                           obscured_vs_cams_regional=sp(D.obscured, D.creg),
                           obscured_vs_local_deficit=sp(D.obscured, D.off),
                           obscured_vs_local_deficit_excl_al_mafraq=sp(D.obscured, D.off_nm))
    # agreement with the sun photometer, where both exist
    try:
        A = json.load(open(os.path.join(RES, "aeronet_dust_check.json")))
        from aeronet_dust_check import read_sda
        a = read_sda()
        full = a.reindex(pd.date_range(a.index.min() - pd.Timedelta("1D"), a.index.max()))
        aw = pd.concat([full.shift(1), full]).groupby(level=0).mean()
        J = D.join(aw, how="inner").dropna(subset=["tau_c"])
        res["aeronet"] = dict(n_days=int(len(J)), obscured_vs_coarse_aod=sp(J.obscured, J.tau_c),
                              coarse_aod_median_dust_report_days=float(J.tau_c[J.dust_report].median()),
                              coarse_aod_median_other_days=float(J.tau_c[~J.dust_report].median()),
                              site=A["_meta"]["site"])
    except FileNotFoundError:
        pass
    groups = {}
    for lab, sel in (("dust_report", D.dust_report), ("no_dust_report", ~D.dust_report),
                     ("clear", D.obscured == 0)):
        x = D[sel]
        groups[lab] = dict(n_days=int(len(x)), regional_obs_median=float(x.reg.median()),
                           regional_bias_median=float((x.creg - x.reg).median()),
                           local_deficit_median=float(x.off.median()),
                           local_deficit_excl_al_mafraq_median=float(x.off_nm.median()),
                           n_days_with_deficit=int(x.off.notna().sum()))
    res["groups"] = groups
    a, b = D.off[D.dust_report].dropna(), D.off[~D.dust_report].dropna()
    res["mann_whitney_local_deficit_dust_vs_other_p"] = float(stats.mannwhitneyu(a, b).pvalue)
    a, b = (D.creg - D.reg)[D.dust_report], (D.creg - D.reg)[~D.dust_report]
    res["mann_whitney_regional_bias_dust_vs_other_p"] = float(stats.mannwhitneyu(a, b).pvalue)
    cl = D[D.obscured == 0]
    res["clear_day_metro_bias_by_monitor"] = {name[s]: float(bias.loc[bias.index.intersection(cl.index), s].mean())
                                              for s in metro}
    json.dump(res, open(os.path.join(RES, "metar_dust_check.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
