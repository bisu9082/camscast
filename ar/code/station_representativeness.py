"""Reanalysis error, sub-grid spread and spatial coherence of PM10 over the Abu Dhabi network.

Reported values are trailing 24-h running means (common_ae.py); the model counterpart of each is
the CAMS mean over the same 24 h. Observation-only statistics use every reported value; every
model-against-observation statistic uses the matched record (values with a model counterpart).

  A. per-monitor statistics and Moran's I of the mean bias
  B. sub-grid spread inside each 0.75 degree CAMS cell
  C. the operator's station typology against mean bias
  D. coherence between monitor pairs, and how often two monitors that share one model value
     fall on opposite sides of the guideline
  E. bias and spread stratified by the regional level of the day

Writes results/station_representativeness.json
Usage:  CAMSCAST_DATA=../data_processed python station_representativeness.py
"""
import itertools
import json
import math
import os

import numpy as np
import pandas as pd
from scipy import stats

from common_ae import (ALAIN, METRO, RES, SEED, SOURCE_TYPES, STATIONS, WHO24, cell_of, daily,
                       haversine, load_all)


def morans_i(z, D, band, rng, nperm=9999):
    z = np.asarray(z, float) - np.mean(z)
    W = ((D <= band) & (D > 0)).astype(float)
    if W.sum() == 0:
        return None

    def stat(v):
        return len(v) / W.sum() * (v @ W @ v) / (v @ v)

    obs = stat(z)
    perm = np.array([stat(rng.permutation(z)) for _ in range(nperm)])
    return dict(I=float(obs), n_pairs=int(W.sum() // 2),
                p_perm=float((np.sum(perm >= obs) + 1) / (nperm + 1)))


def main():
    rng = np.random.default_rng(SEED)
    rec = load_all()
    sids = list(rec)
    res = {"_meta": dict(who_24h_guideline=WHO24, seed=SEED, n_stations=len(rec),
                         note="obs are trailing 24-h running means; cams24 is the matched mean")}

    # ---- A. per-monitor statistics --------------------------------------------------------
    st = []
    for sid, g in rec.items():
        name, lat, lon, typ = STATIONS[sid]
        o = g.obs.dropna()
        m = g[["obs", "cams24"]].dropna()
        err = m.cams24 - m.obs
        st.append(dict(
            sid=sid, name=name, lat=lat, lon=lon, type=typ, cell=list(cell_of(lat, lon)),
            n=int(len(o)), obs_mean=float(o.mean()), obs_median=float(o.median()),
            frac_above_guideline=float((o > WHO24).mean()),
            n_matched=int(len(m)), obs_mean_matched=float(m.obs.mean()),
            cams_mean=float(m.cams24.mean()), bias=float(err.mean()),
            bias_median=float(err.median()),
            mae=float(err.abs().mean()), rmse=float(np.sqrt((err ** 2).mean())),
            pearson_r=float(np.corrcoef(m.cams24, m.obs)[0, 1]),
            spearman_r=float(stats.spearmanr(m.cams24, m.obs).statistic),
            bias_unmatched=float((g.cams - g.obs).mean()),
            pearson_r_unmatched=float(g.obs.corr(g.cams)),
        ))
    S = pd.DataFrame(st)
    res["stations"] = st
    sp = stats.spearmanr(S.obs_mean_matched, S.cams_mean)
    res["network"] = dict(
        n_obs=int(S.n.sum()), n_matched=int(S.n_matched.sum()),
        obs_mean_min=float(S.obs_mean.min()), obs_mean_max=float(S.obs_mean.max()),
        n_station_means_above_guideline=int((S.obs_mean > WHO24).sum()),
        frac_above_min=float(S.frac_above_guideline.min()),
        frac_above_max=float(S.frac_above_guideline.max()),
        cams_mean_min=float(S.cams_mean.min()), cams_mean_max=float(S.cams_mean.max()),
        bias_min=float(S.bias.min()), bias_max=float(S.bias.max()),
        bias_station_mean=float(S.bias.mean()), bias_station_median=float(S.bias.median()),
        mae_median=float(S.mae.median()), rmse_median=float(S.rmse.median()),
        pearson_r_median=float(S.pearson_r.median()),
        pearson_r_min=float(S.pearson_r.min()), pearson_r_max=float(S.pearson_r.max()),
        spearman_r_median=float(S.spearman_r.median()),
        pearson_r_unmatched_median=float(S.pearson_r_unmatched.median()),
        obs_between_station_sd=float(S.obs_mean_matched.std(ddof=1)),
        cams_between_station_sd=float(S.cams_mean.std(ddof=1)),
        spearman_station_means=float(sp.statistic), spearman_station_means_p=float(sp.pvalue),
        bias_vs_unmatched_spearman=float(stats.spearmanr(S.bias, S.bias_unmatched).statistic),
        bias_vs_unmatched_max_abs_diff=float((S.bias - S.bias_unmatched).abs().max()),
    )
    allobs = pd.concat([g.obs.dropna() for g in rec.values()])
    mon = allobs.groupby(allobs.index.to_period("M")).size()
    res["coverage"] = {str(k): int(v) for k, v in mon.items()}
    res["coverage_share_sep_dec_2023"] = float(
        mon[[p for p in mon.index if str(p) in ("2023-09", "2023-10", "2023-11", "2023-12")]]
        .sum() / mon.sum())
    res["nn_distance_median_km"] = float(np.median([
        min(haversine(a.lat, a.lon, b.lat, b.lon) for b in S.itertuples() if b.sid != a.sid)
        for a in S.itertuples()]))

    D = np.array([[haversine(a.lat, a.lon, b.lat, b.lon) for b in S.itertuples()]
                  for a in S.itertuples()])
    res["moran"] = {f"<= {b} km": morans_i(S.bias.values, D, b, rng) for b in (15, 25, 50, 100)}
    metro = S[S.name.isin(METRO)].bias
    alain = S[S.name.isin(ALAIN)].bias
    res["regions"] = dict(metro_min=float(metro.min()), metro_max=float(metro.max()),
                          alain_min=float(alain.min()), alain_max=float(alain.max()),
                          min_distance_km=float(min(
                              haversine(a.lat, a.lon, b.lat, b.lon)
                              for a in S[S.name.isin(METRO)].itertuples()
                              for b in S[S.name.isin(ALAIN)].itertuples())))

    # ---- B. sub-grid spread ----------------------------------------------------------------
    cells = []
    for key, g in S.groupby(S.cell.apply(tuple)):
        cells.append(dict(
            cell_lat=key[0], cell_lon=key[1], n=int(len(g)), stations=list(g.name),
            obs_mean=float(g.obs_mean.mean()), cams_mean=float(g.cams_mean.mean()),
            offset=float(g.cams_mean.mean() - g.obs_mean_matched.mean()),
            obs_range=float(g.obs_mean.max() - g.obs_mean.min()),
            obs_sd=float(g.obs_mean.std(ddof=1)) if len(g) > 1 else None,
            bias_closest=float(g.bias.abs().min()),
        ))
    res["cells"] = sorted(cells, key=lambda c: -c["n"])
    multi = [c for c in cells if c["n"] > 1]
    grand = S.obs_mean.mean()
    ss_tot = float(((S.obs_mean - grand) ** 2).sum())
    ss_within = float(sum(((g.obs_mean - g.obs_mean.mean()) ** 2).sum()
                          for _, g in S.groupby(S.cell.apply(tuple))))
    res["subgrid"] = dict(
        n_cells_occupied=len(cells), n_cells_multi=len(multi),
        n_stations_sharing=int(sum(c["n"] for c in multi)),
        within_cell_sd_median=float(np.median([c["obs_sd"] for c in multi])),
        within_cell_sd_max=float(max(c["obs_sd"] for c in multi)),
        within_cell_range_max=float(max(c["obs_range"] for c in multi)),
        share_of_station_mean_variance_within_cells=ss_within / ss_tot,
    )
    # the same partition on station medians, which a single episode cannot dominate
    gm = S.obs_median.mean()
    res["subgrid"]["share_of_station_median_variance_within_cells"] = float(
        sum(((g.obs_median - g.obs_median.mean()) ** 2).sum()
            for _, g in S.groupby(S.cell.apply(tuple))) / ((S.obs_median - gm) ** 2).sum())
    res["subgrid"]["within_cell_median_range"] = {
        f"{k[0]:.2f},{k[1]:.2f}": float(g.obs_median.max() - g.obs_median.min())
        for k, g in S.groupby(S.cell.apply(tuple)) if len(g) > 1}
    # share of the matched record in September-December 2023
    idx = pd.concat([pd.Series(1, index=g.index[g.obs.notna() & g.cams24.notna()])
                     for g in rec.values()])
    res["coverage_matched_share_sep_dec_2023"] = float(
        ((idx.index >= "2023-09-01") & (idx.index < "2024-01-01")).mean())
    res["coverage_matched"] = {k: int(v) for k, v in
                               idx.groupby(idx.index.strftime("%Y-%m")).size().items()}

    # ---- C. operator typology against bias -------------------------------------------------
    src = S[S.type.isin(SOURCE_TYPES)].bias.values
    bg = S[~S.type.isin(SOURCE_TYPES)].bias.values
    kw = stats.kruskal(*[g.bias.values for _, g in S.groupby("type")])
    mw = stats.mannwhitneyu(src, bg, alternative="less")
    pooled = np.concatenate([src, bg])
    obs_diff = src.mean() - bg.mean()
    perm = np.array([(lambda p: p[:len(src)].mean() - p[len(src):].mean())(
        rng.permutation(pooled)) for _ in range(200_000)])
    spool = math.sqrt(((len(src) - 1) * src.var(ddof=1) + (len(bg) - 1) * bg.var(ddof=1))
                      / (len(src) + len(bg) - 2))
    res["typology"] = dict(
        by_type={t: dict(n=int(len(g)), median=float(g.bias.median()),
                         min=float(g.bias.min()), max=float(g.bias.max()))
                 for t, g in S.groupby("type")},
        kruskal_H=float(kw.statistic), kruskal_p=float(kw.pvalue),
        source_n=int(len(src)), source_median=float(np.median(src)),
        background_n=int(len(bg)), background_median=float(np.median(bg)),
        mannwhitney_p_one_sided=float(mw.pvalue),
        permutation_p_one_sided=float((np.sum(perm <= obs_diff) + 1) / (len(perm) + 1)),
        cohens_d=float(obs_diff / spool),
        note="two-group contrast formed after the five-type omnibus test; report both",
    )

    # ---- D. pairwise coherence --------------------------------------------------------------
    pairs = []
    for a, b in itertools.combinations(sids, 2):
        j = pd.concat([rec[a].obs.rename("a"), rec[b].obs.rename("b")], axis=1, sort=True).dropna()
        if len(j) < 100:
            continue
        na, la, lo, _ = STATIONS[a]
        nb, lb, lob, _ = STATIONS[b]
        dj = pd.concat([daily(rec[a]).obs.rename("a"), daily(rec[b]).obs.rename("b")],
                       axis=1, sort=True).dropna()
        pairs.append(dict(
            a=na, b=nb, dist_km=haversine(la, lo, lb, lob),
            same_cell=bool(cell_of(la, lo) == cell_of(lb, lob)), n_common=int(len(j)),
            spearman_obs=float(stats.spearmanr(j.a, j.b).statistic),
            median_abs_diff=float((j.a - j.b).abs().median()),
            frac_opposite_sides=float(((j.a > WHO24) != (j.b > WHO24)).mean()),
            n_days=int(len(dj)),
            n_days_opposite=int(((dj.a > WHO24) != (dj.b > WHO24)).sum()),
        ))
    P = pd.DataFrame(pairs)
    res["pairs"] = pairs

    def summ(x):
        return dict(n_pairs=int(len(x)), dist_min=float(x.dist_km.min()),
                    dist_median=float(x.dist_km.median()), dist_max=float(x.dist_km.max()),
                    spearman_median=float(x.spearman_obs.median()),
                    spearman_min=float(x.spearman_obs.min()),
                    spearman_max=float(x.spearman_obs.max()),
                    median_abs_diff_median=float(x.median_abs_diff.median()),
                    frac_opposite_median=float(x.frac_opposite_sides.median()),
                    frac_opposite_max=float(x.frac_opposite_sides.max()),
                    n_days_total=int(x.n_days.sum()),
                    frac_days_opposite_pooled=float(x.n_days_opposite.sum()
                                                    / max(x.n_days.sum(), 1)))

    far = P[P.dist_km >= 80]
    res["coherence"] = dict(
        n_pairs=int(len(P)), same_cell=summ(P[P.same_cell]), different_cell=summ(P[~P.same_cell]),
        coherence_ge_80km_median=float(far.spearman_obs.median()), n_ge_80km=int(len(far)),
        cams_obs_spearman_median=float(S.spearman_r.median()),
        cams_obs_spearman_min=float(S.spearman_r.min()),
        cams_obs_spearman_max=float(S.spearman_r.max()),
        closest_pair=min(pairs, key=lambda p: p["dist_km"]),
    )

    # ---- E. regional level of the day ------------------------------------------------------
    obs_d = pd.DataFrame({s: daily(g).obs for s, g in rec.items()})
    cams_d = pd.DataFrame({s: daily(g).cams24 for s, g in rec.items()})
    ok = obs_d.notna().sum(axis=1) >= 10
    reg = obs_d[ok].median(axis=1)
    q1, q2 = reg.quantile([1 / 3, 2 / 3])
    cls = pd.Series(np.where(reg <= q1, "low", np.where(reg <= q2, "mid", "high")),
                    index=reg.index)
    name = {s: STATIONS[s][0] for s in sids}
    cellmap = {}
    for s in sids:
        cellmap.setdefault(cell_of(STATIONS[s][1], STATIONS[s][2]), []).append(s)
    multi_cells = [v for v in cellmap.values() if len(v) > 1]
    classes = {}
    for k in ("low", "mid", "high"):
        days = cls.index[cls == k]
        o, c = obs_d.loc[days], cams_d.loc[days]
        bias = c - o
        wsd = [v for members in multi_cells
               for v in o[members].dropna(thresh=2).std(axis=1, ddof=1).dropna()]
        classes[k] = dict(
            n_days=int(len(days)), regional_median=float(reg[days].median()),
            regional_min=float(reg[days].min()), regional_max=float(reg[days].max()),
            bias_all=float(np.nanmean(bias.values)),
            bias_metro=float(np.nanmean(bias[[s for s in sids if name[s] in METRO]].values)),
            bias_al_ain=float(np.nanmean(bias[[s for s in sids if name[s] in ALAIN]].values)),
            between_monitor_sd_median=float(o.std(axis=1, ddof=1).median()),
            cams_between_monitor_sd_median=float(c.std(axis=1, ddof=1).median()),
            within_cell_sd_median=float(np.median(wsd)),
        )
        classes[k]["metro_minus_al_ain"] = classes[k]["bias_metro"] - classes[k]["bias_al_ain"]
    creg = cams_d[ok].median(axis=1)
    both = pd.concat([reg, creg], axis=1, sort=True).dropna()
    rr = stats.spearmanr(both.iloc[:, 0], both.iloc[:, 1])
    res["regional"] = dict(
        n_days=int(ok.sum()), n_days_with_model=int(len(both)),
        first_day=str(reg.index.min())[:10], last_day=str(reg.index.max())[:10],
        tercile_bounds=[float(q1), float(q2)], classes=classes,
        spearman_regional_obs_vs_cams=float(rr.statistic),
        spearman_regional_obs_vs_cams_p=float(rr.pvalue),
        n_days_obs_above_150=int((both.iloc[:, 0] > 150).sum()),
        n_days_cams_above_150=int((both.iloc[:, 1] > 150).sum()),
        max_obs=float(both.iloc[:, 0].max()), max_cams=float(both.iloc[:, 1].max()),
        ols_slope_cams_on_obs=float(np.polyfit(both.iloc[:, 0], both.iloc[:, 1], 1)[0]),
        series=[dict(day=str(d)[:10], obs=float(o), cams=float(c), cls=str(cls[d]))
                for d, o, c in zip(both.index, both.iloc[:, 0], both.iloc[:, 1])],
    )

    # ---- G. stability across sub-periods (daily values) -----------------------------------
    periods = {"to_2023": (None, "2023-12-31"), "from_2024": ("2024-01-01", None),
               "sep_dec_2023": ("2023-09-01", "2023-12-31"),
               "outside_sep_dec_2023": None}
    full_bias = (cams_d - obs_d).mean()
    sub = {}
    for k, rg in periods.items():
        if rg is None:
            sel = ~((obs_d.index >= "2023-09-01") & (obs_d.index <= "2023-12-31"))
        else:
            sel = np.ones(len(obs_d), bool)
            if rg[0]:
                sel &= obs_d.index >= rg[0]
            if rg[1]:
                sel &= obs_d.index <= rg[1]
        o, c = obs_d[sel], cams_d[sel]
        b = (c - o).mean()
        okd = o.notna().sum(axis=1) >= 10
        rr_ = stats.spearmanr(o[okd].median(axis=1), c[okd].median(axis=1))
        sub[k] = dict(
            n_days=int(sel.sum()), n_monitor_days=int(o.notna().sum().sum()),
            bias_min=float(b.min()), bias_max=float(b.max()),
            bias_metro=float(np.nanmean((c - o)[[s for s in sids if name[s] in METRO]].values)),
            bias_al_ain=float(np.nanmean((c - o)[[s for s in sids if name[s] in ALAIN]].values)),
            spearman_bias_with_full_period=float(stats.spearmanr(b, full_bias).statistic),
            regional_spearman=float(rr_.statistic), regional_n_days=int(okd.sum()))
    res["subperiods"] = sub

    # ---- F. autocorrelation of the matched residual (sets the guard interval) ---------------
    lags = {"1 d": 8, "3 d": 24, "7 d": 56, "10 d": 80}
    acf = {k: [] for k in lags}
    for g in rec.values():
        r = g.obs - g.cams24
        for k, L in lags.items():
            v = r.corr(r.shift(L))
            if np.isfinite(v):
                acf[k].append(v)
    res["residual_acf_median"] = {k: float(np.median(v)) for k, v in acf.items()}

    os.makedirs(RES, exist_ok=True)
    for k, v in res["subperiods"].items():
        print(f"G  {k:22s} days {v['n_days']:3d} bias {v['bias_min']:+.0f}..{v['bias_max']:+.0f} "
              f"metro {v['bias_metro']:+.0f} Al Ain {v['bias_al_ain']:+.0f} | rank agreement with "
              f"full period {v['spearman_bias_with_full_period']:.2f} | regional rho "
              f"{v['regional_spearman']:.2f} ({v['regional_n_days']} d)")
    out = os.path.join(RES, "station_representativeness.json")
    json.dump(res, open(out, "w"), indent=1)

    n = res["network"]
    print(f"A  reported {n['n_obs']}, matched {n['n_matched']} | obs means "
          f"{n['obs_mean_min']:.1f}..{n['obs_mean_max']:.1f} "
          f"({n['n_station_means_above_guideline']}/19 above guideline; share of values above "
          f"{100*n['frac_above_min']:.0f}-{100*n['frac_above_max']:.0f}%)")
    print(f"   CAMS means {n['cams_mean_min']:.1f}..{n['cams_mean_max']:.1f} | bias "
          f"{n['bias_min']:+.1f}..{n['bias_max']:+.1f} (station mean "
          f"{n['bias_station_mean']:+.1f}) | MAE median {n['mae_median']:.1f} RMSE median "
          f"{n['rmse_median']:.1f} | r median {n['pearson_r_median']:.2f} "
          f"({n['pearson_r_min']:.2f}..{n['pearson_r_max']:.2f}); unmatched r median "
          f"{n['pearson_r_unmatched_median']:.2f}")
    print(f"   between-station SD obs {n['obs_between_station_sd']:.1f} CAMS "
          f"{n['cams_between_station_sd']:.1f} | station means obs vs CAMS Spearman "
          f"{n['spearman_station_means']:+.2f} (p={n['spearman_station_means_p']:.2f}) | "
          f"matched vs unmatched bias: Spearman {n['bias_vs_unmatched_spearman']:.3f}, max |diff| "
          f"{n['bias_vs_unmatched_max_abs_diff']:.1f}")
    print("   Moran: " + " | ".join(f"{k}: I={v['I']:.2f} (pairs {v['n_pairs']}, "
                                    f"p={v['p_perm']:.4f})" for k, v in res["moran"].items()))
    r_ = res["regions"]
    print(f"   metro {r_['metro_min']:+.0f}..{r_['metro_max']:+.0f} | Al Ain "
          f"{r_['alain_min']:+.0f}..{r_['alain_max']:+.0f} | nearest metro-Al Ain "
          f"{r_['min_distance_km']:.0f} km | NN distance median "
          f"{res['nn_distance_median_km']:.1f} km | Sep-Dec 2023 share "
          f"{100*res['coverage_share_sep_dec_2023']:.1f}%")
    s = res["subgrid"]
    print(f"B  {s['n_cells_occupied']} cells, {s['n_stations_sharing']}/19 sharing | within-cell "
          f"SD median {s['within_cell_sd_median']:.1f} max {s['within_cell_sd_max']:.1f}, max "
          f"range {s['within_cell_range_max']:.1f} | within-cell variance share "
          f"{100*s['share_of_station_mean_variance_within_cells']:.1f}%")
    for c in res["cells"]:
        if c["n"] > 1:
            print(f"   {c['cell_lat']:.2f},{c['cell_lon']:.2f} n={c['n']} obs {c['obs_mean']:.0f}"
                  f" CAMS {c['cams_mean']:.0f} offset {c['offset']:+.0f} range "
                  f"{c['obs_range']:.0f} SD {c['obs_sd']:.1f} closest |bias| "
                  f"{c['bias_closest']:.0f} | {', '.join(c['stations'])}")
    t = res["typology"]
    print(f"C  KW p={t['kruskal_p']:.3f} | source median {t['source_median']:+.1f} vs background "
          f"{t['background_median']:+.1f} | MW p={t['mannwhitney_p_one_sided']:.4f} perm p="
          f"{t['permutation_p_one_sided']:.4f} d={t['cohens_d']:+.2f} | "
          + ", ".join(f"{k} {v['median']:+.0f}" for k, v in t["by_type"].items()))
    c = res["coherence"]
    for lab in ("same_cell", "different_cell"):
        x = c[lab]
        print(f"D  {lab:14s} {x['n_pairs']:3d} pairs, {x['dist_min']:.1f}-{x['dist_max']:.1f} km"
              f" | Spearman median {x['spearman_median']:.2f} | opposite sides: median "
              f"{100*x['frac_opposite_median']:.1f}% (max {100*x['frac_opposite_max']:.1f}%) of "
              f"values; {100*x['frac_days_opposite_pooled']:.1f}% of {x['n_days_total']} "
              f"pair-days")
    cp = c["closest_pair"]
    print(f"   >=80 km: {c['coherence_ge_80km_median']:.2f} (n={c['n_ge_80km']}) | CAMS vs own "
          f"monitor, matched: {c['cams_obs_spearman_median']:.2f} "
          f"({c['cams_obs_spearman_min']:.2f}..{c['cams_obs_spearman_max']:.2f}) | closest pair "
          f"{cp['a']}-{cp['b']} {cp['dist_km']:.1f} km Spearman {cp['spearman_obs']:.2f}")
    r = res["regional"]
    print(f"E  {r['n_days']} days ({r['first_day']}..{r['last_day']}), terciles "
          f"{r['tercile_bounds'][0]:.0f}/{r['tercile_bounds'][1]:.0f} | regional level obs vs "
          f"CAMS Spearman {r['spearman_regional_obs_vs_cams']:+.2f} "
          f"(n={r['n_days_with_model']})")
    for k, v in r["classes"].items():
        print(f"   {k:4s} {v['n_days']:3d} days, regional {v['regional_median']:.0f} | bias all "
              f"{v['bias_all']:+.0f}, metro {v['bias_metro']:+.0f}, Al Ain "
              f"{v['bias_al_ain']:+.0f} | between-monitor SD {v['between_monitor_sd_median']:.0f}"
              f" (CAMS {v['cams_between_monitor_sd_median']:.0f}) | within-cell SD "
              f"{v['within_cell_sd_median']:.0f}")
    print("F  matched residual ACF (median over monitors): "
          + ", ".join(f"{k} {v:+.2f}" for k, v in res["residual_acf_median"].items()))
    print("wrote", out)


if __name__ == "__main__":
    main()
