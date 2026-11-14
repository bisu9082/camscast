"""Observation-only statistics on the full OpenAQ archive of the network (to 20 May 2024).

The model comparison is confined to the period for which the reanalysis was retrieved
(to 31 December 2023). The statistics that need no model -- level, exceedance, and the spread of
PM10 between monitors inside one model cell -- are recomputed here on everything the archive
holds. One value per monitor and day is used: the value reported at 08 UTC (noon local time),
which is the mean over the preceding 24 h and is the only hour the archive holds for parts of
spring and summer 2023.

Reads  data_raw/abudhabi_openaq_hourly.csv.gz
Writes results/extended_observations.json
"""
import itertools
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

from common_ae import RES, STATIONS, WHO24, cell_of, haversine

RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data_raw",
                   "abudhabi_openaq_hourly.csv.gz")
HOUR, MIN_COMMON = 8, 30


def variance_share(means, cells):
    m = pd.Series(means)
    tot = ((m - m.mean()) ** 2).sum()
    within = sum(((m[[s for s in m.index if cells[s] == c]]
                   - m[[s for s in m.index if cells[s] == c]].mean()) ** 2).sum()
                 for c in set(cells.values()))
    return float(within / tot)


def block(D, cells, label):
    """Statistics for a sub-table of days x monitors."""
    means = D.mean()
    out = dict(label=label, n_days=int(len(D)), n_station_days=int(D.notna().sum().sum()),
               first=str(D.index.min().date()), last=str(D.index.max().date()),
               station_mean_min=float(means.min()), station_mean_max=float(means.max()),
               n_station_means_above_guideline=int((means > WHO24).sum()),
               share_station_days_above=float((D > WHO24).sum().sum() / D.notna().sum().sum()),
               between_station_sd=float(means.std(ddof=1)),
               share_of_station_mean_variance_within_cells=variance_share(means, cells),
               share_of_station_median_variance_within_cells=variance_share(D.median(), cells),
               n_station_medians_above_guideline=int((D.median() > WHO24).sum()))
    rng, sds = {}, {}
    for c in sorted(set(cells.values())):
        mem = [s for s in D.columns if cells[s] == c]
        if len(mem) > 1:
            k = f"{c[0]:.2f},{c[1]:.2f}"
            rng[k] = float(means[mem].max() - means[mem].min())
            sds[k] = float(means[mem].std(ddof=1))
    out["within_cell_range"], out["within_cell_sd"] = rng, sds
    out["within_cell_range_max"] = max(rng.values())
    med = D.median()
    out["within_cell_median_range"] = {
        f"{c[0]:.2f},{c[1]:.2f}": float(med[[s for s in D.columns if cells[s] == c]].max()
                                        - med[[s for s in D.columns if cells[s] == c]].min())
        for c in sorted(set(cells.values())) if sum(cells[s] == c for s in D.columns) > 1}
    out["within_cell_sd_median"] = float(np.median(list(sds.values())))
    sc, dc, opp_n, opp_d = [], [], 0, 0
    for a, b in itertools.combinations(D.columns, 2):
        ok = D[a].notna() & D[b].notna()
        if ok.sum() < MIN_COMMON:
            continue
        rho = float(stats.spearmanr(D[a][ok], D[b][ok])[0])
        if cells[a] == cells[b]:
            sc.append(rho)
            opp_n += int(((D[a][ok] > WHO24) != (D[b][ok] > WHO24)).sum())
            opp_d += int(ok.sum())
        else:
            dc.append(rho)
    out.update(same_cell_pairs=len(sc), same_cell_spearman_median=float(np.median(sc)),
               different_cell_pairs=len(dc), different_cell_spearman_median=float(np.median(dc)),
               same_cell_pair_days=opp_d,
               same_cell_share_pair_days_opposite=float(opp_n / opp_d) if opp_d else None)
    # day-by-day spread
    sd_all = D.std(axis=1, ddof=1)[D.notna().sum(axis=1) >= 12]
    out["daily_between_monitor_sd_median"] = float(sd_all.median())
    wc = []
    for c in sorted(set(cells.values())):
        mem = [s for s in D.columns if cells[s] == c]
        if len(mem) > 1:
            x = D[mem]
            wc.append(x.std(axis=1, ddof=1)[x.notna().sum(axis=1) >= 2])
    out["daily_within_cell_sd_median"] = float(pd.concat(wc).median())
    return out


def main():
    d = pd.read_csv(RAW, parse_dates=["t"])
    d = d[(d.parameter == "pm10") & (d.value >= 0)]
    d["t"] = d.t.dt.tz_convert(None)
    d = d[d.t.dt.hour == HOUR]
    name = {int(k): v[0] for k, v in STATIONS.items()}
    D = d.pivot(index="t", columns="location_id", values="value").rename(columns=name)
    D.index = D.index.floor("D")
    cells = {v[0]: cell_of(v[1], v[2]) for v in STATIONS.values()}

    res = {"_meta": dict(hour_utc=HOUR, guideline=WHO24, min_common_days=MIN_COMMON,
                         note="one value per monitor and day: 24-h mean ending 08 UTC")}
    res["by_month"] = {m: dict(n_days=int(len(g)), n_station_days=int(g.notna().sum().sum()),
                               network_median=float(np.nanmedian(g.median(axis=1))),
                               share_above=float((g > WHO24).sum().sum() / g.notna().sum().sum()))
                       for m, g in D.groupby(D.index.strftime("%Y-%m"))}
    res["by_station"] = {s: dict(n_days=int(D[s].notna().sum()), mean=float(D[s].mean()),
                                 share_above=float((D[s].dropna() > WHO24).mean()))
                         for s in D.columns}
    core = (D.index >= "2023-09-01") & (D.index <= "2023-12-31")
    upto = D.index <= "2023-12-31"
    res["all"] = block(D, cells, "full archive, 29 Oct 2022 to 20 May 2024")
    res["to_2023"] = block(D[upto], cells, "to 31 Dec 2023 (period of the model comparison)")
    res["sep_dec_2023"] = block(D[core], cells, "September to December 2023")
    res["outside_sep_dec_2023"] = block(D[~core], cells, "all days outside September-December 2023")
    res["from_2024"] = block(D[~upto], cells, "1 January to 20 May 2024")
    # ranking of station means, inside vs outside the autumn-2023 core
    a, b = D[core].mean(), D[~core].mean()
    r = stats.spearmanr(a, b)
    res["station_mean_rank_agreement_core_vs_rest"] = dict(spearman=float(r[0]), p=float(r[1]))
    json.dump(res, open(os.path.join(RES, "extended_observations.json"), "w"), indent=1)

    for k in ("all", "to_2023", "sep_dec_2023", "outside_sep_dec_2023", "from_2024"):
        e = res[k]
        print(f"{k:22s} days {e['n_days']:3d} st-days {e['n_station_days']:5d} means "
              f"{e['station_mean_min']:.1f}-{e['station_mean_max']:.1f} "
              f"above {e['n_station_means_above_guideline']}/19 share>45 "
              f"{e['share_station_days_above']:.3f} within-cell var share "
              f"{e['share_of_station_mean_variance_within_cells']:.3f} (medians "
              f"{e['share_of_station_median_variance_within_cells']:.3f}) max range "
              f"{e['within_cell_range_max']:.1f} same-cell rho "
              f"{e['same_cell_spearman_median']:.2f} ({e['same_cell_pairs']}) diff "
              f"{e['different_cell_spearman_median']:.2f} opposite "
              f"{e['same_cell_share_pair_days_opposite']:.3f} of {e['same_cell_pair_days']} "
              f"daily SD all {e['daily_between_monitor_sd_median']:.1f} within "
              f"{e['daily_within_cell_sd_median']:.1f}")
    print(res["station_mean_rank_agreement_core_vs_rest"])
    for m, e in res["by_month"].items():
        print(m, e)


if __name__ == "__main__":
    main()
