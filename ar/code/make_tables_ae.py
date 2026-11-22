"""LaTeX tables of the Atmospheric Environment version, generated from the result artefacts.

Reads  results/station_representativeness.json, audit_matched.json, exceedance_detection.json,
       audit_instantaneous.json,
       averaging_fingerprint.json, feature_sets_loso.json, audit_matched_guard3.json,
       archive_checks.json, extended_observations.json
Writes manuscript/tab_*.tex
"""
import json
import os

from common_ae import RES, bh

OUT = os.path.join(os.path.dirname(RES.rstrip("/")), "manuscript")
J = lambda n: json.load(open(os.path.join(RES, n)))
R, A, E = J("station_representativeness.json"), J("audit_matched.json"), J("exceedance_detection.json")
AI, F, M = (J("audit_instantaneous.json"), J("averaging_fingerprint.json"),
            J("feature_sets_loso.json"))


def sg(x, d=1):
    return f"${x:+.{d}f}$"


def th(n):
    return f"{n:,}".replace(",", "{,}")


def write(name, lines):
    open(os.path.join(OUT, name), "w").write("\n".join(lines) + "\n")


def qvals(a):
    ks = [k for k, v in a.items() if isinstance(v, dict) and "binom_p" in v]
    return dict(zip(ks, bh([a[k]["binom_p"] for k in ks])))


# ---- Table 1: monitors
order = ["Industrial", "Traffic", "Urban", "Suburban", "Rural background"]
rows = sorted(R["stations"], key=lambda s: (order.index(s["type"]), s["bias"]))
L = [r"\begin{tabular}{llrrrrrrr}\toprule",
     r"Monitor & Type & Lat & Lon & $N$ & Observed & CAMS & Bias & $r$ \\ \midrule"]
for s in rows:
    L.append(f"{s['name']} & {s['type']} & {s['lat']:.2f} & {s['lon']:.2f} & {th(s['n'])} & "
             f"{s['obs_mean']:.1f} & {s['cams_mean']:.1f} & "
             f"{sg(s['bias'])} & {s['pearson_r']:.2f} \\\\")
write("tab_stations.tex", L + [r"\bottomrule\end{tabular}"])

# ---- Table 2: cells
cells = sorted([c for c in R["cells"] if c["n"] > 1], key=lambda c: (c["cell_lon"], c["cell_lat"]))
L = [r"\begin{tabular}{ccclrrrrr}\toprule",
     r"\multicolumn{2}{c}{Cell centre} & & & \multicolumn{3}{c}{Mean PM$_{10}$} & "
     r"\multicolumn{2}{c}{Between monitors} \\ \cmidrule(lr){1-2}\cmidrule(lr){5-7}\cmidrule(lr){8-9}",
     r"Lat & Lon & $n$ & Monitors & Observed & CAMS & Offset & Range & SD \\ \midrule"]
for c in cells:
    L.append(f"{c['cell_lat']:.2f} & {c['cell_lon']:.2f} & {c['n']} & "
             f"{', '.join(c['stations'])} & {c['obs_mean']:.0f} & {c['cams_mean']:.0f} & "
             f"{sg(c['offset'], 0)} & {c['obs_range']:.0f} & {c['obs_sd']:.1f} \\\\")
write("tab_cells.tex", L + [r"\bottomrule\end{tabular}"])

# ---- Table 3: correction under matched averaging
stages = [("buffer", "Monitors within $R$ withheld"),
          ("buffer_time", "+ time blocks with ten-day guard"),
          ("buffer_time_cell", "+ monitors sharing the model cell withheld")]


def ladder(a, name, stages=stages):
    q = qvals(a)
    L = [r"\begin{tabular}{rrcrrcrr}\toprule",
         r"& & \multicolumn{3}{c}{With coordinates} & \multicolumn{3}{c}{Coordinate-free} \\"
         r" \cmidrule(lr){3-5}\cmidrule(lr){6-8}",
         r"$R$ (km) & Sources & Median $\Delta$ [95\% CI] & Improved & $q$ & "
         r"Median $\Delta$ [95\% CI] & Improved & $q$ \\ \midrule"]
    for stg, lab in stages:
        L.append(rf"\multicolumn{{8}}{{l}}{{\emph{{{lab}}}}} \\")
        for Rk in (0, 25, 50, 100):
            cs = []
            for spec in ("full", "coordfree"):
                v = a[f"{stg}_{spec}_R{Rk}"]
                lo, hi = v["delta_median_ci_cell"]
                cs.append(f"${v['delta_median']:+.3f}$ [${lo:+.3f}$, ${hi:+.3f}$] & "
                          f"{v['wins']}/19 & {q[f'{stg}_{spec}_R{Rk}']:.3f}")
            L.append(f"{Rk} & {a[f'{stg}_full_R{Rk}']['n_src_mean']:.1f} & "
                     + " & ".join(cs) + r" \\")
    write(name, L + [r"\bottomrule\end{tabular}"])


ladder(A, "tab_transfer.tex")

# ---- Table 4: exceedance
names = [("raw", "CAMS reanalysis", "none"),
         ("level_shifted", "CAMS, level-shifted", "own monitor, preceding 8 days"),
         ("persistence", "Persistence", "own monitor, previous day")]
L = [r"\begin{tabular}{llrrcccc}\toprule",
     r"Classifier & Local information & Days & Exceeding & POD [95\% CI] & POFD & "
     r"ETS [95\% CI] & HSS [95\% CI] \\ \midrule"]
for k, lab, info in names:
    v = E[k]
    L.append(f"{lab} & {info} & {th(v['n'])} & {th(v['events'])} & "
             f"{v['POD']:.3f} [{v['POD_lo']:.3f}, {v['POD_hi']:.3f}] & {v['POFD']:.3f} & "
             f"${v['ETS']:.3f}$ [${v['ETS_lo']:.3f}$, ${v['ETS_hi']:.3f}$] & "
             f"${v['HSS']:.3f}$ [${v['HSS_lo']:.3f}$, ${v['HSS_hi']:.3f}$] \\\\")
write("tab_exceed.tex", L + [r"\bottomrule\end{tabular}"])

# ---- SI: fingerprint
L = [r"\begin{tabular}{lrrrrrrrr}\toprule",
     r"Series & $N$ & ACF(3 h) & Var.\ ratio & $r_\mathrm{inst}$ & $r_\mathrm{trail}$ & "
     r"$r_\mathrm{centre}$ & $r_\mathrm{lead}$ & Diurnal (\%) \\ \midrule"]
name = {s["sid"]: s["name"] for s in R["stations"]}
for sid, v in F["pm10"].items():
    L.append(f"{name[sid]} & {th(v['n'])} & {v['acf_3h']:.3f} & {v['var_ratio']:.3f} & "
             f"{v['r_instant']:.2f} & {v['r_trailing']:.2f} & {v['r_centred']:.2f} & "
             f"{v['r_leading']:.2f} & {100*v['diurnal_amp_rel']:.1f} \\\\")
write("tab_si_fingerprint.tex", L + [r"\bottomrule\end{tabular}"])

# ---- SI: coverage
cov = R["coverage"]
L = [r"\begin{tabular}{lr}\toprule Month & Reported values \\ \midrule"]
L += [f"{k} & {th(v)} \\\\" for k, v in cov.items()]
write("tab_si_coverage.tex", L + [r"\bottomrule\end{tabular}"])

# ---- SI: instantaneous against matched model value, same rows
qa, qi = qvals(A), qvals(AI)
L = [r"\begin{tabular}{llrcc}\toprule",
     r"Stage & Features & $R$ & Instantaneous model value & Matched 24-h mean \\ \midrule"]
for stg, lab in (("buffer", "buffer"), ("buffer_time", "+ time"),
                 ("buffer_time_cell", "+ time + cell")):
    for spec, sl in (("full", "with coordinates"), ("coordfree", "coordinate-free")):
        for Rk in (0, 25, 50, 100):
            k = f"{stg}_{spec}_R{Rk}"
            f = lambda a, q: f"${a[k]['delta_median']:+.3f}$, {a[k]['wins']}/19, {q[k]:.3f}"
            L.append(f"{lab} & {sl} & {Rk} & {f(AI, qi)} & {f(A, qa)} \\\\")
write("tab_si_three.tex", L + [r"\bottomrule\end{tabular}"])

# ---- SI: feature sets
L = [r"\begin{tabular}{lcc}\toprule",
     r"Predictors & Instantaneous CAMS & Matched 24-h mean \\ \midrule"]
lab = {"cams": "CAMS value", "cams+hour": "CAMS value, hour of day",
       "cams+doy": "CAMS value, day of year",
       "cams+hour+doy": "CAMS value, hour of day, day of year"}
for k, t in lab.items():
    a, b = M["instantaneous"][k], M["matched"][k]
    L.append(f"{t} & ${a['median']:+.3f}$, {a['wins']}/19 & ${b['median']:+.3f}$, {b['wins']}/19 \\\\")
write("tab_si_features.tex", L + [r"\bottomrule\end{tabular}"])

# ---- SI: three-day guard
G3, AC, X = J("audit_matched_guard3.json"), J("archive_checks.json"), J("extended_observations.json")
ladder(G3, "tab_si_guard3.tex",
       stages=[("buffer_time", "Time blocks with three-day guard"),
               ("buffer_time_cell", "+ monitors sharing the model cell withheld")])

# ---- SI: hourly signature of the reporting chain, by parameter
pl = {"pm10": r"PM$_{10}$", "no2": r"NO$_2$", "so2": r"SO$_2$", "o3": r"O$_3$", "co": "CO"}
L = [r"\begin{tabular}{lrrrrrrrr}\toprule",
     r"& & & & & \multicolumn{4}{c}{Autocorrelation of hourly increments} \\ \cmidrule(lr){6-9}",
     r"Species & Stations & Hourly values & ACF(1 h) & Unchanged & 1 h & 12 h & 24 h [range] & "
     r"Negative at 24 h \\ \midrule"]
for k, lab in pl.items():
    v = AC["hourly_signature"][k]
    L.append(f"{lab} & {v['n_stations']} & {th(v['n_values'])} & {v['acf_1h_median']:.3f} & "
             f"{100*v['frac_unchanged_median']:.1f}\\% & ${v['increment_acf_1h_median']:+.2f}$ & "
             f"${v['increment_acf_12h_median']:+.2f}$ & ${v['increment_acf_24h_median']:+.2f}$ "
             f"[${v['increment_acf_24h_min']:+.2f}$, ${v['increment_acf_24h_max']:+.2f}$] & "
             f"{v['n_stations_negative_at_24h']}/{v['n_stations']} \\\\")
write("tab_si_hourly.tex", L + [r"\bottomrule\end{tabular}"])

# ---- SI: months held by the analysed record and by the archive
L = [r"\begin{tabular}{lrrrrr}\toprule",
     r"& Analysed & \multicolumn{4}{c}{OpenAQ archive} \\ \cmidrule(lr){3-6}",
     r"Month & 3-hourly values & Hourly values & Days with 08 UTC value & Network median & "
     r"Above 45 (\%) \\ \midrule"]
for m, v in AC["pm10_by_month"].items():
    e = X["by_month"].get(m)
    L.append(f"{m} & {th(cov.get(m, 0)) if cov.get(m) else '--'} & {th(v['n'])} & "
             f"{v['n_days_with_08utc']} & {e['network_median']:.0f} & {100*e['share_above']:.0f} \\\\")
write("tab_si_coverage.tex", L + [r"\bottomrule\end{tabular}"])

# ---- SI: observation-only statistics on the full archive
cols = [("to_2023", "To 31 Dec 2023"), ("sep_dec_2023", "Sep--Dec 2023"),
        ("outside_sep_dec_2023", "Outside Sep--Dec 2023"), ("from_2024", "Jan--May 2024"),
        ("all", "Full archive")]
rowsx = [("Days", lambda e: f"{e['n_days']}"),
         ("Monitor-days", lambda e: th(e["n_station_days"])),
         (r"Monitor means, range (\ugm{})",
          lambda e: f"{e['station_mean_min']:.0f}--{e['station_mean_max']:.0f}"),
         ("Monitor means above 45", lambda e: f"{e['n_station_means_above_guideline']}/19"),
         (r"Monitor-days above 45 (\%)", lambda e: f"{100*e['share_station_days_above']:.1f}"),
         (r"Variance of monitor means within cells (\%)",
          lambda e: f"{100*e['share_of_station_mean_variance_within_cells']:.1f}"),
         (r"Variance of monitor medians within cells (\%)",
          lambda e: f"{100*e['share_of_station_median_variance_within_cells']:.1f}"),
         (r"Largest within-cell range of means (\ugm{})",
          lambda e: f"{e['within_cell_range_max']:.0f}"),
         (r"Day-by-day SD between monitors, median (\ugm{})",
          lambda e: f"{e['daily_between_monitor_sd_median']:.1f}"),
         (r"Day-by-day SD within cells, median (\ugm{})",
          lambda e: f"{e['daily_within_cell_sd_median']:.1f}"),
         ("Spearman, same-cell pairs (median)", lambda e: f"{e['same_cell_spearman_median']:.2f}"),
         ("Spearman, different-cell pairs (median)",
          lambda e: f"{e['different_cell_spearman_median']:.2f}"),
         (r"Same-cell pair-days on opposite sides of 45 (\%)",
          lambda e: f"{100*e['same_cell_share_pair_days_opposite']:.1f}")]
L = [r"\begin{tabular}{l" + "r" * len(cols) + r"}\toprule",
     " & " + " & ".join(c[1] for c in cols) + r" \\ \midrule"]
for lab, f in rowsx:
    L.append(lab + " & " + " & ".join(f(X[c[0]]) for c in cols) + r" \\")
write("tab_si_extended.tex", L + [r"\bottomrule\end{tabular}"])
print("wrote tables")

# ---- SI: model comparison by sub-period (daily record)
sp = R["subperiods"]
lab = {"to_2023": "Oct 2022--Dec 2023", "sep_dec_2023": "Sep--Dec 2023",
       "outside_sep_dec_2023": "Outside Sep--Dec 2023", "from_2024": "Jan--May 2024"}
L = [r"\begin{tabular}{lrrrrrr}\toprule",
     r"Sub-period & Days & Bias range & Metropolitan & Al Ain & Rank agreement & Regional $\rho$ \\ \midrule"]
for k, t in lab.items():
    v = sp[k]
    L.append(f"{t} & {v['n_days']} & ${v['bias_min']:+.0f}$ to ${v['bias_max']:+.0f}$ & "
             f"${v['bias_metro']:+.0f}$ & ${v['bias_al_ain']:+.0f}$ & "
             f"{v['spearman_bias_with_full_period']:.2f} & {v['regional_spearman']:.2f} \\\\")
write("tab_si_subperiods.tex", L + [r"\bottomrule\end{tabular}"])
print("wrote sub-period table")
