"""SI tables of Section S6 from results/revision_diagnostics.json.

Writes manuscript/tab_si_season.tex, tab_si_offset.tex, tab_si_interp.tex and, from
results/aeronet_dust_check.json, tab_si_aeronet.tex.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "..", "results", "revision_diagnostics.json")))
OUT = os.path.join(HERE, "..", "manuscript")
K = ("low", "mid", "high")


def f(x, nd=0, sign=True):
    s = f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"
    return s.replace("-", "$-$").replace("+", "$+$")


def season():
    A = R["A_season"]
    L = [r"\begin{tabular}{llrrrrr}", r"\toprule",
         r"Days & Tercile & Days & Regional level & Metropolitan & Al~Ain & Difference \\", r"\midrule"]
    comp = A["tercile_by_season"]
    for k in K:
        c = comp[k]
        L.append(rf"All days & {k} & {sum(c.values())} & \multicolumn{{4}}{{l}}{{cool {c['cool']}, "
                 rf"transition {c['transition']}, warm {c['warm']}}} \\")
    L.append(r"\midrule")
    for s, lab in (("cool", "Cool (Nov--Feb)"), ("warm", "Warm (May--Sep)")):
        for k in K:
            v = A["within_season_terciles"][s]["classes"][k]
            L.append(rf"{lab if k == 'low' else ''} & {k} & {v['n_days']} & {f(v['regional_median'], 0, False)} & "
                     rf"{f(v['bias_metro'])} & {f(v['bias_al_ain'])} & {f(v['metro_minus_al_ain'])} \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(OUT, "tab_si_season.tex"), "w").write("\n".join(L) + "\n")


def offset():
    B = R["B_offset"]
    lab = dict(pooled_all="Mean over monitor-days (main text)",
               pooled_excl_al_mafraq="As above, without Al~Mafraq",
               monitor_mean="Mean of per-monitor means",
               monitor_median="Median of per-monitor means")
    L = [r"\begin{tabular}{lrrr}", r"\toprule", r"Aggregation & Low & Mid & High \\", r"\midrule"]
    for k, v in B["variants"].items():
        L.append(rf"{lab[k]} & {f(v['low'])} & {f(v['mid'])} & {f(v['high'])} \\")
    lo = B["loo_range"]
    L.append(rf"One metropolitan monitor left out (range) & {f(lo['low'][0])} to {f(lo['low'][1])} & "
             rf"{f(lo['mid'][0])} to {f(lo['mid'][1])} & {f(lo['high'][0])} to {f(lo['high'][1])} \\")
    W = R["A_season"]["within_season_terciles"]
    for s, l in (("cool", "Terciles within the cool season"), ("warm", "Terciles within the warm season")):
        v = W[s]["classes"]
        L.append(rf"{l} & {f(v['low']['metro_minus_al_ain'])} & {f(v['mid']['metro_minus_al_ain'])} & "
                 rf"{f(v['high']['metro_minus_al_ain'])} \\")
    L += [r"\midrule", r"Daily regression on the regional level & Days & Slope & 95\% interval \\", r"\midrule"]
    nm = dict(all="All days", warm="Warm season", cool="Cool season",
              all_excl_al_mafraq="All days, without Al~Mafraq", warm_excl_al_mafraq="Warm season, without Al~Mafraq")
    for k, v in B["daily_regression"]["fits"].items():
        L.append(rf"{nm[k]} & {v['n_days']} & {f(v['slope'], 2)} & {f(v['ci95'][0], 2)} to {f(v['ci95'][1], 2)} \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(OUT, "tab_si_offset.tex"), "w").write("\n".join(L) + "\n")


def interp():
    C = R["C_interpolation"]
    L = [r"\begin{tabular}{lrr}", r"\toprule", r"Statistic (daily record) & Nearest grid point & Bilinear \\", r"\midrule"]
    n, b = C["nearest"], C["bilinear"]
    rows = [("Mean bias, range over monitors", f"{f(n['station_bias']['min'])} to {f(n['station_bias']['max'])}",
             f"{f(b['station_bias']['min'])} to {f(b['station_bias']['max'])}"),
            ("Mean bias, mean over monitors", f(n['station_bias']['mean']), f(b['station_bias']['mean'])),
            ("Model station mean, range", f"{n['cams_station_mean_range'][0]:.0f} to {n['cams_station_mean_range'][1]:.0f}",
             f"{b['cams_station_mean_range'][0]:.0f} to {b['cams_station_mean_range'][1]:.0f}"),
            ("Spearman, observed vs model station means", f(n['spearman_station_means'], 2), f(b['spearman_station_means'], 2)),
            ("Pearson $r$ in time, median over monitors", f(n['pearson_r_median'], 2, False), f(b['pearson_r_median'], 2, False)),
            ("Regional level, Spearman", f(n['regional_spearman'], 2, False), f(b['regional_spearman'], 2, False)),
            ("Model SD between monitors, median over days", f(n['cams_between_monitor_sd_median'], 0, False),
             f(b['cams_between_monitor_sd_median'], 0, False))]
    for k in K:
        rows.append((f"Metropolitan minus Al~Ain, {k} tercile", f(n['metro_minus_al_ain'][k]), f(b['metro_minus_al_ain'][k])))
    for r in rows:
        L.append(" & ".join(r) + r" \\")
    L += [r"\midrule",
          rf"\multicolumn{{3}}{{l}}{{Largest change of a monitor's mean bias: {C['max_abs_change_station_bias']:.0f}; "
          rf"Spearman of the two bias rankings: {C['spearman_bias_rank_nearest_vs_bilinear']:.2f}}} \\",
          r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(OUT, "tab_si_interp.tex"), "w").write("\n".join(L) + "\n")


def aeronet():
    A = json.load(open(os.path.join(HERE, "..", "results", "aeronet_dust_check.json")))
    B = A["by_coarse_aod_tercile"]
    L = [r"\begin{tabular}{lrrr}", r"\toprule",
         r"Tercile of coarse-mode AOD & Low & Mid & High \\", r"\midrule"]
    rows = [("Days", "n_days", 0, False), ("Coarse-mode AOD, median", "coarse_aod_median", 2, False),
            ("Fine-mode fraction, median", "fine_fraction_median", 2, False),
            ("Observed regional level, median", "regional_obs_median", 0, False),
            ("Modelled regional level, median", "regional_cams_median", 0, False),
            ("Model minus observed regional level, median", "regional_bias_median", 0, True),
            ("Bias, metropolitan monitors", "bias_metro", 0, True),
            ("Bias, Al~Ain monitors", "bias_al_ain", 0, True),
            ("Metropolitan minus Al~Ain", "metro_minus_al_ain", 0, True),
            ("Metropolitan minus Al~Ain, without Al~Mafraq", "metro_minus_al_ain_excl_al_mafraq", 0, True)]
    for lab, k, nd, sg in rows:
        L.append(lab + " & " + " & ".join(f(B[t][k], nd, sg) for t in K) + r" \\")
    X = A["regional_tercile_x_dust_tercile"]
    L += [r"\midrule", r"Regional tercile (rows) by coarse-AOD tercile & Low & Mid & High \\", r"\midrule"]
    for r in K:
        L.append(f"Regional {r} & " + " & ".join(str(X[r][t]) for t in K) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(OUT, "tab_si_aeronet.tex"), "w").write("\n".join(L) + "\n")


def dustwind():
    W = json.load(open(os.path.join(HERE, "..", "results", "dust_wind_analysis.json")))
    L = [r"\begin{tabular}{lrrr}", r"\toprule",
         r"Stratum & Days & Median of stratifier & Metropolitan minus Al~Ain (without Al~Mafraq) \\", r"\midrule"]
    for key, lab, nd in (("by_dust_aod_tercile", "Dust AOD", 2), ("by_wind_speed_tercile", "10-m wind speed (m\\,s$^{-1}$)", 1)):
        for k in K:
            v = W[key][k]
            L.append(f"{lab}, {k} tercile & {v['n_days']} & {v['median']:.{nd}f} & {f(v['off_median'])} ({f(v['off_excl_al_mafraq_median'])}) \\\\")
        L.append(r"\midrule")
    for k, lab in (("NW_270_360", "Wind from 270--360$^{\\circ}$"), ("other", "Wind from other directions")):
        v = W["by_sector"][k]
        L.append(f"{lab} & {v['n_days']} & -- & {f(v['off_median'])} ({f(v['off_excl_al_mafraq_median'])}) \\\\")
    S = W["spearman"]; C = W["aeronet_consistency"]
    L += [r"\midrule", r"Spearman correlation & Days & \multicolumn{2}{r}{$\rho$ (without Al~Mafraq)} \\", r"\midrule",
          f"Difference vs dust AOD & {W['_meta']['n_days']} & \\multicolumn{{2}}{{r}}{{{f(S['off_vs_dust_aod']['rho'], 2)} ({f(S['off_excl_vs_dust_aod']['rho'], 2)})}} \\\\",
          f"Difference vs wind speed & {W['_meta']['n_days']} & \\multicolumn{{2}}{{r}}{{{f(S['off_vs_wind_speed']['rho'], 2)} ({f(S['off_excl_vs_wind_speed']['rho'], 2)})}} \\\\",
          f"Regional level vs dust AOD & {W['_meta']['n_days']} & \\multicolumn{{2}}{{r}}{{{f(S['reg_vs_dust_aod']['rho'], 2)}}} \\\\",
          f"CAMS dust AOD vs AERONET coarse-mode AOD & {C['n_days']} & \\multicolumn{{2}}{{r}}{{{f(C['cams_dust_aod_vs_aeronet_coarse_aod']['rho'], 2)}}} \\\\",
          f"CAMS dust share vs AERONET fine-mode fraction & {C['n_days']} & \\multicolumn{{2}}{{r}}{{{f(C['cams_dust_share_vs_aeronet_fine_fraction']['rho'], 2)}}} \\\\",
          r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(OUT, "tab_si_dustwind.tex"), "w").write("\n".join(L) + "\n")


def metar():
    M = json.load(open(os.path.join(HERE, "..", "results", "metar_dust_check.json")))
    G = M["groups"]
    L = [r"\begin{tabular}{lrrrrr}", r"\toprule",
         r"Days & $N$ & Regional level & Model minus observed & Metropolitan minus Al~Ain & (without Al~Mafraq) \\",
         r"\midrule"]
    for k, lab in (("dust_report", "Dust or sand reported"), ("no_dust_report", "No dust reported"),
                   ("clear", "Clear (no obscuration)")):
        v = G[k]
        L.append(f"{lab} & {v['n_days']} & {v['regional_obs_median']:.0f} & {f(v['regional_bias_median'])} & "
                 f"{f(v['local_deficit_median'])} & {f(v['local_deficit_excl_al_mafraq_median'])} \\\\")
    L += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(OUT, "tab_si_metar.tex"), "w").write("\n".join(L) + "\n")


season(); offset(); interp(); aeronet(); dustwind(); metar()
print("written")
