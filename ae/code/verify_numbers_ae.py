"""Check the numbers quoted in the manuscript, the SI and the cover letter against the result files.

Each entry rebuilds a phrase from results/*.json; the phrase must occur in the text. Then every
number of the main text that no entry accounts for is listed for inspection. Exit 1 on failure.
"""
import json
import math
import os
import re
import sys

from common_ae import RES, bh

MS = os.path.join(os.path.dirname(RES.rstrip("/")), "manuscript")
J = lambda n: json.load(open(os.path.join(RES, n)))
R, A, E, F, M = (J("station_representativeness.json"), J("audit_matched.json"),
                 J("exceedance_detection.json"), J("averaging_fingerprint.json"),
                 J("feature_sets_loso.json"))
AI, G3, AC, X, CC, QC = (J("audit_instantaneous.json"), J("audit_matched_guard3.json"),
                         J("archive_checks.json"), J("extended_observations.json"),
                         J("cams_cells_check.json"), J("qc_flatline_full.json"))
V13I = J(os.path.join("v13period", "audit_instantaneous_subsample.json"))
V13M = J(os.path.join("v13period", "audit_matched.json"))
FILES = ("body_ae.tex", "main_ae.tex", "SI_ae.tex", "cover_letter_ae.tex")
text = {f: open(os.path.join(MS, f)).read() for f in FILES}
norm = lambda s: re.sub(r"\s+", " ", s.replace("{,}", ","))
ALLN = norm("\n".join(text.values()))


def q(a):
    ks = [k for k, v in a.items() if isinstance(v, dict) and "binom_p" in v]
    return dict(zip(ks, bh([a[k]["binom_p"] for k in ks])))


QA, QI, QG = q(A), q(AI), q(G3)
N, S, CO, RG, CL = R["network"], R["subgrid"], R["coherence"], R["regional"], R["regional"]["classes"]
st = {s["name"]: s for s in R["stations"]}
cell = {(round(c["cell_lat"], 2), round(c["cell_lon"], 2)): c for c in R["cells"]}
sg = lambda x, d=2: f"${x:+.{d}f}$"
th = lambda n: f"{n:,}"
pc = lambda d: f"{100 * (1 - math.exp(-d)):.0f}"
H = AC["hourly_signature"]
SP = R["subperiods"]
checks = []
chk = lambda ph: checks.append(norm(ph))

# --- data and running mean
chk(f"{th(AC['_meta']['rows_by_parameter']['pm10'])} hourly \\pmten{{}} values")
chk(f"{th(N['n_obs'])} values, between {th(min(s['n'] for s in st.values()))} and {th(max(s['n'] for s in st.values()))} per monitor")
D = X["all"]
chk(f"{th(D['n_station_days'])} monitor-days on {D['n_days']} days")
chk(f"The {th(AC['record_vs_archive']['n_values'])} three-hourly values up to December 2023")
assert AC["record_vs_archive"]["identical_to_archive_value_at_same_time"] == AC["record_vs_archive"]["n_values"]
chk(f"median of {H['pm10']['acf_1h_median']:.3f} across the monitors")
chk(f"{100*H['pm10']['frac_unchanged_median']:.0f}\\% of them equal the value one hour earlier")
gas = [H[k]["acf_1h_median"] for k in ("no2", "so2", "o3", "co")]
chk(f"autocorrelations of {min(gas):.2f} to {max(gas):.2f}")
chk(f"negative at all {H['pm10']['n_stations_negative_at_24h']} monitors (median {sg(H['pm10']['increment_acf_24h_median'])}, range {sg(H['pm10']['increment_acf_24h_min'])} to {sg(H['pm10']['increment_acf_24h_max'])})")
g24 = [H[k]["increment_acf_24h_median"] for k in ("no2", "so2", "o3", "co")]
chk(f"positive ({sg(min(g24))} to {sg(max(g24))})")
chk(f"Pearson $r$ over the monitors {N['pearson_r_median']:.2f} against {N['pearson_r_unmatched_median']:.2f})")
chk(f"the correlation falls at {F['n_pm10_trailing_best']} of 19 monitors")
assert CC["missing_three_hourly_steps"] == 0 and CC["v13_values_compared"] == 20818
assert CC["v13_max_abs_diff"] < 2e-5
chk(f"all {th(CC['v13_values_compared'])} model values used in earlier work")
chk(f"at most {N['bias_vs_unmatched_max_abs_diff']:.0f}~\\ugm{{}}, but lowers")
chk(f"separations from {CO['same_cell']['dist_min']:.1f} to {CO['same_cell']['dist_max']:.0f}~km")
chk(f"({QC['removed']} values, {100*QC['removed_fraction']:.1f}\\%), which leaves {th(QC['kept'])} values")
ac = R["residual_acf_median"]
chk(f"autocorrelation {sg(ac['1 d'])} at one day, {sg(ac['3 d'])} at three days and {sg(ac['10 d'])} at ten days)")
chk(f"matched model counterpart: {th(E['raw']['n'])} monitor-days")
# --- 3.1
chk(f"from {N['obs_mean_min']:.0f} to {N['obs_mean_max']:.0f}~\\ugm{{}} across the network")
chk(f"between {100*N['frac_above_min']:.0f}\\% and {100*N['frac_above_max']:.0f}\\% of the archived 24-h means")
chk(f"Of the {th(D['n_station_days'])} monitor-days, {100*D['share_station_days_above']:.0f}\\% exceed")
bm = X["by_month"]
winter = [100 * bm[m]["share_above"] for m in ("2022-12", "2023-01", "2024-02")]
chk(f"from {min(winter):.0f}--{max(winter):.0f}\\% in the cleanest winter months")
assert all(bm[m]["share_above"] > 0.9 for m in ("2023-05", "2023-06", "2023-07", "2023-08", "2023-09"))
chk(f"between {N['cams_mean_min']:.0f} and {N['cams_mean_max']:.0f}~\\ugm{{}}")
chk(f"from ${N['bias_min']:.0f}$~\\ugm{{}} at Al~Mafraq to ${N['bias_max']:+.0f}$~\\ugm{{}} at Sweihan, with a station mean of ${N['bias_station_mean']:.0f}$~\\ugm{{}}")
chk(f"median absolute error is {N['mae_median']:.0f}~\\ugm{{}} and the median root-mean-square error {N['rmse_median']:.0f}~\\ugm{{}}")
chk(f"is {N['obs_between_station_sd']:.0f}~\\ugm{{}} in the observations and {N['cams_between_station_sd']:.0f}~\\ugm{{}} in the model")
chk(f"(Spearman $\\rho={N['spearman_station_means']:.2f}$, $p={N['spearman_station_means_p']:.2f}$)")
mo = R["moran"]
chk(f"is {mo['<= 15 km']['I']:.2f} for pairs within 15~km ($p={mo['<= 15 km']['p_perm']:.3f}$) and {mo['<= 100 km']['I']:.2f} for pairs within 100~km ($p<0.001$)")
assert mo["<= 100 km"]["p_perm"] < 0.001
rg = R["regions"]
chk(f"by {-rg['metro_max']:.0f} to {-rg['metro_min']:.0f}~\\ugm{{}}")
chk(f"within ${rg['alain_min']:.0f}$ to ${rg['alain_max']:+.0f}$~\\ugm{{}}")
assert rg["min_distance_km"] > 100
ty = R["typology"]
chk(f"by a median {-ty['source_median']:.0f}~\\ugm{{}} and suburban and rural-background monitors by a median {-ty['background_median']:.0f}~\\ugm{{}} (Mann--Whitney $p={ty['mannwhitney_p_one_sided']:.3f}$; permutation $p={ty['permutation_p_one_sided']:.3f}$; Cohen's $d={ty['cohens_d']:.2f}$)")
chk(f"(Kruskal--Wallis $p={ty['kruskal_p']:.2f}$)")
chk(f"Ruwais (industrial, ${st['Ruwais']['bias']:.0f}$~\\ugm{{}}), Al~Ain Street (traffic, ${st['Al Ain St.']['bias']:.0f}$~\\ugm{{}}) and Zakher (urban, ${st['Zakher']['bias']:+.0f}$~\\ugm{{}})")
ra = [SP[k]["spearman_bias_with_full_period"] for k in SP]
chk(f"(Spearman {min(ra):.2f} to {max(ra):.2f}; Section~S4)")
# --- 3.2
rngs = sorted(c["obs_range"] for c in R["cells"] if c["n"] > 1)
chk(f"agree to within {rngs[0]:.0f}--{rngs[3]:.0f}~\\ugm{{}}")
c1, c2 = cell[(24.0, 54.0)], cell[(24.0, 54.75)]
chk(f"differ by {c1['obs_range']:.0f}~\\ugm{{}} in the cell containing E11~Road")
chk(f"by {c2['obs_range']:.0f}~\\ugm{{}} in the cell containing Al~Mafraq")
chk(f"standard deviations of {c1['obs_sd']:.0f} and {c2['obs_sd']:.0f}~\\ugm{{}}")
chk(f"{100*S['share_of_station_mean_variance_within_cells']:.0f}\\% of the between-monitor variance of mean \\pmten{{}} lies within cells")
chk(f"its median ({st['E11 Road']['obs_median']:.0f}~\\ugm{{}})")
chk(f"within-cell share of the variance is {100*S['share_of_station_median_variance_within_cells']:.0f}\\%")
assert S["within_cell_median_range"]["24.00,54.75"] > 70
blocks = ["all", "to_2023", "sep_dec_2023", "outside_sep_dec_2023", "from_2024"]
assert all(X[b]["within_cell_range"]["24.00,54.75"] > 100 for b in blocks)
assert all(X[b]["n_station_means_above_guideline"] == 19 for b in blocks)
km, kd = "share_of_station_mean_variance_within_cells", "share_of_station_median_variance_within_cells"
chk(f"between {min(100*X[b][km] for b in blocks):.0f}\\% and {max(100*X[b][km] for b in blocks):.0f}\\% of the variance of station means, and between {min(100*X[b][kd] for b in blocks):.0f}\\% and {max(100*X[b][kd] for b in blocks):.0f}\\% of that of station medians")
ca, cb = cell[(24.75, 54.0)], cell[(24.75, 54.75)]
assert round(ca["obs_range"]) == round(cb["obs_range"]) == 1
chk(f"(ranges of {ca['obs_range']:.0f}~\\ugm{{}} each), yet the model lies {-ca['offset']:.0f} and {-cb['offset']:.0f}~\\ugm{{}} below them")
chk(f"within {c1['bias_closest']:.0f} and {c2['bias_closest']:.0f}~\\ugm{{}} of it")
sc = CO["same_cell"]
chk(f"in a median {100*sc['frac_opposite_median']:.0f}\\% of common values, and in up to {100*sc['frac_opposite_max']:.0f}\\%")
chk(f"{100*sc['frac_days_opposite_pooled']:.0f}\\% of {th(sc['n_days_total'])} pair-days")
lo_s = min(100 * min(S[km], S[kd]), min(100 * X[b][kd] for b in blocks))
hi_s = max(100 * X[b][km] for b in blocks)
chk(f"between {lo_s:.0f}\\% and {hi_s:.0f}\\% of the spatial variance")
chk(f"Between {lo_s:.0f}\\% and {hi_s:.0f}\\% of the between-monitor variance")
# --- 3.3
chk(f"(Spearman $\\rho={RG['spearman_regional_obs_vs_cams']:.2f}$ over {RG['n_days_with_model']}~days")
rs = [SP[k]["regional_spearman"] for k in SP]
chk(f"in each sub-period ($\\rho={min(rs):.2f}$ to {max(rs):.2f})")
chk(f"is {N['spearman_r_median']:.2f} (range {CO['cams_obs_spearman_min']:.2f} to {CO['cams_obs_spearman_max']:.2f}; median Pearson $r={N['pearson_r_median']:.2f}$)")
chk(f"(median {sc['spearman_median']:.2f} over the 15 pairs)")
chk(f"in different cells ({CO['different_cell']['spearman_median']:.2f};")
cp = CO["closest_pair"]
chk(f"{cp['dist_km']:.1f}~km apart with mean concentrations of {st['Khadeeja School']['obs_mean']:.0f} and {st['Hamdan St.']['obs_mean']:.0f}~\\ugm{{}}, correlate at only {cp['spearman_obs']:.2f}")
lo, mi, hi = CL["low"], CL["mid"], CL["high"]
chk(f"by {-lo['bias_metro']:.0f}~\\ugm{{}} in the lowest tercile (regional level {lo['regional_min']:.0f} to {lo['regional_max']:.0f}~\\ugm{{}}), by {-mi['bias_metro']:.0f}~\\ugm{{}} in the middle one and by {-hi['bias_metro']:.0f}~\\ugm{{}} in the highest ({hi['regional_min']:.0f} to {hi['regional_max']:.0f}~\\ugm{{}})")
assert round(mi["bias_al_ain"]) == 0
chk(f"move from ${lo['bias_al_ain']:+.0f}$ through 0 to ${hi['bias_al_ain']:.0f}$~\\ugm{{}}")
chk(f"(regression slope of the model network median on the observed one {RG['ols_slope_cams_on_obs']:.2f})")
chk(f"exceeds 150~\\ugm{{}} on {RG['n_days_cams_above_150']} of the {RG['n_days_with_model']} days, against {RG['n_days_obs_above_150']} in the observations")
mm = [CL[k]["metro_minus_al_ain"] for k in ("low", "mid", "high")]
chk(f"is ${mm[0]:.0f}$, ${mm[1]:.0f}$ and ${mm[2]:.0f}$~\\ugm{{}} in the three terciles")
k1, k2 = "between_monitor_sd_median", "cams_between_monitor_sd_median"
chk(f"median of {lo[k1]:.0f}, {mi[k1]:.0f} and {hi[k1]:.0f}~\\ugm{{}} in the three terciles in the observations, against {lo[k2]:.0f}, {mi[k2]:.0f} and {hi[k2]:.0f}~\\ugm{{}} in the model")
# --- 3.4
a0, b0 = A["buffer_full_R0"], A["buffer_coordfree_R0"]
chk(f"{a0['wins']} of 19 monitors with coordinates as predictors (median $\\Delta={a0['delta_median']:+.2f}$, cell-bootstrap interval ${a0['delta_median_ci_cell'][0]:+.2f}$ to ${a0['delta_median_ci_cell'][1]:+.2f}$) and at {b0['wins']} of 19 without them (${b0['delta_median']:+.2f}$, ${b0['delta_median_ci_cell'][0]:+.2f}$ to ${b0['delta_median_ci_cell'][1]:+.2f}$;")
chk(f"($p={a0['binom_p']:.2f}$ and {b0['binom_p']:.2f} before correction")
assert abs(A["buffer_full_R25"]["delta_median"]) < 0.005 and abs(A["buffer_coordfree_R25"]["delta_median"]) < 0.005
t0, u0 = A["buffer_time_full_R0"], A["buffer_time_coordfree_R0"]
chk(f"(medians ${t0['delta_median']:+.2f}$ and ${u0['delta_median']:+.2f}$, {t0['wins']} and {u0['wins']} of 19 monitors improved)")
assert min(QA.values()) > 0.99
tb = [G3[k]["delta_median"] for k in G3 if isinstance(G3[k], dict) and k.startswith("buffer_time")]
chk(f"medians between ${min(tb):.2f}$ and ${max(tb):+.2f}$ in the stages")
assert min(QG.values()) > 0.05
assert QC["kept"] == sum(v["n"] for v in A["buffer_full_R0"]["per_station"].values())
i0, j0 = AI["buffer_full_R0"], AI["buffer_coordfree_R0"]
chk(f"the same {th(QC['kept'])} archived values")
chk(f"is ${i0['delta_median']:+.2f}$ with coordinates and ${j0['delta_median']:+.2f}$ without, an apparent reduction of the error by {pc(i0['delta_median'])}\\% and {pc(j0['delta_median'])}\\% at {i0['wins']} and {j0['wins']} of 19 monitors ($q={QI['buffer_full_R0']:.3f}$ and {QI['buffer_coordfree_R0']:.3f})")
nest = sum(v < 0.05 for v in QI.values())
chk(f"{nest} of the 24 comparisons establish it")
k = AI["buffer_time_cell_coordfree_R100"]
chk(f"(${k['delta_median']:+.2f}$, {k['wins']} of 19 monitors, $q={QI['buffer_time_cell_coordfree_R100']:.3f}$")
chk(f"alone gains ${M['instantaneous']['cams']['median']:+.2f}$ under the instantaneous comparison and loses ${M['matched']['cams']['median']:.2f}$ under matched averaging")
assert abs(M["matched"]["cams+hour+doy"]["median"] - b0["delta_median"]) < 1e-9
assert abs(M["instantaneous"]["cams+hour+doy"]["median"] - j0["delta_median"]) < 1e-9
chk(f"reduction of {pc(i0['delta_median'])}--{pc(j0['delta_median'])}\\%, established at {j0['wins']} and {i0['wins']} of 19 monitors")
chk(f"apparent gain of {pc(i0['delta_median'])}--{pc(j0['delta_median'])}\\%")
chk(f"leave-one-station-out gain was not established even before anything else was withheld ({a0['wins']} of 19")
chk(f"Under the instantaneous treatment {nest} of the 24 comparisons establish a gain")
chk(f"gain of ${V13I['buffer_full_R0']['delta_median']:+.2f}$ to ${V13I['buffer_coordfree_R0']['delta_median']:+.2f}$ under the instantaneous treatment on the matched rows, and of ${V13M['buffer_full_R0']['delta_median']:+.2f}$ to ${V13M['buffer_coordfree_R0']['delta_median']:+.2f}$")
# --- 3.5
r, l, p = E["raw"], E["level_shifted"], E["persistence"]
chk(f"Of the {th(r['n'])} monitor-days evaluated, {th(r['events'])} ({100*r['base_rate']:.0f}\\%) exceed")
chk(f"detects {th(r['hit'])} of the {th(r['events'])} exceedances (POD {r['POD']:.2f}, 95\\% interval {r['POD_lo']:.2f} to {r['POD_hi']:.2f}) and reports exceedance on {r['fa']} of the {r['fa']+r['cn']} days that did not exceed (POFD {r['POFD']:.2f})")
chk(f"score is {r['ETS']:.2f} ({r['ETS_lo']:.2f} to {r['ETS_hi']:.2f}) and its Heidke skill score {r['HSS']:.2f} ({r['HSS_lo']:.2f} to {r['HSS_hi']:.2f})")
assert r["ETS_lo"] > 0
chk(f"raises the equitable threat score to {l['ETS']:.2f} ({l['ETS_lo']:.2f} to {l['ETS_hi']:.2f})")
chk(f"previous day reaches {p['ETS']:.2f} ({p['ETS_lo']:.2f} to {p['ETS_hi']:.2f})")
chk(f"exceedance score from {r['ETS']:.2f} to {l['ETS']:.2f}")
# --- abstract, highlights, conclusions
chk(f"({th(N['n_obs'])} three-hourly values; {th(r['n'])} monitor-days)")
chk(f"correlated with its model counterpart at {RG['spearman_regional_obs_vs_cams']:.2f} (Spearman)")
chk(f"(Spearman {RG['spearman_regional_obs_vs_cams']:.2f}) but not its spatial pattern")
chk(f"from ${N['bias_min']:.0f}$ to ${N['bias_max']:+.0f}$~\\ugm{{}} and was spatially organised")
chk(f"{100*S[km]:.0f}\\% lay inside single model cells ({100*S[kd]:.0f}\\% for medians)")
chk(f"(equitable threat score {r['ETS']:.2f})")
chk(f"regional level (Spearman {RG['spearman_regional_obs_vs_cams']:.2f})")
chk(f"from ${N['bias_min']:.0f}$ to ${N['bias_max']:+.0f}$~\\ugm{{}} and was organised in space")
# --- GATE-8 revision: season, aggregation, interpolation (revision_diagnostics.json)
RV = J("revision_diagnostics.json")
TS, WS = RV["A_season"]["tercile_by_season"], RV["A_season"]["within_season_terciles"]
chk(f"77 of the 107 days of the highest tercile")
assert TS["high"]["warm"] == 77 and sum(TS["high"].values()) == 107
assert TS["low"]["cool"] == 94 and sum(TS["low"].values()) == 109
chk("94 of the 109 days of the lowest")
wc, ww = WS["cool"]["classes"], WS["warm"]["classes"]
mn = lambda x: f"{x:+.0f}".replace("+", "$+$") if x > 0 else f"${x:.0f}$"
chk(f"from ${wc['low']['bias_metro']:.0f}$ to ${wc['high']['bias_metro']:.0f}$~\\ugm{{}} in the cool season and from ${ww['low']['bias_metro']:.0f}$ to ${ww['high']['bias_metro']:.0f}$~\\ugm{{}} in the warm one")
chk(f"Al~Ain bias from ${wc['low']['bias_al_ain']:+.0f}$ to ${wc['high']['bias_al_ain']:.0f}$ and from ${ww['low']['bias_al_ain']:+.0f}$ to ${ww['high']['bias_al_ain']:.0f}$~\\ugm{{}}")
B8 = RV["B_offset"]
lows = [v["low"] for v in B8["variants"].values()] + B8["loo_range"]["low"] + [wc["low"]["metro_minus_al_ain"], ww["low"]["metro_minus_al_ain"]]
chk(f"between ${-round(-max(lows)):.0f}$ and ${round(min(lows)):.0f}$~\\ugm{{}} however it is aggregated")
chk(f"${B8['variants']['pooled_excl_al_mafraq']['high']:.0f}$~\\ugm{{}} without Al~Mafraq")
chk(f"${B8['variants']['monitor_median']['high']:.0f}$~\\ugm{{}} with the median over monitors, and ${ww['high']['metro_minus_al_ain']:.0f}$~\\ugm{{}}")
fa, fx = B8["daily_regression"]["fits"]["all"], B8["daily_regression"]["fits"]["all_excl_al_mafraq"]
chk(f"slope of $-{abs(fa['slope']):.2f}$ (95\\% interval ${fa['ci95'][0]:.2f}$ to ${fa['ci95'][1]:+.2f}$) over all monitors and ${fx['slope']:+.2f}$ (${fx['ci95'][0]:.2f}$ to ${fx['ci95'][1]:+.2f}$)")
assert all(f["ci95"][0] < 0 < f["ci95"][1] for f in B8["daily_regression"]["fits"].values())   # "No slope differs from zero"
CI8 = RV["C_interpolation"]
chk(f"at most {CI8['max_abs_change_station_bias']:.0f}~\\ugm{{}}, leaves the ranking of the monitors by bias almost unchanged (Spearman {CI8['spearman_bias_rank_nearest_vs_bilinear']:.2f})")
chk(f"inversely related ($\\rho={CI8['bilinear']['spearman_station_means']:.2f}$)")
assert CI8["bilinear"]["cams_between_monitor_sd_median"] < CI8["nearest"]["cams_between_monitor_sd_median"]   # "narrows"
assert 49.5 <= -max(lows) and -min(lows) <= 62.5     # "about 50 to 60"
AE8 = J("aeronet_dust_check.json")
AB, AS, AX = AE8["by_coarse_aod_tercile"], AE8["spearman"], AE8["regional_tercile_x_dust_tercile"]
chk(f"On the {AE8['_meta']['n_days_matched']} days with a sun-photometer value")
chk(f"$\\rho={AS['coarse_aod_vs_obs_regional']['rho']:.2f}$ (with the modelled one at {AS['coarse_aod_vs_cams_regional']['rho']:.2f})")
chk(f"of the {sum(AX['high'].values())} days of the highest regional tercile among them, {AX['high']['high']} fall in the highest tercile")
assert AX["high"]["low"] == 0
chk("level is " + ", ".join(f"${AB[k]['regional_bias_median']:.0f}$" for k in ("low", "mid")) + f" and ${AB['high']['regional_bias_median']:.0f}$~\\ugm{{}} from the least")
chk("difference is " + ", ".join(f"${AB[k]['metro_minus_al_ain']:.0f}$" for k in ("low", "mid")) + f" and ${AB['high']['metro_minus_al_ain']:.0f}$~\\ugm{{}} (" + ", ".join(f"${AB[k]['metro_minus_al_ain_excl_al_mafraq']:.0f}$" for k in ("low", "mid")) + f" and")
lowb = AE8["low_dust_days_metro_bias_by_monitor"]
assert max(lowb.values()) < 0
chk(f"underestimated, by {-max(lowb.values()):.0f} to {-min(lowb.values()):.0f}~\\ugm{{}}")
chk(f"Its Spearman correlation with the observed regional level is {AS['coarse_aod_vs_obs_regional']['rho']:.2f}")
chk(f"fine-mode fraction correlates with the observed regional level at ${AS['fine_fraction_vs_obs_regional']['rho']:.2f}$")
DW = J("dust_wind_analysis.json")
DS, DC, DB = DW["spearman"], DW["aeronet_consistency"], DW["by_sector"]
chk(f"measured there (Spearman $\\rho={DC['cams_dust_aod_vs_aeronet_coarse_aod']['rho']:.2f}$ over {DC['n_days']} days)")
chk(f"dust optical depth over the metropolitan cells ($\\rho={DS['off_vs_dust_aod']['rho']:.2f}$; ${DS['off_excl_vs_dust_aod']['rho']:+.2f}$ without Al~Mafraq")
chk(f"10-m wind speed ($\\rho={DS['off_vs_wind_speed']['rho']:+.2f}$; ${DS['off_excl_vs_wind_speed']['rho']:+.2f}$ without Al~Mafraq)")
chk(f"prevails on {DW['nw_sector_days']} of {DW['_meta']['n_days']} days")
chk(f"(medians ${DB['NW_270_360']['off_median']:.0f}$ and ${DB['other']['off_median']:.0f}$~\\ugm{{}}, Mann--Whitney $p={DB['mann_whitney_p']:.2f}$)")
chk(f"a median {100*DW['seasalt_share_median']['metro']:.1f}\\% of the modelled optical depth over the metropolitan cells")
chk(f"AERONET coarse-mode optical depth at {DC['cams_dust_aod_vs_aeronet_coarse_aod']['rho']:.2f}, and its dust share of the optical depth with the AERONET fine-mode fraction at ${DC['cams_dust_share_vs_aeronet_fine_fraction']['rho']:.2f}$ ({DC['n_days']} days)")
chk(f"dust accounts for a median {100*DW['dust_share_median']['metro']:.0f}\\% of the modelled optical depth and sea salt for a median {100*DW['seasalt_share_median']['metro']:.1f}\\% (maximum {100*DW['seasalt_share_max']['metro']:.1f}\\%)")
chk(f"gives $p={DB['mann_whitney_p']:.2f}$")
assert DS["off_excl_vs_dust_aod"]["rho"] > 0 and DS["off_excl_vs_wind_speed"]["rho"] > 0   # "deficit smaller", "does not grow"
MT = J("metar_dust_check.json"); MG = MT["groups"]
chk(f"On the {MG['dust_report']['n_days']} days with a reported dust event the regional level has a median of {MG['dust_report']['regional_obs_median']:.0f}~\\ugm{{}}, against {MG['no_dust_report']['regional_obs_median']:.0f} on the other days, and the model falls short of it by {-MG['dust_report']['regional_bias_median']:.0f}~\\ugm{{}}, against {-MG['no_dust_report']['regional_bias_median']:.0f}")
chk(f"difference is ${MG['dust_report']['local_deficit_median']:.0f}$ against ${MG['no_dust_report']['local_deficit_median']:.0f}$~\\ugm{{}} ($p={MT['mann_whitney_local_deficit_dust_vs_other_p']:.2f}$)")
cb = MT["clear_day_metro_bias_by_monitor"]; assert max(cb.values()) < 0
chk(f"On the {MG['clear']['n_days']} clear days the difference is ${MG['clear']['local_deficit_median']:.0f}$~\\ugm{{}} and every metropolitan monitor is underestimated, by {-max(cb.values()):.0f} to {-min(cb.values()):.0f}~\\ugm{{}}")
chk(f"record dust on only {MG['dust_report']['n_days']} days")
chk(f"{MT['_meta']['n_reports']:,} reports")
chk(f"Dust or sand is reported {MT['_meta']['n_reports_dust']} times (DU {MT['_meta']['dust_codes_found']['DU']}, BLDU {MT['_meta']['dust_codes_found']['BLDU']}), on {MG['dust_report']['n_days']} days; haze is reported {MT['_meta']['n_reports_haze']} times")
chk(f"median coarse-mode optical depth at the sun photometer is {MT['aeronet']['coarse_aod_median_dust_report_days']:.2f}, against {MT['aeronet']['coarse_aod_median_other_days']:.2f}")
chk(f"correlates with the observed regional level at {MT['spearman']['obscured_vs_obs_regional']['rho']:.2f} and with the coarse-mode optical depth at {MT['aeronet']['obscured_vs_coarse_aod']['rho']:.2f}")
chk(f"reports has $p={MT['mann_whitney_local_deficit_dust_vs_other_p']:.2f}$ (Mann--Whitney), that in regional bias $p={MT['mann_whitney_regional_bias_dust_vs_other_p']:.2f}$")
chk(f"{MT['_meta']['n_days_with_record']} days of the daily record qualify")
# --- SI
m = AC["_meta"]
chk(f"6,803 daily files holding {th(m['rows'])} rows")
chk(f"({th(m['rows_by_parameter']['no2'])}), sulfur dioxide ({th(m['rows_by_parameter']['so2'])}), ozone ({th(m['rows_by_parameter']['o3'])}) and carbon monoxide ({th(m['rows_by_parameter']['co'])})")
chk(f"extends the same sampling to 20~May 2024 ({th(N['n_obs'])} values")
chk(f"every one of the {th(CC['n_times'])} three-hourly times")
chk(f"within ${CC['v13_max_abs_diff']*1e5:.1f}\\times10^{{-5}}$~\\ugm{{}} (median ${CC['v13_median_abs_diff']*1e6:.1f}\\times10^{{-6}}$)")
chk(f"sulfur dioxide at {19 - H['so2']['n_stations_negative_at_24h']} of 19")
chk(f"exceeds the instantaneous comparison at {sum(v['r_trailing'] > v['r_instant'] for v in F['pm10'].values())}")
chk(f"is {100*F['pm10_median']['diurnal_amp_rel']:.0f}\\% of the mean")
chk(f"observations over the monitors is {N['pearson_r_median']:.2f}; with the reanalysis taken at the same instant it is {N['pearson_r_unmatched_median']:.2f}")
chk(f"at most {N['bias_vs_unmatched_max_abs_diff']:.1f}~\\ugm{{}}")
chk(f"$\\rho={N['bias_vs_unmatched_spearman']:.3f}$")
mb = [SP[k]["bias_metro"] for k in SP]
chk(f"underestimated by {-max(mb):.0f} to {-min(mb):.0f}~\\ugm{{}} in each")
assert max(abs(SP[k]["bias_al_ain"]) for k in SP) <= 11.5

chk(f"understated ({N['pearson_r_unmatched_median']:.2f} instead of {N['pearson_r_median']:.2f})")
chk(f"median nearest-neighbour separation of {R['nn_distance_median_km']:.1f}~km")
bad = [ph for ph in checks if ph not in ALLN]
for ph in bad:
    print("MISSING:", ph)
print(f"{len(checks)} phrases checked, {len(bad)} missing")
body = norm(text["body_ae.tex"] + text["main_ae.tex"].split("\\begin{abstract}")[1].split("\\end{highlights}")[0])
masked = body
for ph in checks:
    masked = masked.replace(ph, " ")
masked = re.sub(r"\\(ref|label|citep|citet|cite)\{[^}]*\}", " ", masked)
left = sorted(set(re.findall(r"(?<![A-Za-z])[-+]?\d[\d,]*\.?\d*", masked)))
print("numbers outside checked phrases:", " ".join(left))
sys.exit(1 if bad else 0)
