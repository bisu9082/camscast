"""Graphical abstract for Atmospheric Research (at least 1328 px wide; 13 x 5 cm at 600 dpi).
Left: CAMS dust AOD against AERONET coarse-mode AOD in the cell of the sun photometer (the column).
Right: by tercile of the coarse-mode AOD, the regional bias of CAMS and the metropolitan minus
Al Ain difference (the surface).  Reads results/column_surface_daily.csv,
column_surface_check.json, aeronet_dust_check.json and dust_wind_analysis.json;
writes figures/graphical_abstract.{png,pdf}.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
D = pd.read_csv(os.path.join(RES, "column_surface_daily.csv"))
CS = json.load(open(os.path.join(RES, "column_surface_check.json")))
AE = json.load(open(os.path.join(RES, "aeronet_dust_check.json")))["by_coarse_aod_tercile"]
DW = json.load(open(os.path.join(RES, "dust_wind_analysis.json")))["aeronet_consistency"]
INK, INK2, MUTED, AXISC = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.edgecolor": AXISC,
                     "axes.linewidth": 0.6, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "pdf.fonttype": 42, "legend.frameon": False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13 / 2.54, 5 / 2.54),
                             gridspec_kw=dict(left=0.09, right=0.97, top=0.74, bottom=0.2,
                                              wspace=0.45, width_ratios=[1, 1.35]))
hi = 0.8
a1.plot([0, hi], [0, hi], color=AXISC, lw=0.8)
a1.scatter(D.tau_c, D.dewa_dust, s=5, c=INK2, edgecolors="none")
a1.set_xlim(0, hi); a1.set_ylim(0, hi); a1.set_aspect("equal")
a1.set_xticks([0, 0.4, 0.8]); a1.set_yticks([0, 0.4, 0.8])
a1.set_xlabel("AERONET coarse-mode AOD"); a1.set_ylabel("CAMS dust AOD")
rho, r = DW["cams_dust_aod_vs_aeronet_coarse_aod"]["rho"], CS["all"]["ratio_dust_to_coarse"]
a1.set_title(f"Column: timed (ρ = {rho:.2f}),\namount ×{r:.2f}", fontsize=7.5, color=INK, loc="left", pad=4)
for sp in ("top", "right"):
    a1.spines[sp].set_visible(False)

ks = ("low", "mid", "high"); x = np.arange(3); w = 0.36
reg = [AE[k]["regional_bias_median"] for k in ks]
loc = [AE[k]["metro_minus_al_ain"] for k in ks]
a2.bar(x - w / 2, reg, w, color=ORANGE, label="regional level, CAMS − observed")
a2.bar(x + w / 2, loc, w, color=BLUE, label="city minus inland (local deficit)")
a2.axhline(0, color=INK2, lw=0.6)
a2.set_xticks(x); a2.set_xticklabels(["least", "middle", "most"])
a2.set_xlabel("Dust tercile (AERONET coarse AOD)")
a2.set_ylim(-100, 5); a2.set_yticks([-60, -40, -20, 0])
a2.set_ylabel("PM$_{10}$ (µg m$^{-3}$)")
a2.set_title("Surface: regional deficit grows with dust;\nlocal deficit on all days",
             fontsize=7.5, color=INK, loc="left", pad=4)
a2.legend(loc="lower left", fontsize=5.5, handlelength=1.2, labelspacing=0.25)
for sp in ("top", "right"):
    a2.spines[sp].set_visible(False)
fig.text(0.01, 0.95, "CAMS reanalysis vs 19 Abu Dhabi PM$_{10}$ monitors and an AERONET sun photometer, 2022–2024",
         fontsize=6.8, color=INK, fontweight="bold")
fig.savefig(os.path.join(FIG, "graphical_abstract.png"), dpi=600)
fig.savefig(os.path.join(FIG, "graphical_abstract.pdf"))
print("graphical abstract written")
