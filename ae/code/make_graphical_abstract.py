"""Graphical abstract (Atmospheric Environment: 531 x 1328 px or proportionally more, readable at
5 x 13 cm). Left: the regional daily level, observed against CAMS (when). Right: mean bias by
regional tercile for the Abu Dhabi metropolitan monitors and Al Ain, and their difference (where).
Reads results/station_representativeness.json; writes figures/graphical_abstract.{png,pdf}.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "..", "results", "station_representativeness.json")))
FIG = os.path.join(HERE, "..", "figures")
INK, INK2, MUTED, GRIDC, AXISC = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.edgecolor": AXISC,
                     "axes.linewidth": 0.6, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "pdf.fonttype": 42, "legend.frameon": False})
reg = R["regional"]
W, H = 13 / 2.54, 5 / 2.54
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W, H),
                             gridspec_kw=dict(left=0.09, right=0.95, top=0.80, bottom=0.2,
                                              wspace=0.55, width_ratios=[1, 1.3]))
o = np.array([d["obs"] for d in reg["series"]])
c = np.array([d["cams"] for d in reg["series"]])
hi = 270
a1.plot([0, hi], [0, hi], color=AXISC, lw=0.8)
a1.scatter(o, c, s=5, c=INK2, edgecolors="none")
a1.set_xlim(0, hi); a1.set_ylim(0, hi); a1.set_aspect("equal")
a1.set_xticks([0, 100, 200]); a1.set_yticks([0, 100, 200])
a1.set_xlabel("Observed daily median (µg m$^{-3}$)"); a1.set_ylabel("CAMS")
a1.set_title(f"When: tracked (ρ = {reg['spearman_regional_obs_vs_cams']:.2f})",
             fontsize=7.5, color=INK, loc="left", pad=4)
for sp in ("top", "right"):
    a1.spines[sp].set_visible(False)
ks = ("low", "mid", "high")
x = np.arange(3)
for key, lab, col, mk, ls in (("bias_metro", "Abu Dhabi conurbation", BLUE, "o", "-"),
                              ("bias_al_ain", "Al Ain (inland)", ORANGE, "s", "-"),
                              ("metro_minus_al_ain", "difference", INK2, "D", (0, (3, 2)))):
    y = [reg["classes"][k][key] for k in ks]
    a2.plot(x, y, color=col, lw=1.6 if ls == "-" else 1.0, ls=ls, marker=mk,
            ms=4 if ls == "-" else 3, label=lab)
a2.axhline(0, color=INK2, lw=0.6)
a2.set_xticks(x); a2.set_xticklabels(["clean", "", "dusty"])
a2.set_xlabel("Regional level of the day")
a2.set_ylim(-100, 45); a2.set_yticks([-80, -40, 0, 40])
a2.set_ylabel("Bias (µg m$^{-3}$)")
d = [reg["classes"][k]["metro_minus_al_ain"] for k in ks]
a2.set_title(f"Where: not ({min(-v for v in d):.0f}–{max(-v for v in d):.0f} µg m$^{{-3}}$ local offset)",
             fontsize=7.5, color=INK, loc="left", pad=4)
a2.legend(loc="upper right", fontsize=5.5, handlelength=1.8, labelspacing=0.25)
for sp in ("top", "right"):
    a2.spines[sp].set_visible(False)
fig.text(0.01, 0.95, "CAMS vs 19 Abu Dhabi PM$_{10}$ monitors, matched 24-h means, Oct 2022–May 2024",
         fontsize=7, color=INK, fontweight="bold")
fig.savefig(os.path.join(FIG, "graphical_abstract.png"), dpi=600)
fig.savefig(os.path.join(FIG, "graphical_abstract.pdf"))
print("graphical abstract written")
