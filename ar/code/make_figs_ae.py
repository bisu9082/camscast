"""Figures 1-4 of the Atmospheric Environment version.

Fig. 1  station map: mean CAMS bias (diverging colour), operator station type (marker shape),
        0.75 degree cell boundaries.
Fig. 2  (a) station-mean observed PM10 inside each multi-monitor cell against the single CAMS
        value; (b) coherence of observed series between monitor pairs against separation.
Fig. 3  what the reanalysis reproduces (the regional level from day to day) and what it does
        not (the local offset), by regional level of the day.
Fig. 4  gain of the statistical correction at a monitor treated as unmonitored, as information
        is removed from its training data (model and observations averaged over the same 24 h).

Reads  results/station_representativeness.json, results/audit_matched.json and
       data_aux/naturalearth_clip.json (Natural Earth 1:10m coastline and borders)
Writes figures/Fig1_map, Fig2_subgrid, Fig3_regional, Fig4_transfer as 600 dpi PNG and vector PDF

Authored at the final two-column width (190 mm = 7.48 in); smallest text 7 pt.
Identity is never carried by colour alone: station type is a marker shape, the two feature sets
differ in marker and colour, and established/unestablished gains differ in fill.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.environ.get("CAMSCAST_RESULTS", os.path.join(ROOT, "results"))
FIG = os.environ.get("CAMSCAST_FIGDIR", os.path.join(ROOT, "figures"))
os.makedirs(FIG, exist_ok=True)

INK, INK2, MUTED, GRIDC, AXISC = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, SURFACE = "#2a78d6", "#eb6834", "#ffffff"
W = 7.48
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": AXISC,
    "axes.linewidth": 0.7, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.labelsize": 8,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "savefig.facecolor": SURFACE,
    "figure.facecolor": SURFACE, "axes.grid": False, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def save(fig, name):
    """600 dpi PNG for review and a vector PDF with embedded fonts for production."""
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=600)
    fig.savefig(os.path.join(FIG, name + ".pdf"))
DIV = LinearSegmentedColormap.from_list(
    "bias", ["#0d366b", "#256abf", "#86b6ef", "#f0efec", "#f2a99c", "#d03b3b", "#7f1717"])
SHAPE = {"Industrial": "s", "Traffic": "^", "Urban": "D", "Suburban": "o",
         "Rural background": "v"}
GRID = 0.75

R = json.load(open(os.path.join(RES, "station_representativeness.json")))
A = json.load(open(os.path.join(RES, "audit_matched.json")))
NE = json.load(open(os.path.join(ROOT, "data_aux", "naturalearth_clip.json")))
ST = R["stations"]


def panel_letter(ax, s, dx=-0.02, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            color=INK, va="bottom", ha="right")


# ------------------------------------------------------------------ Fig. 1
def fig1():
    fig = plt.figure(figsize=(W, 3.75))
    ax = fig.add_axes([0.075, 0.13, 0.545, 0.80])
    axi = fig.add_axes([0.675, 0.455, 0.315, 0.475])
    cax = fig.add_axes([0.690, 0.315, 0.285, 0.028])
    norm = TwoSlopeNorm(vmin=-140, vcenter=0, vmax=140)

    def cells(a, x0, x1, y0, y1):
        for e in np.arange(-0.375 + GRID * np.floor(x0 / GRID), x1 + GRID, GRID):
            a.axvline(e, color=GRIDC, lw=0.8, zorder=1)
        for e in np.arange(-0.375 + GRID * np.floor(y0 / GRID), y1 + GRID, GRID):
            a.axhline(e, color=GRIDC, lw=0.8, zorder=1)

    def base(a):
        for seg in NE["coastline"]:
            xs, ys = zip(*seg)
            a.plot(xs, ys, color=INK2, lw=0.7, zorder=2)
        # national borders are not drawn (Atmospheric Research map policy: study area only)

    def pts(a, size):
        for s in ST:
            a.scatter(s["lon"], s["lat"], s=size, marker=SHAPE[s["type"]],
                      c=[DIV(norm(s["bias"]))], edgecolors=INK, linewidths=0.6, zorder=4)

    x0, x1, y0, y1 = 52.3, 56.15, 22.85, 24.95
    cells(ax, x0, x1, y0, y1)
    base(ax)
    pts(ax, 46)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
    ax.set_aspect(1 / np.cos(np.radians(24.0)))
    lab = {"Liwa": (7, 2), "Ruwais": (-4, 7), "Gayathi": (-10, -11), "Bida Zayed": (6, -10),
           "Habshan S.": (7, -3), "E11 Road": (7, 2), "Sweihan": (7, 2), "Al Qua'a": (7, -3)}
    for s in ST:
        if s["name"] in lab:
            ax.annotate(s["name"], (s["lon"], s["lat"]), xytext=lab[s["name"]],
                        textcoords="offset points", fontsize=7, color=INK2)
    ax.annotate("Al Ain (4 monitors)", (55.72, 24.22), xytext=(-8, -22),
                textcoords="offset points", fontsize=7, color=INK2, ha="center")
    bx0, bx1, by0, by1 = 54.22, 54.93, 24.235, 24.545
    ax.add_patch(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fill=False, ec=INK2,
                           lw=0.8, zorder=3))
    ax.annotate("Abu Dhabi\n(7 monitors, panel b)", (bx0, by1), xytext=(-4, 6),
                textcoords="offset points", fontsize=7, color=INK2, ha="center")
    ax.text(0.985, 0.03, "grey grid: 0.75° CAMS cells · dark line: coast",
            transform=ax.transAxes, fontsize=7, color=MUTED, ha="right",
            bbox=dict(fc=SURFACE, ec="none", pad=1.0, alpha=0.85))
    panel_letter(ax, "a", dx=-0.07)

    cells(axi, bx0, bx1, by0, by1)
    base(axi)
    pts(axi, 60)
    axi.set_xlim(bx0, bx1); axi.set_ylim(by0, by1)
    axi.set_aspect(1 / np.cos(np.radians(24.4)))
    axi.set_xticks([54.3, 54.5, 54.7]); axi.set_yticks([24.3, 24.4, 24.5])
    labi = {"Khadeeja School": (8, -10, "left"), "Hamdan St.": (8, 4, "left"),
            "Mussafah": (-7, -3, "right"), "Khalifa City": (8, -3, "left"),
            "Al Maqtaa": (-7, -3, "right"), "Baniyas School": (7, -1, "left"),
            "Al Mafraq": (-7, -6, "right")}
    for s in ST:
        if s["name"] in labi:
            dx, dy, ha = labi[s["name"]]
            axi.annotate(f"{s['name']} {s['bias']:+.0f}".replace("-", "−")
                         .replace("School", "Sch."),
                         (s["lon"], s["lat"]), xytext=(dx, dy), ha=ha,
                         textcoords="offset points", fontsize=7, color=INK2, zorder=5,
                         bbox=dict(fc=SURFACE, ec="none", pad=0.6, alpha=0.8))
    panel_letter(axi, "b", dx=-0.04)

    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=DIV), cax=cax,
                      orientation="horizontal", ticks=[-120, -60, 0, 60, 120])
    cb.ax.tick_params(labelsize=7)
    cb.outline.set_edgecolor(AXISC); cb.outline.set_linewidth(0.6)
    cb.set_label("Mean bias, CAMS − observed (µg m$^{-3}$)", fontsize=7.5, color=INK2)
    h = [Line2D([], [], marker=m, ls="", mfc="#f0efec", mec=INK, mew=0.6, ms=5.5, label=t)
         for t, m in SHAPE.items()]
    fig.legend(handles=h, loc="lower left", bbox_to_anchor=(0.685, 0.0), ncol=2,
               fontsize=7, handletextpad=0.2, columnspacing=1.2, labelcolor=INK2,
               title="Station type (operator classification)", title_fontsize=7,
               alignment="left")
    save(fig, "Fig1_map")
    plt.close(fig)


# ------------------------------------------------------------------ Fig. 2
def fig2():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W, 3.25),
                                 gridspec_kw=dict(left=0.07, right=0.985, top=0.93,
                                                  bottom=0.215, wspace=0.24))
    cells = sorted([c for c in R["cells"] if c["n"] > 1], key=lambda c: (c["cell_lon"],
                                                                         c["cell_lat"]))
    by = {s["name"]: s for s in ST}
    a1.axhline(R["_meta"]["who_24h_guideline"], color=INK2, lw=0.9, ls=(0, (4, 3)), zorder=2)
    a1.text(5.45, 45, "WHO 24-h\nguideline", fontsize=7, color=INK2, va="center", ha="right",
            bbox=dict(fc=SURFACE, ec="none", pad=1.0))
    for i, c in enumerate(cells):
        ys = sorted(by[n]["obs_mean"] for n in c["stations"])
        a1.plot([i, i], [ys[0], ys[-1]], color=AXISC, lw=1.4, zorder=2,
                solid_capstyle="round")
        off = np.linspace(-0.10, 0.10, len(c["stations"])) if c["n"] > 2 else [-0.07, 0.07]
        for o, n in zip(off, sorted(c["stations"], key=lambda n: by[n]["obs_mean"])):
            a1.scatter(i + o, by[n]["obs_mean"], marker=SHAPE[by[n]["type"]], s=34,
                       c=INK2, edgecolors=SURFACE, linewidths=0.8, zorder=4)
        a1.plot([i - 0.30, i + 0.30], [c["cams_mean"]] * 2, color=BLUE, lw=2.0, zorder=3,
                solid_capstyle="round")
    a1.set_xticks(range(len(cells)))
    a1.set_xticklabels([f"{c['cell_lat']:.2f}°N\n{c['cell_lon']:.2f}°E\n(n = {c['n']})"
                        for c in cells])
    a1.set_xlim(-0.55, 5.55); a1.set_ylim(0, 190)
    a1.set_ylabel("Mean PM$_{10}$ (µg m$^{-3}$)")
    a1.set_xlabel("CAMS cell centre (monitors in cell)")
    a1.yaxis.grid(True, color=GRIDC, lw=0.6); a1.set_axisbelow(True)
    for sp in ("top", "right"):
        a1.spines[sp].set_visible(False)
    for nm in ("Al Mafraq", "E11 Road"):
        i = [k for k, c in enumerate(cells) if nm in c["stations"]][0]
        a1.annotate(nm, (i + 0.10, by[nm]["obs_mean"]), xytext=(8, -2.5),
                    textcoords="offset points", fontsize=7, color=INK2)
    a1.legend(handles=[
        Line2D([], [], marker="o", ls="", mfc=INK2, mec=SURFACE, ms=5.5,
               label="monitor, observed mean"),
        Line2D([], [], color=BLUE, lw=2.0, label="CAMS (one value per cell)")],
        loc="upper left", bbox_to_anchor=(0.0, 1.02), ncol=1, fontsize=7,
        handletextpad=0.5, labelspacing=0.6, labelcolor=INK2)
    panel_letter(a1, "a", dx=-0.08, dy=1.0)

    P = R["pairs"]
    dc = [p for p in P if not p["same_cell"]]
    sc = [p for p in P if p["same_cell"]]
    a2.scatter([p["dist_km"] for p in dc], [p["spearman_obs"] for p in dc], s=13,
               facecolors="none", edgecolors=MUTED, linewidths=0.7, zorder=2,
               label=f"different cells ({len(dc)} pairs)")
    a2.scatter([p["dist_km"] for p in sc], [p["spearman_obs"] for p in sc], s=30, c=BLUE,
               edgecolors=SURFACE, linewidths=0.8, zorder=4,
               label=f"same cell ({len(sc)} pairs)")
    ref = R["coherence"]["cams_obs_spearman_median"]
    a2.axhline(ref, color=INK2, lw=0.9, ls=(0, (4, 3)), zorder=3)
    a2.text(1.25, ref - 0.02, f"CAMS vs. its own monitor (median {ref:.2f})",
            fontsize=7, color=INK2, va="top",
            bbox=dict(fc=SURFACE, ec="none", pad=0.8, alpha=0.85))
    a2.axvline(80, color=AXISC, lw=0.9, zorder=1)
    a2.text(80, 1.012, "grid spacing ≈ 80 km", fontsize=7, color=MUTED, va="bottom",
            ha="center", transform=a2.get_xaxis_transform())
    a2.set_xscale("log"); a2.set_xlim(1.1, 420); a2.set_ylim(-0.1, 1.0)
    a2.set_xticks([2, 5, 10, 20, 50, 100, 200, 400])
    a2.set_xticklabels(["2", "5", "10", "20", "50", "100", "200", "400"])
    a2.minorticks_off()
    a2.set_xlabel("Separation between monitors (km)")
    a2.set_ylabel("Spearman correlation of 24-h mean PM$_{10}$")
    a2.yaxis.grid(True, color=GRIDC, lw=0.6); a2.set_axisbelow(True)
    for sp in ("top", "right"):
        a2.spines[sp].set_visible(False)
    a2.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=1, fontsize=7,
              handletextpad=0.2, labelcolor=INK2)
    panel_letter(a2, "b", dx=-0.09, dy=1.0)
    save(fig, "Fig2_subgrid")
    plt.close(fig)


# ------------------------------------------------------------------ Fig. 3
def bh(p):
    p = np.asarray(p, float)
    n = len(p); o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n); q[o] = np.minimum(r, 1.0)
    return q


def fig3():
    keys = [k for k, v in A.items() if isinstance(v, dict) and "binom_p" in v]
    q = dict(zip(keys, bh([A[k]["binom_p"] for k in keys])))
    LOW = -0.60
    stages = [("buffer", "Distance buffer only"),
              ("buffer_time", "+ time blocks, 10-day guard"),
              ("buffer_time_cell", "+ same-cell monitors withheld")]
    specs = [("full", BLUE, "o", -0.11, "with coordinates"),
             ("coordfree", ORANGE, "s", +0.11, "coordinate-free")]
    radii = [0, 25, 50, 100]
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.95), sharey=True,
                             gridspec_kw=dict(left=0.10, right=0.985, top=0.80,
                                              bottom=0.165, wspace=0.07))
    for ax, (stg, title) in zip(axes, stages):
        ax.axhline(0, color=INK2, lw=0.8, zorder=2)
        for spec, col, mk, off, _ in specs:
            for i, Rk in enumerate(radii):
                k = f"{stg}_{spec}_R{Rk}"
                v = A[k]
                lo, hi = v["delta_median_ci_cell"]
                ax.plot([i + off] * 2, [max(lo, LOW), hi], color=col, lw=1.3, zorder=3,
                        solid_capstyle="round")
                if lo < LOW:
                    ax.annotate("", (i + off, LOW - 0.045), (i + off, LOW + 0.04),
                                arrowprops=dict(arrowstyle="-|>", color=col, lw=1.1,
                                                mutation_scale=7))
                    ax.text(i + off + (0.07 if off > 0 else -0.07), LOW,
                            f"{lo:+.2f}".replace("-", "−"),
                            fontsize=7, color=INK2, va="center",
                            ha="left" if off > 0 else "right")
                est = q[k] < 0.05
                ax.scatter(i + off, v["delta_median"], marker=mk, s=30, zorder=4,
                           facecolors=col if est else SURFACE, edgecolors=col,
                           linewidths=1.2)
        ax.set_title(title, fontsize=8, color=INK, pad=5, loc="center")
        ax.set_xticks(range(4)); ax.set_xticklabels([str(r) for r in radii])
        ax.set_xlim(-0.5, 3.5); ax.set_ylim(LOW - 0.06, 0.42)
        ax.yaxis.grid(True, color=GRIDC, lw=0.6); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_xlabel("Distance buffer $R$ (km)")
    axes[0].set_ylabel("Median $\\Delta$ = log(MAE$_{\\mathrm{raw}}$ / "
                       "MAE$_{\\mathrm{corrected}}$)\n> 0: correction reduces error")
    h = [Line2D([], [], marker=mk, ls="", mfc=col, mec=col, ms=5.5, label=lab)
         for _, col, mk, _, lab in specs]
    h += [Line2D([], [], marker="o", ls="", mfc=INK2, mec=INK2, ms=5.5,
                 label="differs from zero (sign test, BH $q$ < 0.05)"),
          Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=INK2, mew=1.2, ms=5.5,
                 label="does not")]
    fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.08, 1.0), ncol=4, fontsize=7,
               handletextpad=0.2, columnspacing=1.4, labelcolor=INK2)
    for ax, s in zip(axes, "abc"):
        ax.text(0.0, 1.03, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
                color=INK, va="bottom", ha="left")
    save(fig, "Fig4_transfer")
    plt.close(fig)


# ------------------------------------------------------------------ Fig. 4
def fig4():
    reg = R["regional"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W, 3.3),
                                 gridspec_kw=dict(left=0.075, right=0.985, top=0.93,
                                                  bottom=0.235, wspace=0.26,
                                                  width_ratios=[1.0, 1.15]))
    o = np.array([d["obs"] for d in reg["series"]])
    c = np.array([d["cams"] for d in reg["series"]])
    hi = 1.05 * max(o.max(), c.max())
    a1.plot([0, hi], [0, hi], color=AXISC, lw=0.9, zorder=1)
    a1.scatter(o, c, s=14, c=INK2, edgecolors=SURFACE, linewidths=0.5, zorder=3)
    a1.set_xlim(0, hi); a1.set_ylim(0, hi); a1.set_aspect("equal")
    a1.set_xlabel("Observed network median, 24-h mean PM$_{10}$ (µg m$^{-3}$)")
    a1.set_ylabel("CAMS network median (µg m$^{-3}$)")
    a1.text(0.04, 0.95, f"Spearman {reg['spearman_regional_obs_vs_cams']:.2f}\n"
            f"{reg['n_days_with_model']} days", transform=a1.transAxes, fontsize=7.5,
            color=INK, va="top")
    a1.text(hi * 0.98, hi * 0.90, "1:1", fontsize=7, color=MUTED, ha="right")
    a1.grid(True, color=GRIDC, lw=0.6); a1.set_axisbelow(True)
    for sp in ("top", "right"):
        a1.spines[sp].set_visible(False)
    panel_letter(a1, "a", dx=-0.13, dy=1.0)

    ks = ("low", "mid", "high")
    x = np.arange(3)
    series = [("bias_metro", "Abu Dhabi metropolitan area (7 monitors)", BLUE, "o", "-"),
              ("bias_al_ain", "Al Ain (4 monitors)", ORANGE, "s", "-"),
              ("metro_minus_al_ain", "difference, metropolitan minus Al Ain", INK2, "D",
               (0, (3, 2)))]
    a2.axhline(0, color=INK2, lw=0.8, zorder=2)
    for key, lab, col, mk, ls in series:
        y = [reg["classes"][k][key] for k in ks]
        a2.plot(x, y, color=col, lw=2.0 if ls == "-" else 1.3, ls=ls, marker=mk,
                ms=6 if ls == "-" else 4.5, mec=SURFACE, mew=0.9, zorder=4, label=lab)
        for xi, yi in zip(x, y):
            off = {"bias_al_ain": 8, "bias_metro": 8, "metro_minus_al_ain": -12}[key]
            lab_ = "0" if round(yi) == 0 else f"{yi:+.0f}".replace("-", "−")
            a2.annotate(lab_, (xi, yi), xytext=(0, off),
                        textcoords="offset points", ha="center", fontsize=7,
                        color=col if key != "metro_minus_al_ain" else INK2)
    a2.set_xticks(x)
    a2.set_xticklabels([f"{k}\n{reg['classes'][k]['regional_min']:.0f}–"
                        f"{reg['classes'][k]['regional_max']:.0f} µg m$^{{-3}}$\n"
                        f"({reg['classes'][k]['n_days']} days)" for k in ks])
    a2.set_xlim(-0.45, 2.45); a2.set_ylim(-110, 30)
    a2.set_xlabel("Regional level of the day (network median, terciles)")
    a2.set_ylabel("Mean bias, CAMS − observed (µg m$^{-3}$)")
    a2.yaxis.grid(True, color=GRIDC, lw=0.6); a2.set_axisbelow(True)
    for sp in ("top", "right"):
        a2.spines[sp].set_visible(False)
    a2.legend(loc="lower left", fontsize=7, labelcolor=INK2, handlelength=2.2)
    panel_letter(a2, "b", dx=-0.10, dy=1.0)
    save(fig, "Fig3_regional")
    plt.close(fig)


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
    print("wrote Fig1-Fig4 to", FIG)
