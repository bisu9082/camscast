"""Shared definitions for the Atmospheric Environment analysis.

The archived PM10 values of the Abu Dhabi network behave as trailing 24-h running means, not as
hourly concentrations (averaging_fingerprint.py gives the evidence; the national reporting
convention for particulate matter is a 24-h rolling average). Every reported value is therefore
treated as the mean over the preceding 24 h, and the model counterpart of a value reported at
time t is the mean of the eight three-hourly CAMS values in the 24 h ending at t, taken from the
continuous EAC4 series of the monitor's cell (data_processed/cams_eac4_pm10_cells.csv).
"""
import math
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("CAMSCAST_DATA", os.path.join(HERE, "..", "data_processed", "full"))
CAMS_CELLS = os.path.join(HERE, "..", "data_processed", "cams_eac4_pm10_cells.csv")
DAILY_HOUR = 8        # 08 UTC: the only hour the archive holds in every month
RES = os.environ.get("CAMSCAST_RESULTS", os.path.join(HERE, "..", "results"))
WHO24 = 45.0          # WHO 2021 24-h PM10 guideline value, ug/m3
GRID = 0.75           # CAMS EAC4 grid spacing, degrees
WINDOW = 8            # three-hourly steps in 24 h
MIN_STEPS = 7         # CAMS steps required in the trailing window
SEED = 20261113

# OpenAQ sensor id -> (name, lat, lon, operator station type)
# Coordinates copied programmatically from loso_v2.PM10C (repository release v1.3).
# Station type: Statistics Centre - Abu Dhabi, "Air Quality Statistics Methodology".
STATIONS = {
    "369953": ("Bida Zayed", 23.6522, 53.7038, "Suburban"),
    "369954": ("Khadeeja School", 24.4819, 54.3691, "Urban"),
    "369955": ("Liwa", 23.0957, 53.6064, "Rural background"),
    "369957": ("Hamdan St.", 24.4889, 54.3637, "Traffic"),
    "369958": ("Mussafah", 24.3471, 54.5028, "Industrial"),
    "369959": ("Al Ain St.", 24.2258, 55.7658, "Traffic"),
    "369960": ("Sweihan", 24.4665, 55.3428, "Suburban"),
    "369961": ("Khalifa City", 24.4199, 54.5781, "Suburban"),
    "369962": ("Gayathi", 23.8355, 52.8103, "Suburban"),
    "369963": ("Al Maqtaa", 24.4035, 54.516, "Urban"),
    "369965": ("Ruwais", 24.0908, 52.7548, "Industrial"),
    "369966": ("Baniyas School", 24.3213, 54.6359, "Suburban"),
    "369967": ("Habshan S.", 23.7504, 53.7452, "Industrial"),
    "369968": ("Al Ain Islamic", 24.219, 55.7348, "Suburban"),
    "369969": ("Zakher", 24.1634, 55.7021, "Urban"),
    "369971": ("Al Qua'a", 23.5311, 55.4859, "Rural background"),
    "369972": ("Al Tawia", 24.2591, 55.7048, "Suburban"),
    "369973": ("E11 Road", 24.0351, 53.8853, "Traffic"),
    "369974": ("Al Mafraq", 24.287, 54.5861, "Industrial"),
}
SOURCE_TYPES = ("Industrial", "Traffic", "Urban")
METRO = ("Khadeeja School", "Hamdan St.", "Mussafah", "Khalifa City", "Al Maqtaa",
         "Baniyas School", "Al Mafraq")
ALAIN = ("Al Ain St.", "Al Ain Islamic", "Zakher", "Al Tawia")


def haversine(a, b, c, d):
    p1, p2 = math.radians(a), math.radians(c)
    h = (math.sin(math.radians(c - a) / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(math.radians(d - b) / 2) ** 2)
    return 2 * 6371.0 * math.asin(math.sqrt(h))


def cell_of(lat, lon):
    return (round(lat / GRID) * GRID, round(lon / GRID) * GRID)


_CAMS = {}
_DAILY = {}


def cams_cells():
    if "c" not in _CAMS:
        c = pd.read_csv(CAMS_CELLS, index_col=0, parse_dates=True)
        _CAMS["c"] = c.reindex(pd.date_range(c.index.min(), c.index.max(), freq="3h"))
    return _CAMS["c"]


def load(sid, data=None):
    """One monitor on the regular three-hourly grid of the continuous CAMS series.

    obs     reported value (trailing 24-h running mean ending at the time stamp); NaN where the
            archive holds no value at that stamp
    cams    instantaneous three-hourly CAMS value of the monitor's cell
    cams24  mean of the eight CAMS values in the 24 h ending at the time stamp
    """
    d = pd.read_csv(os.path.join(data or DATA, f"aligned_pm10_{sid}.csv"),
                    index_col=0, parse_dates=True)[["obs"]].astype(float)
    c = cams_cells()
    la, lo = STATIONS[sid][1], STATIONS[sid][2]
    k = f"{cell_of(la, lo)[0]:.2f}_{cell_of(la, lo)[1]:.2f}"
    g = pd.DataFrame({"cams": c[k]})
    g["cams24"] = g.cams.rolling(WINDOW, min_periods=WINDOW).mean()
    g = g.join(d, how="left")
    g = g[(g.index >= d.index.min()) & (g.index <= d.index.max())]
    g.attrs["sid"] = sid
    return g


def load_all(data=None):
    return {sid: load(sid, data) for sid in STATIONS}


def daily(g):
    """One value per monitor and day: the reported 24-h mean ending at 08 UTC and the mean of the
    eight CAMS values in the same 24 h (data_processed/full/daily_08utc.csv)."""
    if "d" not in _DAILY:
        D = pd.read_csv(os.path.join(DATA, "daily_08utc.csv"), dtype={"sid": str},
                        parse_dates=["day"])
        _DAILY["d"] = {s: x.set_index("day")[["obs", "cams24"]] for s, x in D.groupby("sid")}
    return _DAILY["d"][g.attrs["sid"]]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n)
    q[o] = np.minimum(r, 1.0)
    return q
