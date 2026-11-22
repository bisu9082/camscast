"""Build the records used for the correction analysis.

From the flat-line-screened three-hourly records (data_processed/full/aligned_qc6), two sets
with identical rows:
  matched        model value = mean of the eight CAMS values in the 24 h ending at the time stamp
  instantaneous  model value = CAMS at the time stamp (how an hourly record would be treated)
The residual obs - model is recomputed for each.

Then, from code/audit/:
  CAMSCAST_DATA=<out>/matched        GUARD_DAYS=10 BLOCKS=quantile python audit_v7.py
  CAMSCAST_DATA=<out>/instantaneous  GUARD_DAYS=10 BLOCKS=quantile python audit_v7.py

Usage: python make_matched_records.py [out_dir]   (default data_processed/full/derived)
"""
import os
import sys

import pandas as pd

from common_ae import DATA, STATIONS, load

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA, "derived")
for k in ("matched", "instantaneous"):
    os.makedirs(os.path.join(out, k), exist_ok=True)
n = 0
for sid in STATIONS:
    q = pd.read_csv(os.path.join(DATA, "aligned_qc6", f"aligned_pm10_{sid}.csv"),
                    index_col=0, parse_dates=True)
    g = load(sid)
    x = q[["obs"]].join(g[["cams", "cams24"]], how="left").dropna()
    m = pd.DataFrame({"cams": x.cams24, "obs": x.obs, "resid": x.obs - x.cams24})
    i = pd.DataFrame({"cams": x.cams, "obs": x.obs, "resid": x.obs - x.cams})
    m.to_csv(os.path.join(out, "matched", f"aligned_pm10_{sid}.csv"))
    i.to_csv(os.path.join(out, "instantaneous", f"aligned_pm10_{sid}.csv"))
    n += len(x)
print("rows per set:", n)
