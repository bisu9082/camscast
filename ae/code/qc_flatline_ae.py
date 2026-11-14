"""Flat-line screen of the three-hourly records (rule of release v1.3, primary threshold).

Removes every value that belongs to a run of six or more identical consecutive values of a
monitor's record. Reads data_processed/full/aligned_pm10_<id>.csv, writes
data_processed/full/aligned_qc6/ and results/qc_flatline_full.json.
"""
import json
import os

import pandas as pd

from common_ae import DATA, RES, STATIONS

K = 6


def main():
    out = os.path.join(DATA, "aligned_qc6")
    os.makedirs(out, exist_ok=True)
    per, rem, kept = {}, 0, 0
    for sid in STATIONS:
        d = pd.read_csv(os.path.join(DATA, f"aligned_pm10_{sid}.csv"), index_col=0,
                        parse_dates=True)
        run = d.obs.groupby((d.obs != d.obs.shift()).cumsum()).transform("size")
        keep = run < K
        d[keep].to_csv(os.path.join(out, f"aligned_pm10_{sid}.csv"))
        per[sid] = dict(before=int(len(d)), removed=int((~keep).sum()))
        rem += int((~keep).sum())
        kept += int(keep.sum())
    res = dict(threshold=K, removed=rem, kept=kept, removed_fraction=rem / (rem + kept),
               stations_affected=int(sum(v["removed"] > 0 for v in per.values())), per_station=per)
    json.dump(res, open(os.path.join(RES, "qc_flatline_full.json"), "w"), indent=1)
    print(f"removed {rem} of {rem + kept} ({100 * rem / (rem + kept):.1f}%), "
          f"{res['stations_affected']} monitors affected")


if __name__ == "__main__":
    main()
