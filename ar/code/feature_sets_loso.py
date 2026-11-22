"""Which predictors carry the leave-one-station-out gain of the correction?

The coordinate-free correction is refitted with four sets of predictors, leaving one monitor out
at a time and withholding nothing else, on two record sets: the screened record with the
instantaneous model value, and the matched record (model averaged over the 24 h ending at the
observation time). delta = log(MAE_raw / MAE_corrected) per monitor. Features are built
as in audit/loso_v2.py, so the last set reproduces the coordinate-free leave-one-station-out
result of the correction analysis.

Usage: python feature_sets_loso.py <dir instantaneous records> <dir matched records>
Writes results/feature_sets_loso.json
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.ensemble import HistGradientBoostingRegressor

from common_ae import RES, STATIONS

SETS = {"cams": ["cams"], "cams+hour": ["cams", "hsin", "hcos"],
        "cams+doy": ["cams", "dsin", "dcos"],
        "cams+hour+doy": ["cams", "hsin", "hcos", "dsin", "dcos"]}


def features(d):
    h, doy = d.index.hour + d.index.minute / 60, d.index.dayofyear
    return pd.DataFrame(dict(cams=d.cams.values,
                             hsin=np.sin(2 * np.pi * h / 24), hcos=np.cos(2 * np.pi * h / 24),
                             dsin=np.sin(2 * np.pi * doy / 365),
                             dcos=np.cos(2 * np.pi * doy / 365)), index=d.index)


def run(folder):
    data = {}
    for sid in STATIONS:
        d = pd.read_csv(os.path.join(folder, f"aligned_pm10_{sid}.csv"), index_col=0,
                        parse_dates=True).dropna()
        data[sid] = (d, features(d))
    out = {}
    for name, cols in SETS.items():
        dl = []
        for t in data:
            src = [j for j in data if j != t]
            X = pd.concat([data[j][1][cols] for j in src])
            y = np.concatenate([data[j][0].resid.values for j in src])
            m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_depth=6,
                                              random_state=0, early_stopping=False)
            p = m.fit(X, y).predict(data[t][1][cols])
            d = data[t][0]
            dl.append(float(np.log(np.mean(np.abs(d.cams - d.obs))
                                   / np.mean(np.abs(d.cams + p - d.obs)))))
        out[name] = dict(median=float(np.median(dl)), wins=int(sum(x > 0 for x in dl)), deltas=dl)
        print(f"{os.path.basename(folder.rstrip('/')):14s} {name:14s} median delta "
              f"{np.median(dl):+.3f} improved {sum(x > 0 for x in dl)}/19", flush=True)
    return out, int(sum(len(v[0]) for v in data.values()))


def main():
    inst, n_i = run(sys.argv[1])
    mat, n_m = run(sys.argv[2])
    res = {"_meta": dict(
        design="leave-one-station-out, nothing else withheld; "
               "delta = log(MAE_raw / MAE_corrected)",
        n_instantaneous=n_i, n_matched=n_m,
        software=f"scikit-learn {sklearn.__version__}, numpy {np.__version__}, "
                 f"pandas {pd.__version__}, scipy {scipy.__version__}"),
        "instantaneous": inst, "matched": mat}
    json.dump(res, open(os.path.join(RES, "feature_sets_loso.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
