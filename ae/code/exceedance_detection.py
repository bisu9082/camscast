"""Does the reanalysis separate days that exceed the WHO 24-h PM10 guideline value from days that
do not?

One value per monitor and day is used: the reported 24-h mean ending at 08 UTC, against the CAMS
mean over the same 24 h (common_ae.daily). Three classifiers of "observed 24-h mean > 45":

  raw            the matched CAMS 24-h mean; needs no monitor at the site
  level_shifted  CAMS plus the monitor's own mean observed-minus-CAMS difference over the
                 preceding eight daily values; needs a monitor at the site
  persistence    the monitor's own value one day earlier

Intervals: 2,000 bootstrap resamples of monitors, the contingency table rebuilt from the drawn
monitors. Writes results/exceedance_detection.json
"""
import json
import os

import numpy as np
import pandas as pd

from common_ae import RES, SEED, WHO24, daily, load_all

NBOOT, LOOKBACK = 2000, 8


def tables(x):
    full = pd.date_range(x.index.min(), x.index.max(), freq="D")
    x = x.reindex(full)
    diff = (x.obs - x.cams24)
    shift = diff.shift(1).rolling(LOOKBACK, min_periods=LOOKBACK // 2).mean()
    cls = {"raw": x.cams24, "level_shifted": x.cams24 + shift, "persistence": x.obs.shift(1)}
    out = {}
    for k, f in cls.items():
        ok = f.notna() & x.obs.notna() & x.cams24.notna()
        fo, ob = f[ok] > WHO24, x.obs[ok] > WHO24
        out[k] = np.array([(fo & ob).sum(), (~fo & ob).sum(), (fo & ~ob).sum(),
                           (~fo & ~ob).sum()], float)
    # all three on the days where every classifier is defined
    ok = x.obs.notna() & x.cams24.notna() & cls["level_shifted"].notna() \
        & cls["persistence"].notna()
    for k, f in cls.items():
        fo, ob = f[ok] > WHO24, x.obs[ok] > WHO24
        out[k + "_common"] = np.array([(fo & ob).sum(), (~fo & ob).sum(), (fo & ~ob).sum(),
                                       (~fo & ~ob).sum()], float)
    return out


def scores(hit, miss, fa, cn):
    n = hit + miss + fa + cn
    ar = (hit + fa) * (hit + miss) / n
    ex = ((hit + miss) * (hit + fa) + (cn + miss) * (cn + fa)) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        return dict(POD=hit / (hit + miss), FAR=fa / (hit + fa),
                    POFD=fa / (fa + cn), CSI=hit / (hit + miss + fa),
                    ETS=(hit - ar) / (hit + miss + fa - ar), HSS=(hit + cn - ex) / (n - ex))


def main():
    rng = np.random.default_rng(SEED)
    rec = load_all()
    T = [tables(daily(g)) for g in rec.values()]
    res = {"_meta": dict(guideline=WHO24, rule="strictly above", n_boot=NBOOT, seed=SEED,
                         lookback_days=LOOKBACK, n_monitors=len(T))}
    for key in T[0]:
        M = np.array([t[key] for t in T])
        tot = M.sum(axis=0)
        pt = scores(*tot)
        draws = {k: [] for k in pt}
        for _ in range(NBOOT):
            s = scores(*M[rng.integers(0, len(M), len(M))].sum(axis=0))
            for k, v in s.items():
                draws[k].append(v)
        r = dict(hit=int(tot[0]), miss=int(tot[1]), fa=int(tot[2]), cn=int(tot[3]),
                 n=int(tot.sum()), events=int(tot[0] + tot[1]),
                 base_rate=float((tot[0] + tot[1]) / tot.sum()))
        for k, v in pt.items():
            a = np.array(draws[k], float)
            a = a[np.isfinite(a)]
            r[k] = float(v)
            r[k + "_lo"], r[k + "_hi"] = (float(np.percentile(a, 2.5)),
                                          float(np.percentile(a, 97.5)))
        res[key] = r
        print(f"{key:22s} N={r['n']:5d} events={r['events']:5d} ({100*r['base_rate']:.1f}%) "
              f"h/m/f/c={r['hit']}/{r['miss']}/{r['fa']}/{r['cn']} POD {r['POD']:.3f} "
              f"[{r['POD_lo']:.3f},{r['POD_hi']:.3f}] POFD {r['POFD']:.3f} FAR {r['FAR']:.3f} "
              f"ETS {r['ETS']:.3f} [{r['ETS_lo']:.3f},{r['ETS_hi']:.3f}] "
              f"HSS {r['HSS']:.3f} [{r['HSS_lo']:.3f},{r['HSS_hi']:.3f}]")
    json.dump(res, open(os.path.join(RES, "exceedance_detection.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
