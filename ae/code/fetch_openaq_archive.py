"""Rebuild data_raw/abudhabi_openaq_hourly.csv.gz from the OpenAQ public archive (no API key).

Lists every object under records/csv.gz/locationid=<id>/ for the 19 Abu Dhabi locations in the
bucket openaq-data-archive, downloads the daily files and concatenates them. Executed on
13 November 2026: 6,803 daily files, 401,693 rows (pm10 98,459; no2 99,270; so2 99,218;
o3 68,171; co 36,575), 29 October 2022 to 20 May 2024.

Usage: python fetch_openaq_archive.py [workdir]
"""
import glob
import gzip
import io
import os
import re
import sys
import urllib.parse
import urllib.request

import pandas as pd

from common_ae import STATIONS

BASE = "https://openaq-data-archive.s3.amazonaws.com/"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data_raw", "abudhabi_openaq_hourly.csv.gz")


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except OSError:
            if k == tries - 1:
                raise


def keys(loc):
    tok = None
    while True:
        q = {"list-type": "2", "prefix": f"records/csv.gz/locationid={loc}/"}
        if tok:
            q["continuation-token"] = tok
        x = get(BASE + "?" + urllib.parse.urlencode(q)).decode()
        yield from re.findall(r"<Key>([^<]+)</Key>", x)
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", x)
        if not m:
            return
        tok = m.group(1)


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "data_raw", "files")
    os.makedirs(work, exist_ok=True)
    for loc in STATIONS:
        for k in keys(loc):
            f = os.path.join(work, os.path.basename(k))
            if not (os.path.exists(f) and os.path.getsize(f) > 0):
                open(f, "wb").write(get(BASE + k))
    parts = [pd.read_csv(io.BytesIO(gzip.open(f).read()))
             for f in sorted(glob.glob(os.path.join(work, "*.csv.gz")))]
    d = pd.concat(parts, ignore_index=True)
    d["t"] = pd.to_datetime(d.datetime, utc=True)
    d.drop(columns=["datetime"]).to_csv(OUT, index=False)
    print(len(parts), "files,", len(d), "rows ->", OUT)


if __name__ == "__main__":
    main()
