"""Download the AERONET Version 3 daily averages used in aeronet_dust_check.py (no account needed).

    python fetch_aeronet.py [output_dir]      (default: data_raw/aeronet)

Writes DEWA_SDA15.txt (Spectral Deconvolution Algorithm, Level 1.5) and
DEWA_ResearchCentre_AOD15.txt (direct-sun AOD, Level 1.5) for DEWA_ResearchCentre,
1 October 2022 to 31 May 2024, from the AERONET web service print_web_data_v3.
Retrieved for the manuscript on 14 November 2026; Level 1.5 data may change.
"""
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "data_raw", "aeronet")
BASE = ("https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site=DEWA_ResearchCentre"
        "&year=2022&month=10&day=1&year2=2024&month2=5&day2=31&AVG=20&if_no_html=1")
os.makedirs(OUT, exist_ok=True)
for prod, name in (("SDA15", "DEWA_SDA15.txt"), ("AOD15", "DEWA_ResearchCentre_AOD15.txt")):
    with urllib.request.urlopen(f"{BASE}&{prod}=1", timeout=180) as r:
        data = r.read()
    open(os.path.join(OUT, name), "wb").write(data)
    print(name, len(data), "bytes")
