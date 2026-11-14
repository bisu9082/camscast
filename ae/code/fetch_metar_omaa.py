"""Download the METAR reports of Abu Dhabi International Airport (ICAO OMAA) from the Iowa
Environmental Mesonet ASOS/METAR archive (no account needed), 30 September 2022 to 31 May 2024.

    python fetch_metar_omaa.py [output_dir]      (default: data_raw/metar)

Writes OMAA_metar_20220930_20240531.csv with the station, time (UTC), visibility (statute miles),
present-weather codes and the raw METAR text. Retrieved for the manuscript on 14 November 2026.
"""
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "data_raw", "metar")
URL = ("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=OMAA"
       "&data=vsby&data=wxcodes&data=metar&year1=2022&month1=9&day1=30&year2=2024&month2=6&day2=1"
       "&tz=Etc/UTC&format=onlycomma&latlon=yes&elev=no&missing=M&trace=T&direct=no&report_type=3")
os.makedirs(OUT, exist_ok=True)
with urllib.request.urlopen(URL, timeout=300) as r:
    data = r.read()
f = os.path.join(OUT, "OMAA_metar_20220930_20240531.csv")
open(f, "wb").write(data)
print(f, len(data), "bytes")
