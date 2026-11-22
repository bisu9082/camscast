"""Download CAMS global reanalysis (EAC4) aerosol optical depths and 10-m wind over Abu Dhabi,
Oct 2022 - May 2024 (companion to fetch_cams_eac4_pm10.py; same area, period and times).

Run on a machine with an ADS account:
    pip install "cdsapi>=0.7"
    ~/.cdsapirc (or ~/.adsapirc) must contain
        url: https://ads.atmosphere.copernicus.eu/api
        key: <your ADS personal access token>
    and the EAC4 licence must have been accepted once on the dataset page.
    python fetch_cams_eac4_aod_wind.py [output_dir]

Writes one NetCDF file per calendar year into output_dir (default: this script's folder):
    cams_eac4_aodwind_abudhabi_2022.nc, ..._2023.nc, ..._2024.nc
Area 25.5N 51.75E 22.5S 56.25E; eight analysis times per day; variables: total, dust and
sea-salt aerosol optical depth at 550 nm (dimensionless) and the 10-m wind components (m s-1).
Variable names were checked against the ADS process description of cams-global-reanalysis-eac4
on 14 November 2026.
The script never prints or stores the token.
"""
import os
import sys

import cdsapi

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)
VARIABLES = ["total_aerosol_optical_depth_550nm", "dust_aerosol_optical_depth_550nm",
             "sea_salt_aerosol_optical_depth_550nm",
             "10m_u_component_of_wind", "10m_v_component_of_wind"]
PERIODS = {"2022": "2022-10-01/2022-12-31",
           "2023": "2023-01-01/2023-12-31",
           "2024": "2024-01-01/2024-05-31"}


def client():
    for rc in ("~/.adsapirc", "~/.cdsapirc"):
        p = os.path.expanduser(rc)
        if os.path.exists(p) and "ads.atmosphere" in open(p).read():
            os.environ["CDSAPI_RC"] = p
            return cdsapi.Client()
    if os.environ.get("ADS_KEY"):
        return cdsapi.Client(url="https://ads.atmosphere.copernicus.eu/api",
                             key=os.environ["ADS_KEY"])
    sys.exit("No ADS credentials found: create ~/.adsapirc with the ADS url and key "
             "(see the docstring), or set ADS_KEY in the environment.")


def main():
    c = client()
    for year, dates in PERIODS.items():
        f = os.path.join(OUT, f"cams_eac4_aodwind_abudhabi_{year}.nc")
        if os.path.exists(f) and os.path.getsize(f) > 0:
            print("exists, skipped:", f)
            continue
        print("requesting", year, dates, flush=True)
        c.retrieve("cams-global-reanalysis-eac4",
                   {"variable": VARIABLES,
                    "date": dates,
                    "time": ["00:00", "03:00", "06:00", "09:00",
                             "12:00", "15:00", "18:00", "21:00"],
                    "area": [25.5, 51.75, 22.5, 56.25],
                    "data_format": "netcdf"}, f)
        print("written", f, os.path.getsize(f), "bytes", flush=True)


if __name__ == "__main__":
    main()
