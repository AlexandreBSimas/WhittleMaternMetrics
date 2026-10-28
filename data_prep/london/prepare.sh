#!/bin/sh
# Rebuilds data/london/ from National Highways WebTRIS (loop detector speeds)
# and OpenStreetMap (road network). Needs curl, python3, R (sf, jsonlite)
# and internet access. The WebTRIS download takes about an hour.
set -e
cd "$(dirname "$0")"
mkdir -p raw
curl -s -o raw/sites.json https://webtris.nationalhighways.co.uk/api/v1.0/sites
python3 select_sites.py           # sites_selected.csv (514 sites)
python3 download.py               # raw/site_<id>.csv, 2025-01-01 .. 2025-06-30
Rscript make_panel.R              # panel_1700.csv (88 days x 256 sites)
Rscript fetch_osm.R               # raw/osm_london.rds (Overpass API)
cp panel_1700.csv sites_selected.csv raw/osm_london.rds ../../data/london/
echo "data/london updated"
