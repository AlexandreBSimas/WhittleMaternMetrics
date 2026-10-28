#!/bin/sh
# Rebuilds data/paris/ from the City of Paris open data (permanent traffic
# counters). Needs curl, python3 with pyproj, and internet access.
# NOTE: the source keeps only about the last 13 months of hourly data, so a
# rerun gives a different period than the files in data/paris/ (downloaded
# on 23 October 2026: 1 September 2025 - 18 October 2026).
set -e
cd "$(dirname "$0")"
mkdir -p raw
B=https://opendata.paris.fr/api/explore/v2.1/catalog/datasets
# Hour 08:00-09:00 local time (time stamp 09:00) on all days: ~50 MB.
curl -s -G "$B/comptages-routiers-permanents/exports/csv" \
  --data-urlencode "select=iu_ac,t_1h,q,k" --data-urlencode "where=hour(t_1h)=9" \
  --data-urlencode "timezone=Europe/Paris" --data-urlencode "delimiter=," \
  -o raw/h09.csv --max-time 1800
# Geometry of the counted arcs, with upstream and downstream node ids.
curl -s -G "$B/referentiel-comptages-routiers/exports/geojson" -o raw/arcs.geojson --max-time 600
python3 make_panel.py 0.9
cp network_edges.csv sites.csv panel_k.csv panel_q.csv ../../data/paris/
echo "data/paris updated"
