# Select WebTRIS (National Highways, England) traffic sites for the London
# application: active MIDAS loop sites in a box around the M25, thinned
# greedily so that no two selected sites are closer than MIN_KM (straight-line
# distance). This keeps at most one carriageway per location.
# Input raw/sites.json (https://webtris.nationalhighways.co.uk/api/v1.0/sites).
import json, math, csv
BOX = (-0.65, 51.22, 0.40, 51.78)
MIN_KM = 0.3
s = json.load(open("raw/sites.json"))["sites"]
a = [x for x in s if x["Status"] == "Active" and x["Description"].startswith("MIDAS")
     and BOX[0] <= x["Longitude"] <= BOX[2] and BOX[1] <= x["Latitude"] <= BOX[3]]
a.sort(key=lambda x: int(x["Id"]))
def km(p, q):
    dx = (p["Longitude"] - q["Longitude"]) * 111.32 * math.cos(math.radians(51.5))
    dy = (p["Latitude"] - q["Latitude"]) * 110.57
    return math.hypot(dx, dy)
keep = []
for x in a:
    if all(km(x, y) >= MIN_KM for y in keep):
        keep.append(x)
with open("sites_selected.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["site_id", "name", "lon", "lat", "description"])
    for x in keep:
        w.writerow([x["Id"], x["Name"], x["Longitude"], x["Latitude"], x["Description"]])
print(len(a), "active MIDAS sites in box;", len(keep), "selected")
