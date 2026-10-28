# Paris permanent traffic counters (opendata.paris.fr, "comptages-routiers-
# permanents"), hour 08:00-09:00 local time (time stamp 09:00), working days
# Monday-Friday, 2025-09-01 .. 2026-10-16.
#
# Inputs (raw/): h09.csv (iu_ac, t_1h, q, k) and arcs.geojson (reference
# geometry of the counted arcs, with upstream/downstream node ids).
# Outputs:
#   network_edges.csv  one row per pair of nodes joined by a counted arc:
#                      edge, node_a, node_b, then the polyline as
#                      "x1 y1;x2 y2;..." in Lambert-93 km (EPSG:2154); the
#                      end points are moved to the node positions (mean of the
#                      arc end points at that node), so arcs meet at vertices.
#   sites.csv          one site per edge (the direction with the best
#                      coverage): iu_ac, edge, libelle, x_km, y_km (arc
#                      midpoint), relative position along the edge.
#   panel_q.csv / panel_k.csv   rows = days, columns = sites (flow in veh/h,
#                      occupancy in %); complete (no missing values).
# Usage: python3 make_panel.py [min_coverage]
import csv, json, sys, math, datetime, collections
from pyproj import Transformer

MIN_COV = float(sys.argv[1]) if len(sys.argv) > 1 else 0.9
tr = Transformer.from_crs(4326, 2154, always_xy=True)
def km(lon, lat):
    x, y = tr.transform(lon, lat); return (x / 1000.0, y / 1000.0)

arcs = json.load(open("raw/arcs.geojson"))["features"]
ends = collections.defaultdict(list); arc = {}
for f in arcs:
    p = f["properties"]
    co = [km(*c) for c in f["geometry"]["coordinates"]]
    a, b = p["iu_nd_amont"], p["iu_nd_aval"]
    arc[p["iu_ac"]] = dict(a=a, b=b, co=co, lib=p["libelle"])
    ends[a].append(co[0]); ends[b].append(co[-1])
node = {v: (sum(x for x, _ in c) / len(c), sum(y for _, y in c) / len(c))
        for v, c in ends.items()}

q = collections.defaultdict(dict); k = collections.defaultdict(dict)
for r in csv.DictReader(open("raw/h09.csv", encoding="utf-8-sig")):
    d = r["t_1h"][:10]
    if datetime.date.fromisoformat(d).weekday() >= 5: continue
    if r["q"] and r["k"] and float(r["q"]) > 0 and float(r["k"]) > 0:
        q[r["iu_ac"]][d] = float(r["q"]); k[r["iu_ac"]][d] = float(r["k"])
days = sorted({d for s in q.values() for d in s})

# Edges: unordered node pairs; the polyline of the best-covered arc.
by_pair = collections.defaultdict(list)
for i, a in arc.items():
    if a["a"] != a["b"]: by_pair[frozenset((a["a"], a["b"]))].append(i)
edges, sites = [], []
for e, (pair, ids) in enumerate(sorted(by_pair.items(), key=lambda t: min(t[1]))):
    ids.sort(key=lambda i: -len(q.get(i, {})))
    best = arc[ids[0]]
    co = [node[best["a"]]] + best["co"][1:-1] + [node[best["b"]]]
    edges.append((e + 1, best["a"], best["b"], co))
    if len(q.get(ids[0], {})) >= MIN_COV * len(days):
        # site at the arc midpoint (half of its length along the polyline)
        L = [math.dist(u, v) for u, v in zip(co, co[1:])]
        sites.append(dict(iu_ac=ids[0], edge=e + 1, libelle=best["lib"],
                          t=0.5, length_km=sum(L)))
print(f"{len(edges)} edges, {len(node)} nodes, {len(days)} working days; "
      f"{len(sites)} sites with coverage >= {MIN_COV}")
# Make the panel complete: repeatedly remove the day or the site with the
# largest share of missing values (relative to the current panel).
keep_day = list(days)
while True:
    ids = [s["iu_ac"] for s in sites]
    md = {d: sum(d not in q[i] for i in ids) / len(ids) for d in keep_day}
    ms = {i: sum(d not in q[i] for d in keep_day) / len(keep_day) for i in ids}
    wd, ws = max(md, key=md.get), max(ms, key=ms.get)
    if md[wd] == 0: break
    if md[wd] >= ms[ws]: keep_day.remove(wd)
    else: sites = [s for s in sites if s["iu_ac"] != ws]
print(f"panel: {len(keep_day)} days x {len(sites)} sites")

with open("network_edges.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["edge", "node_a", "node_b", "coords"])
    for e, a, b, co in edges:
        w.writerow([e, a, b, ";".join(f"{x:.5f} {y:.5f}" for x, y in co)])
with open("sites.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["iu_ac", "edge", "libelle", "t", "length_km"])
    w.writeheader(); w.writerows(sites)
for name, src in (("panel_q.csv", q), ("panel_k.csv", k)):
    with open(name, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["date"] + [s["iu_ac"] for s in sites])
        for d in keep_day: w.writerow([d] + [src[s["iu_ac"]][d] for s in sites])
