# Download 15-minute WebTRIS daily reports (average speed and volume) for the
# selected sites, 2025-01-01 to 2025-06-30, into raw/site_<id>.csv.
# Usage: python3 download.py   (resumes; existing files are skipped)
import csv, json, os, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
URL = ("https://webtris.nationalhighways.co.uk/api/v1.0/reports/Daily?sites={}"
       "&start_date=01012025&end_date=30062025&page=1&page_size=20000")
def get(site):
    out = f"raw/site_{site}.csv"
    if os.path.exists(out): return site, "cached"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(URL.format(site), timeout=180) as r:
                body = r.read()
            rows = json.loads(body)["Rows"] if body.strip() else []
            break
        except Exception as e:
            err = e; time.sleep(5 * (attempt + 1))
    else:
        return site, f"failed: {err}"
    with open(out + ".tmp", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["date", "time_end", "avg_mph", "volume"])
        for x in rows:
            w.writerow([x["Report Date"][:10], x["Time Period Ending"],
                        x["Avg mph"], x["Total Volume"]])
    os.replace(out + ".tmp", out)
    return site, len(rows)
sites = [r["site_id"] for r in csv.DictReader(open("sites_selected.csv"))]
with ThreadPoolExecutor(4) as ex:
    for i, (s, n) in enumerate(ex.map(get, sites)):
        if i % 25 == 0 or not isinstance(n, int): print(i, s, n, flush=True)
print("done")
