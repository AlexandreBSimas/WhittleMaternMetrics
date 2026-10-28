# OpenStreetMap motorways and trunk roads (with link roads) in the London box.
# Downloads raw/osm_london.json from an Overpass server (if absent) and
# converts it to raw/osm_london.rds (an sf object of LINESTRINGs).
# Usage: Rscript fetch_osm.R
suppressPackageStartupMessages({library(sf); library(jsonlite)})
json <- "raw/osm_london.json"
if (!file.exists(json)) {
  q <- paste0('[out:json][timeout:600];',
              'way["highway"~"^(motorway|motorway_link|trunk|trunk_link)$"]',
              '(51.22,-0.65,51.78,0.40);out geom tags;')
  status <- system2("curl", c("-s", "--max-time", "700", "-o", json,
                              "--data-urlencode", shQuote(paste0("data=", q)),
                              "https://overpass.private.coffee/api/interpreter"))
  stopifnot(status == 0)
}
el <- fromJSON(json, simplifyVector = FALSE)$elements
tag <- function(e, k) if (is.null(e$tags[[k]])) NA_character_ else e$tags[[k]]
geoms <- lapply(el, function(e)
  st_linestring(do.call(rbind, lapply(e$geometry, function(p) c(p$lon, p$lat)))))
lines <- st_sf(osm_id = vapply(el, function(e) as.character(e$id), ""),
               highway = vapply(el, tag, "", k = "highway"),
               ref = vapply(el, tag, "", k = "ref"),
               maxspeed = vapply(el, tag, "", k = "maxspeed"),
               geometry = st_sfc(geoms, crs = 4326))
saveRDS(lines, "raw/osm_london.rds")
print(table(lines$highway))
cat("total km:", round(sum(as.numeric(st_length(lines))) / 1000), "\n")
