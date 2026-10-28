#!/usr/bin/env Rscript
# Figure 4 of the paper: road network (black) and observation sites coloured
# by the mean over the daily maps (Esri light grey basemap, viridis scale). Paris: occupancy in %,
# averaged on the log scale; London: speed in mph.
# Needs leaflet, webshot2, Google Chrome and internet access for the tiles.
# Run from the repository root. Writes output/{paris_occupancy,london_speed}.png.
suppressPackageStartupMessages({library(MetricGraph); library(leaflet); library(sf)})
source("R/data.R")
art <- "output"; dir.create(art, showWarnings = FALSE)

edges_sf <- function(g, crs, scale = 1) {
  ls <- lapply(g$edges, function(e) st_linestring(e[, 1:2] * scale))
  st_transform(st_sfc(ls, crs = crs), 4326)
}
make_map <- function(ds, value, unit, edges, out) {
  gd <- ds$graph$get_data()
  ll <- if (ds$name == "london") cbind(gd$.coord_x, gd$.coord_y) else
    st_coordinates(st_transform(st_sfc(lapply(seq_len(nrow(gd)), function(i)
      st_point(c(gd$.coord_x[i], gd$.coord_y[i]) * 1000)), crs = 2154), 4326))
  pal <- colorNumeric(viridis::viridis(100, option = "D"), domain = value)
  m <- leaflet(options = leafletOptions(zoomControl = FALSE)) |>
    addProviderTiles(providers$Esri.WorldGrayCanvas) |>
    addTiles(urlTemplate = paste0("https://server.arcgisonline.com/ArcGIS/rest/",
      "services/Canvas/World_Light_Gray_Reference/MapServer/tile/{z}/{y}/{x}"),
      attribution = "Esri") |>
    addPolylines(data = edges, color = "black", weight = 1.2, opacity = 1) |>
    addCircleMarkers(lng = ll[, 1], lat = ll[, 2], radius = 4, stroke = FALSE,
                     fillOpacity = 0.85, fillColor = pal(value)) |>
    addLegend("topright", pal = pal, values = value, title = unit, opacity = 1) |>
    fitBounds(min(ll[, 1]), min(ll[, 2]), max(ll[, 1]), max(ll[, 2]))
  html <- file.path(tempdir(), paste0(ds$name, "_map.html"))
  htmlwidgets::saveWidget(m, html, selfcontained = FALSE)
  webshot2::webshot(html, out, vwidth = 992, vheight = 744, delay = 4)
  message("Wrote ", out)
}

ds <- load_dataset("paris_k")
make_map(ds, exp(rowMeans(ds$Y)), "%", edges_sf(ds$graph, 2154, 1000),
         file.path(art, "paris_occupancy.png"))
ds <- load_dataset("london")
make_map(ds, rowMeans(ds$Y), "mph", edges_sf(ds$graph, 4326),
         file.path(art, "london_speed.png"))
