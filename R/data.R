# Data loaders for the Paris and London applications (Section 5 of the paper).
# Each returns a list with the metric graph (observation locations as
# vertices), the response matrix Y (locations x days), the mean covariates
# (easting, northing), coordinates in km and geodesic distances in km.
# The inputs are in data/ (see README.md and data_prep/).

# London: daily 17:00-18:00 mean speeds (mph) at National Highways MIDAS loop
# sites (WebTRIS) on the M25 and the radial motorways, on the OpenStreetMap
# motorway and trunk network (with link roads) of the box -0.65..0.40 E,
# 51.22..51.78 N; see data_prep/london/. The largest connected
# component of the network is used.
load_london <- function(dir) {
  lines <- readRDS(file.path(dir, "osm_london.rds"))
  lines <- suppressWarnings(sf::st_cast(lines, "LINESTRING", warn = FALSE))
  edges <- lapply(sf::st_geometry(lines), function(g)
    unname(sf::st_coordinates(g)[, 1:2, drop = FALSE]))
  keep <- vapply(edges, nrow, integer(1)) >= 2L
  graph <- MetricGraph::metric_graph$new(edges = edges[keep], longlat = TRUE,
                                         verbose = 0, check_connected = FALSE)
  comps <- graph$get_components()
  graph <- comps[[which.max(vapply(comps, function(x) sum(x$edge_lengths), 1))]]
  raw <- read.csv(file.path(dir, "panel_1700.csv"), check.names = FALSE)
  sites <- read.csv(file.path(dir, "sites_selected.csv"),
                    colClasses = c(site_id = "character"))
  sites <- sites[match(colnames(raw)[-1], sites$site_id), ]
  Y <- t(as.matrix(raw[-1])); storage.mode(Y) <- "double"
  graph$add_observations(
    data = data.frame(y = rowMeans(Y), coord_x = sites$lon, coord_y = sites$lat,
                      site_index = seq_len(nrow(Y))),
    data_coords = "spatial", tolerance = 0.1, verbose = 0)
  graph$observation_to_vertex()
  # Merge the degree-2 vertices without data (OSM shape points) into longer
  # edges: same metric graph (15,325 -> 3,734 vertices, same length), but a
  # much smaller engine (2 GB instead of 24 GB of memory per process).
  graph$prune_vertices(verbose = 0)
  graph$observation_to_vertex()   # refresh PtV after pruning
  gd <- graph$get_data()
  ord <- gd$site_index
  stopifnot(!anyNA(ord), !anyDuplicated(ord))
  Y <- Y[ord, , drop = FALSE]
  stopifnot(max(abs(rowMeans(Y) - gd$y)) < 1e-8)
  xy <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(gd, coords = c(".coord_x", ".coord_y"), crs = 4326),
    27700)) / 1000
  graph$compute_geodist(obs = TRUE, verbose = 0)
  geo <- graph$geo_dist[[1]]
  stopifnot(all(dim(geo) == nrow(Y)), all(is.finite(geo)))
  list(name = "london", graph = graph, Y = unname(Y),
       covariates = cbind(easting = xy[, 1], northing = xy[, 2]),
       xy_km = unname(xy), geo_km = unname(geo), unit = "mph",
       site_id = sites$site_id[ord], dates = as.character(raw[[1]]))
}

# Paris: hourly counts 08:00-09:00 on working days at the permanent loop
# counters of the City of Paris (opendata.paris.fr), on the network of the
# counted arcs themselves (one edge per pair of nodes, Lambert-93 km); see
# data_prep/paris/make_panel.py. variable = "q" (flow, veh/h) or "k"
# (occupancy, %); log = TRUE models the logarithm.
load_paris <- function(dir, variable = "q", log = TRUE) {
  ed <- read.csv(file.path(dir, "network_edges.csv"))
  edges <- lapply(strsplit(ed$coords, ";"), function(p)
    do.call(rbind, lapply(strsplit(p, " "), as.numeric)))
  graph <- MetricGraph::metric_graph$new(edges = edges, longlat = FALSE,
                                         verbose = 0, check_connected = FALSE)
  comps <- graph$get_components()
  graph <- comps[[which.max(vapply(comps, function(x) sum(x$edge_lengths), 1))]]
  raw <- read.csv(file.path(dir, paste0("panel_", variable, ".csv")), check.names = FALSE)
  sites <- read.csv(file.path(dir, "sites.csv"), colClasses = c(iu_ac = "character"))
  sites <- sites[match(colnames(raw)[-1], sites$iu_ac), ]
  # Site position: midpoint of its edge polyline.
  mid <- t(vapply(sites$edge, function(e) {
    co <- edges[[e]]; L <- c(0, cumsum(sqrt(rowSums(diff(co)^2)))); h <- L[length(L)] / 2
    j <- max(which(L <= h)); j <- min(j, nrow(co) - 1L); w <- (h - L[j]) / (L[j + 1] - L[j])
    co[j, ] + w * (co[j + 1, ] - co[j, ])
  }, numeric(2)))
  Y <- t(as.matrix(raw[-1])); storage.mode(Y) <- "double"
  if (log) Y <- log(Y)
  graph$add_observations(
    data = data.frame(y = rowMeans(Y), coord_x = mid[, 1], coord_y = mid[, 2],
                      site_index = seq_len(nrow(Y))),
    data_coords = "spatial", tolerance = 0.02, verbose = 0)
  graph$observation_to_vertex()
  gd <- graph$get_data()
  ord <- gd$site_index
  stopifnot(!anyNA(ord), !anyDuplicated(ord))
  if (length(ord) < nrow(Y)) message(nrow(Y) - length(ord), " sites outside the main component")
  Y <- Y[ord, , drop = FALSE]
  stopifnot(max(abs(rowMeans(Y) - gd$y)) < 1e-8)
  xy <- cbind(gd$.coord_x, gd$.coord_y)
  graph$compute_geodist(obs = TRUE, verbose = 0)
  geo <- graph$geo_dist[[1]]
  stopifnot(all(dim(geo) == nrow(Y)), all(is.finite(geo)))
  list(name = paste0("paris_", variable), graph = graph, Y = unname(Y),
       covariates = cbind(easting = xy[, 1], northing = xy[, 2]),
       xy_km = unname(xy), geo_km = unname(geo),
       unit = if (variable == "q") "log veh/h" else "log %",
       site_id = sites$iu_ac[ord], dates = as.character(raw[[1]]))
}

load_dataset <- function(name, root = ".") {
  switch(name,
         paris_k = load_paris(file.path(root, "data", "paris"), "k"),
         paris_q = load_paris(file.path(root, "data", "paris"), "q"),
         london = load_london(file.path(root, "data", "london")),
         stop("Unknown data set: ", name))
}
