#!/usr/bin/env Rscript
# Cross-validation of the nine models of the Paris and London applications.
#
# Usage (from the repository root):
#   Rscript run_cv.R <dataset> <design> [cores] [repeats] [models] [pcts]
#     dataset  paris_k (log occupancy, paper) | paris_q (log flow) | london
#     design   lpo  (leave 75% and 90% of the sites out, `repeats` draws each)
#              lro  (buffered leave-region-out; not used in the paper)
#     models   all (default) | nofem (skip the slow FEM model; it can be
#              added later to the same splits by rerunning with `all`)
#     pcts     lpo only: comma-separated subset of 75,90 to run. The draws
#              are always generated for both, so split i is the same split
#              whichever subset is run.
# Results are written to output/<dataset>_<design>/. Each split is saved as
# it finishes, so an interrupted run resumes where it stopped.
suppressPackageStartupMessages({
  library(MetricGraph); library(Matrix); library(rSPDE); library(parallel)
})
source("R/data.R")
source("R/models.R")
source("R/design.R")

SETTINGS <- list(
  # Paris is about 10 km across: shorter buffers.
  paris_k = list(K = 8L, buffers = c(0, 0.5, 1, 2), min_train = 30L),
  paris_q = list(K = 8L, buffers = c(0, 0.5, 1, 2), min_train = 30L),
  london = list(K = 8L, buffers = c(0, 2, 5, 10), min_train = 30L))

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1L) args[[1]] else "paris_k"
design <- if (length(args) >= 2L) args[[2]] else "lpo"
cores <- if (length(args) >= 3L) as.integer(args[[3]]) else 1L
repeats <- if (length(args) >= 4L) as.integer(args[[4]]) else 25L
model_set <- if (length(args) >= 5L) args[[5]] else "all"
run_pcts <- if (length(args) >= 6L) as.numeric(strsplit(args[[6]], ",")[[1]]) else c(75, 90)
stopifnot(model_set %in% c("all", "nofem"))
stopifnot(dataset %in% names(SETTINGS), design %in% c("lpo", "lro"))
RUN_MODELS <- MODELS

out_dir <- file.path("output", paste(dataset, design, sep = "_"))
dir.create(file.path(out_dir, "splits"), recursive = TRUE, showWarnings = FALSE)

message("Loading ", dataset)
ds <- load_dataset(dataset)
message("Building metric engine (n = ", nrow(ds$Y), ", T = ", ncol(ds$Y), ")")
MESH_H <- fem_mesh_h(dataset)
RUN_FEM_VERSION <- fem_version(MESH_H)
eng <- make_engine(ds, mesh_h = MESH_H, fem = model_set == "all")
n <- nrow(ds$Y)
cfg <- SETTINGS[[dataset]]
splits <- if (design == "lpo") {
  splits_lpo(n, pcts = c(75, 90), repeats = repeats, seed = 20261022L)
} else {
  splits_lro(ds$xy_km, ds$geo_km, K = cfg$K, buffers = cfg$buffers,
             seed = 20261022L, min_train = cfg$min_train)
}
saveRDS(list(splits = splits, region = attr(splits, "region"),
             site_id = ds$site_id, xy_km = ds$xy_km),
        file.path(out_dir, "splits.rds"))
message(length(splits), " splits")

run_split <- function(i) {
  file <- file.path(out_dir, "splits", sprintf("split_%03d.rds", i))
  sp <- splits[[i]]
  cached <- NULL
  if (file.exists(file)) {
    cached <- readRDS(file)
    if (!identical(cached$test, sp$test) || !identical(cached$train, sp$train))
      cached <- NULL
    else if (model_set == "nofem")
      return(cached$scores[cached$scores$model != "wm_fem", ])
    else if (identical(cached$fem_version, RUN_FEM_VERSION))
      return(cached$scores)
  }
  train <- sp$train; test <- sp$test
  # Responses and covariates are standardised with training data only.
  y_center <- mean(ds$Y[train, ])
  y_scale <- sd(as.vector(ds$Y[train, ]))
  Y <- (ds$Y - y_center) / y_scale
  Z <- ds$covariates
  Z <- sweep(sweep(Z, 2, colMeans(Z[train, , drop = FALSE])), 2,
             pmax(apply(Z[train, , drop = FALSE], 2, sd), 1e-8), "/")
  X <- cbind(1, Z)
  nearest <- apply(ds$geo_km[test, train, drop = FALSE], 1, min)
  # A split cached under an older FEM setting only needs WM1 (its warm
  # start) and FEM refitted; the other models are unchanged.
  res <- if (model_set == "nofem")
    fit_all(eng, Y, X, train, test, RUN_MODELS[RUN_MODELS$name != "wm_fem", ]) else
    if (is.null(cached)) fit_all(eng, Y, X, train, test, RUN_MODELS) else
    fit_all(eng, Y, X, train, test,
            MODELS[MODELS$name %in% c("wm1", "wm_fem"), ])["wm_fem"]
  rows <- lapply(names(res), function(nm) {
    model <- as.list(ALL_MODELS[ALL_MODELS$name == nm, ])
    r <- res[[nm]]
    cbind(data.frame(dataset = dataset, design = sp$design, split = i,
                     pct = sp$pct, repeat_id = sp$repeat_id,
                     region = sp$region, buffer_km = sp$buffer_km,
                     model = nm, n_train = length(train), n_test = length(test),
                     nearest_median_km = median(nearest),
                     nearest_min_km = min(nearest),
                     train_nll = r$fit$nll, convergence = r$fit$convergence,
                     kappa = model_kappa(model, r$fit$par),
                     alpha = model_alpha(model, r$fit$par),
                     seconds = r$seconds),
          score_split(r$pred, ds$Y[test, , drop = FALSE], y_center, y_scale))
  })
  out <- do.call(rbind, rows)
  if (!is.null(cached))
    out <- rbind(cached$scores[cached$scores$model != "wm_fem", ], out)
  saveRDS(list(test = test, train = train, scores = out,
               fem_version = if (model_set == "nofem") NA else RUN_FEM_VERSION),
          file)
  message(sprintf("split %d/%d done (%s, %d train / %d test)", i,
                  length(splits), sp$design, length(train), length(test)))
  out
}

t0 <- Sys.time()
todo <- seq_along(splits)
if (design == "lpo")
  todo <- todo[vapply(splits, `[[`, numeric(1), "pct") %in% run_pcts]
res <- mclapply(todo, function(i)
  tryCatch(run_split(i), error = function(e) {
    message("split ", i, " failed: ", conditionMessage(e)); NULL }),
  mc.cores = cores, mc.preschedule = FALSE)
scores <- do.call(rbind, res)
write.csv(scores, file.path(out_dir, if (model_set == "nofem")
  "scores_nofem.csv" else "scores.csv"), row.names = FALSE)
message("Finished ", sum(!vapply(res, is.null, logical(1))), "/",
        length(todo), " splits in ", format(Sys.time() - t0))
