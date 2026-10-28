#!/usr/bin/env Rscript
# Fit the nine models to all sites (Table 1 of the paper), with the same code,
# starting values and FEM mesh as the cross-validation. The FEM model is also
# fitted on the finer h = 0.1 km mesh to show the effect of the mesh size.
# Usage: Rscript fit_fulldata.R [paris_k | london]
suppressPackageStartupMessages({library(MetricGraph); library(Matrix); library(rSPDE)})
source("R/data.R"); source("R/models.R")
dataset <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "paris_k"
ds <- load_dataset(dataset)
eng <- make_engine(ds, mesh_h = fem_mesh_h(dataset))
n <- nrow(ds$Y); all <- seq_len(n)
# Standardised as in the CV; the log-likelihood is mapped back to the
# original scale by adding N log(sd).
yc <- mean(ds$Y); ys <- sd(as.vector(ds$Y)); Y <- (ds$Y - yc) / ys
X <- cbind(1, scale(ds$covariates))
n_par <- c(exp_d0 = 5, cauchy_d0 = 7, matern_d0 = 6, exp_dkappa = 6,
           cauchy_dkappa = 8, matern_dkappa = 7, wm1 = 5, wm_fem = 6,
           exp_dres = 5) + (ncol(X) - 2)      # Table 1 counts two mean terms
fits <- list(); rows <- list()
fit_one <- function(nm, engine, warm, exp_warm = NULL) {
  model <- as.list(MODELS[MODELS$name == nm, ])
  t0 <- proc.time()[[3]]
  fit <- fit_model(model, engine, Y, X, all,
                   make_starts(model, engine, Y, X, all, warm, exp_warm))
  nll <- fit$nll + length(Y) * log(ys)
  list(fit = fit, row = data.frame(
    model = nm, p = n_par[[nm]], neg_loglik = nll,
    kappa = model_kappa(model, fit$par), alpha = model_alpha(model, fit$par),
    convergence = fit$convergence, seconds = proc.time()[[3]] - t0))
}
for (nm in c("exp_d0", "exp_dres", "wm1", "exp_dkappa", "cauchy_d0",
             "matern_d0", "cauchy_dkappa", "matern_dkappa", "wm_fem")) {
  warm <- switch(nm, exp_dkappa = fits$exp_d0$par,
                 cauchy_dkappa = fits$cauchy_d0$par,
                 matern_dkappa = fits$matern_d0$par, wm_fem = fits$wm1$par, NULL)
  exp_warm <- switch(nm, cauchy_d0 = , matern_d0 = fits$exp_d0$par,
                     cauchy_dkappa = , matern_dkappa = fits$exp_dkappa$par, NULL)
  r <- fit_one(nm, eng, warm, exp_warm)
  fits[[nm]] <- r$fit; rows[[nm]] <- r$row
  print(r$row)
}
fine <- fit_one("wm_fem", make_engine(ds, mesh_h = 0.1), fits$wm1$par)$row
fine$model <- "wm_fem_h0.1"
out <- rbind(do.call(rbind, rows), fine)
out <- out[order(match(out$model, c(MODELS$name, "wm_fem_h0.1"))), ]
dir.create("output", showWarnings = FALSE)
write.csv(out, file.path("output", paste0(dataset, "_fulldata_fit.csv")),
          row.names = FALSE)
print(out, digits = 7, row.names = FALSE)
