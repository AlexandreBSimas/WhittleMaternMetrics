#!/usr/bin/env Rscript
# Figures 1 and 3 of the paper. Run from the repository root; writes
#   output/metric_multigraph.png   Figure 1: d_kappa(s0, .) for kappa in {0, 1, 5}
#                                  and the geodesic distance d(s0, .) on a
#                                  multigraph with three parallel edges and a
#                                  self-loop, with a common colour bar.
#   output/samples_isotropic.png   Figure 3: two realisations each of a centred
#                                  Gaussian field with Matern covariance in d_0,
#                                  for smoothness 1/4 and 1/2.
# Geodesic distances come from MetricGraph's compute_geodist_mesh and the
# Matern covariance from rSPDE::matern.covariance.
# Needs MetricGraph, rSPDE, Matrix, MASS, ggplot2, viridis and patchwork.
library(MetricGraph)
library(rSPDE)
library(Matrix)
library(MASS)
library(ggplot2)
library(viridis)
library(patchwork)

set.seed(42)

out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE)
message("Writing figures to: ", out_dir)

## ----------------------------------------------------------------------------
## Helpers shared by both figures
## ----------------------------------------------------------------------------

## Generalized graph Laplacian Q (eq:Q0 for kappa = 0; eq:expprec for kappa > 0)
## evaluated on the mesh of a metric_graph object.
build_Q <- function(graph, kappa = 0) {
  endpts  <- graph$mesh$E    # (n_seg x 2) endpoint indices
  seg_len <- graph$mesh$h_e  # length of each mesh segment
  M       <- nrow(graph$mesh$V)

  ii <- jj <- xx <- numeric(0)
  for (e in seq_len(nrow(endpts))) {
    a <- endpts[e, 1]; b <- endpts[e, 2]; le <- seg_len[e]
    if (kappa == 0) {
      if (a == b) next   # loop contributes nothing
      ii <- c(ii, a, b, a, b)
      jj <- c(jj, a, b, b, a)
      xx <- c(xx, 1/le, 1/le, -1/le, -1/le)
    } else {
      if (a == b) {
        ii <- c(ii, a); jj <- c(jj, a)
        xx <- c(xx, 2 * kappa * tanh(kappa * le / 2))
      } else {
        d_term <- 2 * kappa * (0.5 + exp(-2 * kappa * le) /
                                       (1 - exp(-2 * kappa * le)))
        o_term <- -2 * kappa * exp(-kappa * le) / (1 - exp(-2 * kappa * le))
        ii <- c(ii, a, b, a, b)
        jj <- c(jj, a, b, b, a)
        xx <- c(xx, d_term, d_term, o_term, o_term)
      }
    }
  }
  Q <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(M, M))
  (Q + t(Q)) / 2
}

## Variogram matrix d_kappa[i,j] = Sigma[i,i] + Sigma[j,j] - 2*Sigma[i,j],
## where Sigma = Q^{-1} (kappa>0) or Q^{+} (kappa=0).
variogram_from_Q <- function(Q, kappa) {
  if (kappa == 0) {
    Sigma <- as.matrix(MASS::ginv(as.matrix(Q)))
  } else {
    Sigma <- as.matrix(solve(Q))
  }
  diag_S <- diag(Sigma)
  outer(diag_S, diag_S, `+`) - 2 * Sigma
}

## ============================================================================
## Figure 1: multigraph with parallel edges and a loop -- d_kappa and geodesic
## ============================================================================

edge_straight <- cbind(seq(0, 1, length.out = 80), 0)

a <- seq(pi, 0, length.out = 80)
edge_upper <- cbind(0.5 + 0.5 * cos(a),  0.5 * sin(a))
edge_lower <- cbind(0.5 + 0.5 * cos(a), -0.75 * sin(a) - 0.05 * sin(2*a))

b <- seq(0, 2 * pi, length.out = 81)
edge_loop2 <- cbind(-0.25 + 0.25 * cos(b), 0.25 * sin(b))
edge_loop2[1, ]                  <- c(0, 0)
edge_loop2[nrow(edge_loop2), ]   <- c(0, 0)

theta_graph <- metric_graph$new(
  edges = list(edge_straight, edge_upper, edge_lower, edge_loop2)
)
theta_graph$build_mesh(h = 0.02)

s0_idx <- 1L
kappas <- c(0, 1, 5)
metric_vals <- list()
for (k in kappas) {
  metric_vals[[as.character(k)]] <-
    variogram_from_Q(build_Q(theta_graph, kappa = k), kappa = k)[s0_idx, ]
}
## Geodesic distance on the mesh from the source vertex s0_idx, computed by
## MetricGraph's built-in routine.  `compute_geodist_mesh` populates the
## mesh$geo_dist element of the graph object with the full mesh-to-mesh
## geodesic-distance matrix; we extract the row corresponding to s0.
theta_graph$compute_geodist_mesh()
metric_vals[["geo"]] <- as.numeric(theta_graph$mesh$geo_dist[s0_idx, ])

## A common colour-scale across all four panels.  Because d_5 is roughly an
## order of magnitude smaller than the geodesic, a linear common scale would
## make the d_5 panel uniformly dark.  We therefore use a square-root colour
## transformation (`trans = "sqrt"`), which compresses the high end and lets
## the structure of all four panels remain visible while keeping a single
## shared colour bar.
metric_vmax <- max(unlist(metric_vals))

plot_metric_panel <- function(X, title) {
  theta_graph$plot_function(X = X, plotly = FALSE) +
    scale_colour_viridis_c(name = "metric value",
                           limits = c(0, metric_vmax),
                           trans  = "sqrt") +
    coord_equal() +
    theme_minimal(base_size = 14) +
    labs(title = title) +
    theme(axis.text  = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())
}

p0 <- plot_metric_panel(metric_vals[["0"]],   expression(d[0](s[0], "."  )))
p1 <- plot_metric_panel(metric_vals[["1"]],   expression(d[1](s[0], "."  )))
p5 <- plot_metric_panel(metric_vals[["5"]],   expression(d[5](s[0], "."  )))
pg <- plot_metric_panel(metric_vals[["geo"]], expression(d  (s[0], "."  )))

## Single-row 1x4 layout with a common colour bar on the right.  The bar
## height is set to span the panels' plotting area without overflowing.
fig1 <- (p0 | p1 | p5 | pg) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(colour = guide_colourbar(barheight = unit(4.0, "cm")))
ggsave(file.path(out_dir, "metric_multigraph.png"),
       fig1, width = 14, height = 4.5, dpi = 300)


## ============================================================================
## Figure 3: sample paths of an isotropic Gaussian field with Matern covariance
## ============================================================================

D0_th <- variogram_from_Q(build_Q(theta_graph, kappa = 0), kappa = 0)
D0_th <- (D0_th + t(D0_th)) / 2

settings <- list(
  list(nu = 0.25, beta = 2, label = expression(nu == 1/4 * "," ~ beta == 2)),
  list(nu = 0.5,  beta = 2, label = expression(nu == 1/2 * "," ~ beta == 2))
)

set.seed(7)
plots <- list()
for (col in seq_along(settings)) {
  s <- settings[[col]]
  ## Matern covariance via rSPDE.  rSPDE uses the convention
  ##   r(h) = sigma^2 * 2^{1-nu}/Gamma(nu) * (kappa h)^nu K_nu(kappa h),
  ## so its `kappa` plays the role of our `beta` and `sigma` is the marginal
  ## standard deviation (sigma^2 is the marginal variance).
  Sigma <- rSPDE::matern.covariance(D0_th, kappa = s$beta,
                                    nu = s$nu, sigma = 1)
  Sigma <- as.matrix(Sigma) + 1e-8 * diag(nrow(D0_th))
  L <- chol(Sigma)
  for (row in 1:2) {
    z <- rnorm(nrow(Sigma))
    u <- as.numeric(t(L) %*% z)
    p <- theta_graph$plot_function(X = u, plotly = FALSE) +
      scale_colour_gradient2(low = "navy", mid = "white",
                             high = "firebrick", midpoint = 0,
                             limits = c(-2.5, 2.5), name = "u(s)") +
      coord_equal() +
      theme_minimal(base_size = 14) +
      labs(title = s$label) +
      theme(axis.text  = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank())
    plots[[length(plots) + 1]] <- p
  }
}

## Single-row 1x4 layout: panels 1-2 are two realisations with nu = 1/4,
## panels 3-4 are two realisations with nu = 1/2; titles repeat above each
## panel for clarity.
fig2 <- (plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(colour = guide_colourbar(barheight = unit(4.0, "cm")))
ggsave(file.path(out_dir, "samples_isotropic.png"),
       fig2, width = 14, height = 4.5, dpi = 300)

cat("Done. Wrote two figures to", out_dir, "\n")
