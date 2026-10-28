# The nine covariance models of the Paris and London applications (Table 1
# of the paper), implemented so that every parameter can be re-estimated from an
# arbitrary training subset and predictions at held-out sites use the exact
# joint Gaussian conditional distribution.
#
# All models share the structure
#   y_t(s) = x(s)^T beta + u_t(s) + e_t(s),  t = 1, ..., T,
# with independent replicates u_t of the latent field, a mean shared across
# replicates, and Gaussian noise with variance sigma_e^2. beta is profiled out
# by generalised least squares, so the optimiser only sees covariance
# parameters. Parameters are estimated by maximum likelihood.

MODELS <- data.frame(
  name = c("exp_d0", "cauchy_d0", "matern_d0",
           "exp_dkappa", "cauchy_dkappa", "matern_dkappa",
           "wm1", "wm_fem", "exp_dres"),
  family = c("exp", "cauchy", "matern", "exp", "cauchy", "matern",
             "wm1", "fem", "exp"),
  metric = c("d0", "d0", "d0", "dkappa", "dkappa", "dkappa",
             NA, NA, "dres"),
  label = c("$r(d_0)$, exponential", "$r(d_0)$, generalised Cauchy",
            "$r(d_0)$, Mat\\'ern",
            "$r(d_\\kappa)$, exponential", "$r(d_\\kappa)$, generalised Cauchy",
            "$r(d_\\kappa)$, Mat\\'ern",
            "Whittle--Mat\\'ern $\\alpha=1$",
            "Whittle--Mat\\'ern FEM, estimated $\\alpha$",
            "$r(d_{\\mathrm{res}})$, exponential"),
  stringsAsFactors = FALSE)

# Extra model, not used in the paper: the exponential in d_0 plus a constant
# covariance c, i.e. an independent random level per map.
EXTRA_MODELS <- data.frame(name = "exp_d0c", family = "expc", metric = "d0",
                           label = "$r(d_0)+c$, exponential",
                           stringsAsFactors = FALSE)
ALL_MODELS <- rbind(MODELS, EXTRA_MODELS)

# kappa (km^-1) is estimated without bounds. KAPPA_START = 1e-3 (screening
# length 1000 km) is numerically d_0; it is the first starting value of the
# d_kappa models and is used in the sanity check of d_kappa -> d_0.
KAPPA_START <- 1e-3
NU_FEM_MAX <- 2      # alpha = nu + 1/2 in (1/2, 5/2)
FEM_MESH_H <- 0.25   # km; FEM mesh, used for both the full-data fit and the CV
FEM_VERSION <- "single-start-h0.25"   # identifies the FEM fitting settings
# FEM mesh per data set (0.25 km for both data sets).
fem_mesh_h <- function(dataset) FEM_MESH_H
fem_version <- function(h) if (h == FEM_MESH_H) FEM_VERSION else paste0("single-start-h", h)

# ---------------------------------------------------------------------------
# Metric engine: everything that depends only on the graph and the locations.
# ---------------------------------------------------------------------------
make_engine <- function(ds, mesh_h = FEM_MESH_H, fem = TRUE) {
  g <- ds$graph
  n <- length(g$PtV)
  S <- Matrix::sparseMatrix(i = g$PtV, j = seq_len(n), x = 1,
                            dims = c(g$nV, n))
  variogram <- function(Sigma) {
    D <- outer(diag(Sigma), diag(Sigma), "+") - 2 * Sigma
    D <- pmax((D + t(D)) / 2, 0)
    diag(D) <- 0
    D
  }
  # d_0 from the generalised graph Laplacian with conductances 1/l_e,
  # grounded at a vertex without an observation (Section 3 of the paper).
  uv <- g$E
  keep <- uv[, 1] != uv[, 2]
  u <- uv[keep, 1]; v <- uv[keep, 2]; w <- 1 / g$edge_lengths[keep]
  Q0 <- Matrix::sparseMatrix(i = c(u, u, v, v), j = c(u, v, u, v),
                             x = c(w, -w, -w, w), dims = c(g$nV, g$nV))
  ground <- setdiff(seq_len(g$nV), g$PtV)[1]
  idx <- setdiff(seq_len(g$nV), ground)
  Sidx <- S[idx, , drop = FALSE]
  D0 <- variogram(as.matrix(Matrix::crossprod(Sidx,
                    Matrix::solve(Q0[idx, idx], Sidx))))
  # Q^kappa (alpha = 1) assembled from the edge lengths; identical to
  # MetricGraph::spde_precision(kappa, tau = 1, alpha = 1, BC) (checked to
  # relative error 1e-11), but faster, and the symbolic Cholesky factor is
  # reused across kappa. BC = 1 adds kappa on the diagonal at the vertices of
  # degree one.
  lp <- uv[, 1] == uv[, 2]
  ue <- uv[!lp, 1]; ve <- uv[!lp, 2]; le <- g$edge_lengths[!lp]
  ul <- uv[lp, 1]; ll <- g$edge_lengths[lp]
  deg1 <- which(tabulate(c(uv[, 1], uv[, 2]), g$nV) == 1L)
  Qk <- function(k, BC = 0) {
    a <- k / tanh(k * le); b <- -k / sinh(k * le)
    d1 <- if (BC == 1) deg1 else integer(0)
    Matrix::sparseMatrix(i = c(ue, ve, ue, ve, ul, d1),
                         j = c(ue, ve, ve, ue, ul, d1),
                         x = c(a, a, b, b, 2 * k * tanh(k * ll / 2),
                               rep(k, length(d1))),
                         dims = c(g$nV, g$nV))
  }
  chol0 <- Matrix::Cholesky(Matrix::forceSymmetric(Qk(0.1, BC = 1)), LDL = FALSE)
  sigma_wm <- function(kappa, BC) {
    ch <- Matrix::update(chol0, Matrix::forceSymmetric(Qk(kappa, BC)))
    as.matrix(Matrix::crossprod(S, Matrix::solve(ch, S, system = "A")))
  }
  D_kappa <- function(kappa) variogram(sigma_wm(kappa, BC = 0))
  # Resistance distance computed directly on the supplied representation.
  gr <- g$clone()
  gr$compute_resdist(verbose = 0)
  Dres <- unname(as.matrix(gr$res_dist[[1]]))
  stopifnot(all(dim(Dres) == n))
  eng <- list(n = n, D0 = D0, Dres = Dres, D_kappa = D_kappa,
              Sigma_wm1 = function(kappa) sigma_wm(kappa, BC = 1))
  # Sanity checks: d_kappa decreases in kappa and tends to d_0.
  Dsmall <- D_kappa(KAPPA_START)
  if (max(abs(Dsmall - D0)) > 1e-2 * max(D0))
    stop("d_kappa at kappa = 1e-3 does not approximate d_0")
  if (max(D_kappa(0.5) - D_kappa(0.05)) > 1e-8)
    stop("d_kappa is not monotone in kappa")
  if (fem) {
    gm <- g$clone()
    gm$build_mesh(h = mesh_h)
    op <- rSPDE::matern.operators(alpha = 1.5, kappa = 1, tau = 1, m = 2,
                                  graph = gm, parameterization = "spde")
    A0 <- gm$fem_basis(gm$get_PtE())
    A <- kronecker(matrix(1, 1, nrow(op$Q) / ncol(A0)), A0)
    eng$Sigma_fem <- function(nu, kappa) {
      opk <- update(op, alpha = nu + 0.5, kappa = kappa, tau = 1,
                    check_stationarity = FALSE)
      as.matrix(A %*% Matrix::solve(opk$Q, Matrix::t(A)))
    }
  }
  eng
}

# ---------------------------------------------------------------------------
# Latent covariance at all n sites for a model and covariance parameters.
# theta[1] is always log sigma_e; theta[2] log of the marginal variance.
# ---------------------------------------------------------------------------
radial <- function(family, D, th) {
  b <- exp(th[1])
  switch(family,
         exp = exp(-b * D),
         cauchy = {
           a <- th[2]; xi <- exp(th[3])          # a in (0, 1], xi > 0
           exp(-(xi / a) * log1p(b * D^a))     # (1 + b D^a)^(-xi / a)
         },
         matern = {
           nu <- th[2]                           # nu in (0, 1/2]
           C <- matrix(rSPDE::matern.covariance(as.vector(D), b, nu, 1),
                       nrow = nrow(D))
           (C + t(C)) / 2
         })
}

n_radial <- c(exp = 1L, cauchy = 3L, matern = 2L)

latent_cov <- function(model, theta, eng) {
  s2 <- exp(theta[2])
  fam <- model$family
  if (fam == "wm1") {
    Sig <- eng$Sigma_wm1(exp(theta[3]))
    return(s2 * Sig / median(diag(Sig)))
  }
  if (fam == "fem") {
    Sig <- eng$Sigma_fem(NU_FEM_MAX * plogis(theta[4]), exp(theta[3]))
    Sig <- (Sig + t(Sig)) / 2
    return(s2 * Sig / median(diag(Sig)))
  }
  if (fam == "expc")
    return(s2 * exp(-exp(theta[3]) * eng$D0) + exp(theta[4]))
  k <- n_radial[[fam]]
  th <- theta[2 + seq_len(k)]
  D <- switch(model$metric,
              d0 = eng$D0, dres = eng$Dres,
              dkappa = eng$D_kappa(exp(theta[3 + k])))
  s2 * radial(fam, D, th)
}

bounds <- function(model) {
  fam <- model$family
  lo <- c(log(1e-3), log(1e-3)); hi <- c(log(3), log(50))
  if (fam %in% names(n_radial)) {
    lo <- c(lo, -10); hi <- c(hi, 6)                    # log decay
    # The upper bounds include the exponential model: nu = 1/2 exactly for
    # the Matern, and a = 1 with xi -> infinity for the generalised Cauchy.
    # For the Cauchy, the exponential limit needs xi large and b = decay / xi
    # small, so the log-decay bound is widened accordingly.
    if (fam == "cauchy") { lo[3] <- -30; lo <- c(lo, 1e-3, -5); hi <- c(hi, 1, 16) }
    if (fam == "matern") { lo <- c(lo, 1e-3); hi <- c(hi, 0.5) }
    if (model$metric == "dkappa") {
      lo <- c(lo, -Inf); hi <- c(hi, Inf)              # log kappa, unbounded
    }
  } else if (fam == "expc") {
    lo <- c(lo, -10, log(1e-4)); hi <- c(hi, 6, log(50))
  } else if (fam == "wm1") {
    lo <- c(lo, log(1e-2)); hi <- c(hi, log(50))
  } else if (fam == "fem") {
    lo <- c(lo, log(1e-2), -6); hi <- c(hi, log(50), 6)
  }
  list(lower = lo, upper = hi)
}

# ---------------------------------------------------------------------------
# Profile likelihood with a mean shared across the T replicates.
# ---------------------------------------------------------------------------
profile_nll <- function(K, Y, X, details = FALSE) {
  U <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(U)) return(if (details) NULL else 1e12)
  Wx <- forwardsolve(t(U), X)
  Wy <- forwardsolve(t(U), Y)
  beta <- tryCatch(solve(crossprod(Wx), crossprod(Wx, rowMeans(Wy))),
                   error = function(e) NULL)
  if (is.null(beta)) return(if (details) NULL else 1e12)
  R <- Wy - as.vector(Wx %*% beta)
  val <- 0.5 * (length(Y) * log(2 * pi) +
                  ncol(Y) * 2 * sum(log(diag(U))) + sum(R^2))
  if (!is.finite(val)) return(if (details) NULL else 1e12)
  if (details) list(value = val, beta = as.vector(beta), U = U) else val
}

fit_model <- function(model, eng, Y, X, train, starts) {
  b <- bounds(model)
  Ytr <- Y[train, , drop = FALSE]; Xtr <- X[train, , drop = FALSE]
  obj <- function(theta) {
    C <- tryCatch(latent_cov(model, theta, eng), error = function(e) NULL)
    if (is.null(C) || any(!is.finite(C))) return(1e12)
    K <- C[train, train, drop = FALSE]
    diag(K) <- diag(K) + exp(2 * theta[1])
    profile_nll(K, Ytr, Xtr)
  }
  run <- function(s, maxit) tryCatch(
    optim(pmin(pmax(s, b$lower), b$upper), obj, method = "L-BFGS-B",
          lower = b$lower, upper = b$upper, control = list(maxit = maxit)),
    error = function(e) NULL)
  opts <- Filter(function(o) !is.null(o) && o$value < 1e11,
                 lapply(starts, run, maxit = 300))
  if (!length(opts)) return(NULL)
  best <- opts[[which.min(vapply(opts, `[[`, numeric(1), "value"))]]
  polish <- run(best$par, 300)                  # restart from the optimum
  if (!is.null(polish) && polish$value <= best$value) best <- polish
  list(par = best$par, nll = best$value, convergence = best$convergence)
}

# Exponential parameters (log sigma_e, log sigma^2, log decay[, log kappa])
# mapped into the Matern (nu = 1/2, which is exactly exponential) or the
# generalised Cauchy family (a = 1, xi = e^14, decay / xi), whose covariance
# then differs from the exponential by a relative amount of order 1e-6.
from_exp <- function(family, p) {
  k <- if (length(p) == 4L) p[4] else NULL
  switch(family,
         matern = c(p[1:3], 0.5, k),
         cauchy = c(p[1:2], p[3] - 14, 1, 14, k))
}

# Starting values. d_kappa models start from the fitted d_0 model of the same
# family, including kappa = KAPPA_START, so their training likelihood can never
# fall below that of d_0 (the models are nested as kappa -> 0). The Matern and
# generalised Cauchy models also start from the fitted exponential model in
# the same metric (exp_warm), so they cannot fall below it either.
make_starts <- function(model, eng, Y, X, train, warm = NULL, exp_warm = NULL) {
  Ytr <- Y[train, , drop = FALSE]; Xtr <- X[train, , drop = FALSE]
  res <- Ytr - as.vector(Xtr %*% qr.coef(qr(Xtr), rowMeans(Ytr)))
  v <- max(var(as.vector(res)), 1e-2)
  base <- c(log(sqrt(0.2 * v)), log(0.8 * v))
  fam <- model$family
  if (fam %in% names(n_radial)) {
    from_e <- if (!is.null(exp_warm) && fam != "exp")
      list(from_exp(fam, exp_warm)) else list()
    if (model$metric == "dkappa") {
      stopifnot(!is.null(warm))
      return(c(lapply(log(c(KAPPA_START, 0.02, 0.1, 0.5)),
                      function(lk) c(warm, lk)), from_e))
    }
    D <- if (model$metric == "dres") eng$Dres else eng$D0
    med <- median(D[train, train][upper.tri(D[train, train])])
    extra <- switch(fam, exp = numeric(0), cauchy = c(0.5, 0), matern = 0.25)
    return(c(lapply(c(1, 4) / med, function(dec) c(base, log(dec), extra)),
             from_e))
  }
  if (fam == "expc") {
    stopifnot(!is.null(warm))            # warm = fitted exp_d0 parameters
    return(lapply(log(c(1e-3, 0.1, 0.5)), function(lc) c(warm, lc)))
  }
  if (fam == "wm1") return(lapply(log(c(0.3, 2)), function(lk) c(base, lk)))
  if (fam == "fem") {
    stopifnot(!is.null(warm))            # warm = fitted WM1 parameters
    return(list(c(warm, qlogis(0.5 / NU_FEM_MAX))))     # alpha = 1
  }
}

# Exact Gaussian predictive distribution at the test sites (plug-in).
predict_model <- function(model, fit, eng, Y, X, train, test) {
  C <- latent_cov(model, fit$par, eng)
  K <- C[train, train, drop = FALSE]
  diag(K) <- diag(K) + exp(2 * fit$par[1])
  pd <- profile_nll(K, Y[train, , drop = FALSE], X[train, , drop = FALSE],
                    details = TRUE)
  U <- pd$U
  Ctt <- C[test, train, drop = FALSE]
  R <- Y[train, , drop = FALSE] - as.vector(X[train, , drop = FALSE] %*% pd$beta)
  mean <- as.vector(X[test, , drop = FALSE] %*% pd$beta) +
    Ctt %*% backsolve(U, forwardsolve(t(U), R))
  V <- C[test, test, drop = FALSE] - crossprod(forwardsolve(t(U), t(Ctt)))
  diag(V) <- diag(V) + exp(2 * fit$par[1])
  list(mean = mean, V = (V + t(V)) / 2, beta = pd$beta)
}

# Fit all nine models on one split and return predictions.
fit_all <- function(eng, Y, X, train, test, models = MODELS) {
  out <- list()
  fits <- list()
  fit_order <- c("exp_d0", "exp_d0c", "exp_dres", "wm1", "exp_dkappa", "cauchy_d0",
                 "matern_d0", "cauchy_dkappa", "matern_dkappa", "wm_fem")
  fit_order <- intersect(fit_order, models$name)
  for (nm in fit_order) {
    model <- as.list(models[models$name == nm, ])
    warm <- switch(nm,
                   exp_dkappa = , exp_d0c = fits$exp_d0$par,
                   cauchy_dkappa = fits$cauchy_d0$par,
                   matern_dkappa = fits$matern_d0$par,
                   wm_fem = fits$wm1$par, NULL)
    if (nm %in% c("exp_dkappa", "exp_d0c", "cauchy_dkappa", "matern_dkappa",
                  "wm_fem") &&
        is.null(warm)) { out[[nm]] <- NULL; next }
    exp_warm <- switch(nm, cauchy_d0 = , matern_d0 = fits$exp_d0$par,
                       cauchy_dkappa = , matern_dkappa = fits$exp_dkappa$par,
                       NULL)
    t0 <- proc.time()[["elapsed"]]
    fit <- fit_model(model, eng, Y, X, train,
                     make_starts(model, eng, Y, X, train, warm, exp_warm))
    if (is.null(fit)) next
    fits[[nm]] <- fit
    pred <- predict_model(model, fit, eng, Y, X, train, test)
    out[[nm]] <- list(fit = fit, pred = pred,
                      seconds = proc.time()[["elapsed"]] - t0)
  }
  out
}

# Interpretable summaries of fitted parameters.
model_kappa <- function(model, par) {
  fam <- model$family
  if (fam %in% names(n_radial) && model$metric == "dkappa")
    return(exp(par[3 + n_radial[[fam]]]))
  if (fam %in% c("wm1", "fem")) return(exp(par[3]))
  0
}
model_alpha <- function(model, par) {
  if (model$family == "fem") return(NU_FEM_MAX * plogis(par[4]) + 0.5)
  if (model$family == "wm1") return(1)
  NA_real_
}
