# Cross-validation designs and predictive scores.

# Leave-percentage-out: in each repeat a fixed percentage of the sites is
# drawn at random and predicted jointly from the remaining sites.
splits_lpo <- function(n, pcts = c(75, 90), repeats = 50L, seed = 1L) {
  set.seed(seed)
  out <- list()
  for (pct in pcts) for (r in seq_len(repeats)) {
    test <- sort(sample.int(n, round(pct / 100 * n)))
    out[[length(out) + 1L]] <- list(design = "lpo", pct = pct, repeat_id = r,
                                    region = NA_integer_, buffer_km = NA_real_,
                                    test = test, train = setdiff(seq_len(n), test))
  }
  out
}

# Buffered leave-region-out: the sites are partitioned into K contiguous
# regions by k-means on their coordinates (locations only, no responses).
# Each region is predicted jointly from the sites whose road-network distance
# to every site of the region is at least buffer_km, so the training data are
# removed from a neighbourhood of the prediction region.
splits_lro <- function(xy, geo, K = 8L, buffers = c(0, 2, 4, 6),
                       seed = 1L, min_train = 30L) {
  set.seed(seed)
  km <- kmeans(xy, centers = K, nstart = 100, iter.max = 100)
  # Number regions west to east for readability.
  relabel <- match(seq_len(K), order(km$centers[, 1]))
  region <- relabel[km$cluster]
  out <- list()
  for (h in buffers) for (k in seq_len(K)) {
    test <- which(region == k)
    others <- which(region != k)
    near <- apply(geo[others, test, drop = FALSE], 1, min)
    train <- others[near >= h]
    if (length(train) < min_train) {
      message(sprintf("  region %d, buffer %g km: %d training sites < %d, skipped",
                      k, h, length(train), min_train))
      next
    }
    out[[length(out) + 1L]] <- list(design = "lro", pct = NA_real_,
                                    repeat_id = NA_integer_, region = k,
                                    buffer_km = h, test = test, train = train)
  }
  attr(out, "region") <- region
  out
}

normal_crps <- function(y, m, s) {
  z <- (y - m) / s
  s * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))
}

# Scores on the original response scale; lower is better except coverage.
#   rmse, crps         pointwise accuracy of the predictive mean/distribution
#   log_score          mean marginal negative log predictive density
#   joint_log_score    negative log joint predictive density of all held-out
#                      values of a replicate, divided by their number
#   cover95            empirical coverage of 95% predictive intervals
score_split <- function(pred, Ytest, y_center, y_scale) {
  M <- pred$mean * y_scale + y_center
  V <- pred$V * y_scale^2
  E <- Ytest - M
  sd <- sqrt(pmax(diag(V), 1e-12))
  Z <- sweep(E, 1, sd, "/")
  U <- chol(V)
  W <- forwardsolve(t(U), E)
  joint <- 0.5 * (length(E) * log(2 * pi) +
                    ncol(E) * 2 * sum(log(diag(U))) + sum(W^2)) / length(E)
  data.frame(rmse = sqrt(mean(E^2)),
             crps = mean(normal_crps(Ytest, M, sd)),
             log_score = mean(0.5 * log(2 * pi) + log(sd) + 0.5 * Z^2),
             joint_log_score = joint,
             cover95 = mean(abs(Z) <= qnorm(0.975)))
}
