#!/usr/bin/env Rscript
# Summaries of the cross-validation runs. Run from the repository root after
# run_cv.R. Writes CSV summaries to output/ and LaTeX rows to output/tex/.
source("R/models.R")             # MODELS (names and labels)

read_scores <- function(dataset, design) {
  f <- file.path("output", paste(dataset, design, sep = "_"), "scores.csv")
  if (!file.exists(f)) return(NULL)
  read.csv(f)
}

SCORES <- c("rmse", "crps", "log_score", "joint_log_score", "cover95")
PAIRS <- c(exp_dkappa = "exp_d0", cauchy_dkappa = "cauchy_d0",
           matern_dkappa = "matern_d0")

# Leave-percentage-out: equal-weight means over the random draws, paired
# differences against the exponential d_0 model, and paired win counts of
# each d_kappa model against the d_0 model of the same family.
summarise_lpo <- function(x) {
  x <- x[order(x$pct, x$split, x$model), ]
  complete <- names(which(table(x$split) == nrow(MODELS)))
  x <- x[x$split %in% as.integer(complete), ]
  out <- list()
  for (p in sort(unique(x$pct))) {
    xp <- x[x$pct == p, ]
    ref <- xp[xp$model == "exp_d0", ]
    for (m in MODELS$name) {
      xm <- xp[xp$model == m, ]
      stopifnot(identical(xm$split, ref$split))
      row <- data.frame(pct = p, model = m, draws = nrow(xm),
                        n_train = mean(xm$n_train),
                        nearest_median_km = mean(xm$nearest_median_km),
                        kappa_median = median(xm$kappa),
                        alpha_median = median(xm$alpha))
      for (s in SCORES) {
        row[[s]] <- mean(xm[[s]])
        row[[paste0(s, "_median")]] <- median(xm[[s]])
        d <- xm[[s]] - ref[[s]]
        row[[paste0(s, "_diff")]] <- mean(d)
        row[[paste0(s, "_diff_se")]] <- sd(d) / sqrt(length(d))
      }
      if (m %in% names(PAIRS)) {
        base <- xp[xp$model == PAIRS[[m]], ]
        for (s in c("rmse", "crps", "log_score", "joint_log_score"))
          row[[paste0(s, "_wins_vs_d0")]] <- mean(xm[[s]] < base[[s]])
      }
      out[[length(out) + 1L]] <- row
    }
  }
  do.call(rbind, lapply(out, function(r) {
    miss <- setdiff(unique(unlist(lapply(out, names))), names(r))
    r[miss] <- NA; r }))
}

# Leave-region-out: scores pooled over the K regions of each buffer, weighting
# each region by its number of held-out values (RMSE pooled on the squared
# scale), so every held-out site counts once per buffer.
summarise_lro <- function(x) {
  out <- list()
  for (h in sort(unique(x$buffer_km))) {
    xh <- x[x$buffer_km == h, ]
    regions <- names(which(table(xh$region) == nrow(MODELS)))
    xh <- xh[xh$region %in% as.integer(regions), ]
    ref <- xh[xh$model == "exp_d0", ]
    ref <- ref[order(ref$region), ]
    for (m in MODELS$name) {
      xm <- xh[xh$model == m, ]
      xm <- xm[order(xm$region), ]
      w <- xm$n_test / sum(xm$n_test)
      row <- data.frame(buffer_km = h, model = m, regions = nrow(xm),
                        n_test = sum(xm$n_test), n_train = mean(xm$n_train),
                        nearest_median_km = sum(w * xm$nearest_median_km),
                        kappa_median = median(xm$kappa),
                        rmse = sqrt(sum(w * xm$rmse^2)),
                        crps = sum(w * xm$crps),
                        log_score = sum(w * xm$log_score),
                        joint_log_score = sum(w * xm$joint_log_score),
                        cover95 = sum(w * xm$cover95),
                        regions_crps_better_than_exp_d0 =
                          sum(xm$crps < ref$crps))
      if (m %in% names(PAIRS)) {
        base <- xh[xh$model == PAIRS[[m]], ]
        base <- base[order(base$region), ]
        row$regions_crps_better_than_d0 <- sum(xm$crps < base$crps)
        row$regions_rmse_better_than_d0 <- sum(xm$rmse < base$rmse)
      } else {
        row$regions_crps_better_than_d0 <- NA
        row$regions_rmse_better_than_d0 <- NA
      }
      out[[length(out) + 1L]] <- row
    }
  }
  do.call(rbind, out)
}

fmt <- function(v, d) formatC(v, format = "f", digits = d)
bold_min <- function(v, d) {
  s <- fmt(v, d)
  best <- which(fmt(v, d) == fmt(min(v), d))
  s[best] <- paste0("$\\mathbf{", s[best], "}$")
  s[-best] <- paste0("$", s[-best], "$")
  s
}

# LaTeX body rows: one row per model, columns given as list(table, score, digits).
tex_rows <- function(cols) {
  cells <- lapply(cols, function(cl) {
    v <- cl$tab[[cl$score]][match(MODELS$name, cl$tab$model)]
    bold_min(v, cl$digits)
  })
  rows <- vapply(seq_len(nrow(MODELS)), function(i)
    paste(c(MODELS$label[i], vapply(cells, `[`, character(1), i)),
          collapse = " & "), character(1))
  rows <- paste0(rows, "\\\\")
  # Separate the d_0, d_kappa and reference blocks as in Table 1.
  c(rows[1:3], "\\hline", rows[4:6], "\\hline", rows[7:9])
}

dir.create(file.path("output", "tex"), recursive = TRUE, showWarnings = FALSE)
lpo <- lro <- list()
DATASETS <- if (length(commandArgs(TRUE))) commandArgs(TRUE) else c("paris_k", "london")
for (d in DATASETS) {
  x <- read_scores(d, "lpo")
  if (!is.null(x)) {
    lpo[[d]] <- summarise_lpo(x)
    write.csv(lpo[[d]], file.path("output", paste0(d, "_lpo_summary.csv")),
              row.names = FALSE)
  }
  x <- read_scores(d, "lro")
  if (!is.null(x)) {
    lro[[d]] <- summarise_lro(x)
    write.csv(lro[[d]], file.path("output", paste0(d, "_lro_summary.csv")),
              row.names = FALSE)
  }
}

DIG <- function(d) 3
for (d in names(lpo)) {
  t75 <- lpo[[d]][lpo[[d]]$pct == 75, ]
  t90 <- lpo[[d]][lpo[[d]]$pct == 90, ]
  writeLines(tex_rows(list(
    list(tab = t75, score = "rmse", digits = DIG(d)),
    list(tab = t75, score = "crps", digits = DIG(d)),
    list(tab = t90, score = "rmse", digits = DIG(d)),
    list(tab = t90, score = "crps", digits = DIG(d)))),
    file.path("output", "tex", paste0(d, "_lpo_rows.tex")))
}
for (d in names(lro)) {
  hs <- sort(unique(lro[[d]]$buffer_km))
  cols <- lapply(hs, function(h)
    list(tab = lro[[d]][lro[[d]]$buffer_km == h, ], score = "crps", digits = DIG(d)))
  writeLines(tex_rows(cols),
             file.path("output", "tex", paste0(d, "_lro_crps_rows.tex")))
  cols <- lapply(hs, function(h)
    list(tab = lro[[d]][lro[[d]]$buffer_km == h, ], score = "rmse", digits = DIG(d)))
  writeLines(tex_rows(cols),
             file.path("output", "tex", paste0(d, "_lro_rmse_rows.tex")))
}

options(width = 200)
for (d in names(lpo)) {
  cat("\n=====", d, "leave-percentage-out =====\n")
  print(lpo[[d]][, c("pct", "model", "draws", "n_train", "kappa_median",
                     "rmse", "crps", "log_score", "joint_log_score", "cover95",
                     "rmse_diff", "rmse_diff_se", "crps_diff", "crps_diff_se",
                     "crps_wins_vs_d0", "joint_log_score_wins_vs_d0")],
        digits = 4, row.names = FALSE)
}
for (d in names(lro)) {
  cat("\n=====", d, "leave-region-out =====\n")
  print(lro[[d]][, c("buffer_km", "model", "regions", "n_train",
                     "nearest_median_km", "kappa_median", "rmse", "crps",
                     "log_score", "joint_log_score", "cover95",
                     "regions_crps_better_than_exp_d0",
                     "regions_crps_better_than_d0")],
        digits = 4, row.names = FALSE)
}

# Table 1 rows: number of parameters and maximised negative log-likelihood
# from fit_fulldata.R (Paris, London).
fits <- lapply(c("paris_k", "london"), function(d) {
  f <- file.path("output", paste0(d, "_fulldata_fit.csv"))
  if (file.exists(f)) read.csv(f) else NULL
})
if (!any(vapply(fits, is.null, logical(1)))) {
  v <- lapply(fits, function(f) bold_min(f$neg_loglik[match(MODELS$name, f$model)], 3))
  p <- fits[[1]]$p[match(MODELS$name, fits[[1]]$model)]
  rows <- paste0(MODELS$label, " & ", p, " & ", v[[1]], " & ", v[[2]], "\\\\")
  writeLines(c(rows[1:3], "\\hline", rows[4:6], "\\hline", rows[7:9]),
             file.path("output", "tex", "fulldata_rows.tex"))
}
