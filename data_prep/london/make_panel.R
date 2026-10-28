# Hourly speed maps 17:00-18:00 (mean of the four 15-minute averages) for
# every day of 2025-01-01..06-30, from the WebTRIS downloads in raw/.
# A site-hour needs at least 3 of the 4 intervals. Sites with a value on at
# least 95% of the days are kept, and then the days on which all of them have
# a value, so the panel is complete.
# Writes panel_1700.csv (rows = days, first column = date, then one column per
# site id).
# Usage: Rscript make_panel.R
sites <- read.csv("sites_selected.csv", colClasses = c(site_id = "character"))
days <- format(seq(as.Date("2025-01-01"), as.Date("2025-06-30"), by = "day"))
Y <- matrix(NA_real_, length(days), nrow(sites),
            dimnames = list(days, sites$site_id))
for (j in seq_len(nrow(sites))) {
  f <- file.path("raw", paste0("site_", sites$site_id[j], ".csv"))
  if (!file.exists(f) || file.size(f) < 100) next
  x <- read.csv(f)
  x <- x[x$time_end %in% c("17:14:00", "17:29:00", "17:44:00", "17:59:00") &
           is.finite(x$avg_mph) & x$avg_mph > 0, ]
  if (!nrow(x)) next
  n <- tapply(x$avg_mph, x$date, length)
  m <- tapply(x$avg_mph, x$date, mean)
  m <- m[n >= 3 & names(m) %in% days]
  Y[names(m), j] <- m
}
site_cov <- colMeans(is.finite(Y))
cat("sites with any data:", sum(site_cov > 0), "of", ncol(Y), "\n")
cat("median day coverage over those sites:",
    round(median(rowMeans(is.finite(Y[, site_cov > 0]))), 3), "\n")
keep_site <- site_cov >= 0.95
keep_day <- rowSums(!is.finite(Y[, keep_site, drop = FALSE])) == 0
P <- Y[keep_day, keep_site, drop = FALSE]
cat("panel:", nrow(P), "days x", ncol(P), "sites\n")
write.csv(data.frame(date = rownames(P), P, check.names = FALSE),
          "panel_1700.csv", row.names = FALSE)
