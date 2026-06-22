# dev/region_diag.R
# Epoch-by-epoch divergence of R final refiner vs Python refined_offwrist
# (valid space) around the 60254-60484 block, plus the final period lists.
# Run:  source("dev/region_diag.R")  then paste output.

devtools::load_all(quiet = TRUE)
options(zeitR.debug_offwrist = TRUE)

rec  <- read_acttrust(system.file("extdata", "input1.txt", package = "zeitR"), tz = "UTC")
prep <- detect_offwrist_bimodal(prepare_actigraphy(rec))
dbg  <- get(".zeitR_dbg", envir = .GlobalEnv)

valid_temp <- prep$int_temp > 0
r_off  <- as.integer(prep$state[valid_temp] == 4)            # 1 = offwrist (valid space)
pin    <- read.csv(system.file("extdata", "python_intermediate.csv", package = "zeitR"))
py_off <- as.integer(round(pin$refined_offwrist) == 0)        # 1 = offwrist (valid space)

runs <- function(flag) {
  r <- rle(flag); e <- cumsum(r$lengths); s <- e - r$lengths
  k <- r$values; data.frame(start = s[k], end = e[k], len = r$lengths[k])
}

cat("== FINAL off-wrist periods in valid epochs [59400, 60700] ==\n")
rr <- runs(r_off == 1L);  rr <- rr[rr$start < 60700 & rr$end > 59400, ]
pp <- runs(py_off == 1L); pp <- pp[pp$start < 60700 & pp$end > 59400, ]
cat("R final periods:\n");      print(rr, row.names = FALSE)
cat("Python final periods:\n"); print(pp, row.names = FALSE)

cat("\n== after_initial_periods (R, post feature-filters) in [59400,60700] ==\n")
aip <- dbg$after_initial_periods
print(aip[aip$start < 60700 & aip$end > 59400, c("start","end","length","sleep_prop")],
      row.names = FALSE)

cat("\n== forbidden_zone coverage in [60240,60620] (R) ==\n")
fz <- dbg$forbidden_zone
fzr <- runs(fz == 1L); print(fzr[fzr$start < 60620 & fzr$end > 60240, ], row.names = FALSE)

cat("\n== ALL 9 FP locations: are they whole-period or border? (valid space) ==\n")
dis <- runs(r_off != py_off)
dis$r_says_off <- vapply(seq_len(nrow(dis)), function(i)
  as.integer(any(r_off[(dis$start[i]+1L):dis$end[i]] == 1L)), integer(1))
print(dis[dis$len >= 1, ][order(-dis$len), ][seq_len(min(12, nrow(dis))), ], row.names = FALSE)

# For the biggest disagreement block, show which R/PY period each border belongs to
big <- dis[order(-dis$len), ][1, ]
cat(sprintf("\n== biggest divergence block valid [%d,%d) len=%d ==\n",
            big$start, big$end, big$len))
cat("R period(s) covering this block:\n")
print(rr[rr$start < big$end & rr$end > big$start, ], row.names = FALSE)
cat("Python period(s) covering this block:\n")
print(pp[pp$start < big$end & pp$end > big$start, ], row.names = FALSE)

options(zeitR.debug_offwrist = FALSE)
