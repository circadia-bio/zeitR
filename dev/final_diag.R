# dev/final_diag.R
# Locates the contiguous false-positive / false-negative runs in the FINAL
# off-wrist output and attributes each to Stage 2 vs a later stage.
# Run:  source("dev/final_diag.R")   then paste the output back.

devtools::load_all(quiet = TRUE)

rec  <- read_acttrust(system.file("extdata", "input1.txt", package = "zeitR"), tz = "UTC")
prep <- detect_offwrist_bimodal(prepare_actigraphy(rec))

r_ow  <- as.integer(prep$state == 4)
py    <- read.csv(system.file("extdata", "python_output.csv", package = "zeitR"))
py_ow <- as.integer(py$state == 4)

n <- length(r_ow)
cat(sprintf("Correlation: %.6f   n=%d\n", cor(r_ow, py_ow), n))
print(table(R = r_ow, Python = py_ow))

# Build python_output index space. python_output.csv is the FULL recording
# (incl. temp<=0 epochs forced off-wrist). r_ow is also full length.
stopifnot(length(r_ow) == length(py_ow))

# ── helper: contiguous runs of a logical vector, as 0-based [start,end) ───────
runs <- function(flag) {
  r <- rle(flag); ends <- cumsum(r$lengths); strt <- ends - r$lengths
  keep <- r$values
  data.frame(start = strt[keep], end = ends[keep], len = r$lengths[keep])
}

fp <- runs(r_ow == 1 & py_ow == 0)   # R off, Python on  (over-detection)
fn <- runs(r_ow == 0 & py_ow == 1)   # R on,  Python off (under-detection)

cat(sprintf("\nFP runs: %d  (total %d epochs)\n", nrow(fp), sum(fp$len)))
print(fp[order(-fp$len), ][seq_len(min(15, nrow(fp))), ], row.names = FALSE)
cat(sprintf("\nFN runs: %d  (total %d epochs)\n", nrow(fn), sum(fn$len)))
print(fn[order(-fn$len), ][seq_len(min(15, nrow(fn))), ], row.names = FALSE)

# ── Attribute the biggest FP runs: are they off-wrist in Python Stage 2? ──────
pys2 <- read.csv(system.file("extdata", "python_stage2_periods.csv", package = "zeitR"))
# valid-temp index mapping: stage2 is in valid-temp space; final is full space.
valid_temp <- prep$int_temp > 0
# epoch (full) -> valid index: cumulative count of valid up to that epoch
valid_cumidx <- cumsum(valid_temp)            # full-space epoch i (1-based) -> # valid <= i
in_py_stage2 <- function(full_epoch0) {
  # full_epoch0 is 0-based full-space; convert to 0-based valid-space
  fi <- full_epoch0 + 1L
  if (fi < 1L || fi > n || !valid_temp[fi]) return(NA)
  v0 <- valid_cumidx[fi] - 1L                 # 0-based valid index
  any(pys2$start <= v0 & pys2$end > v0)
}

cat("\n-- biggest FP runs: midpoint membership in Python Stage 2 (valid-space) --\n")
big <- fp[order(-fp$len), ][seq_len(min(8, nrow(fp))), ]
for (k in seq_len(nrow(big))) {
  mid0 <- big$start[k] + big$len[k] %/% 2L
  cat(sprintf("FP [%d,%d) len=%d  mid-epoch=%d  in_python_stage2=%s\n",
              big$start[k], big$end[k], big$len[k], mid0, in_py_stage2(mid0)))
}
