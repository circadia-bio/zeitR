# dev/stage2_diag.R
# Isolates Stage 2 (.refine_borders) and diffs it against python_stage2_periods.csv.
# Run from package root:  source("dev/stage2_diag.R")
# Paste the printed output back.

devtools::load_all(quiet = TRUE)

# ── Reproduce the exact inputs detect_offwrist_bimodal feeds the refiner ──────
hws <- 10L
rec  <- read_acttrust(system.file("extdata", "input1.txt", package = "zeitR"), tz = "UTC")
x    <- prepare_actigraphy(rec)

valid_temp <- x$int_temp > 0
activity <- as.double(x$activity[valid_temp])
int_temp <- as.double(x$int_temp[valid_temp])
ext_temp <- as.double(x$ext_temp[valid_temp])

act_median   <- zeitR:::median_filter(activity, hws)
norm_act_med <- zeitR:::norm_01(act_median)
zp           <- zeitR:::zero_prop(norm_act_med)
low_q        <- zp + 0.15 * (1 - zp)
low_act_thr  <- max(stats::quantile(norm_act_med, low_q, names = FALSE, type = 1), 0.015)
is_low_act   <- norm_act_med < low_act_thr

norm_temp     <- zeitR:::norm_01(int_temp)
low_act_temps <- norm_temp[is_low_act]
gmm <- zeitR:::.fit_gmm_threshold(low_act_temps, nbins = 100L, min_temp_threshold = 0.35)
temp_threshold <- gmm$threshold; ashman <- gmm$ashman_d
temp_min <- min(int_temp); temp_max <- max(int_temp)
temp_threshold_orig <- temp_min + temp_threshold * (temp_max - temp_min)
med_high <- stats::median(int_temp[!is_low_act], na.rm = TRUE)
if (!is.na(med_high) && temp_threshold_orig > med_high) temp_threshold_orig <- med_high
is_low_temp <- int_temp < temp_threshold_orig

offwrist_raw <- 1L - as.integer(is_low_act & is_low_temp)

norm_tv  <- zeitR:::norm_01(zeitR:::var_filter(int_temp, hws))
n        <- length(offwrist_raw)

# ── Replicate the constants + Stage 1 from .bimodal_refine_acttrust ───────────
minimum_onwrist_length          <- 20L
minimum_low_activity_proportion <- 0.8
max_low_tv_border               <- 0.3
minimum_preceding_onwrist_length<- 40L
act_zero_prop <- zeitR:::zero_prop(activity)
act_thr_q     <- act_zero_prop + (1 - act_zero_prop) * 0.24
activity_thr  <- as.double(stats::quantile(activity, act_thr_q, names = FALSE, type = 1))
if (activity_thr < 200) activity_thr <- 200
temp_var_thr  <- as.double(stats::quantile(norm_tv, 0.9, names = FALSE))

on_periods <- zeitR:::.rle_periods(offwrist_raw == 1L)
on_periods <- on_periods[(on_periods$end - on_periods$start) > minimum_onwrist_length, ]
on_periods <- on_periods[order(on_periods$start), ]; row.names(on_periods) <- NULL
stage1 <- integer(n)
for (i in seq_len(nrow(on_periods))) {
  s <- on_periods$start[i]; e <- on_periods$end[i]; if (s < e) stage1[(s + 1L):e] <- 1L
}
ow_s1 <- zeitR:::.rle_periods(stage1 == 0L)

# ── Run Stage 2 only ──────────────────────────────────────────────────────────
half_filter_hws <- as.integer(hws / 2)
r_s2 <- zeitR:::.refine_borders(
  ow_s1, norm_tv, act_median, as.integer(is_low_act),
  int_temp, temp_threshold_orig, activity_thr, minimum_low_activity_proportion,
  minimum_preceding_onwrist_length, hws, half_filter_hws,
  max_low_tv_border, temp_var_thr, activity, n
)

# ── Load Python Stage 2 ground truth ──────────────────────────────────────────
py <- read.csv(system.file("extdata", "python_stage2_periods.csv", package = "zeitR"))
# normalise column names
names(py) <- tolower(names(py))
scol <- grep("start", names(py), value = TRUE)[1]
ecol <- grep("end",   names(py), value = TRUE)[1]
py <- data.frame(start = py[[scol]], end = py[[ecol]])
py <- py[order(py$start), ]; row.names(py) <- NULL

cat(sprintf("R stage2 periods: %d   |   Python stage2 periods: %d\n",
            nrow(r_s2), nrow(py)))

# ── Align by nearest start and report disagreements ───────────────────────────
match_to_py <- function(rs) which.min(abs(py$start - rs))
diffs <- data.frame(r_start = r_s2$start, r_end = r_s2$end,
                    py_start = NA_integer_, py_end = NA_integer_,
                    ds = NA_integer_, de = NA_integer_)
for (i in seq_len(nrow(r_s2))) {
  j <- match_to_py(r_s2$start[i])
  diffs$py_start[i] <- py$start[j]; diffs$py_end[i] <- py$end[j]
  diffs$ds[i] <- r_s2$start[i] - py$start[j]
  diffs$de[i] <- r_s2$end[i]   - py$end[j]
}

bad <- diffs[abs(diffs$ds) > 1 | abs(diffs$de) > 1, ]
cat(sprintf("\nPeriods with |start diff|>1 or |end diff|>1: %d of %d\n",
            nrow(bad), nrow(diffs)))
cat("start-diff distribution:\n"); print(table(diffs$ds))
cat("end-diff distribution:\n");   print(table(diffs$de))

cat("\n-- worst 25 by |start diff|+|end diff| --\n")
ord <- order(-(abs(bad$ds) + abs(bad$de)))
print(utils::head(bad[ord, ], 25), row.names = FALSE)

# ── Focused trace of the ~60255 region ────────────────────────────────────────
cat("\n-- R stage2 periods overlapping epochs 60000..60800 --\n")
print(r_s2[r_s2$start < 60800 & r_s2$end > 60000, ], row.names = FALSE)
cat("-- Python stage2 periods overlapping epochs 60000..60800 --\n")
print(py[py$start < 60800 & py$end > 60000, ], row.names = FALSE)

cat("\n-- ow_s1 (stage1) periods overlapping 60000..60800 --\n")
print(ow_s1[ow_s1$start < 60800 & ow_s1$end > 60000, ], row.names = FALSE)

cat("\n-- ntv around 60250..60360 (R index = epoch+1) --\n")
rng <- 60250:60360
cat("temp_var_thr =", round(temp_var_thr, 5),
    " activity_thr =", round(activity_thr, 2),
    " temp_threshold =", round(temp_threshold_orig, 4), "\n")
print(data.frame(
  epoch = rng,
  ntv   = round(norm_tv[rng + 1L], 4),
  ge_thr= norm_tv[rng + 1L] >= temp_var_thr,
  is_low_act = as.integer(is_low_act)[rng + 1L],
  low_temp   = (int_temp[rng + 1L] < temp_threshold_orig)
), row.names = FALSE)
