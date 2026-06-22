testthat::test_that(".identify_peaks_and_valleys returns correct region structure", {
  ipv <- getFromNamespace(".identify_peaks_and_valleys", "zeitR")

  # signal: peak(2), valley(3), peak(1)  -> 0-indexed exclusive-end regions
  df <- ipv(c(1, 1, 0, 0, 0, 1))

  expect_identical(df$class,  c("p", "v", "p"))
  expect_identical(df$start,  c(0L, 2L, 5L))
  expect_identical(df$end,    c(2L, 5L, 6L))
  expect_identical(df$length, c(2L, 3L, 1L))
  # features are not computed unless `activity` is supplied
  expect_false("mean" %in% names(df))
})

testthat::test_that(".identify_peaks_and_valleys computes region features when activity supplied", {
  ipv <- getFromNamespace(".identify_peaks_and_valleys", "zeitR")

  activity <- c(10, 20, 0, 0, 0, 30)
  signal   <- c(1,  1,  0, 0, 0, 1)
  df <- ipv(signal, activity = activity, threshold = 0.5)

  expect_identical(df$length,                     c(2L, 3L, 1L))
  expect_identical(df$mean,                       c(15, 0, 30))
  expect_identical(df$median,                     c(15, 0, 30))
  expect_identical(df$zero_proportion,            c(0, 1, 0))
  expect_identical(df$above_threshold_proportion, c(1, 0, 1))
})

testthat::test_that(".peak_valley_length_filter merges a short mid valley+peak (golden case c1)", {
  pvlf <- getFromNamespace(".peak_valley_length_filter", "zeitR")

  # p3, v2, p1, v4, p3 ; min_length = 2 -> the v2 and following p1 are absorbed
  input    <- c(1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1)
  expected <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1)

  expect_identical(pvlf(input, 2L), expected)
})

testthat::test_that(".peak_valley_length_filter special-cases a short first region (golden case c2)", {
  pvlf <- getFromNamespace(".peak_valley_length_filter", "zeitR")

  # v2 (first), p4, v3, p2 ; min_length = 2 -> first valley merges forward,
  # trailing short peak merges back into the valley
  input    <- c(0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1)
  expected <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0)

  expect_identical(pvlf(input, 2L), expected)
})

testthat::test_that(".peak_valley_length_filter merges a short trailing region (golden case c3)", {
  pvlf <- getFromNamespace(".peak_valley_length_filter", "zeitR")

  # p4, v5, p4, v1 (last) ; min_length = 2 -> trailing v1 absorbed into the peak
  input    <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0)
  expected <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)

  expect_identical(pvlf(input, 2L), expected)
})

testthat::test_that(".peak_valley_length_filter handles single-region (degenerate) inputs", {
  pvlf <- getFromNamespace(".peak_valley_length_filter", "zeitR")

  # all sleep -> stays all sleep; all wake -> stays all wake (no neighbour to merge)
  expect_identical(pvlf(c(0, 0, 0, 0), 2L), c(0L, 0L, 0L, 0L))
  expect_identical(pvlf(c(1, 1, 1),    2L), c(1L, 1L, 1L))
})

testthat::test_that(".peak_valley_length_filter preserves length and 0/1 domain", {
  pvlf <- getFromNamespace(".peak_valley_length_filter", "zeitR")

  set.seed(1)
  x <- as.integer(stats::runif(500) > 0.5)
  out <- pvlf(x, 11L)

  expect_length(out, 500L)
  expect_true(all(out %in% c(0L, 1L)))
})

# ── Stage 2: sleep-gap separation ─────────────────────────────────────────────

testthat::test_that(".datetime_diff returns 0-prefixed per-epoch second gaps", {
  dtd <- getFromNamespace(".datetime_diff", "zeitR")

  stamps <- as.POSIXct(c("2024-01-01 00:00:00",
                         "2024-01-01 00:01:00",
                         "2024-01-01 00:02:30"), tz = "UTC")
  expect_identical(dtd(stamps), c(0, 60, 90))
})

testthat::test_that(".np_slice_assign matches NumPy slice semantics", {
  nsa <- getFromNamespace(".np_slice_assign", "zeitR")

  # ordinary interior slice: x[2:5] <- 1 (0-indexed, exclusive end)
  expect_identical(nsa(rep(0L, 10), 2L, 5L, 1L),
                   c(0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L))

  # negative start wraps (NumPy): x[-6:15] on length 20 -> only index 14 set
  expect_identical(nsa(rep(0L, 20), -6L, 15L, 1L),
                   c(rep(0L, 14), 1L, rep(0L, 5)))

  # fully out-of-range -> no-op
  expect_identical(nsa(rep(0L, 5), 10L, 20L, 1L), rep(0L, 5))
})

testthat::test_that(".sleep_gap_separation splits a sleep block at an interior gap", {
  sgs <- getFromNamespace(".sleep_gap_separation", "zeitR")

  # 30-epoch sleep block; a >1h gap lands at epoch 18 (1-indexed) inside it
  detection    <- c(1L, 1L, rep(0L, 30), 1L, 1L)
  datetime_diff <- rep(60, 34); datetime_diff[18] <- 4000
  expected <- c(1L, 1L, 0L, 0L, 0L, 0L, 0L,
                rep(1L, 21),
                0L, 0L, 0L, 0L, 1L, 1L)

  expect_identical(sgs(detection, datetime_diff, 3600), expected)
})

testthat::test_that(".sleep_gap_separation reproduces NumPy boundary-wrap behaviour", {
  sgs <- getFromNamespace(".sleep_gap_separation", "zeitR")

  # gap at epoch 5 (1-indexed) -> window start wraps negative -> only 1 epoch set
  detection     <- c(rep(0L, 15), rep(1L, 5))
  datetime_diff <- rep(60, 20); datetime_diff[5] <- 4000
  expected      <- c(rep(0L, 14), 1L, rep(1L, 5))

  expect_identical(sgs(detection, datetime_diff, 3600), expected)
})

testthat::test_that(".sleep_gap_separation is a no-op without large gaps", {
  sgs <- getFromNamespace(".sleep_gap_separation", "zeitR")

  detection     <- c(1L, 0L, 0L, 0L, 1L)
  datetime_diff <- rep(60, 5)
  expect_identical(sgs(detection, datetime_diff, 3600), detection)
})

# ── Stage 3: shared refiner helpers ───────────────────────────────────────────

testthat::test_that(".below_prop computes proportion <= threshold", {
  bp <- getFromNamespace(".below_prop", "zeitR")
  expect_identical(bp(c(1, 2, 3, 4, 5), 3), 0.6)
  expect_identical(bp(numeric(0), 3), 0)
})

testthat::test_that(".get_peak returns the local extremum value (0-indexed centre)", {
  gp <- getFromNamespace(".get_peak", "zeitR")
  sig <- c(5, 1, 9, 2, 7)
  expect_identical(gp(sig, 2L, FALSE, 1L), 9)   # window c(1, 9, 2) -> max
  expect_identical(gp(sig, 2L, TRUE,  1L), 1)   # window c(1, 9, 2) -> min
  expect_identical(gp(sig, 0L, FALSE, 1L), 5)   # clamped window c(5, 1) -> max
})

testthat::test_that(".cspd_median_filter matches functions.py median_filter", {
  mf <- getFromNamespace(".cspd_median_filter", "zeitR")
  # padding = NULL, hws = 1, monotonic signal -> unchanged
  expect_identical(mf(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 1L), as.numeric(1:10))
  # padding = "padded", hws = 1 -> trims 4*hws, length 6
  expect_identical(mf(c(0, 0, 1, 2, 3, 4, 5, 6, 0, 0), 1L, padding = "padded"),
                   c(1, 2, 3, 4, 5, 5))
})

testthat::test_that(".remove_peak_valley merges first region and recomputes features", {
  rpv <- getFromNamespace(".remove_peak_valley", "zeitR")

  pv <- data.frame(
    class                      = c("v", "p", "v", "p"),
    start                      = c(0L, 2L, 5L, 8L),
    end                        = c(2L, 5L, 8L, 10L),
    length                     = c(2L, 3L, 3L, 2L),
    mean                       = c(0, 9, 0, 9),
    median                     = c(0, 9, 0, 9),
    zero_proportion            = c(1, 0, 1, 0),
    above_threshold_proportion = c(0, 1, 0, 1),
    stringsAsFactors = FALSE
  )
  activity <- c(0, 0, 9, 9, 9, 0, 0, 0, 9, 9)

  out <- rpv(pv, 0L, activity, 3)

  expect_identical(out$class,  c("p", "v", "p"))
  expect_identical(out$start,  c(0L, 5L, 8L))
  expect_identical(out$end,    c(5L, 8L, 10L))
  expect_identical(out$length, c(5L, 3L, 2L))
  expect_equal(out$mean,                       c(5.4, 0, 9))
  expect_equal(out$median,                     c(9,   0, 9))
  expect_equal(out$zero_proportion,            c(0.4, 1, 0))
  expect_equal(out$above_threshold_proportion, c(0.6, 1, 0))
})

# ── Stage 3: bed-time refiner — Unit A (window finding) ────────────────────────

testthat::test_that(".bt_datetime_gap_check validates backward/forward gaps", {
  gc <- getFromNamespace(".bt_datetime_gap_check", "zeitR")
  self <- new.env(parent = emptyenv())
  self$secs <- c(0, 60, 120, 4000)      # big jump into the last epoch
  self$maximum_allowed_gap <- 600

  expect_true(gc(self, 1L))                         # backward: 60s gap, valid
  expect_false(gc(self, 3L))                        # backward: 3880s gap, invalid
  expect_true(gc(self, 0L, direction = "forward"))  # forward from epoch 0: 60s
  expect_false(gc(self, 2L, direction = "forward")) # forward from epoch 2: 3880s
  expect_identical(gc(self, 3L, return_gap = TRUE), list(valid = FALSE, gap = 3880))
})

testthat::test_that(".bt_compute_metric handles both metric methods", {
  cm <- getFromNamespace(".bt_compute_metric", "zeitR")

  s1 <- new.env(parent = emptyenv()); s1$metric_method <- 1; s1$metric_parameter <- 0.5
  # positive medians c(2,4,6,8) -> mean 5 -> 0.5 * 5 = 2.5 (zeros are ignored)
  expect_identical(cm(s1, c(0, 2, 4, 0, 6, 8)), 2.5)

  s2 <- new.env(parent = emptyenv()); s2$metric_method <- 2; s2$metric_parameter <- 0.5
  # median (type-7 quantile @0.5) of positive medians c(2,4,6,8) = 5
  expect_identical(cm(s2, c(0, 2, 4, 6, 8)), 5)

  s0 <- new.env(parent = emptyenv()); s0$metric_method <- 1; s0$metric_parameter <- 0.5
  expect_identical(cm(s0, c(0, 0, 0)), 0)   # no positive medians -> 0
})

testthat::test_that(".bt_compute_zero_proportion_around_end windows around the end", {
  zpa <- getFromNamespace(".bt_compute_zero_proportion_around_end", "zeitR")
  self <- new.env(parent = emptyenv())
  self$activity <- c(0, 0, 5, 0, 5, 0, 0, 0)
  self$half_window_around_border <- 2
  self$data_length <- 8L

  # end = 3 (0-indexed): window activity[1:5] (0-indexed) = c(0,5,0,5,0) -> 3/5
  expect_identical(zpa(self, 3L), 0.6)
})

# ── Stage 3: bed-time refiner — Unit B (filtering + candidates) ────────────────

# helper to build a peaks/valleys data.frame for tests
.mk_pv <- function(class, start, end, mean = 0, zero_proportion = 0,
                   above_threshold_proportion = 0) {
  n <- length(class)
  data.frame(
    class = class, start = as.integer(start), end = as.integer(end),
    length = as.integer(end - start),
    mean = rep_len(as.numeric(mean), n),
    median = rep_len(as.numeric(mean), n),
    zero_proportion = rep_len(as.numeric(zero_proportion), n),
    above_threshold_proportion = rep_len(as.numeric(above_threshold_proportion), n),
    stringsAsFactors = FALSE
  )
}

testthat::test_that(".bt_remove_after_long_valley keeps up to the first long valley", {
  f <- getFromNamespace(".bt_remove_after_long_valley", "zeitR")
  self <- new.env(parent = emptyenv())
  self$do_remove_after_long_valley <- TRUE
  self$bedtime_high_probability_sleep_valley_length <- 5

  pv <- .mk_pv(c("p", "v", "p", "v", "p"),
               start = c(0, 3, 5, 8, 18), end = c(3, 5, 8, 18, 21))  # 2nd valley length 10 >= 5
  out <- f(self, pv)
  expect_identical(out$class, c("p", "v", "p", "v"))   # keep rows 1..4
})

testthat::test_that(".bt_remove_before_long_peak keeps from the last long peak", {
  f <- getFromNamespace(".bt_remove_before_long_peak", "zeitR")
  self <- new.env(parent = emptyenv())
  self$do_remove_before_long_peak <- TRUE
  self$bedtime_high_probability_awake_peak_length <- 5

  pv <- .mk_pv(c("p", "v", "p", "v", "p"),
               start = c(0, 3, 5, 15, 17), end = c(3, 5, 15, 17, 20))  # peak row3 length 10 >= 5
  out <- f(self, pv)
  expect_identical(out$class, c("p", "v", "p"))         # keep rows 3..5
  expect_identical(out$start, c(5L, 15L, 17L))
})

testthat::test_that(".bt_identify_bedtime_candidates applies the interior-valley depth test", {
  f <- getFromNamespace(".bt_identify_bedtime_candidates", "zeitR")
  self <- new.env(parent = emptyenv())
  self$refinement_window_start <- 100L
  self$refinement_window_levels <- rep(0, 200)

  pv <- .mk_pv(c("p", "v", "p", "v", "p", "v"),
               start = c(0, 5, 10, 15, 20, 25), end = c(5, 10, 15, 20, 25, 30),
               mean  = c(50, 2, 40, 1, 45, 0.4))
  # valley row2 (r0=1, edge) -> start 5; valley row4 (interior) 1 < 0.5*40 -> start 15;
  # valley row6 (r0=5 == cnt-1, edge) -> start 25
  expect_identical(f(self, pv), c(5L, 15L, 25L))
})

testthat::test_that(".bt_bedtime_candidates_crossings_filter trims before last up-crossing", {
  f <- getFromNamespace(".bt_bedtime_candidates_crossings_filter", "zeitR")
  self <- new.env(parent = emptyenv())
  self$short_window_activity_median <- c(0, 0, 10, 10, 0, 0, 10, 0, 0, 0)
  self$refinement_window_start <- 0L
  self$refinement_window_end   <- 9L
  self$metric <- 5

  res <- f(self, c(2L, 4L, 8L))   # last up-crossing at position 6 -> keep >= 6
  expect_identical(res$candidates, 8L)
  expect_identical(res$down, c(4L, 7L))
})

# ── Stage 3: bed-time refiner — Unit C (scoring) ──────────────────────────

testthat::test_that(".bt_choose_best_bedtime_candidate picks the sharpest-drop candidate", {
  f <- getFromNamespace(".bt_choose_best_bedtime_candidate", "zeitR")
  self <- new.env(parent = emptyenv())
  self$refinement_window_start <- 100L
  # deep valley in the smoothed median-difference only around window-pos 20
  smoothed <- rep(0, 60); smoothed[21] <- -5   # R idx 21 == 0-indexed 20
  self$refinement_window_activity_median_difference_smoothed <- smoothed
  self$after_candidate_window <- 5L
  self$bedtime_score_last_candidate <- FALSE
  self$bedtime_last_candidate_score <- 0
  self$bedtime_best_crossing_distance_candidate_score <- 0
  self$bedtime_best_median_difference_candidate_score <- 0.7
  self$bedtime_best_epochs_above_metric_after_score <- 0
  self$bedtime_thresholded_candidate_score_amplitude <- 0.5
  self$bedtime_thresholded_candidate_score_minimum <- 0.1
  self$condition <- 0L
  self$bedtime_maximum_epochs_above_metric_after_candidate <- 2L
  self$consider_second_best_candidate <- FALSE
  self$initial_transition_candidate <- 0L
  self$metric <- 1000          # activity always below -> 0 epochs above, all thresholded
  self$maximum_allowed_gap <- 3600
  self$activity <- rep(0, 200)
  self$data_length <- 200L
  self$secs <- seq(0, by = 60, length.out = 200)   # regular spacing, no invalid gaps

  # candidate at window-pos 20 has the deepest amd valley -> gets median bonus -> wins
  expect_identical(f(self, c(2L, 20L, 40L), integer(0)), 120L)
})

testthat::test_that(".bt_choose_best_bedtime_candidate returns last candidate when condition not 0/2", {
  f <- getFromNamespace(".bt_choose_best_bedtime_candidate", "zeitR")
  self <- new.env(parent = emptyenv())
  self$refinement_window_start <- 100L
  self$refinement_window_activity_median_difference_smoothed <- rep(0, 60)
  self$after_candidate_window <- 5L
  self$bedtime_score_last_candidate <- FALSE
  self$bedtime_last_candidate_score <- 0
  self$bedtime_best_crossing_distance_candidate_score <- 0
  self$bedtime_best_median_difference_candidate_score <- 0.7
  self$condition <- 1L          # not in {0, 2} -> refined = rws + last original candidate
  self$initial_transition_candidate <- 0L

  expect_identical(f(self, c(5L, 15L, 25L), integer(0)), 125L)
})

testthat::test_that(".bt_choose_best_bedtime_candidate returns initial candidate when empty", {
  f <- getFromNamespace(".bt_choose_best_bedtime_candidate", "zeitR")
  self <- new.env(parent = emptyenv())
  self$initial_transition_candidate <- 4242L
  expect_identical(f(self, integer(0), integer(0)), 4242L)
})

# ── Stage 3: get-up-time refiner (Unit D) ─────────────────────────────

testthat::test_that(".gt_remove_after_long_tall_peak keeps up to the first long tall peak", {
  f <- getFromNamespace(".gt_remove_after_long_tall_peak", "zeitR")
  self <- new.env(parent = emptyenv())
  self$metric <- 10                      # > 5, so tall-peak rule is active
  self$do_remove_after_long_tall_peak <- TRUE
  self$getuptime_high_probability_awake_peak_length <- 5

  # peak row3: length 10 (>= 5) and mean 200 (> 10*metric = 100) -> first invalid
  pv <- .mk_pv(c("p", "v", "p", "v"),
               start = c(0, 3, 7, 17), end = c(3, 7, 17, 22),
               mean = c(50, 1, 200, 1))
  out <- f(self, pv)
  expect_identical(out$class, c("p", "v", "p"))   # keep rows 1..3
})

testthat::test_that(".gt_remove_before_long_valley keeps from the last long valley", {
  f <- getFromNamespace(".gt_remove_before_long_valley", "zeitR")
  self <- new.env(parent = emptyenv())
  self$do_remove_before_long_valley <- TRUE
  self$getuptime_high_probability_sleep_valley_length <- 5

  # valley row2 length 10 (>= 5) is the last invalid; (2-1) < (5-2) so it triggers
  pv <- .mk_pv(c("p", "v", "p", "v", "p"),
               start = c(0, 3, 13, 16, 18), end = c(3, 13, 16, 18, 21))
  out <- f(self, pv)
  expect_identical(out$class, c("v", "p", "v", "p"))   # keep rows 2..5
})

testthat::test_that(".gt_identify_getuptime_candidates returns peak starts and trailing valley end", {
  f <- getFromNamespace(".gt_identify_getuptime_candidates", "zeitR")
  self <- new.env(parent = emptyenv())
  self$refinement_window_start <- 100L
  self$refinement_window_levels <- rep(0, 200)
  self$data_length <- 200L

  pv <- .mk_pv(c("v", "p", "v", "p", "v"),
               start = c(0, 5, 10, 15, 20), end = c(5, 10, 15, 20, 25))
  # peaks at t0=1 (start 5) and t0=3 (start 15); trailing valley end 25
  expect_identical(f(self, pv), c(5L, 15L, 25L))
})

testthat::test_that(".gt_choose_best_getuptime_candidate picks the sharpest-rise candidate", {
  f <- getFromNamespace(".gt_choose_best_getuptime_candidate", "zeitR")
  self <- new.env(parent = emptyenv())
  self$refinement_window_start <- 100L
  smoothed <- rep(0, 60); smoothed[21] <- 5    # rise (peak) at window-pos 20
  self$refinement_window_activity_median_difference_smoothed <- smoothed
  self$after_candidate_window <- 5L
  self$getuptime_score_first_candidate <- FALSE
  self$getuptime_first_candidate_score <- 0
  self$getuptime_best_median_difference_candidate_score <- 0.7
  self$getuptime_best_crossing_distance_candidate_score <- 0
  self$getuptime_best_epochs_above_metric_after_score <- 0
  self$getuptime_thresholded_candidate_score_amplitude <- 0.5
  self$getuptime_thresholded_candidate_score_minimum <- 0.1
  self$getuptime_maximum_epochs_above_metric_after_candidate <- 2L
  self$metric <- 1
  self$activity <- rep(1000, 200)              # all above metric -> 0 epochs below -> thresholded
  self$data_length <- 200L
  self$maximum_allowed_gap <- 3600
  self$secs <- seq(0, by = 60, length.out = 200)
  self$initial_transition_candidate <- 0L

  # window-pos 20 has the sharpest rise -> median bonus -> wins
  expect_identical(f(self, c(2L, 20L, 40L), integer(0)), 120L)
})

testthat::test_that(".gt_choose_best_getuptime_candidate returns initial candidate when empty", {
  f <- getFromNamespace(".gt_choose_best_getuptime_candidate", "zeitR")
  self <- new.env(parent = emptyenv())
  self$initial_transition_candidate <- 777L
  expect_identical(f(self, integer(0), integer(0)), 777L)
})

# ── Stage 4: boolean_length_filter + ActTrust params ─────────────────────────

testthat::test_that(".boolean_length_filter removes short valleys and keeps long ones", {
  blf <- getFromNamespace(".boolean_length_filter", "zeitR")

  # short valley (length 2) between long peaks -> merged away (all wake)
  sig1 <- c(rep(1, 10), rep(0, 2), rep(1, 10))
  expect_true(all(blf(5, sig1) == 1))

  # long valley (length 8) -> kept unchanged (also exercises the no-op 2nd pass)
  sig2 <- c(rep(1, 10), rep(0, 8), rep(1, 10))
  expect_equal(blf(5, sig2), sig2)

  # single region (all wake) -> returned unchanged
  sig3 <- rep(1, 15)
  expect_equal(blf(5, sig3), sig3)

  # class_to_filter = "p": short peak (length 2) merged away (all sleep)
  sig4 <- c(rep(0, 10), rep(1, 2), rep(0, 10))
  expect_true(all(blf(5, sig4, class_to_filter = "p") == 0))
})

testthat::test_that(".cspd_acttrust_params matches the cspd_wrapper param_set", {
  pf <- getFromNamespace(".cspd_acttrust_params", "zeitR")
  p  <- pf()
  expect_equal(p$length_thresholds, c(8, 20, 17, 16))
  expect_equal(p$candidate_thresholds,
               c(0.4328571428571429, 0.37244897959183676, 0.5619047619047619))
  expect_equal(p$peak_valley_minimum_length, 11)
  expect_equal(p$median_filter_short_window, 41)
  expect_equal(p$sleep_minimum_length, 120)
  expect_true(p$bedtime_do_remove_before_long_peak)
  expect_false(p$bedtime_do_remove_before_tall_peak)
  expect_true(p$bedtime_consider_second_best_candidate)
  expect_true(p$getuptime_do_remove_after_long_tall_peak)
  expect_equal(p$bedtime_high_probability_awake_peak_length, 45)
  expect_equal(p$getuptime_high_probability_sleep_valley_length, 45)
})
