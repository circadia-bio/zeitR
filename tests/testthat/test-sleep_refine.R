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
