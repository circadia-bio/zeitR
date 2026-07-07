# Parity tests for Rcpp rolling filter implementations.
#
# Each test compares the C++ function against a reference R implementation
# using the exact same padding semantics as the original rolling_apply()-based
# functions in utils.R:
#   - Border replication by default (repeat first/last value)
#   - Constant padding for zero_prop_filter (pad_value = 1)
#
# Vectors tested:
#   - Small vector (n=9, hws=2) -- exercises all boundary cases
#   - Constant vector -- trivial case; verifies no division/NaN issues
#   - Zeros vector  -- important for zero_prop

ref_border <- function(x, hws, FUN) {
  vapply(seq_along(x), function(i) {
    pad <- c(rep(x[1L], hws), x, rep(x[length(x)], hws))
    FUN(pad[i:(i + 2L * hws)])
  }, numeric(1L))
}

ref_const_pad <- function(x, hws, pad_value, FUN) {
  vapply(seq_along(x), function(i) {
    pad <- c(rep(pad_value, hws), x, rep(pad_value, hws))
    FUN(pad[i:(i + 2L * hws)])
  }, numeric(1L))
}

x_small    <- c(1, 3, 2, 5, 4, 6, 8, 7, 9)
x_constant <- rep(4.0, 12)
x_zeros    <- rep(0.0, 8)
hws        <- 2L

# ── rolling_median_cpp ────────────────────────────────────────────────────────

test_that("rolling_median_cpp matches R reference on small vector", {
  expect_equal(
    rolling_median_cpp(x_small, hws),
    ref_border(x_small, hws, stats::median)
  )
})

test_that("rolling_median_cpp returns constant for constant input", {
  expect_equal(rolling_median_cpp(x_constant, hws), x_constant)
})

test_that("rolling_median_cpp handles length-1 input", {
  expect_equal(rolling_median_cpp(42, 3L), 42)
})

# ── rolling_mean_cpp ──────────────────────────────────────────────────────────

test_that("rolling_mean_cpp matches R reference on small vector", {
  expect_equal(
    rolling_mean_cpp(x_small, hws),
    ref_border(x_small, hws, mean)
  )
})

test_that("rolling_mean_cpp returns constant for constant input", {
  expect_equal(rolling_mean_cpp(x_constant, hws), x_constant)
})

# ── rolling_var_cpp ───────────────────────────────────────────────────────────

test_that("rolling_var_cpp matches R reference on small vector", {
  expect_equal(
    rolling_var_cpp(x_small, hws),
    ref_border(x_small, hws, stats::var)
  )
})

test_that("rolling_var_cpp returns zero for constant input", {
  expect_equal(rolling_var_cpp(x_constant, hws), rep(0, length(x_constant)))
})

# ── rolling_zero_prop_cpp ─────────────────────────────────────────────────────

test_that("rolling_zero_prop_cpp matches R reference with pad_value=1", {
  x_mixed <- c(0, 1, 0, 0, 2, 0, 1, 3)
  expect_equal(
    rolling_zero_prop_cpp(x_mixed, hws, 1.0),
    ref_const_pad(x_mixed, hws, 1.0, function(w) sum(w == 0) / length(w))
  )
})

test_that("rolling_zero_prop_cpp returns 1 for all-zero input with pad_value=0", {
  expect_equal(rolling_zero_prop_cpp(x_zeros, hws, 0.0), rep(1, length(x_zeros)))
})

test_that("rolling_zero_prop_cpp returns 0 for all-zero input with pad_value=1", {
  # boundary pads are 1 (non-zero), so proportion depends on window position
  result <- rolling_zero_prop_cpp(x_zeros, hws, 1.0)
  expect_true(all(result >= 0 & result <= 1))
  # interior epochs (no boundary pads): all zeros -> proportion 1
  win <- 2L * hws + 1L
  expect_equal(result[(hws + 1L):(length(x_zeros) - hws)], rep(1, length(x_zeros) - 2L * hws))
})

# ── rolling_quantile_cpp ──────────────────────────────────────────────────────

test_that("rolling_quantile_cpp matches stats::quantile type 7 on small vector", {
  expect_equal(
    rolling_quantile_cpp(x_small, hws, q = 0.6),
    ref_border(x_small, hws, function(w) stats::quantile(w, 0.6, names = FALSE))
  )
})

test_that("rolling_quantile_cpp q=0.5 matches rolling_median_cpp", {
  expect_equal(
    rolling_quantile_cpp(x_small, hws, q = 0.5),
    rolling_median_cpp(x_small, hws)
  )
})

test_that("rolling_quantile_cpp returns constant for constant input", {
  expect_equal(rolling_quantile_cpp(x_constant, hws, q = 0.6), x_constant)
})

# ── Edge cases ────────────────────────────────────────────────────────────────

test_that("all functions return empty numeric for length-0 input", {
  expect_equal(rolling_median_cpp(numeric(0), 2L),   numeric(0))
  expect_equal(rolling_mean_cpp(numeric(0), 2L),     numeric(0))
  expect_equal(rolling_var_cpp(numeric(0), 2L),      numeric(0))
  expect_equal(rolling_zero_prop_cpp(numeric(0), 2L), numeric(0))
  expect_equal(rolling_quantile_cpp(numeric(0), 2L), numeric(0))
})

test_that("hws=0 returns the input unchanged for mean and median", {
  expect_equal(rolling_mean_cpp(x_small, 0L),   as.double(x_small))
  expect_equal(rolling_median_cpp(x_small, 0L), as.double(x_small))
})

# ── diff5_cpp ─────────────────────────────────────────────────────────────────

# R reference implementation (original loop-based diff5)
ref_diff5 <- function(x, delta = 1) {
  n <- length(x)
  d <- numeric(n)
  for (i in seq(3L, n - 2L))
    d[i] <- (1 / (12 * delta)) * (-x[i+2] + 8*x[i+1] - 8*x[i-1] + x[i-2])
  d[1]   <- (1/(12*delta)) * (-25*x[1] + 48*x[2] - 36*x[3] + 16*x[4] -  3*x[5])
  d[2]   <- (1/(12*delta)) * (-25*x[2] + 48*x[3] - 36*x[4] + 16*x[5] -  3*x[6])
  d[n-1] <- (1/(12*delta)) * ( 25*x[n-1] - 48*x[n-2] + 36*x[n-3] - 16*x[n-4] + 3*x[n-5])
  d[n]   <- (1/(12*delta)) * ( 25*x[n]   - 48*x[n-1] + 36*x[n-2] - 16*x[n-3] + 3*x[n-4])
  d
}

test_that("diff5_cpp matches R reference on smooth input (x^2)", {
  x <- as.double(seq(1, 10)^2)
  expect_equal(diff5_cpp(x, 1.0), ref_diff5(x, 1))
})

test_that("diff5_cpp matches R reference on noisy input", {
  set.seed(42)
  x <- cumsum(rnorm(40))
  expect_equal(diff5_cpp(x, 1.0), ref_diff5(x, 1))
})

test_that("diff5_cpp respects delta scaling", {
  x <- as.double(seq(1, 10)^2)
  expect_equal(diff5_cpp(x, 30.0), ref_diff5(x, 30))
})

test_that("diff5_cpp errors on input shorter than 5", {
  expect_error(diff5_cpp(1:4, 1.0))
})

# ── score_epochs_cole_kripke_cpp ──────────────────────────────────────────────

# R reference (original vectorised implementation)
ref_cole_kripke <- function(zcm, P = 0.000464,
                             wb = c(34.5,133,529,375,408,400.5,1074,2048.5,2424.5),
                             wa = c(1920,149.5,257.5,125,111.5,120,69,40.5)) {
  n <- length(zcm); nb <- length(wb); s <- numeric(n)
  for (i in seq_along(wb)) {
    off <- nb - i + 1L
    if (off < n) s[seq(off+1,n)] <- s[seq(off+1,n)] + wb[i]*zcm[seq(1,n-off)]
  }
  for (i in seq_along(wa)) {
    off <- i
    if (off < n) s[seq(1,n-off)] <- s[seq(1,n-off)] + wa[i]*zcm[seq(off+1,n)]
  }
  as.integer(s * P >= 1.0)
}

test_that("score_epochs_cole_kripke_cpp matches R reference on random input", {
  set.seed(42)
  zcm <- as.double(rpois(500L, lambda = 50))
  expect_identical(score_epochs_cole_kripke_cpp(zcm), ref_cole_kripke(zcm))
})

test_that("score_epochs_cole_kripke_cpp handles all-zero input (all sleep)", {
  zcm <- as.double(rep(0, 100))
  expect_identical(score_epochs_cole_kripke_cpp(zcm), integer(100))
})

test_that("score_epochs_cole_kripke_cpp handles short vector (n < nb)", {
  zcm <- as.double(c(1000, 2000, 3000))
  expect_identical(
    score_epochs_cole_kripke_cpp(zcm),
    ref_cole_kripke(zcm)
  )
})

test_that("score_epochs_cole_kripke_cpp respects custom weights and P", {
  set.seed(7)
  zcm <- as.double(rpois(200L, lambda = 30))
  wb  <- rep(100, 9); wa <- rep(50, 8)
  expect_identical(
    score_epochs_cole_kripke_cpp(zcm, P = 0.001, wb, wa),
    ref_cole_kripke(zcm, P = 0.001, wb, wa)
  )
})
