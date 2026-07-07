// rolling_filters.cpp
// Fast rolling window filters that replace the vapply-based R implementations
// in utils.R. Semantics exactly match the R originals:
//   - Window of size (2*hws + 1) centred at each epoch position
//   - Border-replication padding by default (repeat first/last value)
//   - Constant padding when replicate = false (use pad_value)
//
// All functions are called from R wrappers and are not exported to the
// user-facing namespace.

#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>

using namespace Rcpp;

// ── Internal: build padded vector ─────────────────────────────────────────────

static inline std::vector<double> make_padded(
    const NumericVector& x, int hws, bool replicate, double pad_value
) {
  int n = x.size();
  std::vector<double> p(n + 2 * hws);

  double left  = replicate ? x[0]     : pad_value;
  double right = replicate ? x[n - 1] : pad_value;

  for (int i = 0; i < hws; i++)     p[i]           = left;
  for (int i = 0; i < n;   i++)     p[hws + i]     = x[i];
  for (int i = 0; i < hws; i++)     p[hws + n + i] = right;
  return p;
}

// ── rolling_median_cpp ────────────────────────────────────────────────────────
// Uses nth_element (O(win) per epoch) -- fast for the small windows used in
// zeitR (hws 5-10, so win 11-21).

// [[Rcpp::export]]
NumericVector rolling_median_cpp(
    NumericVector x, int hws,
    bool replicate = true, double pad_value = 0.0
) {
  int n = x.size();
  if (n == 0) return NumericVector(0);

  int win = 2 * hws + 1;
  std::vector<double> p = make_padded(x, hws, replicate, pad_value);
  NumericVector result(n);
  std::vector<double> w(win);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < win; j++) w[j] = p[i + j];
    std::nth_element(w.begin(), w.begin() + hws, w.end());
    result[i] = w[hws];
  }
  return result;
}

// ── rolling_mean_cpp ──────────────────────────────────────────────────────────
// Sliding window accumulator -- O(n) total.

// [[Rcpp::export]]
NumericVector rolling_mean_cpp(
    NumericVector x, int hws,
    bool replicate = true, double pad_value = 0.0
) {
  int n = x.size();
  if (n == 0) return NumericVector(0);

  int win = 2 * hws + 1;
  std::vector<double> p = make_padded(x, hws, replicate, pad_value);
  NumericVector result(n);

  double sum = 0.0;
  for (int j = 0; j < win; j++) sum += p[j];
  result[0] = sum / win;

  for (int i = 1; i < n; i++) {
    sum += p[i + win - 1] - p[i - 1];
    result[i] = sum / win;
  }
  return result;
}

// ── rolling_var_cpp ───────────────────────────────────────────────────────────
// Two-pass (mean then SS) for numerical stability; Bessel's correction (n-1)
// matches R's var(). Small windows make O(n*win) cost negligible.

// [[Rcpp::export]]
NumericVector rolling_var_cpp(
    NumericVector x, int hws,
    bool replicate = true, double pad_value = 0.0
) {
  int n = x.size();
  if (n == 0) return NumericVector(0);

  int win = 2 * hws + 1;
  std::vector<double> p = make_padded(x, hws, replicate, pad_value);
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < win; j++) sum += p[i + j];
    double mean = sum / win;
    double ss   = 0.0;
    for (int j = 0; j < win; j++) {
      double d = p[i + j] - mean;
      ss += d * d;
    }
    result[i] = ss / (win - 1);
  }
  return result;
}

// ── rolling_zero_prop_cpp ─────────────────────────────────────────────────────
// Sliding window counter -- O(n) total. Always uses constant padding
// (pad_value = 1 matches zero_prop_filter's default of padding with 1,
// so boundary windows are not artificially inflated with zero counts).

// [[Rcpp::export]]
NumericVector rolling_zero_prop_cpp(
    NumericVector x, int hws, double pad_value = 1.0
) {
  int n = x.size();
  if (n == 0) return NumericVector(0);

  int win = 2 * hws + 1;
  std::vector<double> p = make_padded(x, hws, false, pad_value);
  NumericVector result(n);

  int zeros = 0;
  for (int j = 0; j < win; j++) if (p[j] == 0.0) zeros++;
  result[0] = static_cast<double>(zeros) / win;

  for (int i = 1; i < n; i++) {
    if (p[i - 1]       == 0.0) zeros--;
    if (p[i + win - 1] == 0.0) zeros++;
    result[i] = static_cast<double>(zeros) / win;
  }
  return result;
}

// ── rolling_quantile_cpp ──────────────────────────────────────────────────────
// Type-7 interpolation (linear), matching stats::quantile's default.
// Formula: h = (win-1)*q; result = sorted[lo] + frac*(sorted[hi]-sorted[lo]).

// [[Rcpp::export]]
NumericVector rolling_quantile_cpp(
    NumericVector x, int hws, double q = 0.6,
    bool replicate = true, double pad_value = 0.0
) {
  int n = x.size();
  if (n == 0) return NumericVector(0);

  int win = 2 * hws + 1;
  std::vector<double> p = make_padded(x, hws, replicate, pad_value);
  NumericVector result(n);
  std::vector<double> w(win);

  double h    = (win - 1) * q;
  int    lo   = static_cast<int>(h);
  double frac = h - lo;
  int    hi   = (lo + 1 < win) ? lo + 1 : lo;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < win; j++) w[j] = p[i + j];
    std::sort(w.begin(), w.end());
    result[i] = w[lo] + frac * (w[hi] - w[lo]);
  }
  return result;
}

// ── diff5_cpp ─────────────────────────────────────────────────────────────────
// Five-point stencil first-derivative estimate.
// Interior points use the central formula; four boundary points use
// one-sided forward/backward stencils. Matches R's diff5() exactly.

// [[Rcpp::export]]
NumericVector diff5_cpp(NumericVector x, double delta = 1.0) {
  int n = x.size();
  if (n < 5) Rcpp::stop("diff5_cpp requires length >= 5");

  NumericVector d(n);
  double c = 1.0 / (12.0 * delta);

  // Interior: central five-point stencil (0-based i = 2 .. n-3)
  for (int i = 2; i <= n - 3; i++) {
    d[i] = c * (-x[i + 2] + 8.0 * x[i + 1] - 8.0 * x[i - 1] + x[i - 2]);
  }

  // Boundary: forward stencil at i=0 and i=1
  d[0] = c * (-25.0*x[0] + 48.0*x[1] - 36.0*x[2] + 16.0*x[3] -  3.0*x[4]);
  d[1] = c * (-25.0*x[1] + 48.0*x[2] - 36.0*x[3] + 16.0*x[4] -  3.0*x[5]);

  // Boundary: backward stencil at i=n-2 and i=n-1
  d[n - 2] = c * (25.0*x[n-2] - 48.0*x[n-3] + 36.0*x[n-4] - 16.0*x[n-5] + 3.0*x[n-6]);
  d[n - 1] = c * (25.0*x[n-1] - 48.0*x[n-2] + 36.0*x[n-3] - 16.0*x[n-4] + 3.0*x[n-5]);

  return d;
}

// ── rolling_max_cpp ──────────────────────────────────────────────────────────
// Replaces rolling_apply(x, hws, max, pad_value) in .morphological_open_close().
// Uses constant-0 padding (replicate=false, pad_value=0) for morphological ops.

// [[Rcpp::export]]
NumericVector rolling_max_cpp(
    NumericVector x, int hws,
    bool replicate = true, double pad_value = 0.0
) {
  int n = x.size();
  if (n == 0) return NumericVector(0);

  int win = 2 * hws + 1;
  std::vector<double> p = make_padded(x, hws, replicate, pad_value);
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double mx = p[i];
    for (int j = 1; j < win; j++) if (p[i + j] > mx) mx = p[i + j];
    result[i] = mx;
  }
  return result;
}

// ── rolling_min_cpp ──────────────────────────────────────────────────────────
// Replaces rolling_apply(x, hws, min, pad_value) in .morphological_open_close().

// [[Rcpp::export]]
NumericVector rolling_min_cpp(
    NumericVector x, int hws,
    bool replicate = true, double pad_value = 0.0
) {
  int n = x.size();
  if (n == 0) return NumericVector(0);

  int win = 2 * hws + 1;
  std::vector<double> p = make_padded(x, hws, replicate, pad_value);
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double mn = p[i];
    for (int j = 1; j < win; j++) if (p[i + j] < mn) mn = p[i + j];
    result[i] = mn;
  }
  return result;
}

// ── zero_mitigation_cpp ───────────────────────────────────────────────────────
// Pass 1 of .crespo_msp(): add mitigation_level to runs of zeros longer than
// consec_zeros_thr. Trailing zero runs are NOT mitigated (matches R/Python).

// [[Rcpp::export]]
NumericVector zero_mitigation_cpp(
    NumericVector activity, int consec_zeros_thr, double mitigation_level
) {
  int n = activity.size();
  NumericVector result = clone(activity);
  int run = 0;

  for (int i = 0; i < n; i++) {
    if (activity[i] == 0.0) {
      run++;
    } else {
      if (run > consec_zeros_thr) {
        for (int j = i - run; j < i; j++) result[j] += mitigation_level;
      }
      run = 0;
    }
  }
  // Trailing zero run: NOT processed -- matches R/Python loop behaviour.
  return result;
}

// ── mark_invalid_zeros_cpp ───────────────────────────────────────────────────
// Pass 2 of .crespo_msp(): scan activity against the morphological detection
// and return 0-indexed positions of invalid zeros. The caller shifts by
// pad_size to get padded-vector positions and sets those to NA.
// Returns a sorted, deduplicated integer vector of 0-indexed signal positions.

// [[Rcpp::export]]
IntegerVector mark_invalid_zeros_cpp(
    NumericVector activity,
    IntegerVector morph_detection,   // 1 = awake, 0 = sleep
    int awake_zeros_thr,
    int sleep_zeros_thr
) {
  int n = activity.size();
  std::vector<int> invalid;
  int awake_run = 0, sleep_run = 0;

  for (int i = 0; i < n; i++) {
    if (morph_detection[i] == 1) {  // awake
      if (sleep_run > sleep_zeros_thr) {
        for (int j = i - sleep_run; j < i; j++) invalid.push_back(j);
      }
      sleep_run = 0;

      if (activity[i] == 0.0) {
        awake_run++;
      } else {
        if (awake_run > awake_zeros_thr) {
          for (int j = i - awake_run; j < i; j++) invalid.push_back(j);
        }
        awake_run = 0;
      }
    } else {  // sleep
      if (awake_run > awake_zeros_thr) {
        for (int j = i - awake_run; j < i; j++) invalid.push_back(j);
      }
      awake_run = 0;

      if (activity[i] == 0.0) {
        sleep_run++;
      } else {
        if (sleep_run > sleep_zeros_thr) {
          for (int j = i - sleep_run; j < i; j++) invalid.push_back(j);
        }
        sleep_run = 0;
      }
    }
  }
  // Trailing runs: NOT processed -- matches R/Python loop behaviour.

  std::sort(invalid.begin(), invalid.end());
  invalid.erase(std::unique(invalid.begin(), invalid.end()), invalid.end());
  return IntegerVector(invalid.begin(), invalid.end());
}

// ── adaptive_median_filter_cpp ───────────────────────────────────────────────
// Variable half-window adaptive median filter used in .crespo_msp() and
// .crespo_nap_msp(). The window grows from pad_size up to max_hws then
// shrinks back. NaN/NA values in padded_activity are skipped (na.rm=TRUE);
// if the whole window is NA the previous output value is carried forward.
//
// padded_activity : already-padded vector (length n + 2*pad_size, may contain
//                   NA_real_ at invalid-zero positions)
// n               : length of the original (un-padded) signal
// pad_size        : initial (minimum) half-window size; also the padding size
// max_hws         : maximum half-window size
//
// Window growth condition (mirrors Python): i < n - max_hws + pad_size - 1

// [[Rcpp::export]]
NumericVector adaptive_median_filter_cpp(
    NumericVector padded_activity, int n, int pad_size, int max_hws
) {
  int N = padded_activity.size();
  NumericVector result(n);
  int hws = pad_size;

  std::vector<double> window;
  window.reserve(2 * max_hws + 1);

  for (int i = 0; i < n; i++) {
    int center = i + pad_size;          // 0-indexed position in padded_activity
    int lo     = std::max(0, center - hws);
    int hi     = std::min(N - 1, center + hws);

    window.clear();
    for (int j = lo; j <= hi; j++) {
      double v = padded_activity[j];
      if (!ISNAN(v)) window.push_back(v);
    }

    double val;
    if (window.empty()) {
      // All NA: carry forward previous output (or 0 at start)
      val = (i > 0) ? result[i - 1] : 0.0;
    } else {
      int m = window.size();
      if (m % 2 == 1) {
        // Odd: exact middle element
        std::nth_element(window.begin(), window.begin() + m / 2, window.end());
        val = window[m / 2];
      } else {
        // Even: average of two middle elements (matches R's stats::median)
        std::nth_element(window.begin(), window.begin() + m / 2, window.end());
        double hi_val = window[m / 2];
        std::nth_element(window.begin(), window.begin() + m / 2 - 1, window.begin() + m / 2);
        val = (window[m / 2 - 1] + hi_val) / 2.0;
      }
    }
    result[i] = val;

    // Variable window: grow toward max_hws, then shrink near end
    if (i < n - max_hws + pad_size - 1) {
      if (hws < max_hws) hws++;
    } else {
      if (hws > pad_size) hws--;
    }
  }
  return result;
}

// ── score_epochs_cole_kripke_cpp ──────────────────────────────────────────────
// Single-pass Cole-Kripke scorer. For each epoch i, accumulates:
//   before: weights_before[nb - lag] * zcm[i - lag]  for lag = 1..nb
//   after:  weights_after[j]         * zcm[i + j+1]  for j   = 0..na-1
// Multiplies by P and thresholds at 1.0. Matches score_epochs_cole_kripke()
// exactly; replaces 17 vectorised R additions with one O(n) C++ loop.

// [[Rcpp::export]]
IntegerVector score_epochs_cole_kripke_cpp(
    NumericVector zcm,
    double P = 0.000464,
    NumericVector weights_before = NumericVector::create(
        34.5, 133.0, 529.0, 375.0, 408.0, 400.5, 1074.0, 2048.5, 2424.5),
    NumericVector weights_after = NumericVector::create(
        1920.0, 149.5, 257.5, 125.0, 111.5, 120.0, 69.0, 40.5)
) {
  int n  = zcm.size();
  int nb = weights_before.size();
  int na = weights_after.size();

  IntegerVector result(n, 0);
  if (n < 2) return result;

  for (int i = 0; i < n; i++) {
    double score = 0.0;

    // Before: weight for lag is weights_before[nb - lag] (matches R indexing)
    for (int lag = 1; lag <= nb; lag++) {
      int src = i - lag;
      if (src >= 0) score += weights_before[nb - lag] * zcm[src];
    }

    // After: weights_after[j] applied at lead j+1
    for (int j = 0; j < na; j++) {
      int src = i + j + 1;
      if (src < n) score += weights_after[j] * zcm[src];
    }

    result[i] = (score * P >= 1.0) ? 1 : 0;
  }

  return result;
}
