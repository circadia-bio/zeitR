// rolling_filters.cpp
// Fast rolling window filters that replace the vapply-based R implementations
// in utils.R. Semantics exactly match the R originals:
//   - Window of size (2*hws + 1) centred at each epoch position
//   - Border-replication padding by default (repeat first/last value)
//   - Constant padding when replicate = false (use pad_value)
//
// All five functions are called from R wrappers (see R/utils.R) and are not
// exported to the user-facing namespace.

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
// Uses nth_element (O(win) per epoch) — fast for the small windows used in
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
