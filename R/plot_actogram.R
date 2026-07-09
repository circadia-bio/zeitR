# ── Actogram visualisations ───────────────────────────────────────────────────
#
# Three display modes for actigraphy epoch-level state data:
#   plot_actogram()          -- single-column raster (standard)
#   plot_actogram_double()   -- double-plotted raster (classic chronobiology)
#   plot_actogram_activity() -- double-plotted with ZCMn activity bars
#
# All three accept a `zeitr_result` list or a bare tibble with `datetime` +
# `state` columns.  ggplot2 is a Suggested dependency; a clear error is thrown
# if it is not installed.

utils::globalVariables(c(
  "mins_since_midnight", "date_fct",
  "x_mins", "row_date_fct",
  "ymin", "ymax",
  "state_label"
))

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @noRd
.check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    zeitr_abort(
      "Package {.pkg ggplot2} is required for actogram plots. Install with: {.code install.packages(\"ggplot2\")}."
    )
}

#' Default state colour palette for actogram plots
#'
#' Returns a named character vector of hex colours keyed on the four state
#' labels produced by [label_states()].  Pass the output to the `colours`
#' argument of any actogram function to override individual colours.
#'
#' @return Named character vector with elements `"wake"`, `"sleep"`, `"nap"`,
#'   and `"off-wrist"`.
#'
#' @export
#'
#' @examples
#' actogram_colours()
#' # wake        sleep        nap    off-wrist
#' # "#C25E2A" "#3B2F6B" "#F0A500" "#D9C8A0"
actogram_colours <- function() {
  c(
    "wake"      = "#C25E2A",
    "sleep"     = "#3B2F6B",
    "nap"       = "#F0A500",
    "off-wrist" = "#D9C8A0"
  )
}

# Resolve a zeitr_result or bare tibble into an epoch data frame, and attach
# computed columns needed by all three actogram variants.
#' @noRd
.actogram_prep <- function(result, tz) {
  data <- if (is.list(result) && !is.data.frame(result)) {
    if (is.null(result[["data"]]))
      zeitr_abort(
        "{.arg result} must be a `zeitr_result` list with a {.code $data} element, or a tibble with {.code datetime} and {.code state} columns."
      )
    result[["data"]]
  } else {
    result
  }

  missing_cols <- setdiff(c("datetime", "state"), names(data))
  if (length(missing_cols) > 0L)
    zeitr_abort("{.arg result} is missing column(s): {.val {missing_cols}}.")

  tz_use <- if (!is.null(tz)) {
    tz
  } else {
    tz_attr <- attr(data[["datetime"]], "tzone")
    if (!is.null(tz_attr) && nchar(tz_attr) > 0L) tz_attr else "UTC"
  }

  dt <- as.POSIXct(data[["datetime"]])
  data[["date"]] <- as.Date(format(dt, tz = tz_use, format = "%Y-%m-%d"))
  data[["mins_since_midnight"]] <-
    as.integer(format(dt, tz = tz_use, format = "%H")) * 60L +
    as.integer(format(dt, tz = tz_use, format = "%M"))
  data[["state_label"]] <- label_states(data[["state"]])

  list(data = data, tz = tz_use)
}

# Build a default title for an actogram function, optionally including the
# subject ID from a zeitr_result.
#' @noRd
.actogram_title <- function(result, prefix) {
  subj <- if (is.list(result) && !is.null(result[["subject_id"]])) {
    result[["subject_id"]]
  } else {
    NULL
  }
  if (!is.null(subj)) paste0(prefix, " \u2014 ", subj) else prefix
}

# Shared minimal theme used by all three variants.
#' @noRd
.actogram_theme <- function(base_size) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      legend.position    = "bottom",
      legend.title       = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold"),
      axis.text.y        = ggplot2::element_text(size = base_size - 3L)
    )
}

# Duplicate each day into a left (x = 0-1439) and right (x = 1440-2879) column
# keyed on a `row_date` variable.  The left column holds the day itself; the
# right column holds the same day advanced one step so it appears in the row
# for the previous date.  Rows outside the recording span are dropped.
#' @noRd
.build_double_data <- function(d) {
  dates   <- sort(unique(d[["date"]]))
  d_left  <- d
  d_right <- d

  d_left[["row_date"]]  <- d_left[["date"]]
  d_left[["x_mins"]]    <- d_left[["mins_since_midnight"]]

  d_right[["row_date"]] <- d_right[["date"]] - 1L
  d_right[["x_mins"]]   <- d_right[["mins_since_midnight"]] + 1440L

  d_double <- rbind(d_left, d_right)
  d_double[d_double[["row_date"]] >= min(dates) &
             d_double[["row_date"]] <= max(dates), , drop = FALSE]
}

# Shared x-axis scale for the double-plot variants (48 h window).
# Limits are widened by half an epoch on each side: geom_tile() centers each
# tile on its x value, so the epoch at x = 0 (midnight) would otherwise have
# its left half clipped by a hard limit of exactly 0 (and ggplot2 drops the
# whole tile with a "missing values" warning). Symmetric padding on the right
# guards against the same issue at the 48 h boundary.
#' @noRd
.scale_x_double <- function(epoch_min) {
  ggplot2::scale_x_continuous(
    breaks = seq(0L, 42L * 60L, 6L * 60L),
    labels = function(x) sprintf("%02d:00", (x %/% 60L) %% 24L),
    limits = c(-epoch_min / 2, 48L * 60L + epoch_min / 2),
    expand = c(0, 0)
  )
}

# ── plot_actogram ─────────────────────────────────────────────────────────────

#' Plot a single-column actogram
#'
#' Renders epoch-level sleep/wake states as a raster with one row per calendar
#' day and time-of-day (00:00 to 24:00) on the x-axis.  The oldest day is at
#' the top, following standard chronobiology convention.
#'
#' @param result A `zeitr_result` list (from [run_pipeline()] or
#'   [run_pipeline_native()]), or a tibble with at least `datetime` and
#'   `state` columns.
#' @param tz `character(1)` or `NULL`. Timezone for date and time-of-day
#'   extraction.  `NULL` (default) auto-detects from the timezone attribute
#'   embedded in the `datetime` POSIXct column; falls back to `"UTC"` when
#'   absent.
#' @param title `character(1)` or `NULL`. Plot title.  `NULL` constructs
#'   `"Actogram \u2014 <subject_id>"` from `result$subject_id` when available.
#' @param colours Named character vector mapping state labels (`"wake"`,
#'   `"sleep"`, `"nap"`, `"off-wrist"`) to hex colours.  `NULL` uses
#'   [actogram_colours()].
#' @param date_label_every `integer(1)`. Label every Nth row on the y-axis.
#'   Default `7` (weekly ticks).
#' @param epoch_min `numeric(1)`. Epoch duration in minutes.  Used as the tile
#'   width in [ggplot2::geom_tile()].  Default `1` (ActTrust standard).
#' @param base_size `numeric(1)`. Base font size passed to
#'   [ggplot2::theme_minimal()].  Default `13`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_actogram_double()], [plot_actogram_activity()],
#'   [actogram_colours()], [label_states()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")
#'
#' # Default colours
#' plot_actogram(result)
#'
#' # Custom colours
#' plot_actogram(result, colours = c(wake = "#E8D5B0", sleep = "#1C1A2E",
#'                                   nap = "#F0A500", "off-wrist" = "#C25E2A"))
#'
#' # Label every 14 days
#' plot_actogram(result, date_label_every = 14L)
#' }
plot_actogram <- function(result,
                          tz               = NULL,
                          title            = NULL,
                          colours          = NULL,
                          date_label_every = 7L,
                          epoch_min        = 1,
                          base_size        = 13) {
  .check_ggplot2()

  inp <- .actogram_prep(result, tz)
  d   <- inp[["data"]]

  col_use   <- if (is.null(colours)) actogram_colours() else colours
  title_use <- if (!is.null(title)) title else .actogram_title(result, "Actogram")

  dates <- sort(unique(d[["date"]]))
  d[["date_fct"]] <- factor(
    as.character(d[["date"]]),
    levels = rev(as.character(dates))
  )

  ggplot2::ggplot(d, ggplot2::aes(
    x    = mins_since_midnight,
    y    = date_fct,
    fill = state_label
  )) +
    ggplot2::geom_tile(height = 0.9, width = epoch_min) +
    ggplot2::scale_fill_manual(values = col_use, name = "State", drop = TRUE) +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, 23L * 60L, 4L * 60L),
      labels = function(x) sprintf("%02d:00", x %/% 60L),
      expand = c(0, 0),
      # Widened by half an epoch on each side: geom_tile() centers each tile
      # on its x value, so the midnight epoch (x = 0) would otherwise have
      # its left half clipped by a hard limit of exactly 0.
      limits = c(-epoch_min / 2, 24L * 60L + epoch_min / 2)
    ) +
    ggplot2::scale_y_discrete(
      breaks = function(x) x[seq(1L, length(x), by = as.integer(date_label_every))],
      labels = function(x) format(as.Date(x), "%d %b"),
      expand = c(0.01, 0.01)
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title_use) +
    .actogram_theme(base_size)
}

# ── plot_actogram_double ──────────────────────────────────────────────────────

#' Plot a double-plotted actogram
#'
#' The classic chronobiology double-plot format.  Each recording day is drawn
#' twice: in the left column of its own row (x = 00:00 to 24:00) and in the
#' right column of the row above (x = 24:00 to 48:00).  This means consecutive
#' pairs of days share a row, making circadian phase drift visible as a
#' diagonal band across rows.
#'
#' @details
#' Row \eqn{i} shows day \eqn{i} on the left and day \eqn{i+1} on the right.
#' Every day therefore appears twice in the plot (except the first, which has
#' no left-column predecessor, and the last, which has no right-column
#' successor).  A dashed vertical line marks the 24 h boundary between the
#' two columns.
#'
#' @inheritParams plot_actogram
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_actogram()], [plot_actogram_activity()], [actogram_colours()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")
#' plot_actogram_double(result)
#' }
plot_actogram_double <- function(result,
                                 tz               = NULL,
                                 title            = NULL,
                                 colours          = NULL,
                                 date_label_every = 7L,
                                 epoch_min        = 1,
                                 base_size        = 13) {
  .check_ggplot2()

  inp <- .actogram_prep(result, tz)
  d   <- inp[["data"]]

  col_use   <- if (is.null(colours)) actogram_colours() else colours
  title_use <- if (!is.null(title)) title else
    .actogram_title(result, "Double actogram")

  d_double    <- .build_double_data(d)
  row_dates   <- sort(unique(d_double[["row_date"]]))

  # Factor: newest = first level = bottom; oldest = last level = top
  d_double[["row_date_fct"]] <- factor(
    as.character(d_double[["row_date"]]),
    levels = rev(as.character(row_dates))
  )

  ggplot2::ggplot(d_double, ggplot2::aes(
    x    = x_mins,
    y    = row_date_fct,
    fill = state_label
  )) +
    ggplot2::geom_tile(height = 0.9, width = epoch_min) +
    ggplot2::geom_vline(
      xintercept = 1440L,
      colour     = "grey45",
      linewidth  = 0.4,
      linetype   = "dashed"
    ) +
    ggplot2::scale_fill_manual(values = col_use, name = "State", drop = TRUE) +
    .scale_x_double(epoch_min) +
    ggplot2::scale_y_discrete(
      breaks = function(x) x[seq(1L, length(x), by = as.integer(date_label_every))],
      labels = function(x) format(as.Date(x), "%d %b"),
      expand = c(0.01, 0.01)
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title_use) +
    .actogram_theme(base_size)
}

# ── plot_actogram_activity ────────────────────────────────────────────────────

#' Plot a double-plotted actogram with activity bars
#'
#' Same double-plot layout as [plot_actogram_double()] — each day appears in
#' both the left and right column of adjacent rows — but each epoch is rendered
#' as a vertical bar whose height is proportional to the raw activity count
#' (default: `ZCMn`, the zero-crossing mean from ActTrust).  Bars are coloured
#' by sleep/wake state, so the plot simultaneously conveys activity intensity
#' and state classification.
#'
#' @details
#' Bar heights are capped at the `activity_cap_quantile` quantile of
#' non-zero epochs to prevent outlier activity bursts from compressing the
#' visible range for the rest of the recording. A thin baseline stub
#' (2 % of row height) is drawn for zero-activity epochs so that sleep periods
#' and off-wrist blocks remain faintly visible.
#'
#' Actigraphy activity counts are typically right-skewed, with occasional
#' bursts far above the typical waking level. On the default linear scale
#' this compresses most of the meaningful variation among low-to-moderate
#' activity epochs into a thin sliver near the baseline. Set `log_scale =
#' TRUE` to apply a `log1p()` transform (`log(1 + x)`, so zero-activity
#' epochs map to `0` rather than `-Inf`) before capping and normalising --
#' this expands the low-activity range at the cost of visually compressing
#' the difference between already-high activity bursts.
#'
#' @inheritParams plot_actogram
#' @param activity_col `character(1)`. Name of the column in `result$data` to
#'   use as the activity signal. Default `"ZCMn"` (zero-crossing mean;
#'   present in all ActTrust recordings processed by `zeitR`).
#' @param activity_cap_quantile `numeric(1)` in (0, 1]. Quantile of non-zero
#'   activity values used to cap bar heights before normalising. Default
#'   `0.99` (top 1 % of active epochs are clipped to full bar height).
#' @param log_scale `logical(1)`. If `TRUE`, applies a `log1p()` transform to
#'   the activity signal before capping and normalising, compressing the
#'   dynamic range so lower-activity variation is easier to see. Default
#'   `FALSE` (linear scale).
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_actogram()], [plot_actogram_double()], [actogram_colours()]
#'
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")
#' plot_actogram_activity(result)
#'
#' # Cap at 95th percentile to highlight moderate activity
#' plot_actogram_activity(result, activity_cap_quantile = 0.95)
#'
#' # Log scale: reveals structure among low-activity epochs that a linear
#' # scale would otherwise compress near the baseline
#' plot_actogram_activity(result, log_scale = TRUE)
#' }
plot_actogram_activity <- function(result,
                                   tz                    = NULL,
                                   title                 = NULL,
                                   colours               = NULL,
                                   activity_col          = "ZCMn",
                                   activity_cap_quantile = 0.99,
                                   log_scale             = FALSE,
                                   date_label_every      = 7L,
                                   epoch_min             = 1,
                                   base_size             = 13) {
  .check_ggplot2()

  inp <- .actogram_prep(result, tz)
  d   <- inp[["data"]]

  if (!activity_col %in% names(d))
    zeitr_abort(
      "Activity column {.val {activity_col}} not found in the data. Check that the pipeline was run on a device that provides this signal, or supply a different {.arg activity_col}."
    )

  col_use   <- if (is.null(colours)) actogram_colours() else colours
  title_use <- if (!is.null(title)) title else
    .actogram_title(result, "Activity actogram")

  # ── Build double-plot data ────────────────────────────────────────────────
  d_double   <- .build_double_data(d)
  row_dates  <- sort(unique(d_double[["row_date"]]))  # ascending (oldest first)
  n_rows     <- length(row_dates)

  # Numeric y: 1 = oldest row, n_rows = newest row.
  # scale_y_reverse() puts 1 at the top so oldest is at the top.
  d_double[["y_row"]] <- match(d_double[["row_date"]], row_dates)

  # ── Normalise activity ────────────────────────────────────────────────────
  # log_scale applies log1p() (never -Inf at zero activity) before capping
  # and normalising, so low-activity structure is easier to see against a
  # right-skewed raw signal.
  act_transform <- if (isTRUE(log_scale)) function(x) log1p(pmax(x, 0)) else identity

  act_raw <- act_transform(d[[activity_col]])
  pos_act <- act_raw[is.finite(act_raw) & act_raw > 0]
  cap     <- if (length(pos_act) > 0L) {
    quantile(pos_act, activity_cap_quantile, na.rm = TRUE)
  } else {
    1
  }
  if (!is.finite(cap) || cap <= 0) cap <- max(act_raw, na.rm = TRUE)
  if (!is.finite(cap) || cap <= 0) cap <- 1

  d_double[["act_norm"]] <- pmin(
    act_transform(d_double[[activity_col]]) / cap,
    1.0
  )
  d_double[["act_norm"]][is.na(d_double[["act_norm"]])] <- 0

  # ── Bar geometry (anchored at bottom of each row, growing upward) ─────────
  # With scale_y_reverse(), increasing y goes downward on screen.
  # "Bottom of row" = y_row + 0.46 (higher numeric y = lower on screen).
  # Bars grow toward lower y (upward on screen).
  bar_available <- 0.90   # fraction of row height available for bars
  min_stub      <- 0.02   # minimum bar as fraction of bar_available

  d_double[["ymax"]] <- d_double[["y_row"]] + 0.46
  d_double[["ymin"]] <- d_double[["ymax"]] -
    pmax(d_double[["act_norm"]], min_stub) * bar_available

  # ── Y-axis ticks ─────────────────────────────────────────────────────────
  tick_idx  <- seq(1L, n_rows, by = as.integer(date_label_every))
  tick_pos  <- tick_idx
  tick_labs <- format(row_dates[tick_idx], "%d %b")

  ggplot2::ggplot(d_double) +
    ggplot2::geom_rect(ggplot2::aes(
      xmin = x_mins,
      xmax = x_mins + epoch_min,
      ymin = ymin,
      ymax = ymax,
      fill = state_label
    )) +
    ggplot2::geom_vline(
      xintercept = 1440L,
      colour     = "grey45",
      linewidth  = 0.4,
      linetype   = "dashed"
    ) +
    ggplot2::scale_fill_manual(values = col_use, name = "State", drop = TRUE) +
    .scale_x_double(epoch_min) +
    ggplot2::scale_y_reverse(
      breaks = tick_pos,
      labels = tick_labs,
      expand = c(0.01, 0.01)
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title_use) +
    .actogram_theme(base_size)
}
