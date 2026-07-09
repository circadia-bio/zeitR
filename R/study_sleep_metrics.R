# ── Internal: expected column names for compute_sleep_metrics()/compute_cpd_metrics() ──
# Generated the same way .group_metrics() itself builds names (base names
# suffixed with "_wd"/"_fd"), so a participant whose metrics computation
# fails (or silently returns an empty list, as compute_sleep_metrics() does
# when zero nights pass the TIB filter) can still be given a full-width
# NA row -- keeping every participant's row the same shape for rbind().

.sm_metric_names <- c(
  "sleep_onset_h", "sleep_offset_h", "fpr_tib_h", "fps_h",
  "tst_h", "latencia_min", "inertia_min", "waso_min",
  "sleep_eff_pct", "tst_24h_h", "dp_midsleep_min", "dp_tst_min"
)

.sm_all_names <- c(
  "n_overall", .sm_metric_names,
  "n_wd", paste0(.sm_metric_names, "_wd"),
  "n_fd", paste0(.sm_metric_names, "_fd")
)

.cpd_all_names <- c(
  "n_nights_cpd", "n_free_days", "n_workdays",
  "msw_h", "msw_hms", "msf_h", "msf_hms", "msfsc_h", "msfsc_hms",
  "sjl_h", "sjl_min", "sjla_h", "sjla_min",
  "cpd_s", "cpd_min", "cpd_h"
)

.cpd_char_names <- c("msw_hms", "msf_hms", "msfsc_hms")

#' @noRd
.sm_na_row <- function() {
  stats::setNames(as.list(rep(NA_real_, length(.sm_all_names))), .sm_all_names)
}

#' @noRd
.cpd_na_row <- function() {
  vals <- as.list(rep(NA_real_, length(.cpd_all_names)))
  vals[.cpd_all_names %in% .cpd_char_names] <- list(NA_character_)
  stats::setNames(vals, .cpd_all_names)
}

#' Batch sleep-timing and chronotype metrics across a study
#'
#' Computes [compute_sleep_metrics()] and [compute_cpd_metrics()] for every
#' participant in a batch of pipeline results and stacks them into a single
#' tibble with one row per participant -- the sleep-timing/chronotype
#' counterpart to [study_summary()], which covers NPCRA activity-rhythm
#' variables instead.
#'
#' @details
#' [compute_sleep_metrics()] and [compute_cpd_metrics()] each return a single
#' named list per participant with no participant identifier, and there was
#' previously no batch wrapper analogous to [study_summary()] for them. This
#' is that wrapper -- intended to make these chronobiological phenotyping
#' metrics database-ready for tools like `syncR::sync()`, which expect one
#' row per participant with a shared `participant_id` column across sources.
#'
#' If either metric computation fails for a participant (e.g. no free days
#' found, or no nights pass the `min_tib_h` filter), that participant's row
#' is filled with `NA` for the affected metrics and a warning is emitted --
#' the rest of the study is unaffected.
#'
#' @param results A named list of `zeitr_result` objects, as returned by
#'   [run_pipeline_batch()] or [run_pipeline_native_batch()]. `participant_id`
#'   is taken from each result's own `$subject_id`, falling back to the list
#'   name if that is unavailable.
#' @param min_tib_h `numeric(1)`. Minimum total in-bed time (hours) for a
#'   night to be included in [compute_sleep_metrics()]. Default `5.0`.
#' @param min_tib_eve_h `numeric(1)`. Minimum TBT (hours) for a night to
#'   qualify as a free-day-eve night in [compute_cpd_metrics()]. Default `3.0`.
#' @param tz `character(1)`. Time zone for extracting clock hours. Default
#'   `"UTC"`.
#' @param holidays,free_days Forwarded to both [compute_sleep_metrics()] and
#'   [compute_cpd_metrics()] for every participant. Default `NULL` for both,
#'   in which case each participant's own `result$holidays`/`result$free_days`
#'   (set when the pipeline was run) are used instead. Supplying either here
#'   overrides that per-participant default for the whole study.
#'
#' @return A tibble with one row per participant: `participant_id`, all
#'   [compute_sleep_metrics()] columns (`n_overall`/`n_wd`/`n_fd` and the
#'   twelve sleep-timing metrics with `_wd`/`_fd` suffixes), and all
#'   [compute_cpd_metrics()] columns (`n_nights_cpd`, `n_free_days`,
#'   `n_workdays`, `msw_h`/`msf_h`/`msfsc_h` and their `_hms` forms, `sjl_h`,
#'   `sjla_h`, `cpd_s`/`cpd_min`/`cpd_h`).
#'
#' @seealso [study_summary()] for the NPCRA (activity-rhythm) analogue,
#'   [compute_sleep_metrics()], [compute_cpd_metrics()],
#'   [run_pipeline_native_batch()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- run_pipeline_native_batch("recordings/", tz = "America/Sao_Paulo")
#' study_sleep_metrics(results)
#'
#' # Feed straight into syncR::sync()
#' sync(zeit = study_sleep_metrics(results))
#' }
study_sleep_metrics <- function(results,
                                min_tib_h     = 5.0,
                                min_tib_eve_h = 3.0,
                                tz            = "UTC",
                                holidays      = NULL,
                                free_days     = NULL) {
  if (!is.list(results) || length(results) == 0L) {
    zeitr_abort(
      "{.arg results} must be a non-empty list of {.cls zeitr_result} objects.
       Did you run {.fn run_pipeline_batch} or {.fn run_pipeline_native_batch}?"
    )
  }

  list_names <- names(results)
  if (is.null(list_names)) list_names <- rep(NA_character_, length(results))

  rows <- lapply(seq_along(results), function(i) {
    result <- results[[i]]

    if (!inherits(result, "zeitr_result")) {
      pid <- list_names[i]
      zeitr_warn("Skipping {.val {pid}}: not a {.cls zeitr_result}.")
      return(NULL)
    }

    pid <- result$subject_id %||% list_names[i]

    sm_args <- list(x = result, min_tib_h = min_tib_h, tz = tz)
    if (!is.null(holidays))  sm_args$holidays  <- holidays
    if (!is.null(free_days)) sm_args$free_days <- free_days

    sm <- tryCatch(
      do.call(compute_sleep_metrics, sm_args),
      error = function(e) {
        zeitr_warn("compute_sleep_metrics() failed for {.val {pid}}: {conditionMessage(e)}")
        NULL
      }
    )
    if (is.null(sm) || length(sm) == 0L) sm <- .sm_na_row()

    cpd_args <- list(x = result, min_tib_h = min_tib_h,
                     min_tib_eve_h = min_tib_eve_h, tz = tz)
    if (!is.null(holidays))  cpd_args$holidays  <- holidays
    if (!is.null(free_days)) cpd_args$free_days <- free_days

    cpd <- tryCatch(
      do.call(compute_cpd_metrics, cpd_args),
      error = function(e) {
        zeitr_warn("compute_cpd_metrics() failed for {.val {pid}}: {conditionMessage(e)}")
        NULL
      }
    )
    if (is.null(cpd) || length(cpd) == 0L) cpd <- .cpd_na_row()

    tibble::as_tibble(c(list(participant_id = pid), sm, cpd))
  })

  rows <- Filter(Negate(is.null), rows)

  if (length(rows) == 0L) {
    zeitr_warn("No valid results found.")
    empty_sm  <- stats::setNames(rep(list(double()), length(.sm_all_names)), .sm_all_names)
    empty_cpd <- stats::setNames(
      lapply(.cpd_all_names, function(nm) if (nm %in% .cpd_char_names) character() else double()),
      .cpd_all_names
    )
    return(tibble::as_tibble(c(list(participant_id = character()), empty_sm, empty_cpd)))
  }

  do.call(rbind, rows)
}
