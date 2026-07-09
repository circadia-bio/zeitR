# Classifying free days: weekends and public holidays

Free-day classification sits at the heart of every chronotype metric
that zeitR computes. MSF (mid-sleep on free days), MSW (mid-sleep on
workdays), SJL (social jet lag), and CPD all depend on correctly
labelling which nights precede a free day and which do not.

By default, only **Saturdays and Sundays** are treated as free days. For
studies conducted in countries with public holidays this is usually
wrong: participants sleep differently on a Monday that follows a bank
holiday compared to a regular Monday. This vignette shows how to supply
holidays, the three supported input formats, and the mechanics of how
holidays flow through the pipeline.

------------------------------------------------------------------------

## Holiday formats

The `holidays` parameter accepted by
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md),
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md),
and
[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
understands three forms, which can be freely combined in the same
vector.

**Specific dates** — `Date` objects or `"YYYY-MM-DD"` character strings.
Use these for holidays whose calendar date changes from year to year
(Carnival, Good Friday, Easter Monday).

``` r

specific <- as.Date(c(
  "2019-03-04",  # Carnaval
  "2019-03-05",  # Carnaval
  "2019-04-19"   # Sexta-feira Santa
))
# character strings work too:
specific <- c("2019-03-04", "2019-03-05", "2019-04-19")
```

**Recurring holidays** — `"DD-MM"` character strings. Use these for
fixed-date holidays that fall on the same calendar day every year
(Christmas, New Year, Independence Day). The year is ignored; any date
whose day and month match is treated as a free day.

``` r

recurring <- c(
  "07-09",  # Independencia do Brasil
  "12-10",  # N. Sra. Aparecida
  "02-11",  # Finados
  "15-11",  # Proclamacao da Republica
  "20-11",  # Consciencia Negra (SP municipal)
  "25-12",  # Natal
  "01-01",  # Confraternizacao Universal
  "25-01",  # Aniversario de Sao Paulo (SP municipal)
  "21-04",  # Tiradentes
  "01-05"   # Dia do Trabalho
)
```

**Mixed** — both forms in one [`c()`](https://rdrr.io/r/base/c.html)
call. This is the typical setup for a national holiday calendar:
fixed-date holidays as `"DD-MM"` and variable-date ones as
`"YYYY-MM-DD"`.

``` r

sp_holidays <- c(
  # ── Recurring fixed-date holidays ────────────────────────────────────────
  "07-09", "12-10", "02-11", "15-11", "20-11",
  "25-12", "01-01", "25-01", "21-04", "01-05",
  # ── Year-specific variable-date holidays ─────────────────────────────────
  "2019-03-04", "2019-03-05",   # Carnaval 2019
  "2019-04-19"                  # Sexta-feira Santa 2019
)
```

------------------------------------------------------------------------

## The no-holidays warning

Calling
[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
or
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
without supplying any holidays emits a warning to alert you that only
weekends are being used.

``` r

FILE   <- system.file("extdata", "input1.txt", package = "zeitR")
TZ     <- "America/Sao_Paulo"
result <- run_pipeline_native(FILE, tz = TZ, quiet = TRUE)
#> ℹ Reading input1.txt ...
#> ✔ [input1] Done. 52 main night(s), 0 secondary episode(s).

# No holidays supplied — warning fires
cpd_wd_only <- compute_cpd_metrics(result$nights, tz = TZ)
#> Warning: compute_cpd_metrics(): no `holidays` supplied.
#> ℹ Only the days in `free_days` are treated as free days.
#> ℹ If you ran `run_pipeline_native()` with holidays, pass `holidays =
#>   result$holidays`.
#> ℹ Suppress with `options(zeitR.no_holidays_warn = FALSE)`.
```

This is intentional: a silent default that misclassifies bank holidays
as workdays is hard to catch. Once you have confirmed that weekends-only
is correct for your data (e.g. a study where participants had no common
holidays), suppress the warning globally for the session:

``` r

options(zeitR.no_holidays_warn = FALSE)
```

The rest of this vignette runs with the warning suppressed.

------------------------------------------------------------------------

## Passing holidays to the pipeline

Supply `holidays` when calling
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md).
The value is stored in `result$holidays` and does **not** affect the
sleep-period detection itself (the Vallim rules are day-of-week
agnostic); it is carried through so that downstream metric functions can
pick it up automatically.

``` r

sp_holidays <- c(
  # Recurring
  "07-09", "12-10", "02-11", "15-11", "20-11",
  "25-12", "01-01", "25-01", "21-04", "01-05",
  # Carnaval + Sexta-feira Santa 2019 (dates change each year)
  "2019-03-04", "2019-03-05", "2019-04-19"
)

result <- run_pipeline_native(
  FILE,
  tz       = TZ,
  holidays = sp_holidays,
  quiet    = TRUE
)
#> ℹ Reading input1.txt ...
#> ✔ [input1] Done. 52 main night(s), 0 secondary episode(s).

# holidays are stored on the result object
length(result$holidays)
#> [1] 13
head(result$holidays, 5)
#> [1] "07-09" "12-10" "02-11" "15-11" "20-11"
```

------------------------------------------------------------------------

## Auto-forwarding with S3 dispatch

[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
and
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
are S3 generics. Passing a `zeitr_result` directly — rather than
`result$nights` — picks up `result$holidays` automatically as the
default for the `holidays` argument.

``` r

# S3 dispatch: holidays forwarded automatically
cpd <- compute_cpd_metrics(result, tz = TZ)

# Equivalent explicit call
cpd_explicit <- compute_cpd_metrics(result$nights,
                                    tz       = TZ,
                                    holidays = result$holidays)

identical(cpd$msf_h, cpd_explicit$msf_h)
#> [1] TRUE
```

The same applies to
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md):

``` r

sm <- compute_sleep_metrics(result, tz = TZ)
```

If you need to override the stored holidays (e.g. to test a different
calendar), pass `holidays` explicitly and the explicit value wins:

``` r

cpd_alt <- compute_cpd_metrics(result,
                                tz       = TZ,
                                holidays = c("25-12", "01-01"))
```

------------------------------------------------------------------------

## Effect on metrics

Adding public holidays to the weekend pool changes how nights are split
between the free-day and workday groups, which in turn affects MSF, MSW,
and SJL.

``` r

cpd_we  <- compute_cpd_metrics(result$nights, tz = TZ)              # weekends only
cpd_all <- compute_cpd_metrics(result, tz = TZ)                     # + SP holidays

cat(sprintf("                  Weekends only   With holidays\n"))
#>                   Weekends only   With holidays
cat(sprintf("Free days         %13d   %d\n",
            cpd_we$n_free_days, cpd_all$n_free_days))
#> Free days                    15   15
cat(sprintf("Workdays          %13d   %d\n",
            cpd_we$n_workdays, cpd_all$n_workdays))
#> Workdays                     37   37
cat(sprintf("MSW               %13s   %s\n",
            cpd_we$msw_hms, cpd_all$msw_hms))
#> MSW                       03:54   03:54
cat(sprintf("MSF               %13s   %s\n",
            cpd_we$msf_hms, cpd_all$msf_hms))
#> MSF                       03:26   03:26
cat(sprintf("SJL               %10.1f min   %.1f min\n",
            cpd_we$sjl_min, cpd_all$sjl_min))
#> SJL                     28.1 min   28.1 min
cat(sprintf("CPD               %10.1f min   %.1f min\n",
            cpd_we$cpd_min, cpd_all$cpd_min))
#> CPD                     58.0 min   58.0 min
```

Reclassifying public holidays as free days moves nights from the workday
pool into the free-day pool. MSW and MSF shift because the composition
of each group changes; SJL reflects the true contrast between
unconstrained and constrained sleep across the whole study rather than
across weekends alone.

The
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
output shows the same split:

``` r

sm_we  <- compute_sleep_metrics(result$nights, tz = TZ)
sm_all <- compute_sleep_metrics(result, tz = TZ)

data.frame(
  group        = c("Overall", "Workday", "Free day"),
  n_we         = c(sm_we$n_overall,  sm_we$n_wd,  sm_we$n_fd),
  n_all        = c(sm_all$n_overall, sm_all$n_wd, sm_all$n_fd),
  tst_wd_we    = c(NA, round(sm_we$tst_h_wd, 2), NA),
  tst_fd_we    = c(NA, NA, round(sm_we$tst_h_fd, 2)),
  tst_wd_all   = c(NA, round(sm_all$tst_h_wd, 2), NA),
  tst_fd_all   = c(NA, NA, round(sm_all$tst_h_fd, 2))
) |>
  knitr::kable(
    col.names = c("Group",
                  "N (WE only)", "N (+ holidays)",
                  "TST wd", "TST fd",
                  "TST wd", "TST fd"),
    align = "lrrrrrr",
    caption = "Night counts and mean TST (h) by group."
  )
```

| Group    | N (WE only) | N (+ holidays) | TST wd | TST fd | TST wd | TST fd |
|:---------|------------:|---------------:|-------:|-------:|-------:|-------:|
| Overall  |          52 |             52 |     NA |     NA |     NA |     NA |
| Workday  |          37 |             37 |   7.42 |     NA |   7.42 |     NA |
| Free day |          15 |             15 |     NA |   7.73 |     NA |   7.73 |

Night counts and mean TST (h) by group. {.table}

------------------------------------------------------------------------

## Batch studies spanning multiple years

For a multi-year longitudinal study, use `"DD-MM"` recurring strings for
all fixed-date national holidays and add year-specific entries for each
year that falls within the recording window.

``` r

# Brazilian national holidays (recurring every year)
br_national <- c(
  "01-01",  # Confraternizacao Universal
  "21-04",  # Tiradentes
  "01-05",  # Dia do Trabalho
  "07-09",  # Independencia do Brasil
  "12-10",  # N. Sra. Aparecida
  "02-11",  # Finados
  "15-11",  # Proclamacao da Republica
  "25-12"   # Natal
)

# Variable-date holidays for each study year
variable <- c(
  # Carnaval
  "2018-02-12", "2018-02-13",
  "2019-03-04", "2019-03-05",
  "2020-02-24", "2020-02-25",
  # Sexta-feira Santa
  "2018-03-30",
  "2019-04-19",
  "2020-04-10"
)

study_holidays <- c(br_national, variable)

results <- run_pipeline_native_batch(
  "recordings/",
  tz       = TZ,
  holidays = study_holidays,
  quiet    = TRUE
)

# compute_cpd_metrics(result, tz = TZ) automatically uses each result's holidays
cpd_list <- lapply(results, compute_cpd_metrics, tz = TZ)
```

Because the `"DD-MM"` entries match any year, you only need to update
`variable` when a new study year is added — the recurring national
holidays require no changes.

------------------------------------------------------------------------

## Non-standard schedules with `free_days`

Not every participant works Monday-to-Friday. The `free_days` parameter
lets you specify which days of the week are unconditionally treated as
free days, replacing the default Saturday + Sunday.

Accepted inputs are English day names (case-insensitive) or ISO integers
(1 = Monday, …, 6 = Saturday, 7 = Sunday).

``` r

# Compressed work week: Monday off, works Saturday
compute_cpd_metrics(result$nights,
                    tz        = TZ,
                    holidays  = sp_holidays,
                    free_days = c("Sunday", "Monday"))

# Middle-East schedule: Friday + Saturday as weekend
compute_cpd_metrics(result$nights,
                    tz        = TZ,
                    free_days = c("Friday", "Saturday"))

# ISO integers work too (5 = Friday, 6 = Saturday)
compute_cpd_metrics(result$nights,
                    tz        = TZ,
                    free_days = c(5L, 6L))
```

Passing `free_days` to
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
stores it in `result$free_days` so the S3 dispatch picks it up
automatically:

``` r

result_mon <- run_pipeline_native(
  FILE,
  tz        = TZ,
  holidays  = sp_holidays,
  free_days = c("Saturday", "Sunday", "Monday"),
  quiet     = TRUE
)
#> ℹ Reading input1.txt ...
#> ✔ [input1] Done. 52 main night(s), 0 secondary episode(s).

# S3 dispatch forwards both holidays and free_days
cpd_mon <- compute_cpd_metrics(result_mon, tz = TZ)

cat(sprintf("Free days (Sat+Sun):        %d\n", cpd_all$n_free_days))
#> Free days (Sat+Sun):        15
cat(sprintf("Free days (Sat+Sun+Mon):    %d\n", cpd_mon$n_free_days))
#> Free days (Sat+Sun+Mon):    23
cat(sprintf("SJL (Sat+Sun):        %.1f min\n", cpd_all$sjl_min))
#> SJL (Sat+Sun):        28.1 min
cat(sprintf("SJL (Sat+Sun+Mon):    %.1f min\n", cpd_mon$sjl_min))
#> SJL (Sat+Sun+Mon):    27.2 min
```

`free_days` and `holidays` are additive: a night is a free-day night if
the get-up date falls on any day in `free_days` **or** matches any entry
in `holidays`.

------------------------------------------------------------------------

## Summary

| Format | Example | When to use |
|----|----|----|
| `Date` / `"YYYY-MM-DD"` | `as.Date("2019-03-04")` | Holidays whose date changes year to year |
| `"DD-MM"` | `"25-12"` | Fixed-date holidays that repeat annually |
| Mixed | `c("25-12", "2019-03-04")` | Typical national calendar (most common) |

| Parameter | Default | Purpose |
|----|----|----|
| `free_days` | `c("Saturday", "Sunday")` | Days of the week that are always free |
| `holidays` | `NULL` | Specific or recurring dates that are also free |

Passing `holidays` to
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
stores them on the result object; calling
`compute_cpd_metrics(result, tz = ...)` or
`compute_sleep_metrics(result, tz = ...)` with a `zeitr_result` picks
them up automatically. To silence the reminder that fires when no
holidays are supplied, use `options(zeitR.no_holidays_warn = FALSE)`.

------------------------------------------------------------------------

## See also

- [`vignette("vallim-pipeline")`](https://zeitr.circadia-lab.uk/articles/vallim-pipeline.md)
  — the Vallim classification rule set
- [`?compute_cpd_metrics`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  — MSF, MSW, SJL, MSFsc, and CPD reference
- [`?compute_sleep_metrics`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  — per-day-type sleep summary metrics
- [`?run_pipeline_native`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
  — full pipeline parameter reference
