# dev/sleep_diag.R
# Why does Python remove the long off-wrist period at valid-space ~[60350,60598)
# while R keeps it?  Candidates: sleep filter (sleep_prop>=0.4) or forbidden zone.
# Run:  source("dev/sleep_diag.R")  then paste output back.

devtools::load_all(quiet = TRUE)
options(zeitR.debug_offwrist = TRUE)

rec  <- read_acttrust(system.file("extdata", "input1.txt", package = "zeitR"), tz = "UTC")
prep <- detect_offwrist_bimodal(prepare_actigraphy(rec))
dbg  <- get(".zeitR_dbg", envir = .GlobalEnv)

es   <- dbg$estimated_sleep            # R, valid-space, 0=sleep 1=wake
fz   <- dbg$forbidden_zone             # R, valid-space
aip  <- dbg$after_initial_periods      # pre-sleep-filter periods + sleep_prop
sfp  <- dbg$sleep_filtered_post        # post-forbidden survivors
nv   <- length(es)
cat(sprintf("valid-space length: %d   max_offwrist_sleep_prop=%.2f  long_len=%d\n",
            nv, dbg$max_offwrist_sleep_prop, dbg$long_offwrist_length))

# Python intermediates (valid-space). estimated_sleep column is stored INVERTED
# (file = 1 - actual). refined_offwrist is Python's final refiner output.
pin <- read.csv(system.file("extdata", "python_intermediate.csv", package = "zeitR"))
cat(sprintf("python_intermediate rows: %d  (matches valid-space: %s)\n",
            nrow(pin), nrow(pin) == nv))
py_es  <- 1L - as.integer(round(pin$estimated_sleep))   # de-invert -> 0=sleep 1=wake
py_ref <- as.integer(round(pin$refined_offwrist))        # 0=offwrist 1=onwrist (Python)

# ── The period in question (valid-space, overlapping 60350..60598) ────────────
cat("\n-- after_initial_periods overlapping valid 60300..60650 (R) --\n")
sel <- aip[aip$start < 60650 & aip$end > 60300, ]
for (k in seq_len(nrow(sel))) {
  s <- sel$start[k]; e <- sel$end[k]
  fov <- any(fz[(s + 1L):min(e, nv)] == 1L)
  cat(sprintf("  [%d,%d) len=%d sleep_prop=%.3f passes_sleep=%s forbidden_overlap=%s in_final=%s\n",
              s, e, e - s, sel$sleep_prop[k],
              sel$sleep_prop[k] < dbg$max_offwrist_sleep_prop, fov,
              any(sfp$start == s & sfp$end == e)))
}

# ── estimated_sleep agreement, global and local ───────────────────────────────
if (nrow(pin) == nv) {
  cat(sprintf("\nglobal estimated_sleep agreement: %.5f  (mismatches=%d)\n",
              mean(es == py_es), sum(es != py_es)))
  mm <- which(es != py_es)
  if (length(mm)) {
    mmr <- rle(diff(mm) == 1L)
    cat("mismatch runs (valid-space epoch ranges, 0-based):\n")
    # report contiguous mismatch blocks
    starts <- mm[c(TRUE, diff(mm) != 1L)]
    ends   <- mm[c(diff(mm) != 1L, TRUE)]
    blk <- data.frame(start0 = starts - 1L, end0 = ends - 1L, len = ends - starts + 1L)
    print(utils::head(blk[order(-blk$len), ], 15), row.names = FALSE)
  }

  cat("\n-- local view, valid epochs 60345..60600 (sleep: 0=sleep,1=wake) --\n")
  rng <- 60345:60600
  tab <- data.frame(
    epoch  = rng,
    R_sleep  = es[rng + 1L],
    PY_sleep = py_es[rng + 1L],
    R_forbid = fz[rng + 1L],
    PY_final_onwrist = py_ref[rng + 1L]
  )
  # only print rows where something is interesting (disagreement or boundaries)
  show <- which(tab$R_sleep != tab$PY_sleep | tab$R_forbid == 1L |
                c(0, diff(tab$PY_final_onwrist)) != 0 | c(0, diff(tab$R_sleep)) != 0)
  if (length(show) == 0) show <- seq_len(min(40, nrow(tab)))
  print(tab[unique(pmax(1, c(show, show + 1))), ], row.names = FALSE)
}

options(zeitR.debug_offwrist = FALSE)
