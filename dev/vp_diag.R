# dev/vp_diag.R
# Inspect R's valley-peak output (vp_df) around the bogus [60066,60296] FP and
# the FN regions, to see WHY R forms/keeps them. Run: source("dev/vp_diag.R")

devtools::load_all(quiet = TRUE)
options(zeitR.debug_offwrist = TRUE)
rec  <- read_acttrust(system.file("extdata", "input1.txt", package = "zeitR"), tz = "UTC")
prep <- detect_offwrist_bimodal(prepare_actigraphy(rec))
dbg  <- get(".zeitR_dbg", envir = .GlobalEnv)

vp <- dbg$vp_df
cat(sprintf("R vp_df rows (kept VP periods): %d\n", if (is.null(vp)) 0 else nrow(vp)))
if (!is.null(vp) && nrow(vp)) {
  cat("\n-- VP periods overlapping the FP/FN regions --\n")
  show <- vp[(vp$start < 60400 & vp$end > 60000) |
             (vp$start < 53300 & vp$end > 53000) |
             (vp$start < 39700 & vp$end > 39500) |
             (vp$start < 43000 & vp$end > 42900), ]
  print(show, row.names = FALSE)
  cat("\n-- full vp_df (all kept VP periods) --\n")
  print(vp, row.names = FALSE)
}
options(zeitR.debug_offwrist = FALSE)
