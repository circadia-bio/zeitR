# dev/fp_attrib.R  — attribute the current false positives to VP periods.
# Run: source("dev/fp_attrib.R")  then paste output.

devtools::load_all(quiet = TRUE)
options(zeitR.debug_offwrist = TRUE)
rec  <- read_acttrust(system.file("extdata","input1.txt",package="zeitR"), tz="UTC")
prep <- detect_offwrist_bimodal(prepare_actigraphy(rec))
dbg  <- get(".zeitR_dbg", envir = .GlobalEnv)

valid_temp <- prep$int_temp > 0
r_ow  <- as.integer(prep$state == 4)
py    <- read.csv(system.file("extdata","python_output.csv",package="zeitR"))
py_ow <- as.integer(py$state == 4)
cat(sprintf("Correlation: %.6f\n", cor(r_ow, py_ow)))

runs <- function(flag){r<-rle(flag);e<-cumsum(r$lengths);s<-e-r$lengths;k<-r$values;data.frame(start=s[k],end=e[k],len=r$lengths[k])}
fp <- runs(r_ow==1 & py_ow==0)
fn <- runs(r_ow==0 & py_ow==1)
cat(sprintf("\nFP runs: %d (total %d).  FN runs: %d (total %d)\n",
            nrow(fp), sum(fp$len), nrow(fn), sum(fn$len)))
cat("Top FP runs:\n"); print(fp[order(-fp$len),][seq_len(min(15,nrow(fp))),], row.names=FALSE)

vp <- dbg$vp_df
cat(sprintf("\nR vp_df rows (kept VP periods): %d\n", if(is.null(vp)) 0 else nrow(vp)))
if(!is.null(vp) && nrow(vp)){
  print(vp, row.names=FALSE)
  # For each FP run (valid-space ~ full-space since temp<=0 forced offwrist on both),
  # mark if it falls inside a VP period.
  in_vp <- function(s,e){ any(vp$start < e & vp$end > s) }
  fp$in_vp <- vapply(seq_len(nrow(fp)), function(i) in_vp(fp$start[i], fp$end[i]), logical(1))
  cat("\nFP runs that fall inside a VP period:\n")
  print(fp[fp$in_vp, ][order(-fp$len[fp$in_vp]),], row.names=FALSE)
  cat(sprintf("\nFP epochs inside VP periods: %d of %d\n",
              sum(fp$len[fp$in_vp]), sum(fp$len)))
}
options(zeitR.debug_offwrist = FALSE)
