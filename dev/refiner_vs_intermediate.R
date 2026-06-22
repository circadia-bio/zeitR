# dev/refiner_vs_intermediate.R
# Is the remaining gap inside the refiner, or in a POST-refiner pipeline step?
# Compares R refiner output vs Python refined_offwrist (valid space), and
# Python refined_offwrist vs python_output state (the validation target).
# Run:  source("dev/refiner_vs_intermediate.R")  then paste output.

devtools::load_all(quiet = TRUE)

rec  <- read_acttrust(system.file("extdata", "input1.txt", package = "zeitR"), tz = "UTC")
prep <- detect_offwrist_bimodal(prepare_actigraphy(rec))

valid_temp <- prep$int_temp > 0
nv <- sum(valid_temp); nf <- nrow(prep)

pin <- read.csv(system.file("extdata", "python_intermediate.csv", package = "zeitR"))
pout<- read.csv(system.file("extdata", "python_output.csv",       package = "zeitR"))

# R refiner offwrist in valid space: state==4 over valid epochs (1 = offwrist)
r_off_valid  <- as.integer(prep$state[valid_temp] == 4)
# Python refiner offwrist in valid space: refined_offwrist==0 (Python 0=offwrist)
py_off_valid <- as.integer(round(pin$refined_offwrist) == 0)

cat(sprintf("valid n: R=%d  py_intermediate=%d  py_output=%d  full=%d\n",
            length(r_off_valid), nrow(pin), nrow(pout), nf))

cat("\n[A] R refiner  vs  Python refined_offwrist  (valid space)\n")
if (length(r_off_valid) == length(py_off_valid)) {
  cat(sprintf("    correlation=%.6f  mismatches=%d\n",
              suppressWarnings(cor(r_off_valid, py_off_valid)),
              sum(r_off_valid != py_off_valid)))
  print(table(R = r_off_valid, PY_refined = py_off_valid))
} else cat("    length mismatch, cannot compare directly\n")

# Python refined_offwrist mapped to FULL space (offwrist where refined==0),
# invalid-temp epochs are forced offwrist by the pipeline.
py_ref_full <- integer(nf)
py_ref_full[!valid_temp] <- 1L
py_ref_full[valid_temp]  <- py_off_valid
py_out_full <- as.integer(pout$state == 4)

cat("\n[B] Python refined_offwrist(full)  vs  python_output state==4  (validation target)\n")
cat(sprintf("    correlation=%.6f  mismatches=%d\n",
            suppressWarnings(cor(py_ref_full, py_out_full)),
            sum(py_ref_full != py_out_full)))
print(table(PY_refined_full = py_ref_full, PY_output = py_out_full))

# Where do refined_offwrist and python_output disagree? (post-refiner pipeline edits)
runs <- function(flag) {
  r <- rle(flag); e <- cumsum(r$lengths); s <- e - r$lengths
  k <- r$values; data.frame(start = s[k], end = e[k], len = r$lengths[k])
}
d_post <- runs(py_ref_full != py_out_full)
cat(sprintf("\nPost-refiner edit runs (Python refined != Python output): %d (total %d epochs)\n",
            nrow(d_post), sum(d_post$len)))
print(d_post[order(-d_post$len), ][seq_len(min(15, nrow(d_post))), ], row.names = FALSE)

# And R final vs python_output for reference
r_out_full <- as.integer(prep$state == 4)
cat(sprintf("\n[C] R output vs python_output: correlation=%.6f mismatches=%d\n",
            cor(r_out_full, py_out_full), sum(r_out_full != py_out_full)))
