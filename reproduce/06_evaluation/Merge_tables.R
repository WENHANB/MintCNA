###############################################################
# Merge tables across all 9 simulation settings
# Adds coverage and noise columns parsed from folder names
# Default: coverage = 0.1, noise = 0 (when not in folder name)
###############################################################

library(data.table)

TAB_ROOT <- "<MINTCNA_ROOT>/Simulation/results/tables"
OUT_DIR  <- file.path(TAB_ROOT, "merged")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

settings <- c(
  "c3_N100_100kb_l1e4",
  "c3_N100_100kb_l1e4_c0.05",
  "c3_N100_100kb_l1e4_c0.05_n0.05",
  "c3_N100_100kb_l1e4_c0.05_n0.1",
  "c3_N100_100kb_l1e4_c0.2",
  "c3_N100_100kb_l1e4_c0.2_n0.05",
  "c3_N100_100kb_l1e4_c0.2_n0.1",
  "c3_N100_100kb_l1e4_n0.05",
  "c3_N100_100kb_l1e4_n0.1"
)

# Parse coverage and noise from folder name
parse_params <- function(s) {
  # Extract coverage: _c0.05 or _c0.2; default 0.1
  cov_match <- regmatches(s, regexpr("_c(0\\.[0-9]+)", s))
  coverage <- if (length(cov_match) == 1) as.numeric(sub("^_c", "", cov_match)) else 0.1

  # Extract noise: _n0.05 or _n0.1; default 0
  noise_match <- regmatches(s, regexpr("_n(0\\.[0-9]+)", s))
  noise <- if (length(noise_match) == 1) as.numeric(sub("^_n", "", noise_match)) else 0

  list(coverage = coverage, noise = noise)
}

table_names <- c("summary_del", "summary_dup", "summary_mss", "summary_one", "bb_summary", "bb_percell", "RMSE_percell")

for (tbl in table_names) {
  merged <- rbindlist(lapply(settings, function(s) {
    params <- parse_params(s)
    dt <- fread(file.path(TAB_ROOT, s, paste0(tbl, ".csv")))
    dt[, coverage := params$coverage]
    dt[, noise   := params$noise]
    dt
  }), use.names = TRUE, fill = TRUE)

  # Move coverage and noise to front
  front_cols <- c("coverage", "noise")
  other_cols <- setdiff(names(merged), front_cols)
  setcolorder(merged, c(front_cols, other_cols))

  out_path <- file.path(OUT_DIR, paste0(tbl, ".csv"))
  fwrite(merged, out_path)
  message("Saved: ", out_path)
}

message("Done. All merged tables in: ", OUT_DIR)
