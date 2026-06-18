rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(GEOquery)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_raw_import_utils.R"))

cohort_id <- get_arg("--scrna-cohort-id")
database_dir <- get_arg("--database-dir")
raw_dir <- get_arg("--raw-dir")
geo_filter_regex <- get_arg("--geo-filter-regex", required = FALSE, default = "none")
outdir <- get_arg("--outdir")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

ensure_dir(outdir)
ensure_dir(database_dir)
ensure_dir(raw_dir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_geo_supp_download")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_geo_supp_download")
log_message("scRNA cohort ID:", cohort_id)
log_message("Database dir:", database_dir)
log_message("Raw dir:", raw_dir)

if (!grepl("^GSE[0-9]+$", cohort_id, ignore.case = TRUE)) {
  save_note_and_stop("scRNA_cohort_id must be a GEO GSE accession for GEOquery download.")
}

filter_value <- if (is_null_path(geo_filter_regex)) NULL else geo_filter_regex

download_base_dir <- dirname(normalizePath(raw_dir, winslash = "/", mustWork = FALSE))
ensure_dir(download_base_dir)

res <- tryCatch(
  GEOquery::getGEOSuppFiles(
    GEO = cohort_id,
    makeDirectory = TRUE,
    baseDir = download_base_dir,
    fetch_files = TRUE,
    filter_regex = filter_value
  ),
  error = function(e) save_note_and_stop(paste0("GEO supplementary download failed: ", conditionMessage(e)))
)

download_dir <- file.path(download_base_dir, cohort_id)
if (!dir.exists(download_dir)) save_note_and_stop(paste0("Expected GEO download directory not found: ", download_dir))

if (normalizePath(download_dir, winslash = "/", mustWork = FALSE) != normalizePath(raw_dir, winslash = "/", mustWork = FALSE)) {
  ensure_dir(raw_dir)
  downloaded_files <- list.files(download_dir, recursive = TRUE, full.names = TRUE)
  downloaded_files <- downloaded_files[file.exists(downloaded_files) & !dir.exists(downloaded_files)]
  for (src in downloaded_files) {
    rel <- substring(normalizePath(src, winslash = "/", mustWork = FALSE), nchar(normalizePath(download_dir, winslash = "/", mustWork = FALSE)) + 2)
    dst <- file.path(raw_dir, rel)
    ensure_dir(dirname(dst))
    if (!file.exists(dst) || file.info(dst)$size == 0) file.copy(src, dst, overwrite = TRUE)
  }
}

cohort_dir <- raw_dir
files <- list.files(cohort_dir, recursive = TRUE, full.names = TRUE)
files <- files[file.exists(files) & !dir.exists(files)]
if (length(files) == 0) save_note_and_stop("No supplementary files were downloaded or found.")

manifest <- data.frame(
  cohort_id = cohort_id,
  file_name = basename(files),
  file_path = normalizePath(files, winslash = "/", mustWork = FALSE),
  file_size = as.numeric(file.info(files)$size),
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$file_name), , drop = FALSE]
write.csv(manifest, file.path(outdir, "01.geo_supp_files.csv"), row.names = FALSE)

log_message("Downloaded/found files:", nrow(manifest))
log_message("Script completed: scrna_geo_supp_download")

write_summary(c(
  "scRNA GEO supplementary download completed.",
  paste0("scRNA cohort ID: ", cohort_id),
  paste0("Downloaded/found files: ", nrow(manifest)),
  paste0("Raw directory: ", normalizePath(cohort_dir, winslash = "/", mustWork = FALSE)),
  "Supplementary files were downloaded by GEOquery and recorded in 01.geo_supp_files.csv.",
  "",
  "限制与说明",
  "This step requires network access to GEO. File roles are resolved in the next step."
))
