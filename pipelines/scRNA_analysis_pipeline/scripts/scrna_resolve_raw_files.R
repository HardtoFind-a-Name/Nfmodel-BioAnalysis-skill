rm(list = ls())
gc()

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_raw_import_utils.R"))

manifest_file <- get_arg("--manifest")
raw_dir <- get_arg("--raw-dir")
raw_umi_file <- get_arg("--raw-umi-file", required = FALSE, default = "none")
annotation_file <- get_arg("--annotation-file", required = FALSE, default = "none")
raw_umi_pattern <- get_arg("--raw-umi-pattern", required = FALSE, default = "raw_UMI|raw.*matrix|UMI_matrix")
annotation_pattern <- get_arg("--annotation-pattern", required = FALSE, default = "annotation|cell_annotation|metadata")
outdir <- get_arg("--outdir")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_resolve_raw_files")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_resolve_raw_files")
if (!file.exists(manifest_file)) save_note_and_stop(paste0("Manifest not found: ", manifest_file))
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"file_path" %in% colnames(manifest)) save_note_and_stop("Manifest must contain file_path column.")
files <- manifest$file_path[file.exists(manifest$file_path)]
if (length(files) == 0) save_note_and_stop("No existing files in manifest.")

resolve_one <- function(explicit_value, pattern, role) {
  if (!is_null_path(explicit_value)) {
    candidates <- c(explicit_value, file.path(raw_dir, explicit_value))
    hit <- candidates[file.exists(candidates)]
    if (length(hit) == 0) {
      hit <- files[basename(files) == basename(explicit_value)]
    }
    if (length(hit) != 1) {
      save_note_and_stop(paste0("Could not resolve explicit ", role, " file: ", explicit_value))
    }
    return(normalizePath(hit[1], winslash = "/", mustWork = TRUE))
  }

  hits <- safe_basename_match(files, pattern)
  if (length(hits) == 0) {
    save_note_and_stop(paste0("No ", role, " file matched pattern: ", pattern,
                              ". Available files: ", paste(basename(files), collapse = ", ")))
  }
  if (length(hits) > 1) {
    save_note_and_stop(paste0("Multiple ", role, " candidates matched pattern: ",
                              paste(basename(hits), collapse = ", "),
                              ". Specify the file explicitly by CLI."))
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

raw_umi_path <- resolve_one(raw_umi_file, raw_umi_pattern, "raw UMI")
annotation_path <- resolve_one(annotation_file, annotation_pattern, "annotation")

resolved <- data.frame(
  role = c("raw_umi", "annotation"),
  file_path = c(raw_umi_path, annotation_path),
  file_name = basename(c(raw_umi_path, annotation_path)),
  stringsAsFactors = FALSE
)
write.csv(resolved, file.path(outdir, "01.resolved_raw_files.csv"), row.names = FALSE)

log_message("Raw UMI:", raw_umi_path)
log_message("Annotation:", annotation_path)
log_message("Script completed: scrna_resolve_raw_files")

write_summary(c(
  "scRNA raw file resolution completed.",
  paste0("Raw UMI file: ", basename(raw_umi_path)),
  paste0("Annotation file: ", basename(annotation_path)),
  "File roles were resolved from explicit CLI values or filename patterns.",
  "",
  "限制与说明",
  "If filename pattern matching is ambiguous, specify raw UMI and annotation files explicitly."
))
