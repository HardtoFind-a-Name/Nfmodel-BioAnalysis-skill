rm(list = ls())
gc()

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_raw_import_utils.R"))

raw_source_manifest <- get_arg("--raw-source-manifest")
outdir <- get_arg("--outdir")
unpack_root <- get_arg("--unpack-root")
raw_10x_prefer <- get_arg("--raw-10x-prefer", required = FALSE, default = "filtered_feature_bc_matrix")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

ensure_dir(outdir)
ensure_dir(unpack_root)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_unpack_10x")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

has_10x_files <- function(dir_path) {
  if (!dir.exists(dir_path)) return(FALSE)
  files <- basename(list.files(dir_path, all.files = FALSE, full.names = FALSE))
  has_matrix <- any(grepl("^matrix\\.mtx(\\.gz)?$", files, ignore.case = TRUE))
  has_barcodes <- any(grepl("^barcodes\\.tsv(\\.gz)?$", files, ignore.case = TRUE))
  has_features <- any(grepl("^(features|genes)\\.tsv(\\.gz)?$", files, ignore.case = TRUE))
  has_matrix && has_barcodes && has_features
}

find_10x_dirs <- function(root) {
  if (!dir.exists(root)) return(character())
  matrix_files <- list.files(root, pattern = "^matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  dirs <- unique(dirname(matrix_files))
  dirs[vapply(dirs, has_10x_files, logical(1))]
}

choose_10x_dir <- function(dirs, prefer) {
  if (length(dirs) == 0) return(NA_character_)
  if (!is_null_path(prefer)) {
    preferred <- dirs[grepl(prefer, dirs, ignore.case = TRUE, perl = TRUE)]
    if (length(preferred) >= 1) return(normalizePath(sort(preferred)[1], winslash = "/", mustWork = TRUE))
  }
  normalizePath(sort(dirs)[1], winslash = "/", mustWork = TRUE)
}

log_message("Script started: scrna_unpack_10x")
if (!file.exists(raw_source_manifest)) save_note_and_stop(paste0("Raw source manifest not found: ", raw_source_manifest))
manifest <- read.csv(raw_source_manifest, stringsAsFactors = FALSE, check.names = FALSE)
raw_type <- manifest$raw_source_type[manifest$role == "raw_source"][1]
source_path <- manifest$source_path[manifest$role == "raw_source"][1]
if (is.na(raw_type) || is.na(source_path) || !file.exists(source_path)) save_note_and_stop("Raw source manifest is invalid.")

if (!raw_type %in% c("10x_dir", "10x_archive")) {
  log_message("Raw source type is", raw_type, "; 10X unpack step passes manifest through unchanged")
  write.csv(manifest, file.path(outdir, "01.raw_source_ready.csv"), row.names = FALSE)
  write_summary(c(
    "scRNA 10X unpack skipped.",
    paste0("Raw source type: ", raw_type),
    "The manifest was passed through unchanged for mainline construction.",
    "",
    "限制与说明",
    "Only 10x_archive sources are unpacked here; 10x_dir sources are normalized to a matrix directory."
  ))
  quit(save = "no", status = 0)
}

matrix_dir <- ""
if (raw_type == "10x_dir") {
  dirs <- if (has_10x_files(source_path)) source_path else find_10x_dirs(source_path)
  if (length(dirs) == 0) save_note_and_stop(paste0("No 10X matrix directory found under: ", source_path))
  matrix_dir <- choose_10x_dir(dirs, raw_10x_prefer)
} else if (raw_type == "10x_archive") {
  label <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(source_path)))
  dest <- file.path(unpack_root, safe_file_label(label))
  ensure_dir(dest)
  lower <- tolower(source_path)
  log_message("Unpack 10X archive:", source_path, "->", dest)
  if (grepl("\\.zip$", lower)) {
    utils::unzip(source_path, exdir = dest)
  } else if (grepl("\\.(tar|tar\\.gz|tgz)$", lower)) {
    utils::untar(source_path, exdir = dest)
  } else {
    save_note_and_stop(paste0("Unsupported 10X archive extension: ", source_path))
  }
  dirs <- find_10x_dirs(dest)
  if (length(dirs) == 0) save_note_and_stop("Archive unpacked but no 10X matrix directory was found.")
  matrix_dir <- choose_10x_dir(dirs, raw_10x_prefer)
}

manifest$raw_source_type[manifest$role == "raw_source"] <- "10x_dir"
manifest$source_path[manifest$role == "raw_source"] <- matrix_dir
manifest$matrix_dir[manifest$role == "raw_source"] <- matrix_dir
manifest$file_name[manifest$role == "raw_source"] <- basename(matrix_dir)
manifest$note[manifest$role == "raw_source"] <- "10X matrix directory ready for Read10X"
write.csv(manifest, file.path(outdir, "01.raw_source_ready.csv"), row.names = FALSE)

log_message("10X matrix dir:", matrix_dir)
log_message("Script completed: scrna_unpack_10x")
write_summary(c(
  "scRNA 10X unpack completed.",
  paste0("10X matrix dir: ", matrix_dir),
  "The raw source manifest now points to a Read10X-compatible matrix directory.",
  "",
  "限制与说明",
  "If multiple 10X directories exist, the preferred pattern is used first and then lexical order."
))
