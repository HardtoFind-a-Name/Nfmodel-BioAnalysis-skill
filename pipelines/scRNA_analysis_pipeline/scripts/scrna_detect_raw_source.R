rm(list = ls())
gc()

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_raw_import_utils.R"))

resolved_file <- get_arg("--resolved-file")
outdir <- get_arg("--outdir")
raw_source_type <- tolower(get_arg("--raw-source-type", required = FALSE, default = "auto"))
raw_10x_pattern <- get_arg("--raw-10x-pattern", required = FALSE, default = "filtered_feature_bc_matrix|raw_feature_bc_matrix|matrix.mtx|10x")
raw_10x_prefer <- get_arg("--raw-10x-prefer", required = FALSE, default = "filtered_feature_bc_matrix")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_detect_raw_source")
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
    if (length(preferred) == 1) return(normalizePath(preferred, winslash = "/", mustWork = TRUE))
    if (length(preferred) > 1) return(normalizePath(sort(preferred)[1], winslash = "/", mustWork = TRUE))
  }
  normalizePath(sort(dirs)[1], winslash = "/", mustWork = TRUE)
}

archive_entries <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.zip$", lower)) {
    z <- tryCatch(utils::unzip(path, list = TRUE), error = function(e) NULL)
    if (is.null(z)) return(character())
    return(z$Name)
  }
  if (grepl("\\.(tar|tar\\.gz|tgz)$", lower)) {
    return(tryCatch(utils::untar(path, list = TRUE), error = function(e) character()))
  }
  character()
}

archive_has_10x <- function(path) {
  entries <- archive_entries(path)
  if (length(entries) == 0) return(FALSE)
  dirs <- unique(dirname(entries[grepl("(^|/)matrix\\.mtx(\\.gz)?$", entries, ignore.case = TRUE)]))
  for (d in dirs) {
    members <- basename(entries[dirname(entries) == d])
    has_matrix <- any(grepl("^matrix\\.mtx(\\.gz)?$", members, ignore.case = TRUE))
    has_barcodes <- any(grepl("^barcodes\\.tsv(\\.gz)?$", members, ignore.case = TRUE))
    has_features <- any(grepl("^(features|genes)\\.tsv(\\.gz)?$", members, ignore.case = TRUE))
    if (has_matrix && has_barcodes && has_features) return(TRUE)
  }
  FALSE
}

log_message("Script started: scrna_detect_raw_source")
if (!file.exists(resolved_file)) save_note_and_stop(paste0("Resolved file table not found: ", resolved_file))
resolved <- read.csv(resolved_file, stringsAsFactors = FALSE, check.names = FALSE)
raw_source_path <- resolved$file_path[resolved$role == "raw_umi"][1]
annotation_path <- resolved$file_path[resolved$role == "annotation"][1]
if (is.na(raw_source_path) || !file.exists(raw_source_path)) save_note_and_stop("Resolved raw source file/path is missing.")
if (is.na(annotation_path) || !file.exists(annotation_path)) save_note_and_stop("Resolved annotation file is missing.")
raw_source_path <- normalizePath(raw_source_path, winslash = "/", mustWork = TRUE)
annotation_path <- normalizePath(annotation_path, winslash = "/", mustWork = TRUE)

source_type <- raw_source_type
source_path <- raw_source_path
matrix_dir <- ""
note <- ""

if (source_type == "auto") {
  if (dir.exists(raw_source_path)) {
    dirs <- find_10x_dirs(raw_source_path)
    if (length(dirs) > 0) {
      source_type <- "10x_dir"
      matrix_dir <- choose_10x_dir(dirs, raw_10x_prefer)
      source_path <- matrix_dir
    } else {
      save_note_and_stop(paste0("Raw source directory is not a 10X matrix directory: ", raw_source_path))
    }
  } else {
    lower <- tolower(raw_source_path)
    parent_10x <- if (grepl("matrix\\.mtx(\\.gz)?$", basename(lower))) has_10x_files(dirname(raw_source_path)) else FALSE
    if (parent_10x) {
      source_type <- "10x_dir"
      matrix_dir <- normalizePath(dirname(raw_source_path), winslash = "/", mustWork = TRUE)
      source_path <- matrix_dir
    } else if (grepl("\\.h5$", lower)) {
      source_type <- "10x_h5"
    } else if (archive_has_10x(raw_source_path)) {
      source_type <- "10x_archive"
    } else if (grepl("\\.rds(\\.gz)?$", lower) || grepl("raw_UMI|raw.*matrix|UMI_matrix", basename(raw_source_path), ignore.case = TRUE, perl = TRUE)) {
      source_type <- "umi_rds"
    } else {
      save_note_and_stop(paste0("Could not auto-detect raw source type for: ", raw_source_path))
    }
  }
} else if (source_type == "10x_dir") {
  if (file.exists(raw_source_path) && !dir.exists(raw_source_path) && grepl("matrix\\.mtx(\\.gz)?$", basename(raw_source_path), ignore.case = TRUE)) {
    raw_source_path <- dirname(raw_source_path)
  }
  dirs <- if (has_10x_files(raw_source_path)) raw_source_path else find_10x_dirs(raw_source_path)
  if (length(dirs) == 0) save_note_and_stop("raw_source_type=10x_dir but no 10X matrix directory was found.")
  matrix_dir <- choose_10x_dir(dirs, raw_10x_prefer)
  source_path <- matrix_dir
} else if (source_type == "10x_archive") {
  if (!archive_has_10x(raw_source_path)) save_note_and_stop("raw_source_type=10x_archive but archive does not contain a 10X matrix structure.")
} else if (source_type == "umi_rds") {
  note <- "UMI RDS will be validated after recursive gzip unpacking."
} else if (source_type == "10x_h5") {
  if (!grepl("\\.h5$", tolower(raw_source_path))) save_note_and_stop("raw_source_type=10x_h5 requires an .h5 file.")
} else {
  save_note_and_stop(paste0("Unsupported raw_source_type: ", source_type))
}

manifest <- data.frame(
  role = c("raw_source", "annotation"),
  raw_source_type = c(source_type, source_type),
  source_path = c(source_path, annotation_path),
  original_path = c(raw_source_path, annotation_path),
  matrix_dir = c(matrix_dir, ""),
  file_name = basename(c(source_path, annotation_path)),
  note = c(note, ""),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(outdir, "01.raw_source_manifest.csv"), row.names = FALSE)

log_message("Raw source type:", source_type)
log_message("Source path:", source_path)
if (nzchar(matrix_dir)) log_message("10X matrix dir:", matrix_dir)
log_message("Script completed: scrna_detect_raw_source")
write_summary(c(
  "scRNA raw source detection completed.",
  paste0("Raw source type: ", source_type),
  paste0("Source path: ", source_path),
  paste0("Annotation file: ", annotation_path),
  "Raw source type was detected before entering the main scRNA analysis line.",
  "",
  "限制与说明",
  "10X and UMI compatibility is handled before mainline_input.rds construction."
))
