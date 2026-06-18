# ============================================================
# r.00_raw_import_utils.R -- raw scRNA import helpers
# Sources r.00_scrna_utils.R internally.
# ============================================================

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_scrna_utils.R"))

parse_chr_vector <- function(x, default = character()) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x)) || tolower(trimws(x)) %in% c("null", "none", "na")) {
    return(default)
  }
  vals <- trimws(unlist(strsplit(as.character(x), ",")))
  vals[!is.na(vals) & nzchar(vals)]
}

has_gzip_magic <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 2) return(FALSE)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2)
  length(magic) == 2 && all(as.integer(magic) == c(0x1f, 0x8b))
}

strip_gzip_suffix <- function(path, depth = 1) {
  out <- sub("\\.gz$", "", path, ignore.case = TRUE)
  if (identical(out, path)) out <- paste0(path, ".unpacked", depth)
  out
}

safe_basename_match <- function(files, pattern) {
  if (is_null_path(pattern)) return(character())
  files[grepl(pattern, basename(files), ignore.case = TRUE, perl = TRUE)]
}

make_barcode_sample <- function(barcode, sample) {
  paste0(as.character(barcode), "_", as.character(sample))
}
