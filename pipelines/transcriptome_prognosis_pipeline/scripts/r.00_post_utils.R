# ============================================================
# r.00_post_utils.R — Post-modeling shared utilities
# Sources r.00_modeling_utils.R internally.
# Provides: setup_logging, save_figure, read_gene_file,
#           standardize_barcode, auto_read, extract_label, describe_rds
# ============================================================

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_modeling_utils.R"))

# --- Logging setup -------------------------------------------

setup_logging <- function(logdir, summary_dir, sn, module_name, subdir = NULL) {
  if (!is.null(subdir) && nzchar(subdir)) {
    log_path     <- file.path(logdir, paste0(sn, "_", module_name), subdir)
    summary_path <- file.path(summary_dir, paste0(sn, "_", module_name), subdir)
  } else {
    log_path     <- logdir
    summary_path <- summary_dir
  }
  ensure_dir(log_path)
  ensure_dir(summary_path)

  file_id <- if (!is.null(subdir) && nzchar(subdir)) subdir else paste0(sn, "_", module_name)
  log_file     <- file.path(log_path, paste0(file_id, ".log"))
  summary_file <- file.path(summary_path, paste0(file_id, "_summary.txt"))

  cat("", file = log_file)
  cat("", file = summary_file)

  list(
    log_file     = log_file,
    summary_file = summary_file,

    log_message = function(...) {
      msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
      cat(msg, "\n")
      write(msg, file = log_file, append = TRUE)
    },

    save_note_and_stop = function(note) {
      log_message("[致命错误] ", note)
      writeLines(c(
        paste0(module_name, "分析失败。"),
        paste0("时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste0("原因: ", note)
      ), con = summary_file)
      stop(note, call. = FALSE)
    },

    write_summary = function(lines) {
      writeLines(lines, con = summary_file)
    }
  )
}

# --- Plot saving ---------------------------------------------

save_figure <- function(p, prefix, width = 6, height = 6, dpi = 300) {
  ggsave(paste0(prefix, ".pdf"), p, width = width, height = height, device = "pdf")
  ggsave(paste0(prefix, ".png"), p, width = width, height = height, dpi = dpi, device = "png")
}

# --- Gene file reading ---------------------------------------

read_gene_file <- function(path) {
  if (!file.exists(path)) stop(paste("Gene file not found:", path))
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0) stop("Gene file is empty")
  gene_col <- intersect(c("gene", "gene_symbol", "symbol", "Gene", "genes"), colnames(df))[1]
  if (is.na(gene_col)) stop("Gene file must have a 'gene' column")
  genes <- unique(trimws(as.character(df[[gene_col]])))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) stop("No valid genes found in gene file")
  genes
}

# --- TCGA barcode standardization ----------------------------

standardize_barcode <- function(x, len = 12) {
  x <- gsub("\\.", "-", as.character(x))
  if (len > 0) substr(x, 1, len) else x
}

# --- Auto file reader ----------------------------------------

auto_read <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    if (grepl("\\.csv", path, ignore.case = TRUE)) {
      return(read.csv(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE))
    }
    return(read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) {
    return(read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# --- GMT label extraction (shared by GSEA scripts) -----------

extract_label <- function(filename) {
  base <- basename(filename)
  base <- sub("\\.gmt$", "", base, ignore.case = TRUE)
  base <- sub("^c[1-8]\\.(cp\\.)?", "", base)
  base <- sub("^h\\.", "", base)
  base <- sub("\\.v[0-9]+(\\.[0-9]+)?\\.symbols$", "", base)
  base <- gsub("\\.", "_", base)
  tolower(base)
}

# --- RDS object structure description ------------------------

describe_rds <- function(path) {
  if (!file.exists(path)) return("(RDS file not found)")
  obj <- readRDS(path)
  cls <- class(obj)
  info <- character()
  info <- c(info, paste0("  Type: ", paste(cls, collapse = ", ")))

  if (is.list(obj) && !is.data.frame(obj)) {
    info <- c(info, paste0("  Names: ", paste(names(obj), collapse = ", ")))
    info <- c(info, paste0("  Length: ", length(obj)))
    for (nm in names(obj)) {
      elem <- obj[[nm]]
      if (is.data.frame(elem) || is.matrix(elem)) {
        info <- c(info, paste0("  $", nm, ": ", paste(dim(elem), collapse = " x ")))
      } else if (is.vector(elem)) {
        info <- c(info, paste0("  $", nm, ": length ", length(elem)))
      }
    }
  } else if (is.data.frame(obj)) {
    info <- c(info, paste0("  Dim: ", paste(dim(obj), collapse = " x ")))
    info <- c(info, paste0("  Columns: ", paste(colnames(obj), collapse = ", ")))
  } else if (is.matrix(obj)) {
    info <- c(info, paste0("  Dim: ", paste(dim(obj), collapse = " x ")))
  } else if (is.vector(obj)) {
    info <- c(info, paste0("  Length: ", length(obj)))
  }

  paste(info, collapse = "\n")
}
