# ============================================================
# r.00_scrna_utils.R -- scRNA shared utilities
# Sources r.00_post_utils.R internally.
# ============================================================

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

as_bool <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

parse_num_vector <- function(x, default = NULL) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(default)
  vals <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) default else vals
}

parse_int_vector <- function(x, default = NULL) {
  vals <- parse_num_vector(x, default = default)
  if (is.null(vals)) NULL else as.integer(vals)
}

is_null_path <- function(x) {
  is.null(x) || is.na(x) || !nzchar(trimws(x)) || tolower(trimws(x)) %in% c("null", "none", "na")
}

parse_chr_vector <- function(x, default = character()) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x)) || tolower(trimws(x)) %in% c("null", "none", "na")) {
    return(default)
  }
  vals <- trimws(unlist(strsplit(as.character(x), ",")))
  vals[!is.na(vals) & nzchar(vals)]
}

require_package_or_stop <- function(pkg, save_note_and_stop = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Required R package is not installed: ", pkg)
    if (!is.null(save_note_and_stop)) save_note_and_stop(msg)
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

get_assay_data_compat <- function(obj, assay = NULL, layer = "data") {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) {
      Seurat::GetAssayData(obj, assay = assay, slot = layer)
    }
  )
}

detect_input_format <- function(obj, input_format = "auto") {
  fmt <- tolower(trimws(input_format))
  if (fmt != "auto") return(fmt)
  if (inherits(obj, "Seurat")) return("seurat_rds")
  if (is.list(obj) && all(c("expr_mat", "meta") %in% names(obj))) return("mainline_rds")
  "unknown"
}

count_detected_genes <- function(obj, assay = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  counts <- get_assay_data_compat(obj, assay = assay, layer = "counts")
  sum(Matrix::rowSums(counts > 0) > 0)
}

add_group_from_origin <- function(obj, origin_col, group_col, normal_origin, tumor_origin) {
  if (!origin_col %in% colnames(obj@meta.data)) return(obj)
  origin <- as.character(obj@meta.data[[origin_col]])
  grp <- ifelse(origin == normal_origin, "normal",
                ifelse(origin == tumor_origin, "tumor", NA_character_))
  if (!group_col %in% colnames(obj@meta.data) || any(is.na(obj@meta.data[[group_col]]))) {
    obj@meta.data[[group_col]] <- grp
  }
  obj
}

read_scrna_input <- function(input_rds,
                             input_format = "auto",
                             project_name = "scRNA",
                             sample_col = "Sample",
                             origin_col = "Sample_Origin",
                             group_col = "group",
                             normal_origin = "nLung",
                             tumor_origin = "tLung",
                             min_cells_per_gene = 3,
                             min_features_create = 0) {
  if (!file.exists(input_rds)) stop(paste0("Input RDS not found: ", input_rds), call. = FALSE)
  obj <- readRDS(input_rds)
  fmt <- detect_input_format(obj, input_format)

  if (fmt == "seurat_rds") {
    if (!inherits(obj, "Seurat")) stop("Input format is seurat_rds but object is not a Seurat object.", call. = FALSE)
    seu <- obj
  } else if (fmt == "mainline_rds") {
    if (!is.list(obj) || !all(c("expr_mat", "meta") %in% names(obj))) {
      stop("mainline_rds input must be a list with expr_mat and meta.", call. = FALSE)
    }
    expr_mat <- obj$expr_mat
    meta <- as.data.frame(obj$meta, stringsAsFactors = FALSE)
    if (is.null(colnames(expr_mat))) stop("expr_mat must have cell barcodes as column names.", call. = FALSE)
    if ("cell_id" %in% colnames(meta)) rownames(meta) <- as.character(meta$cell_id)
    common_cells <- intersect(colnames(expr_mat), rownames(meta))
    if (length(common_cells) == 0) stop("expr_mat and meta have no overlapping cell IDs.", call. = FALSE)
    expr_mat <- expr_mat[, common_cells, drop = FALSE]
    meta <- meta[common_cells, , drop = FALSE]
    seu <- Seurat::CreateSeuratObject(
      counts = expr_mat,
      min.cells = min_cells_per_gene,
      min.features = min_features_create,
      project = project_name
    )
    seu <- Seurat::AddMetaData(seu, metadata = meta[colnames(seu), , drop = FALSE])
  } else {
    stop("Unsupported input_format. Use auto, seurat_rds, or mainline_rds.", call. = FALSE)
  }

  if (sample_col %in% colnames(seu@meta.data)) {
    seu@meta.data[[sample_col]] <- as.character(seu@meta.data[[sample_col]])
  }
  seu <- add_group_from_origin(seu, origin_col, group_col, normal_origin, tumor_origin)
  list(object = seu, format = fmt)
}

read_target_genes <- function(target_gene_file, target_gene_col = "gene") {
  if (is_null_path(target_gene_file)) return(character())
  if (!file.exists(target_gene_file)) stop(paste0("Target gene file not found: ", target_gene_file), call. = FALSE)
  df <- auto_read(target_gene_file)
  if (!target_gene_col %in% colnames(df)) {
    stop(paste0("Target gene column not found: ", target_gene_col), call. = FALSE)
  }
  genes <- unique(trimws(as.character(df[[target_gene_col]])))
  genes[!is.na(genes) & nzchar(genes)]
}

read_seurat_rds <- function(path) {
  if (!file.exists(path)) stop(paste0("Seurat RDS not found: ", path), call. = FALSE)
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop("Input object is not a Seurat object.", call. = FALSE)
  obj
}

save_plot_dual <- function(plot_obj, outdir, prefix, width = 8, height = 6) {
  save_figure(plot_obj, file.path(outdir, prefix), width = width, height = height)
}

sort_cluster_ids <- function(x) {
  x_chr <- unique(as.character(x))
  suppressWarnings(x_num <- as.numeric(x_chr))
  if (all(!is.na(x_num))) x_chr[order(x_num)] else sort(x_chr)
}

resolve_cluster_col <- function(meta_cols, annotation_cluster_col = NULL, annotation_resolution = NULL) {
  if (!is.null(annotation_cluster_col) && nzchar(annotation_cluster_col)) return(annotation_cluster_col)
  if (!is.null(annotation_resolution) && nzchar(as.character(annotation_resolution))) {
    return(paste0("RNA_snn_res.", annotation_resolution))
  }
  if ("seurat_clusters" %in% meta_cols) return("seurat_clusters")
  NA_character_
}

split_markers <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character())
  unique(trimws(unlist(strsplit(as.character(x), ","))))
}
