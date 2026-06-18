rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_raw_import_utils.R"))

raw_source_manifest <- get_arg("--raw-source-manifest")
outdir <- get_arg("--outdir")
keep_origin <- parse_chr_vector(get_arg("--keep-origin", required = FALSE, default = "nLung,tLung"))
exclude_origin <- parse_chr_vector(get_arg("--exclude-origin", required = FALSE, default = "tL/B"))
normal_origins <- parse_chr_vector(get_arg("--normal-origins", required = FALSE, default = "nLung"))
tumor_origins <- parse_chr_vector(get_arg("--tumor-origins", required = FALSE, default = "tLung"))
origin_col <- get_arg("--origin-col", required = FALSE, default = "Sample_Origin")
sample_col <- get_arg("--sample-col", required = FALSE, default = "Sample")
annotation_cell_id_candidates <- parse_chr_vector(get_arg("--annotation-cell-id-candidates", required = FALSE, default = "Index,barcode_sample,Barcode"))
strip_10x_barcode_suffix <- as_bool(get_arg("--strip-10x-barcode-suffix", required = FALSE, default = "false"))
tenx_assay <- get_arg("--tenx-assay", required = FALSE, default = "Gene Expression")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_build_mainline_input")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

strip_barcode <- function(x) sub("-1$", "", as.character(x))

read_expr_source <- function(source_type, source_path, tenx_assay) {
  if (source_type == "umi_rds") {
    expr_obj <- readRDS(source_path)
    if (is.data.frame(expr_obj)) return(as.matrix(expr_obj))
    if (inherits(expr_obj, c("matrix", "dgCMatrix"))) return(expr_obj)
    stop(paste0("Unsupported UMI RDS object class: ", paste(class(expr_obj), collapse = ", ")), call. = FALSE)
  }
  if (source_type == "10x_dir") {
    expr_obj <- Seurat::Read10X(data.dir = source_path)
    if (is.list(expr_obj)) {
      if (tenx_assay %in% names(expr_obj)) return(expr_obj[[tenx_assay]])
      return(expr_obj[[1]])
    }
    return(expr_obj)
  }
  if (source_type == "10x_h5") {
    expr_obj <- Seurat::Read10X_h5(filename = source_path)
    if (is.list(expr_obj)) {
      if (tenx_assay %in% names(expr_obj)) return(expr_obj[[tenx_assay]])
      return(expr_obj[[1]])
    }
    return(expr_obj)
  }
  stop(paste0("Unsupported raw source type for mainline build: ", source_type), call. = FALSE)
}

match_cells <- function(expr_cells, anno, candidates, strip_suffix = FALSE) {
  for (candidate in candidates) {
    if (!candidate %in% colnames(anno)) next
    vals <- as.character(anno[[candidate]])
    names(vals) <- rownames(anno)
    common <- intersect(expr_cells, vals)
    if (length(common) > 0) {
      anno_rows <- names(vals)[match(common, vals)]
      return(list(expr_cells = common, anno_rows = anno_rows, match_mode = candidate, stripped = FALSE))
    }
    if (strip_suffix) {
      expr_key <- strip_barcode(expr_cells)
      anno_key <- strip_barcode(vals)
      keep_expr <- !duplicated(expr_key)
      keep_anno <- !duplicated(anno_key)
      expr_map <- expr_cells[keep_expr]
      names(expr_map) <- expr_key[keep_expr]
      anno_map <- rownames(anno)[keep_anno]
      names(anno_map) <- anno_key[keep_anno]
      common_key <- intersect(names(expr_map), names(anno_map))
      if (length(common_key) > 0) {
        return(list(expr_cells = unname(expr_map[common_key]), anno_rows = unname(anno_map[common_key]),
                    match_mode = paste0(candidate, "_strip_10x_suffix"), stripped = TRUE))
      }
    }
  }
  NULL
}

log_message("Script started: scrna_build_mainline_input")
if (!file.exists(raw_source_manifest)) save_note_and_stop(paste0("Raw source manifest not found: ", raw_source_manifest))
manifest <- read.csv(raw_source_manifest, stringsAsFactors = FALSE, check.names = FALSE)
source_type <- manifest$raw_source_type[manifest$role == "raw_source"][1]
source_path <- manifest$source_path[manifest$role == "raw_source"][1]
annotation_path <- manifest$source_path[manifest$role == "annotation"][1]
if (is.na(source_type) || is.na(source_path) || !file.exists(source_path)) save_note_and_stop("Raw source manifest is invalid.")
if (is.na(annotation_path) || !file.exists(annotation_path)) save_note_and_stop("Annotation file not found.")

log_message("Read raw source:", source_path)
log_message("Raw source type:", source_type)
expr_mat <- tryCatch(read_expr_source(source_type, source_path, tenx_assay),
                     error = function(e) save_note_and_stop(paste0("Expression source read failed: ", conditionMessage(e))))
if (!inherits(expr_mat, c("matrix", "dgCMatrix"))) expr_mat <- as.matrix(expr_mat)
if (is.null(rownames(expr_mat)) || is.null(colnames(expr_mat))) save_note_and_stop("Expression matrix must have rownames and colnames.")

log_message("Read annotation:", annotation_path)
anno <- tryCatch(auto_read(annotation_path), error = function(e) save_note_and_stop(paste0("Annotation read failed: ", conditionMessage(e))))
required_cols <- c(sample_col, origin_col)
missing_required <- setdiff(required_cols, colnames(anno))
if (length(missing_required) > 0) save_note_and_stop(paste0("Annotation missing required columns: ", paste(missing_required, collapse = ", ")))

if ("Barcode" %in% colnames(anno) && sample_col %in% colnames(anno)) {
  anno$barcode_sample <- make_barcode_sample(anno$Barcode, anno[[sample_col]])
}
if (is.null(rownames(anno)) || any(!nzchar(rownames(anno))) || anyDuplicated(rownames(anno)) > 0) {
  rownames(anno) <- paste0("anno_row_", seq_len(nrow(anno)))
}

match_res <- match_cells(colnames(expr_mat), anno, annotation_cell_id_candidates, strip_10x_barcode_suffix)
if (is.null(match_res)) {
  save_note_and_stop(paste0("Could not match expression colnames to annotation candidates: ",
                            paste(annotation_cell_id_candidates, collapse = ", ")))
}
expr_mat <- expr_mat[, match_res$expr_cells, drop = FALSE]
anno2 <- anno[match_res$anno_rows, , drop = FALSE]
rownames(anno2) <- match_res$expr_cells

origin_values <- as.character(anno2[[origin_col]])
if (length(keep_origin) > 0) {
  keep_cells <- rownames(anno2)[origin_values %in% keep_origin]
} else {
  keep_cells <- rownames(anno2)
}
if (length(exclude_origin) > 0) {
  keep_cells <- setdiff(keep_cells, rownames(anno2)[origin_values %in% exclude_origin])
}
if (length(keep_cells) == 0) save_note_and_stop("No cells retained after origin filtering.")

expr_mat_main <- expr_mat[, keep_cells, drop = FALSE]
meta_main <- anno2[keep_cells, , drop = FALSE]
meta_origin <- as.character(meta_main[[origin_col]])
meta_main$group <- ifelse(meta_origin %in% normal_origins, "normal",
                          ifelse(meta_origin %in% tumor_origins, "tumor", NA_character_))
if (any(is.na(meta_main$group))) {
  bad <- unique(meta_origin[is.na(meta_main$group)])
  save_note_and_stop(paste0("Some retained origins cannot be mapped to group: ", paste(bad, collapse = ", ")))
}
meta_main$group <- factor(meta_main$group, levels = c("normal", "tumor"))

origin_summary <- meta_main %>%
  dplyr::count(.data[[origin_col]], group, name = "cell_count") %>%
  dplyr::arrange(group, .data[[origin_col]])
sample_summary <- meta_main %>%
  dplyr::count(.data[[sample_col]], .data[[origin_col]], group, name = "cell_count") %>%
  dplyr::arrange(group, .data[[origin_col]], .data[[sample_col]])

write.csv(origin_summary, file.path(outdir, "01.mainline_origin_summary.csv"), row.names = FALSE)
write.csv(sample_summary, file.path(outdir, "02.mainline_sample_summary.csv"), row.names = FALSE)
meta_out <- meta_main
meta_out$cell_id <- rownames(meta_out)
write.csv(meta_out, file.path(outdir, "03.mainline_metadata.csv"), row.names = FALSE)

main_input <- list(
  expr_mat = expr_mat_main,
  meta = meta_main,
  raw_source_type = source_type,
  keep_origin = keep_origin,
  exclude_origin = exclude_origin,
  normal_origins = normal_origins,
  tumor_origins = tumor_origins,
  match_mode = match_res$match_mode,
  strip_10x_barcode_suffix = strip_10x_barcode_suffix,
  source_expr_file = source_path,
  source_anno_file = annotation_path
)
saveRDS(main_input, file.path(outdir, "01.mainline_input.rds"))

log_message("Matched cells:", length(match_res$expr_cells))
log_message("Retained cells:", ncol(expr_mat_main))
log_message("Genes:", nrow(expr_mat_main))
log_message("Match mode:", match_res$match_mode)
log_message("Script completed: scrna_build_mainline_input")

write_summary(c(
  "scRNA mainline input construction completed.",
  paste0("Raw source type: ", source_type),
  paste0("Matched cells: ", length(match_res$expr_cells)),
  paste0("Retained cells: ", ncol(expr_mat_main)),
  paste0("Genes: ", nrow(expr_mat_main)),
  paste0("Match mode: ", match_res$match_mode),
  paste0("Retained origins: ", paste(unique(meta_origin), collapse = ", ")),
  "A standard list(expr_mat, meta) mainline_input.rds was created for downstream analysis.",
  paste0("Object file: 01.mainline_input.rds\n", describe_rds(file.path(outdir, "01.mainline_input.rds"))),
  "",
  "限制与说明",
  "Raw 10X/UMI compatibility is handled before this standard mainline object enters the analysis pipeline."
))
