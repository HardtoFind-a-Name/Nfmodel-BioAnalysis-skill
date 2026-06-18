rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_scrna_utils.R"))

input_rds <- get_arg("--input-rds")
input_format <- get_arg("--input-format", required = FALSE, default = "auto")
outdir <- get_arg("--outdir")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")
project_name <- get_arg("--project-name", required = FALSE, default = "scRNA")
sample_col <- get_arg("--sample-col", required = FALSE, default = "Sample")
origin_col <- get_arg("--origin-col", required = FALSE, default = "Sample_Origin")
group_col <- get_arg("--group-col", required = FALSE, default = "group")
normal_origin <- get_arg("--normal-origin", required = FALSE, default = "nLung")
tumor_origin <- get_arg("--tumor-origin", required = FALSE, default = "tLung")
min_cells_per_gene <- as.integer(get_arg("--min-cells-per-gene", required = FALSE, default = "3"))
min_features_create <- as.integer(get_arg("--min-features-create", required = FALSE, default = "0"))

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_prepare_input")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_prepare_input")
log_message("Input RDS:", input_rds)
log_message("Input format:", input_format)

res <- tryCatch(
  read_scrna_input(
    input_rds = input_rds,
    input_format = input_format,
    project_name = project_name,
    sample_col = sample_col,
    origin_col = origin_col,
    group_col = group_col,
    normal_origin = normal_origin,
    tumor_origin = tumor_origin,
    min_cells_per_gene = min_cells_per_gene,
    min_features_create = min_features_create
  ),
  error = function(e) save_note_and_stop(conditionMessage(e))
)

obj <- res$object
obj@misc$scrna_pipeline <- list(input_format = res$format, project_name = project_name)
out_rds <- file.path(outdir, "01.input_seurat.rds")
saveRDS(obj, out_rds)

meta_df <- obj@meta.data
meta_df$cell_id <- rownames(meta_df)
write.csv(meta_df, file.path(outdir, "02.input_metadata.csv"), row.names = FALSE)
dim_df <- data.frame(metric = c("cells", "genes"), value = c(ncol(obj), nrow(obj)))
write.csv(dim_df, file.path(outdir, "03.input_dimensions.csv"), row.names = FALSE)

log_message("Cells:", ncol(obj), "Genes:", nrow(obj))
log_message("Script completed: scrna_prepare_input")

group_note <- if (group_col %in% colnames(obj@meta.data)) {
  paste0("Group levels: ", paste(unique(as.character(obj@meta.data[[group_col]])), collapse = ", "))
} else {
  paste0("Group column not available: ", group_col)
}

write_summary(c(
  "scRNA input preparation completed.",
  paste0("Input format: ", res$format),
  paste0("Final cells: ", ncol(obj)),
  paste0("Final genes: ", nrow(obj)),
  group_note,
  "The input object was standardized as a Seurat RDS for downstream scRNA modules.",
  "",
  "限制与说明",
  "This step validates and standardizes the input object only; biological QC is performed downstream."
))
