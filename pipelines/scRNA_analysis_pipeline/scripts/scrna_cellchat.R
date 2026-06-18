rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_scrna_utils.R"))

input_rds <- get_arg("--input-rds")
outdir <- get_arg("--outdir")
sn <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
logdir <- get_arg("--logdir")
group_col <- get_arg("--group-col", required = FALSE, default = "group")
celltype_col <- get_arg("--celltype-col", required = FALSE, default = "celltype_manual")
species <- tolower(get_arg("--species", required = FALSE, default = "human"))
min_cells <- as.integer(get_arg("--min-cells", required = FALSE, default = "10"))

ensure_dir(outdir)
logsetup <- setup_logging(logdir, summary_dir, sn, "scrna_cellchat")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: scrna_cellchat")
if (!requireNamespace("CellChat", quietly = TRUE)) save_note_and_stop("CellChat package is not installed.")
obj <- tryCatch(read_seurat_rds(input_rds), error = function(e) save_note_and_stop(conditionMessage(e)))
missing_cols <- setdiff(c(group_col, celltype_col), colnames(obj@meta.data))
if (length(missing_cols) > 0) save_note_and_stop(paste0("Missing metadata columns: ", paste(missing_cols, collapse = ", ")))

groups <- sort(unique(as.character(obj@meta.data[[group_col]])))
groups <- groups[!is.na(groups) & nzchar(groups)]
if (length(groups) == 0) save_note_and_stop("No valid groups for CellChat.")

run_one <- function(group_name) {
  group_dir <- file.path(outdir, group_name)
  ensure_dir(group_dir)
  cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[group_col]]) == group_name]
  sub <- subset(obj, cells = cells)
  ct_counts <- table(sub@meta.data[[celltype_col]])
  keep_types <- names(ct_counts)[ct_counts >= min_cells]
  sub <- subset(sub, cells = rownames(sub@meta.data)[sub@meta.data[[celltype_col]] %in% keep_types])
  if (ncol(sub) < min_cells || length(unique(sub@meta.data[[celltype_col]])) < 2) {
    stop("Too few cells or cell types after min_cells filtering.")
  }

  data_input <- get_assay_data_compat(sub, assay = DefaultAssay(sub), layer = "data")
  meta <- sub@meta.data
  meta[[celltype_col]] <- factor(meta[[celltype_col]])
  cellchat <- CellChat::createCellChat(object = data_input, meta = meta, group.by = celltype_col)
  if (species == "mouse") {
    cellchat@DB <- CellChat::CellChatDB.mouse
  } else {
    cellchat@DB <- CellChat::CellChatDB.human
  }
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  saveRDS(cellchat, file.path(group_dir, paste0(group_name, "_cellchat.rds")))
  comm <- CellChat::subsetCommunication(cellchat)
  write.csv(comm, file.path(group_dir, paste0(group_name, "_cellchat_communications.csv")), row.names = FALSE)
  data.frame(group = group_name, cells = ncol(sub), celltypes = length(unique(sub@meta.data[[celltype_col]])),
             interactions = nrow(comm), status = "completed", note = "", stringsAsFactors = FALSE)
}

rows <- list()
for (grp in groups) {
  log_message("Running CellChat group:", grp)
  rows[[grp]] <- tryCatch(
    run_one(grp),
    error = function(e) {
      log_message("[警告] CellChat failed for", grp, ":", conditionMessage(e))
      data.frame(group = grp, cells = NA_integer_, celltypes = NA_integer_, interactions = NA_integer_,
                 status = "failed", note = conditionMessage(e), stringsAsFactors = FALSE)
    }
  )
}
status <- dplyr::bind_rows(rows)
write.csv(status, file.path(outdir, "01_cellchat_group_status.csv"), row.names = FALSE)

completed <- sum(status$status == "completed")
log_message("Script completed: scrna_cellchat")
write_summary(c(
  "scRNA CellChat analysis completed.",
  paste0("Groups requested: ", paste(groups, collapse = ", ")),
  paste0("Groups completed: ", completed),
  paste0("Minimum cells per cell type: ", min_cells),
  "CellChat was run separately within each group using manual cell type labels.",
  "",
  "限制与说明",
  "Groups with too few cells or too few cell types are skipped and recorded in the status table."
))
