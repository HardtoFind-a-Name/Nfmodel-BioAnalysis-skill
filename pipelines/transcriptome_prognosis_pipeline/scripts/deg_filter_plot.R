args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(grid)
})

deg_all_file <- get_arg("--deg-all")
vst_file <- get_arg("--vst")
group_file <- get_arg("--group")
outdir <- get_arg("--outdir")
padj_cutoff <- as.numeric(get_arg("--padj"))
logfc_cutoff <- as.numeric(get_arg("--logfc"))
use_fallback <- as.logical(get_arg("--use-fallback", required = FALSE, default = "FALSE"))
logdir <- get_arg("--logdir")
sn          <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

heatmap_note_file <- file.path(outdir, "05.deg_heatmap.note.txt")

logsetup <- setup_logging(logdir, summary_dir, sn, "deg_filter")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

save_pdf_plot <- function(plot_obj, path, width, height) {
  ok <- tryCatch({
    ggsave(path, plot = plot_obj, width = width, height = height, device = cairo_pdf)
    TRUE
  }, error = function(e) {
    log_message("cairo_pdf unavailable, fallback to default pdf device: ", conditionMessage(e))
    FALSE
  })

  if (!ok) {
    ggsave(path, plot = plot_obj, width = width, height = height, device = "pdf")
  }
}

log_message("Script started: deg_filter_plot.R")

if (!file.exists(deg_all_file)) {
  save_note_and_stop(paste0("DEG all file does not exist: ", deg_all_file))
}
if (!file.exists(vst_file)) {
  save_note_and_stop(paste0("VST file does not exist: ", vst_file))
}
if (!file.exists(group_file)) {
  save_note_and_stop(paste0("Group file does not exist: ", group_file))
}

res_df <- tryCatch(
  read.csv(deg_all_file, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) save_note_and_stop(paste0("Failed to read DEG all file: ", conditionMessage(e)))
)
vsd_df <- tryCatch(
  read.csv(vst_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) save_note_and_stop(paste0("Failed to read VST file: ", conditionMessage(e)))
)
group_df <- tryCatch(
  read.csv(group_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) save_note_and_stop(paste0("Failed to read group file: ", conditionMessage(e)))
)

selected_p_col    <- if (use_fallback) "pvalue" else "padj"
selected_p_cutoff <- padj_cutoff
selected_p_label  <- if (use_fallback) "p value" else "adjusted p value"

required_deg_cols <- c("gene_symbol", "log2FoldChange", selected_p_col)
missing_deg_cols <- setdiff(required_deg_cols, colnames(res_df))
if (length(missing_deg_cols) > 0) {
  save_note_and_stop(paste0("DEG all file missing required columns: ", paste(missing_deg_cols, collapse = ", ")))
}

res_df$gene_symbol <- trimws(as.character(res_df$gene_symbol))
res_df$log2FoldChange <- suppressWarnings(as.numeric(res_df$log2FoldChange))

colnames(group_df) <- trimws(colnames(group_df))
rownames(group_df) <- trimws(rownames(group_df))

if ("Type" %in% colnames(group_df) && !("group" %in% colnames(group_df))) {
  group_df$group <- group_df$Type
}
if (!("group" %in% colnames(group_df))) {
  save_note_and_stop("Group table must contain either 'group' or 'Type' column.")
}

group_df$group <- trimws(as.character(group_df$group))
group_df$group[group_df$group %in% c("normal", "NORMAL")] <- "Normal"
group_df$group[group_df$group %in% c("tumor", "TUMOR", "Tumour", "tumour")] <- "Tumor"

annotation_df <- data.frame(
  Group = factor(as.character(group_df$group), levels = c("Normal", "Tumor")),
  row.names = rownames(group_df),
  stringsAsFactors = FALSE
)

deg_sig <- subset(res_df, !is.na(res_df[[selected_p_col]]) & res_df[[selected_p_col]] < selected_p_cutoff & abs(log2FoldChange) >= logfc_cutoff)
deg_threshold_used <- sprintf("%s < %.3f & |log2FC| >= %.2f", selected_p_col, selected_p_cutoff, logfc_cutoff)

write.csv(
  deg_sig,
  file.path(outdir, "03.deg_sig.csv"),
  quote = FALSE,
  row.names = FALSE
)

log_message("Significant DEG count: ", nrow(deg_sig))
log_message("Threshold used: ", deg_threshold_used)

volcano_df <- res_df
volcano_df$change <- ifelse(
  !is.na(volcano_df[[selected_p_col]]) & volcano_df[[selected_p_col]] < selected_p_cutoff & volcano_df$log2FoldChange >= logfc_cutoff,
  "UP",
  ifelse(
    !is.na(volcano_df[[selected_p_col]]) & volcano_df[[selected_p_col]] < selected_p_cutoff & volcano_df$log2FoldChange <= -logfc_cutoff,
    "DOWN",
    "NOT"
  )
)

volcano_df$y_value <- -log10(volcano_df$padj_plot)
volcano_df$change <- factor(volcano_df$change, levels = c("DOWN", "NOT", "UP"))

top_up <- head(volcano_df$gene_symbol[volcano_df$change == "UP"], 10)
top_down <- head(volcano_df$gene_symbol[volcano_df$change == "DOWN"], 10)
label_genes <- unique(c(top_up, top_down))

dat_rep <- subset(volcano_df, gene_symbol %in% label_genes & !is.na(y_value) & is.finite(y_value))
rownames(dat_rep) <- dat_rep$gene_symbol

volcano_plot <- ggplot(
  data = subset(volcano_df, !is.na(y_value) & is.finite(y_value)),
  aes(x = log2FoldChange, y = y_value, color = change)
) +
  scale_color_manual(values = c("DOWN" = "#0000FF", "NOT" = "darkgray", "UP" = "#FF0000")) +
  scale_x_continuous(breaks = c(-logfc_cutoff, 0, logfc_cutoff)) +
  geom_point(size = 1.5, alpha = 0.4, na.rm = TRUE) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), lty = 4, col = "darkgray", linewidth = 0.6) +
  geom_hline(yintercept = -log10(padj_cutoff), lty = 4, col = "darkgray", linewidth = 0.6) +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold", color = "black", size = 15),
    legend.text = element_text(face = "bold", color = "black", family = "Times", size = 13),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face = "bold", color = "black", size = 15),
    axis.text.y = element_text(face = "bold", color = "black", size = 15),
    axis.title.x = element_text(face = "bold", color = "black", size = 15),
    axis.title.y = element_text(face = "bold", color = "black", size = 15),
    text = element_text(family = "Times")
  ) +
  geom_label_repel(
    data = dat_rep,
    aes(label = gene_symbol),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    show.legend = FALSE,
    family = "Times"
  ) +
  labs(x = "log2 (Fold Change)", y = "-log10 (adjusted p value)", color = "Change")

ggsave(file.path(outdir, "04.deg_volcano.png"), plot = volcano_plot, width = 8, height = 6, dpi = 600)
save_pdf_plot(volcano_plot, file.path(outdir, "04.deg_volcano.pdf"), width = 8, height = 6)

vsd_mat <- as.matrix(vsd_df)

cat("", file = heatmap_note_file)

if (nrow(deg_sig) > 1) {
  up_df <- deg_sig[deg_sig$log2FoldChange > 0, , drop = FALSE]
  down_df <- deg_sig[deg_sig$log2FoldChange < 0, , drop = FALSE]

  up_df <- up_df[order(up_df[[selected_p_col]], -up_df$log2FoldChange), , drop = FALSE]
  down_df <- down_df[order(down_df$padj, down_df$log2FoldChange), , drop = FALSE]

  top_up_heatmap <- head(up_df$gene_symbol, 10)
  top_down_heatmap <- head(down_df$gene_symbol, 10)
  heatmap_genes <- unique(c(top_up_heatmap, top_down_heatmap))
  heatmap_genes <- intersect(heatmap_genes, rownames(vsd_mat))

  if (length(heatmap_genes) >= 2) {
    heatmap_fc <- res_df$log2FoldChange[match(heatmap_genes, res_df$gene_symbol)]
    heatmap_genes <- heatmap_genes[order(heatmap_fc, decreasing = TRUE, na.last = NA)]

    common_samples <- intersect(rownames(annotation_df), colnames(vsd_mat))
    heatmap_mat <- vsd_mat[heatmap_genes, common_samples, drop = FALSE]
    annotation_use <- annotation_df[common_samples, , drop = FALSE]
    heatmap_mat <- t(scale(t(heatmap_mat)))
    heatmap_mat[is.na(heatmap_mat)] <- 0
    heatmap_mat[heatmap_mat < -2] <- -2
    heatmap_mat[heatmap_mat > 2] <- 2

    heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
    ann_colors <- list(Group = c("Tumor" = "#E64B35", "Normal" = "#4DBBD5"))

    png(file.path(outdir, "05.deg_heatmap.png"), width = 7, height = 9, units = "in", res = 600)
    pheatmap(
      heatmap_mat,
      color = heat_colors,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      scale = "none",
      show_colnames = FALSE,
      show_rownames = TRUE,
      annotation_col = annotation_use,
      annotation_colors = ann_colors,
      border_color = NA,
      fontsize = 9,
      cellheight = 8,
      angle_col = 0
    )
    dev.off()

    pdf(file.path(outdir, "05.deg_heatmap.pdf"), width = 7, height = 9, family = "Times")
    pheatmap(
      heatmap_mat,
      color = heat_colors,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "none",
      show_colnames = FALSE,
      show_rownames = FALSE,
      annotation_col = annotation_use,
      annotation_colors = ann_colors,
      border_color = NA,
      fontsize = 9,
      cellheight = 8,
      angle_col = 0
    )
    dev.off()

    writeLines("Heatmap generated successfully.", heatmap_note_file)
  } else {
    note <- "Skipped heatmap: fewer than 2 usable significant genes."
    log_message(note)
    writeLines(note, heatmap_note_file)
  }
} else {
  note <- "Skipped heatmap: fewer than 2 significant genes."
  log_message(note)
  writeLines(note, heatmap_note_file)
}

up_genes <- volcano_df$gene_symbol[volcano_df$change == "UP"]
down_genes <- volcano_df$gene_symbol[volcano_df$change == "DOWN"]
up_genes <- up_genes[!is.na(up_genes) & up_genes != ""]
down_genes <- down_genes[!is.na(down_genes) & down_genes != ""]

write_summary(c(
  "=== DEG filter/plot ===",
  paste0("Threshold: ", deg_threshold_used),
  paste0("Significant DEGs: ", nrow(deg_sig)),
  paste0("Upregulated: ", length(up_genes), " | Downregulated: ", length(down_genes)),
  paste0("Top 10 upregulated: ", if (length(top_up) > 0) paste(top_up, collapse = ", ") else "None"),
  paste0("Top 10 downregulated: ", if (length(top_down) > 0) paste(top_down, collapse = ", ") else "None"),
  "",
  "=== 限制与说明 ===",
  "DEG thresholds (padj and log2FC) are user-configurable.",
  "Volcano plot labels show top 10 up/down genes by significance.",
  "Heatmap uses top 10 DEGs per direction, z-score normalized.",
  "Cairo PDF is attempted first; fallback to default PDF device if unavailable."
))

log_message("DEG filter/plot step finished.")
message("DEG filter/plot step finished.")
