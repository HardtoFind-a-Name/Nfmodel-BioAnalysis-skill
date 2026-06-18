set.seed(123)

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(patchwork)
})

# --- CLI args ---
args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

risk_file    <- get_arg("--risk-file")
count_file   <- get_arg("--count-matrix")
train_id     <- get_arg("--train-id")
outdir       <- get_arg("--outdir")
logdir       <- get_arg("--logdir")
sn           <- get_arg("--sn")
summary_dir  <- get_arg("--summary-dir")
gmt_dir      <- get_arg("--gmt-dir")
gmt_list_str <- get_arg("--gmt-gsea_risk")
sort_by      <- get_arg("--gsea-sort-by", required = FALSE, default = "p.adjust")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

top_n_pathway <- 5
nes_cutoff    <- 1
p_cutoff_gsea <- 0.05
fdr_cutoff    <- 0.25

# --- Logging ---
logsetup <- setup_logging(logdir, summary_dir, sn, "gsea_risk_nes")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

# --- Plot helpers ---
save_plot_pdf_png <- function(plot_obj, out_pdf, out_png, width = 10.5, height = 8.2, dpi = 320) {
  grDevices::pdf(out_pdf, width = width, height = height, onefile = FALSE)
  print(plot_obj)
  grDevices::dev.off()
  grDevices::png(out_png, width = width, height = height, units = "in", res = dpi)
  print(plot_obj)
  grDevices::dev.off()
}

filter_gsea_result <- function(gsea_obj) {
  gsea_df <- as.data.frame(gsea_obj)
  if (is.null(gsea_df) || nrow(gsea_df) == 0) return(data.frame())
  gsea_df %>%
    filter(!is.na(NES), !is.na(pvalue), !is.na(p.adjust),
           abs(NES) > nes_cutoff, pvalue < p_cutoff_gsea, p.adjust < fdr_cutoff) %>%
    arrange(pvalue, desc(abs(NES)))
}

plot_gsea_multiline <- function(gsea_obj, gsea_df, title_text, out_pdf, out_png, top_n = 5, sort_by = "p.adjust") {
  if (is.null(gsea_obj) || is.null(gsea_df) || nrow(gsea_df) == 0) return(NULL)

  top_df <- if (sort_by == "NES") {
    gsea_df %>% dplyr::arrange(desc(NES)) %>% dplyr::slice(1:min(top_n, n()))
  } else if (sort_by == "pvalue") {
    gsea_df %>% dplyr::arrange(pvalue) %>% dplyr::slice(1:min(top_n, n()))
  } else {
    gsea_df %>% dplyr::arrange(p.adjust) %>% dplyr::slice(1:min(top_n, n()))
  }

  top_ids <- top_df$ID
  top_desc <- top_df$Description

  gsdata <- do.call(rbind, lapply(seq_along(top_ids), function(i) {
    tmp <- enrichplot:::gsInfo(gsea_obj, geneSetID = top_ids[i])
    tmp$ID <- top_ids[i]
    tmp$Description <- top_desc[i]
    tmp
  }))

  gsdata$Description <- factor(gsdata$Description, levels = top_desc)
  x_min <- min(gsdata$x, na.rm = TRUE)
  x_max <- max(gsdata$x, na.rm = TRUE)

  p1 <- ggplot(gsdata, aes(x = x, y = runningScore, color = Description)) +
    geom_line(linewidth = 0.9, alpha = 0.95) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    labs(title = title_text, x = NULL, y = "Running Enrichment Score") +
    guides(color = guide_legend(byrow = TRUE)) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(0.3, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      legend.spacing.y = grid::unit(0.03, "cm"),
      legend.margin = margin(4, 4, 4, 4),
      legend.background = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.72), color = "grey80"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 0, 5.5)
    )

  hit_df <- gsdata %>% dplyr::filter(position == 1)
  if (nrow(hit_df) > 0) {
    hit_df$Description <- factor(hit_df$Description, levels = rev(top_desc))
    p2 <- ggplot(hit_df, aes(x = x, y = Description, color = Description)) +
      geom_linerange(
        aes(ymin = as.numeric(Description) - 0.35,
            ymax = as.numeric(Description) + 0.35),
        linewidth = 0.45
      ) +
      scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
      labs(x = NULL, y = NULL) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0, 5.5, 0, 5.5)
      )
  } else {
    p2 <- ggplot() + theme_void()
  }

  rank_df <- gsdata %>%
    dplyr::distinct(x, geneList) %>%
    dplyr::arrange(x)

  p3 <- ggplot(rank_df, aes(x = x, y = geneList)) +
    geom_area(fill = "grey70") +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    labs(x = "Rank in Ordered Dataset", y = "Ranked List Metric") +
    theme_classic(base_size = 12) +
    theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

  p <- p1 / p2 / p3 + patchwork::plot_layout(ncol = 1, heights = c(4.2, 1.2, 1.8))
  save_plot_pdf_png(p, out_pdf, out_png)
}

# --- Load GMTs ---
log_message("Loading GMT files")
gmt_files <- trimws(unlist(strsplit(gmt_list_str, ",")))
gmt_files <- gmt_files[nzchar(gmt_files)]
gmts <- lapply(gmt_files, function(fn) {
  full <- file.path(gmt_dir, trimws(fn))
  if (!file.exists(full)) stop(paste("GMT file not found:", full))
  log_message("  Loading: ", full)
  list(label = extract_label(fn), gmt = read.gmt(full), file = full)
})
names(gmts) <- sapply(gmts, `[[`, "label")

# --- Input validation ---
if (!file.exists(count_file)) stop(paste("Count matrix not found:", count_file))
if (!file.exists(risk_file)) stop(paste("Risk file not found:", risk_file))

# --- Read count matrix ---
log_message("Reading count matrix")
count_df <- fread(count_file, data.table = FALSE, check.names = FALSE)
gene_col <- colnames(count_df)[1]
count_df[[gene_col]] <- toupper(trimws(count_df[[gene_col]]))
count_df <- count_df[count_df[[gene_col]] != "" & !is.na(count_df[[gene_col]]), ]
count_df <- count_df[!duplicated(count_df[[gene_col]]), ]
rownames(count_df) <- count_df[[gene_col]]
count_mat <- as.matrix(count_df[, -1, drop = FALSE])
storage.mode(count_mat) <- "numeric"
if (any(count_mat < 0, na.rm = TRUE)) stop("Count matrix contains negative values")

# --- Read risk data ---
log_message("Reading risk file")
risk_df <- read.csv(risk_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample" %in% colnames(risk_df)) colnames(risk_df)[1] <- "sample"
risk_df$sample <- trimws(as.character(risk_df$sample))
risk_df$riskScore <- suppressWarnings(as.numeric(risk_df$riskScore))
risk_df$riskGroup <- trimws(as.character(risk_df$riskGroup))
risk_df <- risk_df[!is.na(risk_df$sample) & nzchar(risk_df$sample), ]
risk_df <- risk_df[risk_df$riskGroup %in% c("Low risk", "High risk"), ]
risk_df <- risk_df[!duplicated(risk_df$sample), ]
rownames(risk_df) <- risk_df$sample

# --- Align samples ---
common_samples <- intersect(colnames(count_mat), rownames(risk_df))
if (length(common_samples) < 10) stop("Too few common samples")
count_mat <- count_mat[, common_samples, drop = FALSE]
risk_df <- risk_df[match(common_samples, rownames(risk_df)), , drop = FALSE]
risk_df$riskGroup <- factor(risk_df$riskGroup, levels = c("Low risk", "High risk"))
rownames(risk_df) <- risk_df$sample
grp_n <- table(risk_df$riskGroup)
if (any(grp_n < 3)) stop("Too few samples in risk groups")

# --- Low-expression filter ---
log_message("Filtering low-expression genes")
keep_genes <- rowSums(count_mat >= 10) >= max(3, floor(ncol(count_mat) * 0.1))
count_mat <- count_mat[keep_genes, , drop = FALSE]
if (nrow(count_mat) < 1000) stop("Too few genes after filtering")

# --- DESeq2 ---
log_message("Running DESeq2")
col_data <- data.frame(row.names = rownames(risk_df), riskGroup = risk_df$riskGroup)
dds <- DESeqDataSetFromMatrix(countData = round(count_mat), colData = col_data, design = ~ riskGroup)
dds <- DESeq(dds)
res <- results(dds, contrast = c("riskGroup", "High risk", "Low risk"))
res_df <- as.data.frame(res)
res_df$GeneSymbol <- rownames(res_df)
res_df <- res_df %>% filter(!is.na(log2FoldChange), !is.na(pvalue)) %>% arrange(desc(log2FoldChange))
fwrite(res_df, file = file.path(outdir, "01.deseq2_results.csv"))

# --- GSEA ranking ---
geneList <- res_df$log2FoldChange
names(geneList) <- res_df$GeneSymbol
geneList <- sort(geneList, decreasing = TRUE)
fwrite(data.frame(GeneSymbol = names(geneList), log2FoldChange = as.numeric(geneList)),
       file = file.path(outdir, "02.gsea_geneList_log2fc.csv"))

# --- Per-GMT GSEA ---
for (gmt_item in gmts) {
  label <- gmt_item$label
  log_message("GSEA with ", label)

  gmt_df <- gmt_item$gmt
  gmt_df$gene <- toupper(trimws(gmt_df$gene))
  gmt_df <- gmt_df[!is.na(gmt_df$gene) & nzchar(gmt_df$gene), ]

  gsea_obj <- GSEA(geneList = geneList, TERM2GENE = gmt_df, pvalueCutoff = 1,
                   pAdjustMethod = "BH", verbose = FALSE)
  gsea_all_df <- as.data.frame(gsea_obj)
  if (is.null(gsea_all_df)) gsea_all_df <- data.frame()
  fwrite(gsea_all_df, file = file.path(outdir, paste0("03.", label, ".gsea_all.csv")))

  gsea_sig_df <- filter_gsea_result(gsea_obj)
  fwrite(gsea_sig_df, file = file.path(outdir, paste0("04.", label, ".gsea_significant.csv")))

  if (nrow(gsea_sig_df) > 0) {
    plot_gsea_multiline(gsea_obj, gsea_sig_df,
      paste0(toupper(label), " risk-based GSEA top ", top_n_pathway, " pathways"),
      file.path(outdir, paste0("05.", label, ".multiline_gsea.pdf")),
      file.path(outdir, paste0("05.", label, ".multiline_gsea.png")),
      top_n_pathway, sort_by)
  }
}

# --- Summary ---
deg_up <- sum(res_df$log2FoldChange > 0, na.rm = TRUE)
deg_down <- sum(res_df$log2FoldChange < 0, na.rm = TRUE)
sig_summary <- sapply(names(gmts), function(label) {
  sig_file <- file.path(outdir, paste0("04.", label, ".gsea_significant.csv"))
  sig_count <- if (file.exists(sig_file)) nrow(read.csv(sig_file)) else 0
  paste0("  ", label, ": ", sig_count, " significant pathways")
})
write_summary(c(
  "GSEA risk-based NES analysis summary",
  paste0("Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("SN: ", sn),
  paste0("Count matrix: ", count_file),
  paste0("Risk file: ", risk_file),
  paste0("GMT dir: ", gmt_dir),
  paste0("GMT files: ", gmt_list_str),
  paste0("Labels: ", paste(names(gmts), collapse = ", ")),
  paste0("Common samples: ", ncol(count_mat)),
  paste0("Low risk N: ", grp_n["Low risk"]),
  paste0("High risk N: ", grp_n["High risk"]),
  paste0("Genes after filter: ", nrow(count_mat)),
  paste0("DEG up-regulated: ", deg_up),
  paste0("DEG down-regulated: ", deg_down),
  "",
  "Significant pathways per GMT:",
  sig_summary,
  "",
  "限制与说明：",
  "1. 本分析基于DESeq2差异分析和GSEA富集分析，高低风险组间比较。",
  "2. log2FC正值表示High risk组上调，NES阈值设为1，FDR阈值设为0.25。",
  "3. 不同阈值设置、多重检验校正方法和GMT数据库选择可能导致结果差异。",
  "4. 本结果为生物信息学预测，需进一步实验验证。"
))

log_message("GSEA risk-based NES analysis complete")
