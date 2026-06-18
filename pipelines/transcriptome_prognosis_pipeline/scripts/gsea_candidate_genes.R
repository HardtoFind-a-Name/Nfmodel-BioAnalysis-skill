args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggvenn)
  library(patchwork)
})

deg_all_file <- get_arg("--deg-all")
deg_sig_file <- get_arg("--deg-sig")
metascape_file <- get_arg("--metascape")
outdir <- get_arg("--outdir")
logdir <- get_arg("--logdir")
sn          <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

note_file <- file.path(logdir, "02.NGRs.summary.txt")

logsetup <- setup_logging(logdir, summary_dir, sn, "gsea_candidate_genes")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

write_note <- function(lines) {
  writeLines(as.character(lines), note_file)
}

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

count_leading_edge_genes <- function(x) {
  if (is.na(x) || x == "") {
    return(NA_real_)
  }
  length(unique(strsplit(x, "/", fixed = TRUE)[[1]]))
}

log_message("Script started: gsea_candidate_genes.R")

if (!file.exists(deg_all_file)) {
  save_note_and_stop(paste0("DEG all file does not exist: ", deg_all_file))
}
if (!file.exists(deg_sig_file)) {
  save_note_and_stop(paste0("DEG significant file does not exist: ", deg_sig_file))
}
if (!file.exists(metascape_file)) {
  save_note_and_stop(paste0("Metascape file does not exist: ", metascape_file))
}

deg_df <- read.csv(deg_all_file, stringsAsFactors = FALSE, check.names = FALSE)
deg_sig_df <- read.csv(deg_sig_file, stringsAsFactors = FALSE, check.names = FALSE)
meta_enrich_df <- read.csv(metascape_file, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(deg_df) == 0) {
  save_note_and_stop("DEG all file is empty; cannot perform GSEA.")
}

required_deg_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene_symbol", "padj_plot")
missing_deg_cols <- setdiff(required_deg_cols, colnames(deg_df))
if (length(missing_deg_cols) > 0) {
  save_note_and_stop(paste0("DEG all file is missing required columns: ", paste(missing_deg_cols, collapse = ", ")))
}

deg_df$gene_symbol <- trimws(as.character(deg_df$gene_symbol))
deg_df$log2FoldChange <- suppressWarnings(as.numeric(deg_df$log2FoldChange))
deg_df$padj <- suppressWarnings(as.numeric(deg_df$padj))
deg_df$pvalue <- suppressWarnings(as.numeric(deg_df$pvalue))
deg_df$stat <- suppressWarnings(as.numeric(deg_df$stat))
deg_df$baseMean <- suppressWarnings(as.numeric(deg_df$baseMean))
deg_df$lfcSE <- suppressWarnings(as.numeric(deg_df$lfcSE))
deg_df$padj_plot <- suppressWarnings(as.numeric(deg_df$padj_plot))

deg_df <- deg_df %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "", !is.na(log2FoldChange))

if (nrow(deg_df) < 20) {
  save_note_and_stop("Too few ranked genes (<20) for stable GSEA.")
}

rank_df <- deg_df %>%
  dplyr::select(gene = gene_symbol, rank_metric = log2FoldChange) %>%
  dplyr::filter(!is.na(rank_metric))

gene_map <- bitr(
  rank_df$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

if (is.null(gene_map) || nrow(gene_map) == 0) {
  save_note_and_stop("Failed to convert gene symbols to ENTREZ IDs for GSEA.")
}

gene_rank_df <- rank_df %>%
  dplyr::inner_join(gene_map, by = c("gene" = "SYMBOL")) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

if (nrow(gene_rank_df) < 20) {
  save_note_and_stop("Too few mapped ENTREZ genes (<20) for stable GSEA.")
}

gene_list <- gene_rank_df$rank_metric
names(gene_list) <- gene_rank_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_res <- as.data.frame(gsea_go@result)
if (is.null(gsea_res) || nrow(gsea_res) == 0) {
  gsea_res <- data.frame()
}

gsea_sig <- if (nrow(gsea_res) > 0) {
  gsea_res %>%
    dplyr::filter(abs(NES) > 1, p.adjust < 0.05)
} else {
  data.frame()
}

write.csv(gsea_res, file.path(outdir, "01.gsea_go_bp_all.csv"), row.names = FALSE)
write.csv(gsea_sig, file.path(outdir, "02.gsea_go_bp_sig.csv"), row.names = FALSE)

top_plot_df <- data.frame()

if (nrow(gsea_sig) > 0) {
  top_up_df <- gsea_sig %>%
    dplyr::arrange(desc(NES), p.adjust) %>%
    dplyr::slice_head(n = 5)

  top_down_df <- gsea_sig %>%
    dplyr::arrange(NES, p.adjust) %>%
    dplyr::slice_head(n = 5)

  top_plot_df <- dplyr::bind_rows(top_up_df, top_down_df) %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::arrange(desc(NES), p.adjust) %>%
    dplyr::mutate(
      leading_edge_gene_count = sapply(core_enrichment, count_leading_edge_genes),
      Description_wrapped = stringr::str_wrap(Description, width = 35)
    )

  top_plot_df$Description_wrapped <- factor(
    top_plot_df$Description_wrapped,
    levels = rev(top_plot_df$Description_wrapped)
  )

  blue_fun <- grDevices::colorRamp(c("#08306B", "#FFFFFF"))
  red_fun <- grDevices::colorRamp(c("#A50F15", "#FFFFFF"))

  fdr_min <- 0
  fdr_max <- 0.5
  top_plot_df$p.adjust_for_plot <- pmin(pmax(top_plot_df$p.adjust, fdr_min), fdr_max)
  top_plot_df$fdr_scaled <- (top_plot_df$p.adjust_for_plot - fdr_min) / (fdr_max - fdr_min)
  top_plot_df$fdr_scaled[top_plot_df$fdr_scaled < 0] <- 0
  top_plot_df$fdr_scaled[top_plot_df$fdr_scaled > 1] <- 1

  top_plot_df$point_fill <- vapply(
    seq_len(nrow(top_plot_df)),
    function(i) {
      rgb_mat <- if (top_plot_df$NES[i] < 0) {
        blue_fun(top_plot_df$fdr_scaled[i])
      } else {
        red_fun(top_plot_df$fdr_scaled[i])
      }
      grDevices::rgb(rgb_mat[1], rgb_mat[2], rgb_mat[3], maxColorValue = 255)
    },
    character(1)
  )

  x_lim <- max(abs(top_plot_df$NES), na.rm = TRUE) * 1.12

  legend_n <- 100
  blue_legend_pal <- colorRampPalette(c("#FFFFFF", "#08306B"))
  red_legend_pal <- colorRampPalette(c("#A50F15", "#FFFFFF"))

  legend_left_df <- data.frame(
    xmin = seq(0.08, 0.40, length.out = legend_n),
    xmax = seq(0.08, 0.40, length.out = legend_n) + (0.32 / legend_n),
    ymin = 0.08,
    ymax = 0.12,
    fill_col = blue_legend_pal(legend_n)
  )

  legend_right_df <- data.frame(
    xmin = seq(0.60, 0.92, length.out = legend_n),
    xmax = seq(0.60, 0.92, length.out = legend_n) + (0.32 / legend_n),
    ymin = 0.08,
    ymax = 0.12,
    fill_col = red_legend_pal(legend_n)
  )

  p_main <- ggplot(top_plot_df, aes(x = NES, y = Description_wrapped)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_point(
      aes(size = leading_edge_gene_count, fill = point_fill),
      shape = 21,
      color = "black",
      stroke = 1.0,
      alpha = 1
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      limits = c(-x_lim, x_lim),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_size_continuous(name = "Leading-edge gene count", range = c(3, 9)) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(face = "plain", lineheight = 0.95, vjust = 0.5, size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 10, r = 15, b = 5, l = 10)
    ) +
    guides(
      size = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        nrow = 1,
        byrow = TRUE,
        order = 1,
        override.aes = list(shape = 21, fill = "white", color = "black", stroke = 1)
      )
    ) +
    labs(
      title = "Top 5 Downregulated and Top 5 Upregulated GSEA Pathways",
      x = "NES",
      y = NULL
    )

  p_fdr_legend <- ggplot() +
    geom_rect(
      data = legend_left_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = legend_left_df$fill_col,
      color = NA
    ) +
    geom_rect(
      data = legend_right_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = legend_right_df$fill_col,
      color = NA
    ) +
    annotate("text", x = 0.50, y = 0.18, label = "FDR", size = 3.8, fontface = "bold") +
    annotate("text", x = 0.08, y = 0.02, label = "0.50", hjust = 0, size = 3.2) +
    annotate("text", x = 0.50, y = 0.02, label = "0.00", hjust = 0.5, size = 3.2) +
    annotate("text", x = 0.92, y = 0.02, label = "0.50", hjust = 1, size = 3.2) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.22), clip = "off") +
    theme_void()

  p_final <- p_main / p_fdr_legend + plot_layout(heights = c(12, 1.8))

  save_pdf_plot(p_final, file.path(outdir, "03.gsea_go_bp_sig_dotplot.pdf"), width = 6.5, height = 7)
  ggsave(file.path(outdir, "04.gsea_go_bp_sig_dotplot.png"), p_final, width = 6.5, height = 7, dpi = 300)

  write.csv(
    top_plot_df %>%
      dplyr::select(
        ID,
        Description,
        Description_wrapped,
        NES,
        p.adjust,
        setSize,
        core_enrichment,
        leading_edge_gene_count,
        point_fill
      ),
    file.path(outdir, "03.gsea_go_bp_top5_up_top5_down_for_plot.csv"),
    row.names = FALSE
  )
}

required_meta_cols <- c("Description", "LogP", "Hits")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_enrich_df))
if (length(missing_meta_cols) > 0) {
  save_note_and_stop(paste0("Metascape file is missing required columns: ", paste(missing_meta_cols, collapse = ", ")))
}

neuroimmune_keywords <- c(
  "neuro", "neuron", "synap", "axon", "dendrite", "glia", "myelin",
  "immune", "immun", "inflamm", "cytokine", "chemokine", "interferon",
  "leukocyte", "lymphocyte", "macrophage", "antigen", "innate immune", "adaptive immune"
)
pattern_kw <- paste(neuroimmune_keywords, collapse = "|")
has_logq_col <- "Log(q-value)" %in% colnames(meta_enrich_df)

meta_selected_df <- meta_enrich_df %>%
  dplyr::mutate(
    Description = as.character(Description),
    LogP = suppressWarnings(as.numeric(LogP)),
    `Log(q-value)` = if (has_logq_col) suppressWarnings(as.numeric(`Log(q-value)`)) else NA_real_,
    Hits = as.character(Hits)
  ) %>%
  dplyr::filter(!is.na(Description), grepl(pattern_kw, Description, ignore.case = TRUE)) %>%
  dplyr::select(dplyr::any_of(c("Description", "LogP", "Log(q-value)", "Hits")))

metascape_genes <- meta_selected_df %>%
  tidyr::separate_rows(Hits, sep = "\\|") %>%
  dplyr::mutate(Hits = trimws(Hits)) %>%
  dplyr::filter(!is.na(Hits), Hits != "") %>%
  dplyr::pull(Hits) %>%
  unique()

write.csv(
  meta_selected_df,
  file.path(outdir, "05.metascape_neuroimmune_selected_terms.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene_symbol = metascape_genes),
  file.path(outdir, "06.metascape_neuroimmune_gene_list.csv"),
  row.names = FALSE
)

if (nrow(meta_selected_df) > 0) {
  meta_plot_df <- meta_selected_df %>%
    dplyr::filter(!is.na(LogP)) %>%
    dplyr::mutate(minus_LogP = -LogP) %>%
    dplyr::arrange(desc(minus_LogP)) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(Description, width = 40),
      Description = factor(Description, levels = rev(unique(Description)))
    )

  p_meta_bar <- ggplot(meta_plot_df, aes(x = minus_LogP, y = Description)) +
    geom_col(fill = "#4C78A8", width = 0.75) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10, lineheight = 0.95),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(face = "bold", size = 11),
      axis.title.y = element_blank()
    ) +
    labs(
      title = "Metascape neuroimmune-related enrichment",
      x = "-LogP",
      y = NULL
    )

  save_pdf_plot(p_meta_bar, file.path(outdir, "07.metascape_neuroimmune_barplot.pdf"), width = 8, height = 5.5)
  ggsave(file.path(outdir, "08.metascape_neuroimmune_barplot.png"), p_meta_bar, width = 8, height = 5.5, dpi = 300)
}

if (nrow(gsea_sig) > 0) {
  gsea_leading_edge_df <- gsea_sig %>%
    dplyr::filter(!is.na(core_enrichment), core_enrichment != "") %>%
    dplyr::select(ID, Description, NES, p.adjust, core_enrichment) %>%
    tidyr::separate_rows(core_enrichment, sep = "/") %>%
    dplyr::rename(ENTREZID = core_enrichment)
} else {
  gsea_leading_edge_df <- data.frame()
}

gsea_leading_edge_genes <- character()

if (nrow(gsea_leading_edge_df) > 0) {
  entrez_to_symbol <- bitr(
    unique(gsea_leading_edge_df$ENTREZID),
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db
  )

  if (!is.null(entrez_to_symbol) && nrow(entrez_to_symbol) > 0) {
    gsea_leading_edge_symbol_df <- gsea_leading_edge_df %>%
      dplyr::inner_join(entrez_to_symbol, by = "ENTREZID")

    gsea_leading_edge_genes <- gsea_leading_edge_symbol_df %>%
      dplyr::pull(SYMBOL) %>%
      unique()
  }
}

write.csv(
  data.frame(gene_symbol = gsea_leading_edge_genes),
  file.path(outdir, "09.gsea_leading_edge_gene_set.csv"),
  row.names = FALSE
)

deg_sig_genes <- if ("gene_symbol" %in% colnames(deg_sig_df)) {
  unique(trimws(as.character(deg_sig_df$gene_symbol)))
} else {
  character()
}
deg_sig_genes <- deg_sig_genes[!is.na(deg_sig_genes) & deg_sig_genes != ""]

candidate_genes <- intersect(intersect(gsea_leading_edge_genes, metascape_genes), deg_sig_genes)

write.csv(
  data.frame(gene_symbol = candidate_genes),
  file.path(outdir, "10.NRGs_from_GSEA_leading_edge_and_Metascape.csv"),
  row.names = FALSE
)

venn_list <- list(
  Metascape = metascape_genes,
  GSEA_leading_edge = gsea_leading_edge_genes,
  DEGs = deg_sig_genes
)

if (sum(lengths(venn_list) > 0) >= 2) {
  p_venn <- ggvenn(
    venn_list,
    names(venn_list),
    fill_color = c("#E64B35", "#4DBBD5", "#ff00ff"),
    show_percentage = TRUE,
    stroke_alpha = 0.5,
    stroke_size = 0.3,
    text_size = 4,
    stroke_color = "white",
    stroke_linetype = "solid",
    set_name_color = c("#E64B35", "#4DBBD5", "#ff00ff"),
    set_name_size = 5,
    text_color = "black"
  )

  save_pdf_plot(p_venn, file.path(outdir, "11.NRGs_venn.pdf"), width = 6, height = 6)
  ggsave(file.path(outdir, "12.NRGs_venn.png"), plot = p_venn, width = 5, height = 5, dpi = 600)
}

note_lines <- c(
  "=== GSEA candidate genes ===",
  paste0("DEGs in ranking: ", nrow(rank_df)),
  paste0("Genes mapped to ENTREZ: ", nrow(gene_rank_df)),
  paste0("Significant GSEA pathways (|NES|>1, FDR<0.05): ", nrow(gsea_sig)),
  paste0("Metascape neuroimmune terms: ", nrow(meta_selected_df)),
  paste0("Metascape genes: ", length(metascape_genes)),
  paste0("GSEA leading-edge genes: ", length(gsea_leading_edge_genes)),
  paste0("DEG significant genes: ", length(deg_sig_genes)),
  paste0("Candidate genes (intersection): ", length(candidate_genes)),
  paste0("Candidate gene list: ", if (length(candidate_genes) > 0) paste(candidate_genes, collapse = ", ") else "None"),
  "",
  "=== 限制与说明 ===",
  "GSEA uses log2FC as ranking metric (pre-ranked mode).",
  "Leading-edge genes are extracted from core_enrichment of significant pathways.",
  "Metascape terms filtered using neuro/immune keyword pattern.",
  "Candidates = intersect(GSEA_leading_edge, Metascape, DEG_sig).",
  "Venn diagram shown for 2+ gene sets with overlap."
)

write_note(note_lines)
write_summary(note_lines)
log_message("GSEA candidate gene step finished.")
