set.seed(123)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
})

# --- CLI args ---
args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

risk_file    <- get_arg("--risk-file", required = FALSE, default = NULL)  # kept for compatibility, unused
gene_file    <- get_arg("--gene-file")
expr_file    <- get_arg("--expr")
train_id     <- get_arg("--train-id")
outdir       <- get_arg("--outdir")
logdir       <- get_arg("--logdir")
sn           <- get_arg("--sn")
summary_dir  <- get_arg("--summary-dir")
gmt_dir      <- get_arg("--gmt-dir")
gmt_list_str <- get_arg("--gmt-gsea_gene")
sort_by      <- get_arg("--gsea-sort-by", required = FALSE, default = "p.adjust")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

top_n_pathway <- 10

# --- Logging ---
logsetup <- setup_logging(logdir, summary_dir, sn, "gsea_gene")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

# --- Read expression matrix ---
log_message("Reading expression matrix")
if (!file.exists(expr_file)) stop(paste("Expression matrix not found:", expr_file))
expr_df <- fread(expr_file, data.table = FALSE, check.names = FALSE)
gene_col <- colnames(expr_df)[1]
expr_df[[gene_col]] <- toupper(trimws(expr_df[[gene_col]]))
expr_df <- expr_df[expr_df[[gene_col]] != "" & !is.na(expr_df[[gene_col]]), ]
expr_df <- expr_df[!duplicated(expr_df[[gene_col]]), ]
rownames(expr_df) <- expr_df[[gene_col]]
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
mode(expr_mat) <- "numeric"
keep_genes <- rowSums(expr_mat > 0) >= max(3, floor(ncol(expr_mat) * 0.1))
expr_mat <- expr_mat[keep_genes, , drop = FALSE]
if (nrow(expr_mat) < 1000) stop("Too few genes after low-expression filter")

# --- Extract target genes from gene selection file ---
log_message("Extracting target genes from: ", gene_file)
if (!file.exists(gene_file)) stop(paste("Gene file not found:", gene_file))
gene_df <- read.csv(gene_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene" %in% colnames(gene_df)) stop("Gene file must have a 'gene' column")
target_genes <- unique(trimws(as.character(gene_df$gene)))
target_genes <- target_genes[!is.na(target_genes) & nzchar(target_genes)]
target_genes <- toupper(target_genes)
available_targets <- intersect(target_genes, rownames(expr_mat))
if (length(available_targets) == 0) stop("None of the target genes found in expression matrix")
log_message("Available target genes: ", paste(available_targets, collapse = ", "))

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

# --- Plot functions ---
calc_spearman_rank <- function(target_gene, expr_mat) {
  target_exp <- as.numeric(expr_mat[target_gene, ])
  cor_df <- lapply(rownames(expr_mat), function(g) {
    gexp <- as.numeric(expr_mat[g, ])
    ct <- suppressWarnings(cor.test(target_exp, gexp, method = "spearman", exact = FALSE))
    data.frame(targetGene = target_gene, GeneSymbol = g,
               rho = unname(ct$estimate), pvalue = ct$p.value, stringsAsFactors = FALSE)
  }) %>% bind_rows()
  cor_df <- cor_df %>% filter(!is.na(rho)) %>% arrange(desc(rho))
  geneList <- cor_df$rho; names(geneList) <- cor_df$GeneSymbol
  geneList <- sort(geneList, decreasing = TRUE)
  list(cor_df = cor_df, geneList = geneList)
}

plot_ridge <- function(gsea_obj, title_text, out_pdf, out_png, top_n = 10) {
  if (is.null(gsea_obj)) return(NULL)
  gsea_df <- as.data.frame(gsea_obj)
  if (is.null(gsea_df) || nrow(gsea_df) == 0) return(NULL)
  gsea_df <- gsea_df %>% filter(!is.na(p.adjust)) %>% arrange(p.adjust)
  if (nrow(gsea_df) == 0) return(NULL)
  n_show <- min(top_n, nrow(gsea_df))
  p <- enrichplot::ridgeplot(gsea_obj, showCategory = n_show, fill = "NES") +
    labs(title = title_text, x = "log2 Fold Change", y = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          axis.text.y = element_text(size = 10))
  ggsave(out_pdf, p, width = 10, height = 7.5)
  ggsave(out_png, p, width = 10, height = 7.5, dpi = 320)
}

plot_gsea_multiline <- function(gsea_obj, gsea_df, title_text, out_pdf, out_png, top_n = 10, sort_by = "p.adjust") {
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
  ggsave(out_pdf, p, width = 10.5, height = 8.2)
  ggsave(out_png, p, width = 10.5, height = 8.2, dpi = 320)
}

# =========================================================
# Per-gene loop
# =========================================================
all_cor_list    <- list()
all_sig_lists   <- setNames(vector("list", length(gmts)), names(gmts))
all_gsea_objs   <- list()   # list of lists: all_gsea_objs[[tg]][[label]] = gsea_obj
summary_lines   <- c()

for (tg in available_targets) {
  log_message("Processing target gene: ", tg)
  rank_res <- calc_spearman_rank(tg, expr_mat)
  cor_df   <- rank_res$cor_df
  geneList <- rank_res$geneList

  fwrite(cor_df, file = file.path(outdir, paste0("01.", tg, ".spearman_all_genes.csv")))
  fwrite(data.frame(targetGene = tg, GeneSymbol = names(geneList), rho = as.numeric(geneList)),
         file = file.path(outdir, paste0("02.", tg, ".ranked_geneList_rho.csv")))

  all_cor_list[[tg]] <- cor_df
  all_gsea_objs[[tg]] <- list()

  for (gmt_item in gmts) {
    label <- gmt_item$label
    log_message("  GSEA: ", tg, " vs ", label)

    gsea_obj <- GSEA(geneList = geneList, TERM2GENE = gmt_item$gmt,
                     pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
    gsea_res <- as.data.frame(gsea_obj)
    all_gsea_objs[[tg]][[label]] <- gsea_obj

    if (nrow(gsea_res) > 0) {
      gsea_res <- gsea_res %>% filter(p.adjust < 0.05)
      if (nrow(gsea_res) > 0) {
        gsea_res$targetGene <- tg
        all_sig_lists[[label]][[tg]] <- gsea_res

        fwrite(gsea_res, file = file.path(outdir,
          paste0("03.", tg, ".", label, ".gsea_significant.csv")))

        plot_gsea_multiline(gsea_obj, gsea_res,
          paste0(tg, " - ", toupper(label), " GSEA top ", top_n_pathway),
          file.path(outdir, paste0("04.", tg, ".", label, ".multiline_gsea.pdf")),
          file.path(outdir, paste0("04.", tg, ".", label, ".multiline_gsea.png")),
          top_n_pathway, sort_by)

        plot_ridge(gsea_obj,
          paste0(tg, " - ", toupper(label), " ridgeplot top ", top_n_pathway),
          file.path(outdir, paste0("05.", tg, ".", label, ".ridgeplot.pdf")),
          file.path(outdir, paste0("05.", tg, ".", label, ".ridgeplot.png")),
          top_n_pathway)
      }
    }

    n_sig <- if (is.null(all_sig_lists[[label]][[tg]])) 0 else nrow(all_sig_lists[[label]][[tg]])
    summary_lines <- c(summary_lines,
      paste0("Target: ", tg, " | ", label, " significant: ", n_sig))
  }
}

# =========================================================
# Merged outputs
# =========================================================
all_cor_df <- bind_rows(all_cor_list)
fwrite(all_cor_df, file = file.path(outdir, "06.all_genes.spearman_all_genes.csv"))

all_merged_list <- list()
for (label in names(gmts)) {
  if (length(all_sig_lists[[label]]) > 0) {
    sig_df <- bind_rows(all_sig_lists[[label]])
    fwrite(sig_df, file = file.path(outdir,
      paste0("07.all_genes.", label, ".gsea_significant.csv")))
    all_merged_list[[label]] <- sig_df
  } else {
    fwrite(data.frame(), file = file.path(outdir,
      paste0("07.all_genes.", label, ".gsea_significant.csv")))
  }
}

if (length(all_merged_list) > 0) {
  all_merged <- bind_rows(all_merged_list)
  fwrite(all_merged, file = file.path(outdir, "08.all_genes.gsea_merged.csv"))
} else {
  fwrite(data.frame(), file = file.path(outdir, "08.all_genes.gsea_merged.csv"))
}

# =========================================================
# Summary multi-line + ridge plots (from merged table)
# =========================================================
for (label in names(gmts)) {
  if (length(all_sig_lists[[label]]) == 0) next
  merged_sig <- bind_rows(all_sig_lists[[label]])
  if (nrow(merged_sig) == 0) next

  # Best per pathway (by chosen sort criteria)
  best_per_pathway <- merged_sig %>%
    group_by(ID) %>%
    summarise(
      best_p_adjust = min(p.adjust, na.rm = TRUE),
      best_pvalue   = min(pvalue, na.rm = TRUE),
      best_NES      = NES[which.max(abs(NES))],
      best_gene = targetGene[if (sort_by == "NES") which.max(abs(NES)) else which.min(p.adjust)],
      Description = dplyr::first(Description),
      .groups = "drop") %>%
    arrange(if (sort_by == "NES") desc(abs(best_NES))
            else if (sort_by == "pvalue") best_pvalue
            else best_p_adjust) %>%
    slice(1:min(top_n_pathway, n()))

  if (nrow(best_per_pathway) == 0) next

  # Use the GSEA object from the best gene for each pathway
  # Build a gsdata frame using gsInfo from each pathway's best gene
  gsdata_list <- lapply(seq_len(nrow(best_per_pathway)), function(i) {
    pid   <- best_per_pathway$ID[i]
    pdesc <- best_per_pathway$Description[i]
    bgene <- best_per_pathway$best_gene[i]
    gobj  <- all_gsea_objs[[bgene]][[label]]
    tmp   <- enrichplot:::gsInfo(gobj, geneSetID = pid)
    tmp$ID <- pid; tmp$Description <- pdesc; tmp
  })
  gsdata <- do.call(rbind, gsdata_list)
  gsdata$Description <- factor(gsdata$Description, levels = best_per_pathway$Description)
  x_min <- min(gsdata$x, na.rm = TRUE); x_max <- max(gsdata$x, na.rm = TRUE)

  p1 <- ggplot(gsdata, aes(x = x, y = runningScore, color = Description)) +
    geom_line(linewidth = 0.9, alpha = 0.95) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    labs(title = paste0("Consensus ", toupper(label), " GSEA top ", top_n_pathway, " pathways"),
         x = NULL, y = "Running Enrichment Score") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "right", legend.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), plot.margin = margin(5.5, 5.5, 0, 5.5))

  hit_df <- gsdata %>% filter(position == 1)
  if (nrow(hit_df) > 0) {
    hit_df$Description <- factor(hit_df$Description, levels = rev(best_per_pathway$Description))
    p2 <- ggplot(hit_df, aes(x = x, y = Description, color = Description)) +
      geom_linerange(aes(ymin = as.numeric(Description) - 0.35,
                         ymax = as.numeric(Description) + 0.35), linewidth = 0.45) +
      scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
      labs(x = NULL, y = NULL) + theme_classic(base_size = 12) +
      theme(legend.position = "none", axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), axis.line.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            plot.margin = margin(0, 5.5, 0, 5.5))
  } else {
    p2 <- ggplot() + theme_void()
  }

  rank_df <- gsdata %>% distinct(x, geneList)
  p3 <- ggplot(rank_df, aes(x = x, y = geneList)) +
    geom_area(fill = "grey70") +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    labs(x = "Rank in Ordered Dataset", y = "Ranked List Metric") +
    theme_classic(base_size = 12) + theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

  p <- p1 / p2 / p3 + patchwork::plot_layout(ncol = 1, heights = c(4.2, 1.2, 1.8))
  ggsave(file.path(outdir, paste0("09.", label, ".multiline_gsea.pdf")), p, width = 10.5, height = 8.2)
  ggsave(file.path(outdir, paste0("09.", label, ".multiline_gsea.png")), p, width = 10.5, height = 8.2, dpi = 320)

  # Summary ridge: use gsea_obj from the gene with most significant pathways
  best_gene_for_label <- merged_sig %>%
    group_by(targetGene) %>% summarise(n = n(), min_p = min(p.adjust), .groups = "drop") %>%
    arrange(min_p) %>% slice(1) %>% pull(targetGene)
  gobj_summary <- all_gsea_objs[[best_gene_for_label]][[label]]
  if (!is.null(gobj_summary)) {
    plot_ridge(gobj_summary,
      paste0("Consensus ", toupper(label), " ridgeplot top ", top_n_pathway),
      file.path(outdir, paste0("10.", label, ".ridgeplot.pdf")),
      file.path(outdir, paste0("10.", label, ".ridgeplot.png")),
      top_n_pathway)
  }
}

# =========================================================
# Intersection NES heatmap
# =========================================================
if (length(all_merged_list) > 0) {
  all_sig <- bind_rows(all_merged_list)
  if (nrow(all_sig) > 0) {
    pathway_counts <- table(all_sig$ID)
    shared_ids <- names(pathway_counts[pathway_counts >= 2])
    if (length(shared_ids) > 0) {
      top_shared <- names(sort(pathway_counts[shared_ids], decreasing = TRUE))
      top_shared <- head(top_shared, 10)

      nes_mat <- all_sig %>%
        filter(ID %in% top_shared) %>%
        select(ID, targetGene, NES) %>%
        pivot_wider(names_from = targetGene, values_from = NES, values_fill = 0)
      nes_mat <- as.data.frame(nes_mat)
      rownames(nes_mat) <- nes_mat$ID; nes_mat$ID <- NULL
      nes_mat <- as.matrix(nes_mat)

      if (nrow(nes_mat) > 0 && ncol(nes_mat) > 0) {
        nes_limit <- max(2, ceiling(max(abs(nes_mat), na.rm = TRUE)))
        pdf(file.path(outdir, "11.intersection_nes_heatmap.pdf"), width = 10, height = 8)
        pheatmap(nes_mat, cluster_rows = TRUE, cluster_cols = TRUE,
          color = colorRampPalette(c("blue", "white", "red"))(100),
          breaks = seq(-nes_limit, nes_limit, length.out = 101),
          main = "Intersection Pathway NES Heatmap", fontsize = 10, display_numbers = FALSE)
        dev.off()
        png(file.path(outdir, "11.intersection_nes_heatmap.png"), width = 10, height = 8,
            units = "in", res = 320)
        pheatmap(nes_mat, cluster_rows = TRUE, cluster_cols = TRUE,
          color = colorRampPalette(c("blue", "white", "red"))(100),
          breaks = seq(-nes_limit, nes_limit, length.out = 101),
          main = "Intersection Pathway NES Heatmap", fontsize = 10, display_numbers = FALSE)
        dev.off()
        fwrite(all_sig %>% filter(ID %in% top_shared),
               file = file.path(outdir, "12.gene_pathway_intersection.csv"))
      }
    }
  }
}

# --- Summary ---
sig_summary <- sapply(names(gmts), function(label) {
  sig_count <- if (!is.null(all_sig_lists[[label]])) {
    sum(sapply(all_sig_lists[[label]], nrow))
  } else 0
  paste0("  ", label, ": ", sig_count, " significant pathways")
})
top_pathways <- sapply(names(gmts), function(label) {
  if (!is.null(all_sig_lists[[label]]) && length(all_sig_lists[[label]]) > 0) {
    merged <- bind_rows(all_sig_lists[[label]])
    top <- merged %>% arrange(p.adjust) %>% slice_head(n = 3)
    paste0("  ", label, ": ", paste(top$Description, collapse = "; "))
  } else ""
})
write_summary(c(
  "GSEA gene-based analysis summary",
  paste0("Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("SN: ", sn),
  paste0("Gene file: ", gene_file),
  paste0("Expression: ", expr_file),
  paste0("GMT dir: ", gmt_dir),
  paste0("GMT files: ", gmt_list_str),
  paste0("Labels: ", paste(names(gmts), collapse = ", ")),
  paste0("Total target genes: ", length(available_targets)),
  paste0("Target genes: ", paste(available_targets, collapse = ", ")),
  paste0("Genes in expression matrix: ", nrow(expr_mat)),
  paste0("Samples: ", ncol(expr_mat)),
  "",
  "Significant pathways per GMT:",
  sig_summary,
  "",
  "Top pathway descriptions:",
  top_pathways,
  "",
  "限制与说明：",
  "1. 本分析基于候选基因与MSigDB通路数据库，结果受样本量和表达矩阵过滤参数影响。",
  "2. GSEA分析依赖于Spearman相关系数排序，不同排序策略可能导致结果差异。",
  "3. 通路显著性经BH多重检验校正，建议结合NES与p.adjust综合判断。",
  "4. 本结果为生物信息学预测，需进一步的实验验证。"
))

log_message("GSEA gene-based analysis complete")
