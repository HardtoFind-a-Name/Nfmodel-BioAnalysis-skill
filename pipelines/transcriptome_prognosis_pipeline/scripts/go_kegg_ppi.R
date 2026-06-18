args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(igraph)
  library(ggraph)
})

candidate_gene_file <- get_arg("--candidate-genes")
string_tsv <- get_arg("--string-tsv")
outdir <- get_arg("--outdir")
logdir <- get_arg("--logdir")
sn          <- get_arg("--sn")
summary_dir <- get_arg("--summary-dir")
go_top_n_each <- as.integer(get_arg("--go-top-n-each", required = FALSE, default = "10"))
kegg_top_n <- as.integer(get_arg("--kegg-top-n", required = FALSE, default = "10"))
adj_cutoff <- as.numeric(get_arg("--adj-cutoff", required = FALSE, default = "0.05"))
fallback_p <- as.numeric(get_arg("--fallback-p", required = FALSE, default = "0.05"))
kegg_gmt <- get_arg("--kegg-gmt", required = FALSE, default = NULL)
kegg_gmt_id_type <- toupper(get_arg("--kegg-gmt-id-type", required = FALSE, default = "SYMBOL"))
ppi_score_cutoff <- as.numeric(get_arg("--ppi-score-cutoff", required = FALSE, default = "0.4"))
ppi_hub_top_n <- as.integer(get_arg("--ppi-hub-top-n", required = FALSE, default = "10"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

logsetup <- setup_logging(logdir, summary_dir, sn, "go_kegg_ppi")
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

log_message("Script started: go_kegg_ppi.R")

if (!file.exists(candidate_gene_file)) save_note_and_stop(paste0("Candidate gene file does not exist: ", candidate_gene_file))
if (!file.exists(string_tsv)) save_note_and_stop(paste0("STRING TSV file does not exist: ", string_tsv))

gene_df <- read.csv(candidate_gene_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(gene_df) == 0) save_note_and_stop("Candidate gene table is empty.")

if (!("gene_symbol" %in% colnames(gene_df))) {
  colnames(gene_df)[1] <- "gene_symbol"
}

gene_symbols <- unique(trimws(as.character(gene_df$gene_symbol)))
gene_symbols <- gene_symbols[gene_symbols != "" & !is.na(gene_symbols)]
if (length(gene_symbols) == 0) save_note_and_stop("No usable gene symbols found in candidate gene file.")

gene_convert <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = c("ENTREZID", "ENSEMBL"),
  OrgDb = org.Hs.eg.db
)
gene_convert <- gene_convert[!duplicated(gene_convert$SYMBOL), , drop = FALSE]
entrez_ids <- unique(gene_convert$ENTREZID)

write.csv(gene_convert, file = file.path(outdir, "01_gene_id_conversion.csv"), row.names = FALSE)

if (length(entrez_ids) == 0) {
  save_note_and_stop("No ENTREZ IDs were converted successfully from candidate genes.")
}

log_message("Input genes: ", length(gene_symbols))
log_message("Converted ENTREZ IDs: ", length(entrez_ids))

go_colors <- c(BP = "#4DBBD5", CC = "#E64B35", MF = "#00A087")

ego <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

go_res <- as.data.frame(ego)
if (is.null(go_res) || nrow(go_res) == 0) {
  go_res <- data.frame()
  go_sig <- data.frame()
  go_used_cutoff <- "no GO terms returned"
} else {
  go_res$gene_count_numeric <- go_res$Count
  go_res$neg_log10_padj <- -log10(go_res$p.adjust)
  go_sig <- go_res %>% dplyr::filter(p.adjust < adj_cutoff)
  go_used_cutoff <- paste0("p.adjust < ", adj_cutoff)
  if (nrow(go_sig) == 0) {
    go_sig <- go_res %>% dplyr::filter(pvalue < fallback_p)
    go_used_cutoff <- paste0("pvalue < ", fallback_p, " (fallback)")
  }
}

write.csv(go_res, file = file.path(outdir, "02_GO_all_results.csv"), row.names = FALSE)
write.csv(go_sig, file = file.path(outdir, "03_GO_significant_results.csv"), row.names = FALSE)

if (nrow(go_sig) > 0) {
  go_plot_df <- go_sig %>%
    dplyr::filter(ONTOLOGY %in% c("BP", "CC", "MF")) %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::arrange(p.adjust, pvalue, .by_group = TRUE) %>%
    dplyr::slice_head(n = go_top_n_each) %>%
    dplyr::ungroup()

  if (nrow(go_plot_df) > 0) {
    go_plot_df$Description <- stringr::str_wrap(go_plot_df$Description, width = 45)
    go_plot_df$ONTOLOGY <- factor(go_plot_df$ONTOLOGY, levels = c("BP", "CC", "MF"))
    go_plot_df <- go_plot_df %>% dplyr::arrange(ONTOLOGY, Count)
    go_plot_df$Description <- factor(go_plot_df$Description, levels = go_plot_df$Description)

    p_go <- ggplot(go_plot_df, aes(x = Count, y = Description, fill = ONTOLOGY)) +
      geom_bar(stat = "identity", width = 0.75) +
      scale_fill_manual(
        values = go_colors,
        labels = c(BP = "Biological Process", CC = "Cellular Component", MF = "Molecular Function")
      ) +
      labs(title = "GO enrichment analysis", x = "Gene count", y = NULL, fill = "Ontology") +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = "top",
        panel.grid.major.y = element_blank()
      )

    save_pdf_plot(p_go, file.path(outdir, "04.GO_barplot.pdf"), width = 11, height = 9)
    ggsave(file.path(outdir, "04.GO_barplot.png"), plot = p_go, width = 11, height = 9, dpi = 600)
  }
}

kegg_ok <- TRUE
kegg_source <- "online_KEGG"
ekk <- NULL

if (!is.null(kegg_gmt) && nzchar(kegg_gmt) && file.exists(kegg_gmt)) {
  kegg_source <- paste0("local_GMT: ", kegg_gmt)
  gmt_df <- tryCatch(
    read.gmt(kegg_gmt),
    error = function(e) save_note_and_stop(paste0("Failed to read KEGG GMT file: ", conditionMessage(e)))
  )

  if (nrow(gmt_df) == 0) {
    save_note_and_stop(paste0("KEGG GMT file is empty: ", kegg_gmt))
  }

  if (!all(c("term", "gene") %in% colnames(gmt_df))) {
    save_note_and_stop("KEGG GMT parse result missing required columns 'term' and 'gene'.")
  }

  if (kegg_gmt_id_type == "SYMBOL") {
    kegg_input_genes <- gene_symbols
  } else if (kegg_gmt_id_type == "ENTREZID") {
    kegg_input_genes <- as.character(entrez_ids)
  } else {
    save_note_and_stop(paste0("Unsupported --kegg-gmt-id-type: ", kegg_gmt_id_type, ". Use SYMBOL or ENTREZID."))
  }

  ekk <- tryCatch(
    enricher(
      gene = unique(kegg_input_genes),
      TERM2GENE = gmt_df,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) {
      kegg_ok <<- FALSE
      log_message("KEGG GMT enrichment skipped: ", conditionMessage(e))
      NULL
    }
  )
} else {
  ekk <- tryCatch(
    enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) {
      kegg_ok <<- FALSE
      log_message("KEGG enrichment skipped: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0) {
    ekk <- tryCatch(
      setReadable(ekk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) {
        log_message("setReadable for KEGG failed, keeping original KEGG result: ", conditionMessage(e))
        ekk
      }
    )
  }
}

kegg_res <- if (!is.null(ekk)) as.data.frame(ekk) else data.frame()
if (is.null(kegg_res) || nrow(kegg_res) == 0) {
  kegg_res <- data.frame()
  kegg_sig <- data.frame()
  kegg_used_cutoff <- if (kegg_ok) "no KEGG terms returned" else "KEGG unavailable (network or remote service failure)"
} else {
  kegg_res$neg_log10_padj <- -log10(kegg_res$p.adjust)
  kegg_sig <- kegg_res %>% dplyr::filter(p.adjust < adj_cutoff)
  kegg_used_cutoff <- paste0("p.adjust < ", adj_cutoff)
  if (nrow(kegg_sig) < 3) {
    kegg_sig <- kegg_res %>% dplyr::filter(pvalue < fallback_p)
    kegg_used_cutoff <- paste0("pvalue < ", fallback_p, " (fallback)")
  }
}

write.csv(kegg_res, file = file.path(outdir, "06_KEGG_all_results.csv"), row.names = FALSE)
write.csv(kegg_sig, file = file.path(outdir, "07_KEGG_significant_results.csv"), row.names = FALSE)

if (nrow(kegg_sig) > 0) {
  kegg_plot_df <- kegg_sig %>%
    dplyr::arrange(desc(neg_log10_padj), desc(Count)) %>%
    dplyr::slice_head(n = kegg_top_n)

  if (nrow(kegg_plot_df) > 0) {
    kegg_plot_df$Description <- stringr::str_wrap(kegg_plot_df$Description, width = 45)
    kegg_plot_df <- kegg_plot_df %>% dplyr::arrange(neg_log10_padj)
    kegg_plot_df$Description <- factor(kegg_plot_df$Description, levels = kegg_plot_df$Description)

    p_kegg <- ggplot(kegg_plot_df, aes(x = neg_log10_padj, y = Description)) +
      geom_bar(stat = "identity", width = 0.75, fill = "#3C5488") +
      labs(title = "KEGG enrichment analysis", x = "-log10(adj.P value)", y = NULL) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        panel.grid.major.y = element_blank()
      )

    save_pdf_plot(p_kegg, file.path(outdir, "08.KEGG_barplot.pdf"), width = 10, height = 7)
    ggsave(file.path(outdir, "08.KEGG_barplot.png"), plot = p_kegg, width = 10, height = 7, dpi = 600)
  }
}

ppi_raw <- fread(string_tsv, sep = "\t", header = TRUE, data.table = FALSE)
if (nrow(ppi_raw) == 0) save_note_and_stop("STRING interaction table is empty.")
colnames(ppi_raw)[1] <- "node1"

required_ppi_cols <- c("node1", "node2", "combined_score")
missing_ppi_cols <- setdiff(required_ppi_cols, colnames(ppi_raw))
if (length(missing_ppi_cols) > 0) {
  save_note_and_stop(paste0("STRING interaction table missing required columns: ", paste(missing_ppi_cols, collapse = ", ")))
}

ppi_df <- ppi_raw %>%
  dplyr::select(node1, node2, combined_score) %>%
  dplyr::filter(node1 != node2) %>%
  dplyr::filter(combined_score >= ppi_score_cutoff) %>%
  distinct()

if (nrow(ppi_df) == 0) {
  save_note_and_stop("No PPI edges retained after score filtering.")
}

g <- igraph::graph_from_data_frame(ppi_df, directed = FALSE)
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

node_metrics <- data.frame(
  gene = igraph::V(g)$name,
  degree = igraph::degree(g, mode = "all"),
  betweenness = igraph::betweenness(g, directed = FALSE, normalized = TRUE),
  closeness = igraph::closeness(g, normalized = TRUE),
  eigenvector = igraph::eigen_centrality(g)$vector,
  stringsAsFactors = FALSE
) %>% dplyr::arrange(desc(degree), desc(betweenness), desc(closeness))

write.csv(node_metrics, file = file.path(outdir, "10_PPI_node_metrics.csv"), row.names = FALSE)

hub_genes <- node_metrics %>% dplyr::slice_head(n = ppi_hub_top_n)
igraph::V(g)$degree <- node_metrics$degree[match(igraph::V(g)$name, node_metrics$gene)]
igraph::V(g)$hub_type <- ifelse(igraph::V(g)$name %in% hub_genes$gene, "Hub gene", "Other gene")

edge_df_graph <- igraph::as_data_frame(g, what = "edges")
colnames(edge_df_graph)[1:2] <- c("from", "to")
ppi_key1 <- paste(ppi_df$node1, ppi_df$node2, sep = "||")
ppi_key2 <- paste(ppi_df$node2, ppi_df$node1, sep = "||")
edge_key <- paste(edge_df_graph$from, edge_df_graph$to, sep = "||")
edge_df_graph$combined_score <- ppi_df$combined_score[match(edge_key, ppi_key1)]
edge_df_graph$combined_score[is.na(edge_df_graph$combined_score)] <-
  ppi_df$combined_score[match(edge_key[is.na(edge_df_graph$combined_score)], ppi_key2)]
igraph::E(g)$combined_score <- edge_df_graph$combined_score

write.csv(edge_df_graph, file = file.path(outdir, "11_PPI_edges_used_for_plot.csv"), row.names = FALSE)

set.seed(123)
p_ppi <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(alpha = combined_score), colour = "grey70", show.legend = FALSE) +
  scale_edge_alpha(range = c(0.2, 0.9)) +
  geom_node_point(aes(size = degree, color = hub_type)) +
  scale_color_manual(values = c("Hub gene" = "#D73027", "Other gene" = "#4575B4")) +
  geom_node_text(aes(label = ifelse(hub_type == "Hub gene", name, "")), repel = TRUE, size = 4) +
  scale_size(range = c(3, 10)) +
  labs(title = "Protein-protein interaction network", color = NULL, size = "Degree") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "right")

save_pdf_plot(p_ppi, file.path(outdir, "12.PPI_network.pdf"), width = 10, height = 8)
ggsave(file.path(outdir, "12.PPI_network.png"), plot = p_ppi, width = 10, height = 8, dpi = 600)

network_summary <- data.frame(
  metric = c("nodes", "edges", "density", "average_degree"),
  value = c(igraph::vcount(g), igraph::ecount(g), igraph::edge_density(g), mean(igraph::degree(g)))
)
write.csv(network_summary, file = file.path(outdir, "13_PPI_network_summary.csv"), row.names = FALSE)

write_summary(c(
  "=== GO/KEGG/PPI ===",
  paste0("Input genes: ", length(gene_symbols)),
  paste0("Converted ENTREZ IDs: ", length(entrez_ids)),
  paste0("GO terms (significant / total): ", nrow(go_sig), " / ", nrow(go_res)),
  paste0("GO cutoff: ", go_used_cutoff),
  paste0("KEGG source: ", kegg_source),
  paste0("KEGG terms (significant / total): ", nrow(kegg_sig), " / ", nrow(kegg_res)),
  paste0("KEGG cutoff: ", kegg_used_cutoff),
  paste0("PPI nodes: ", igraph::vcount(g), " | edges: ", igraph::ecount(g)),
  paste0("PPI score cutoff: ", ppi_score_cutoff),
  paste0("Hub genes (top ", ppi_hub_top_n, "): selected by degree/betweenness centrality"),
  "",
  "=== 限制与说明 ===",
  "GO enrichment covers BP, CC, MF ontologies (clusterProfiler).",
  "KEGG enrichment may use online DB or local GMT file.",
  "PPI network built from STRING DB with combined_score cutoff.",
  "Hub genes identified by degree/betweenness centrality."
))

log_message("GO/KEGG/PPI step finished.")
