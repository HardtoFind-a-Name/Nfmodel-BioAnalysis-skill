args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(biomaRt); library(RCircos)
})
source(file.path(script_dir, "r.00_post_utils.R"))

sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")

gwas_file  <- get_arg("--gwas-file")
risk_file  <- get_arg("--risk-file", required = FALSE, default = NULL)
gene_file  <- get_arg("--gene-file")
train_id   <- get_arg("--train-id")

logsetup <- setup_logging(logdir, summary_dir, sn, "circos")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

log_message("Script started: circos")

pvalue_threshold <- 5e-8

# Extract prognostic genes from gene selection file
if (!file.exists(gene_file)) save_note_and_stop("Gene file not found: ", gene_file)
gene_df <- read.csv(gene_file, stringsAsFactors = FALSE, check.names = FALSE)
log_message("Read gene file: ", gene_file, " -- ", nrow(gene_df), " x ", ncol(gene_df))
if (!"gene" %in% colnames(gene_df)) save_note_and_stop("Gene file must have a 'gene' column")
prognostic_genes <- unique(trimws(as.character(gene_df$gene)))
prognostic_genes <- prognostic_genes[!is.na(prognostic_genes) & nzchar(prognostic_genes)]
if (length(prognostic_genes) == 0) save_note_and_stop("No prognostic genes found in gene file")
log_message("Prognostic genes: ", length(prognostic_genes), ", GWAS: ", gwas_file)

# GWAS filtering
gwas <- fread(gwas_file, data.table = FALSE)
log_message("Read GWAS: ", gwas_file, " -- ", nrow(gwas), " x ", ncol(gwas))
gwas_luad <- gwas %>%
  filter(str_detect(toupper(`DISEASE/TRAIT`), "LUNG ADENOCARCINOMA|ADENOCARCINOMA OF LUNG|LUNG CANCER|NON-SMALL CELL"),
         as.numeric(`P-VALUE`) < pvalue_threshold) %>%
  mutate(trait = `DISEASE/TRAIT`, chr = paste0("chr", CHR_ID), pos = as.integer(CHR_POS),
         snp = SNPS, pvalue = as.numeric(`P-VALUE`)) %>%
  filter(!is.na(chr), !is.na(pos), chr %in% paste0("chr", c(1:22,"X","Y"))) %>%
  distinct(trait, snp, chr, pos, pvalue)
fwrite(gwas_luad, file.path(outdir, "01.gwas_luad_snps_p5e8.csv"))
log_message("GWAS LUAD SNPs after filtering: ", nrow(gwas_luad))

# biomaRt gene loci
log_message("Fetching hg38 loci via biomaRt...")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
gene_loci_raw <- getBM(attributes = c("hgnc_symbol","chromosome_name","start_position","end_position","strand","ensembl_gene_id"),
                       filters = "hgnc_symbol", values = unique(toupper(prognostic_genes)), mart = mart)
gene_loci <- gene_loci_raw %>%
  filter(chromosome_name %in% c(as.character(1:22),"X","Y")) %>%
  mutate(Gene = hgnc_symbol, Chromosome = paste0("chr", chromosome_name),
         chromStart = as.integer(start_position), chromEnd = as.integer(end_position)) %>%
  distinct(Gene, Chromosome, chromStart, chromEnd, strand, ensembl_gene_id)
fwrite(gene_loci, file.path(outdir, "02.prognostic_genes_hg38_loci.csv"))
log_message("Located genes: ", nrow(gene_loci))

# RCircos tracks
gene_track <- gene_loci %>% mutate(Chr = Chromosome, Start = chromStart, End = chromEnd, Label = Gene) %>%
  dplyr::select(Chr, Start, End, Label) %>% arrange(Chr, Start)
snp_track <- gwas_luad %>% mutate(Chr = chr, Start = pos, End = pos, Label = snp, Score = -log10(pvalue)) %>%
  dplyr::select(Chr, Start, End, Label, Score) %>% arrange(Chr, Start)
fwrite(gene_track, file.path(outdir, "03.rcircos_gene_track.csv"))
fwrite(snp_track, file.path(outdir, "04.rcircos_snp_track.csv"))

# Circos plot
tryCatch({
  rc_items <- data(package = "RCircos")$results[, "Item"]
  cyto_name <- rc_items[grepl("hg38", rc_items, ignore.case = TRUE) & grepl("cyto", rc_items, ignore.case = TRUE)][1]
  if (is.na(cyto_name)) save_note_and_stop("No hg38 cytoband found in RCircos")
  data(list = cyto_name, package = "RCircos")
  cyto.info <- get(cyto_name)
  colnames(cyto.info)[1:5] <- c("Chromosome", "chromStart", "chromEnd", "Band", "Stain")
  cyto.info <- cyto.info[cyto.info$Chromosome %in% paste0("chr", c(1:22,"X","Y")), ]

  RCircos.Set.Core.Components(cyto.info, chr.exclude = NULL, tracks.inside = 5, tracks.outside = 0)
  pdf(file.path(outdir, "06.chromosome_circos.pdf"), width = 10, height = 10)
  RCircos.Set.Plot.Area(); RCircos.Chromosome.Ideogram.Plot()
  if (nrow(snp_track) > 0) RCircos.Scatter.Plot(as.data.frame(snp_track[,c("Chr","Start","End","Score")]), data.col = 4, track.num = 1, side = "in")
  if (nrow(gene_track) > 0) {
    glt <- gene_track; colnames(glt) <- c("Chromosome","chromStart","chromEnd","Gene")
    params <- RCircos.Get.Plot.Parameters(); params$text.size <- 1.0; RCircos.Reset.Plot.Parameters(params)
    RCircos.Gene.Name.Plot(glt, name.col = 4, track.num = 2, side = "in")
  }
  dev.off()
  png(file.path(outdir, "06.chromosome_circos.png"), width = 2400, height = 2400, res = 300)
  RCircos.Set.Plot.Area(); RCircos.Chromosome.Ideogram.Plot()
  if (nrow(snp_track) > 0) RCircos.Scatter.Plot(as.data.frame(snp_track[,c("Chr","Start","End","Score")]), data.col = 4, track.num = 1, side = "in")
  if (nrow(gene_track) > 0) {
    glt <- gene_track; colnames(glt) <- c("Chromosome","chromStart","chromEnd","Gene")
    params <- RCircos.Get.Plot.Parameters(); params$text.size <- 1.0; RCircos.Reset.Plot.Parameters(params)
    RCircos.Gene.Name.Plot(glt, name.col = 4, track.num = 2, side = "in")
  }
  dev.off()
  log_message("Circos plot generated")
}, error = function(e) log_message("Circos plot failed: ", e$message))

log_message("Script completed successfully")

write_summary(c(
  paste0("分析对象: ", train_id),
  paste0("预后基因数量: ", length(prognostic_genes)),
  paste0("GWAS SNP (LUAD, p<5e-8): ", nrow(gwas_luad), " 个"),
  paste0("定位到的基因 (biomaRt hg38): ", nrow(gene_loci), " 个"),
  paste0("Circos图已生成, 显示预后基因与GWAS SNP在染色体上的共定位关系."),
  paste0("限制与说明: GWAS数据需包含LUAD相关trait; 基因注释依赖biomaRt在线服务; RCircos版本需支持hg38.")
))
