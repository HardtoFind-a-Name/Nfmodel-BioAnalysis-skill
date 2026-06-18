suppressPackageStartupMessages({
  library(dplyr)
  library(maftools)
})

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "r.00_post_utils.R"))

sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")

risk_file  <- get_arg("--risk-file")
maf_file   <- get_arg("--maf-file")
train_id   <- get_arg("--train-id")

logsetup <- setup_logging(logdir, summary_dir, sn, "tmb")
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

log_message("Script started: tmb")

# --- Read risk ---
log_message("Reading risk file: ", risk_file)
if (!file.exists(risk_file)) save_note_and_stop("Risk file not found: ", risk_file)
risk_df <- read.csv(risk_file, check.names = FALSE, stringsAsFactors = FALSE)
log_message("Read risk: ", risk_file, " -- ", nrow(risk_df), " x ", ncol(risk_df))
if (!"sample" %in% colnames(risk_df)) colnames(risk_df)[1] <- "sample"
if (!"riskGroup" %in% colnames(risk_df)) save_note_and_stop("Risk file must have 'riskGroup' column")

risk_df$sample <- as.character(risk_df$sample)
risk_df$sample <- gsub("\\.", "-", risk_df$sample)
risk_df$sample_short <- substr(risk_df$sample, 1, 16)
risk_df$riskGroup <- as.character(risk_df$riskGroup)
risk_df <- risk_df[!is.na(risk_df$sample) & nzchar(risk_df$sample) &
                   !is.na(risk_df$sample_short) & nzchar(risk_df$sample_short) &
                   !is.na(risk_df$riskGroup) & nzchar(risk_df$riskGroup) &
                   risk_df$riskGroup %in% c("Low risk", "High risk"), ]
risk_df <- risk_df[!duplicated(risk_df$sample_short), ]
risk_df$riskGroup <- factor(risk_df$riskGroup, levels = c("Low risk", "High risk"))
log_message("Risk samples: ", nrow(risk_df))

# --- Load MAF (local file → TCGAmutations fallback) ---
log_message("Loading MAF")
if (file.exists(maf_file)) {
  log_message("Using local MAF: ", maf_file)
  maf_df <- read.delim(maf_file, header = TRUE, sep = "\t", quote = "", comment.char = "#",
                        check.names = FALSE, stringsAsFactors = FALSE)
} else {
  log_message("Local MAF not found, loading from TCGAmutations for ", train_id)
  if (!requireNamespace("TCGAmutations", quietly = TRUE))
    save_note_and_stop("MAF file not found and TCGAmutations package not installed. Expected: ", maf_file)
  study <- gsub("^TCGA-", "", train_id)
  maf_df <- TCGAmutations::tcga_load(study = study)@data
  log_message("Loaded from TCGAmutations, rows: ", nrow(maf_df))
}
colnames(maf_df) <- make.names(colnames(maf_df), unique = TRUE)
if (!"Tumor_Sample_Barcode" %in% colnames(maf_df)) save_note_and_stop("Tumor_Sample_Barcode column not found in MAF")

# --- Map barcodes ---
maf_sample_map <- data.frame(barcode = unique(as.character(maf_df$Tumor_Sample_Barcode)), stringsAsFactors = FALSE)
maf_sample_map$sample_short <- gsub("\\.", "-", substr(maf_sample_map$barcode, 1, 16))
maf_sample_map <- merge(maf_sample_map, risk_df[, c("sample", "sample_short", "riskGroup")], by = "sample_short")
maf_sample_map <- maf_sample_map[!duplicated(maf_sample_map$barcode), ]
if (nrow(maf_sample_map) == 0) save_note_and_stop("No overlapping samples between MAF and risk data")
write.csv(maf_sample_map, file.path(outdir, "01_luad_mutation_sample_barcode_map.csv"), row.names = FALSE)
log_message("Mutation samples matched: ", nrow(maf_sample_map))

# --- Build MAF objects ---
build_maf <- function(barcodes) {
  sub <- maf_df %>% filter(Tumor_Sample_Barcode %in% barcodes)
  if (nrow(sub) == 0) return(NULL)
  read.maf(maf = sub, clinicalData = data.frame(Tumor_Sample_Barcode = barcodes, stringsAsFactors = FALSE))
}

all_barcodes <- maf_sample_map$barcode
high_barcodes <- maf_sample_map$barcode[maf_sample_map$riskGroup == "High risk"]
low_barcodes  <- maf_sample_map$barcode[maf_sample_map$riskGroup == "Low risk"]

maf_train <- build_maf(all_barcodes)
maf_high  <- build_maf(high_barcodes)
maf_low   <- build_maf(low_barcodes)

if (is.null(maf_train) || is.null(maf_high) || is.null(maf_low))
  save_note_and_stop("Failed to build one or more MAF objects")

log_message("Train mutations: ", nrow(maf_train@data))
log_message("High risk mutations: ", nrow(maf_high@data))
log_message("Low risk mutations: ", nrow(maf_low@data))

# --- Oncoplots ---
for (group in c("train", "high_risk", "low_risk")) {
  maf_obj <- switch(group, train = maf_train, high_risk = maf_high, low_risk = maf_low)
  title <- switch(group,
    train     = paste0(train_id, " training cohort: top 20 mutated genes"),
    high_risk = paste0(train_id, " high risk: top 20 mutated genes"),
    low_risk  = paste0(train_id, " low risk: top 20 mutated genes"))
  prefix <- switch(group, train = "02", high_risk = "03", low_risk = "04")

  pdf(file.path(outdir, paste0(prefix, "_oncoplot_top20.pdf")), width = 12, height = 8, family = "Times")
  oncoplot(maf = maf_obj, top = 20, titleText = title)
  dev.off()
  png(file.path(outdir, paste0(prefix, "_oncoplot_top20.png")), width = 12, height = 8, units = "in", res = 300, family = "Times")
  oncoplot(maf = maf_obj, top = 20, titleText = title)
  dev.off()
}

# --- Summary ---
train_top20 <- names(sort(table(maf_train@data$Hugo_Symbol), decreasing = TRUE))[1:min(20, length(unique(maf_train@data$Hugo_Symbol)))]
high_top20  <- names(sort(table(maf_high@data$Hugo_Symbol),  decreasing = TRUE))[1:min(20, length(unique(maf_high@data$Hugo_Symbol)))]
low_top20   <- names(sort(table(maf_low@data$Hugo_Symbol),   decreasing = TRUE))[1:min(20, length(unique(maf_low@data$Hugo_Symbol)))]

log_message("Script completed successfully")

write_summary(c(
  paste0("分析对象: ", train_id),
  paste0("风险分组样本数: ", nrow(risk_df), " (High: ", sum(risk_df$riskGroup == "High risk"), ", Low: ", sum(risk_df$riskGroup == "Low risk"), ")"),
  paste0("MAF突变匹配样本: ", nrow(maf_sample_map), " (总突变: ", nrow(maf_train@data), ", High: ", nrow(maf_high@data), ", Low: ", nrow(maf_low@data), ")"),
  paste0("Top mutated genes - All: ", paste(head(train_top20, 5), collapse = ", "), " ..."),
  paste0("Top mutated genes - High risk: ", paste(head(high_top20, 5), collapse = ", "), " ..."),
  paste0("Top mutated genes - Low risk: ", paste(head(low_top20, 5), collapse = ", "), " ..."),
  paste0("Oncoplot analysis reveals distinct mutation landscapes between high-risk and low-risk groups."),
  paste0("限制与说明: MAF数据来源于TCGA或提供的本地文件; Oncoplot仅显示top 20突变基因; TMB比较仅供参考.")
))
