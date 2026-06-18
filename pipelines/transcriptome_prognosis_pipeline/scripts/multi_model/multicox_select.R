args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, required = TRUE, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (required) stop(paste("Missing argument:", flag), call. = FALSE)
    return(default)
  }
  if (idx == length(args)) stop(paste("Missing value for:", flag), call. = FALSE)
  args[idx + 1]
}

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(ggplot2)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

expr_file <- get_arg("--expr")
surv_file <- get_arg("--surv")
gene_list_file <- get_arg("--gene-list")
model_id <- get_arg("--model-id")
sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")

sub_id <- get_arg("--sub-id", required = FALSE, default = NULL)
if (is.null(sub_id)) sub_id <- paste0(model_id, "/multicox")

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

candidate_genes <- read_gene_list(gene_list_file)
if (length(candidate_genes) < 2) save_note_and_stop("multiCox requires at least 2 genes.")

expr_list <- prepare_expr_matrix(expr_file, candidate_genes, "train")
surv_df <- prepare_surv_data(surv_file)
train_data <- merge_expr_surv(expr_list, surv_df)
cleaned <- drop_invalid_gene_cols(train_data, setdiff(colnames(train_data), c("sample", "OS.time", "OS")))
train_data <- cleaned$data
gene_cols <- cleaned$genes
if (length(gene_cols) < 2) save_note_and_stop("After cleaning, fewer than 2 genes for multiCox.")

formula_str <- paste0("Surv(OS.time, OS) ~ ", paste0("`", gene_cols, "`", collapse = " + "))
full_fit <- tryCatch(
  survival::coxph(as.formula(formula_str), data = train_data),
  error = function(e) { save_note_and_stop(paste0("multiCox fit failed: ", e$message)) }
)

fit_sum <- summary(full_fit)
res_df <- data.frame(
  gene = gene_cols,
  coef = fit_sum$coefficients[, "coef"],
  HR = fit_sum$conf.int[, "exp(coef)"],
  HR.95L = fit_sum$conf.int[, "lower .95"],
  HR.95H = fit_sum$conf.int[, "upper .95"],
  z = fit_sum$coefficients[, "z"],
  pvalue = fit_sum$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
) %>% dplyr::arrange(pvalue)

write.csv(res_df, file.path(outdir, "01.multiCox_all_results.csv"), row.names = FALSE, quote = FALSE)

zph <- tryCatch(survival::cox.zph(full_fit), error = function(e) NULL)
if (!is.null(zph)) {
  ph_df <- as.data.frame(zph$table)
  ph_df$term <- rownames(ph_df)
  rownames(ph_df) <- NULL
  write.csv(ph_df, file.path(outdir, "02.multiCox_ph_test.csv"), row.names = FALSE, quote = FALSE)
}

coef_df <- data.frame(gene = res_df$gene, coef = res_df$coef, stringsAsFactors = FALSE)
write.csv(coef_df, file.path(outdir, "03.multiCox_final_genes_coef.csv"), row.names = FALSE, quote = FALSE)
write.csv(data.frame(gene_symbol = res_df$gene), file.path(outdir, "04.multiCox_gene_list.csv"), row.names = FALSE, quote = FALSE)

log_message("Input genes: ", length(gene_cols), ", samples: ", nrow(train_data))
log_message("Multi-Cox significant genes (p<0.05): ", sum(res_df$pvalue < 0.05), " / ", nrow(res_df))

if (nrow(res_df) > 0) {
  plot_cox_forest(res_df, outdir, "05.multiCox_forest", paste0(model_id, " multiCox"))
}

write_summary(c(
  paste0("模型ID: ", model_id),
  paste0("输入基因数: ", length(gene_cols)),
  paste0("样本数: ", nrow(train_data)),
  paste0("多因素Cox回归基因数: ", nrow(res_df)),
  paste0("显著基因数 (p<0.05): ", sum(res_df$pvalue < 0.05)),
  paste0("Top 基因: ", paste(head(res_df$gene, 5), collapse = ", ")),
  "",
  paste0("结论: 对", length(gene_cols), "个基因进行多因素Cox回归分析，共", sum(res_df$pvalue < 0.05), "个基因显著。"),
  "",
  paste0("限制与说明: 多因素Cox模型可能受多重共线性影响。PH假设未进行全局检验。")
))

log_message("Script completed: ", sub_id)
