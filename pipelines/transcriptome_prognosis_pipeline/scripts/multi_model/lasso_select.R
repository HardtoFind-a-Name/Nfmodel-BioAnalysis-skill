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
  library(glmnet)
  library(ggplot2)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

expr_file <- get_arg("--expr")
surv_file <- get_arg("--surv")
gene_list_file <- get_arg("--gene-list")
sn          <- get_arg("--sn")
logdir      <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")

sub_id <- get_arg("--sub-id", required = FALSE, default = "01_lasso_select")

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expr_file)) save_note_and_stop(paste0("Expression file not found: ", expr_file))
if (!file.exists(surv_file)) save_note_and_stop(paste0("Survival file not found: ", surv_file))
if (!file.exists(gene_list_file)) save_note_and_stop(paste0("Gene list file not found: ", gene_list_file))

candidate_genes <- read_gene_list(gene_list_file)
if (length(candidate_genes) < 2) save_note_and_stop("LASSO requires at least 2 candidate genes.")

expr_list <- prepare_expr_matrix(expr_file, candidate_genes, "train")
surv_df <- prepare_surv_data(surv_file)
train_data <- merge_expr_surv(expr_list, surv_df)
cleaned <- drop_invalid_gene_cols(train_data, setdiff(colnames(train_data), c("sample", "OS.time", "OS")))
train_data <- cleaned$data
gene_cols <- cleaned$genes
if (length(gene_cols) < 2) save_note_and_stop("After cleaning, fewer than 2 genes remain for LASSO.")

write.csv(train_data, file.path(outdir, "00.lasso_input_merged.csv"), row.names = FALSE, quote = FALSE)

x <- as.matrix(train_data[, gene_cols, drop = FALSE])
y <- survival::Surv(train_data$OS.time, train_data$OS)

set.seed(20260423)
cv_fit <- glmnet::cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)

coef_min <- as.matrix(coef(cv_fit, s = "lambda.min"))
selected <- data.frame(
  gene = rownames(coef_min), coef = as.numeric(coef_min), stringsAsFactors = FALSE
) %>%
  dplyr::filter(coef != 0) %>%
  dplyr::mutate(abs_coef = abs(coef)) %>%
  dplyr::arrange(dplyr::desc(abs_coef))

if (nrow(selected) == 0) save_note_and_stop("LASSO selected 0 genes at lambda.min.")

write.csv(data.frame(gene_symbol = selected$gene), file.path(outdir, "01.lasso_selected_genes_lambda_min.csv"), row.names = FALSE, quote = FALSE)
write.csv(selected[, c("gene", "coef")], file.path(outdir, "02.lasso_final_genes_coef.csv"), row.names = FALSE, quote = FALSE)

png(file.path(outdir, "03.lasso_cv_curve.png"), width = 7, height = 5, units = "in", res = 600)
plot(cv_fit, main = "LASSO-Cox CV curve")
dev.off()
tryCatch({ pdf(file.path(outdir, "03.lasso_cv_curve.pdf"), width = 7, height = 5); plot(cv_fit, main = "LASSO-Cox CV curve"); dev.off() }, error = function(e) NULL)

png(file.path(outdir, "04.lasso_coef_path.png"), width = 7, height = 5, units = "in", res = 600)
plot(cv_fit$glmnet.fit, xvar = "lambda", main = "LASSO coefficient path")
abline(v = log(cv_fit$lambda.min), lty = 2, col = "red")
abline(v = log(cv_fit$lambda.1se), lty = 2, col = "blue")
dev.off()
tryCatch({ pdf(file.path(outdir, "04.lasso_coef_path.pdf"), width = 7, height = 5); plot(cv_fit$glmnet.fit, xvar = "lambda", main = "LASSO coefficient path"); abline(v = log(cv_fit$lambda.min), lty = 2, col = "red"); abline(v = log(cv_fit$lambda.1se), lty = 2, col = "blue"); dev.off() }, error = function(e) NULL)

top_genes <- head(selected$gene, 5)

write_summary(c(
  paste0("分析队列: ", basename(dirname(expr_file))),
  paste0("输入基因数: ", length(candidate_genes)),
  paste0("cleaned后有效基因数: ", length(gene_cols)),
  paste0("样本数: ", nrow(train_data)),
  paste0("lambda.min: ", cv_fit$lambda.min),
  paste0("lambda.1se: ", cv_fit$lambda.1se),
  paste0("LASSO选择基因数 (lambda.min): ", nrow(selected)),
  paste0("Top genes by |coef|: ", paste(top_genes, collapse = ", ")),
  "",
  paste0("结论: 在", length(gene_cols), "个输入基因中，LASSO-Cox通过10折交叉验证在lambda.min处筛选出", nrow(selected), "个非零系数基因，作为后续建模的输入。"),
  "",
  paste0("限制与说明: LASSO选择依赖于lambda.min参数，可能过于保守。若需更宽松的筛选，可考虑lambda.1se。")
))

log_message("Script completed: ", sub_id)
