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

suppressPackageStartupMessages({ library(dplyr) })

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "..", "r.00_post_utils.R"))

input_dir <- get_arg("--input-dir")
sn        <- get_arg("--sn")
logdir    <- get_arg("--logdir")
summary_dir <- get_arg("--summary-dir")
outdir      <- get_arg("--outdir")
auc_threshold <- as.numeric(get_arg("--auc-threshold", required = FALSE, default = "0.6"))

sub_id <- get_arg("--sub-id", required = FALSE, default = "99_summary")

logsetup <- setup_logging(logdir, summary_dir, sn, "multi_model", sub_id)
log_message <- logsetup$log_message
write_summary <- logsetup$write_summary
save_note_and_stop <- logsetup$save_note_and_stop

log_message("Script started: ", sub_id)
log_message("Log file: ", logsetup$log_file)
log_message("Summary file: ", logsetup$summary_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

train_files <- list.files(input_dir, pattern = "train_summary\\.csv$", full.names = TRUE, recursive = TRUE)
valid_files <- list.files(input_dir, pattern = "validation_summary\\.csv$", full.names = TRUE, recursive = TRUE)

train_all <- if (length(train_files) > 0) {
  dplyr::bind_rows(lapply(train_files, function(f) tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)))
} else data.frame()

valid_all <- if (length(valid_files) > 0) {
  dplyr::bind_rows(lapply(valid_files, function(f) tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)))
} else data.frame()

comparison <- dplyr::bind_rows(train_all, valid_all)
write.csv(comparison, file.path(outdir, "00.multi_model_comparison.csv"), row.names = FALSE, quote = FALSE)

pass_fail <- data.frame(model_id = character(), time_set = character(), pass = logical(), reason = character(), stringsAsFactors = FALSE)

if (nrow(train_all) > 0 && "time_set" %in% colnames(train_all)) {
  for (ts in unique(train_all$time_set)) {
    t_sub <- train_all[train_all$time_set == ts, , drop = FALSE]
    v_sub <- valid_all[valid_all$time_set == ts, , drop = FALSE]

    auc_cols <- grep("^auc_", names(t_sub), value = TRUE)
    if (length(auc_cols) == 0) next

    for (mid in unique(t_sub$model_id)) {
      t_row <- t_sub[t_sub$model_id == mid, , drop = FALSE]
      if (nrow(t_row) == 0) next
      train_pass <- all(t_row[, auc_cols] > auc_threshold, na.rm = TRUE)

      valid_pass <- FALSE
      pass_cohorts <- ""
      reason <- ""
      if (nrow(v_sub) > 0) {
        v_rows <- v_sub[v_sub$model_id == mid, , drop = FALSE]
        passing_ids <- c()
        for (j in seq_len(nrow(v_rows))) {
          cohort_pass <- all(v_rows[j, auc_cols] > auc_threshold, na.rm = TRUE)
          if (isTRUE(cohort_pass)) { passing_ids <- c(passing_ids, v_rows$cohort_id[j]); valid_pass <- TRUE }
        }
        if (length(passing_ids) > 0) pass_cohorts <- paste(passing_ids, collapse = "; ")
      }

      overall_pass <- isTRUE(train_pass) && isTRUE(valid_pass)
      if (!isTRUE(train_pass)) reason <- paste0("train AUC below ", auc_threshold, " (", ts, ")")
      else if (!isTRUE(valid_pass)) reason <- paste0("no validation cohort passed (", ts, ")")
      else reason <- "pass"

      pass_fail <- rbind(pass_fail, data.frame(model_id = mid, time_set = ts, pass = overall_pass, reason = reason, pass_cohorts = pass_cohorts, stringsAsFactors = FALSE))
    }
  }
}

write.csv(pass_fail, file.path(outdir, "01.pass_fail_report.csv"), row.names = FALSE, quote = FALSE)

log_message("Comparison: ", nrow(comparison), " rows, ", nrow(pass_fail), " models evaluated")

pass_count <- sum(pass_fail$pass, na.rm = TRUE)
fail_count <- sum(!pass_fail$pass, na.rm = TRUE)
best_models <- pass_fail$model_id[pass_fail$pass]

write_summary(c(
  paste0("输入目录: ", input_dir),
  paste0("比较的模型总数: ", nrow(pass_fail)),
  paste0("通过模型数: ", pass_count),
  paste0("未通过模型数: ", fail_count),
  paste0("AUC阈值: ", auc_threshold),
  paste0("最佳模型: ", ifelse(length(best_models) > 0, paste(best_models, collapse = "; "), "无")),
  "",
  paste0("结论: 共比较", nrow(pass_fail), "个模型，其中", pass_count, "个通过AUC>", auc_threshold, "的验证标准。"),
  "",
  paste0("限制与说明: pass/fail判断基于所有时间点AUC均>", auc_threshold, "的标准。未考虑临床适用性。详细结果见00.multi_model_comparison.csv和01.pass_fail_report.csv。")
))

log_message("Script completed: ", sub_id)
