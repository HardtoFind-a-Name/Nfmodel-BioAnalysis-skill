options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

get_arg <- function(flag, required = TRUE, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (required) stop(paste("Missing argument:", flag), call. = FALSE)
    return(default)
  }
  if (idx == length(args)) stop(paste("Missing value for:", flag), call. = FALSE)
  args[idx + 1]
}

write_note <- function(lines, path) {
  ensure_dir(dirname(path))
  writeLines(as.character(lines), path)
  invisible(path)
}

stop_with_note <- function(msg, note_file) {
  write_note(msg, note_file)
  stop(msg, call. = FALSE)
}

safe_file_label <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "model", x)
}

first_existing_file <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

get_validation_files <- function(input_dir, cohort_id) {
  cohort_lower <- tolower(cohort_id)
  nested_expr_file <- file.path(input_dir, cohort_id, paste0("07.", cohort_lower, "_expr.csv"))
  nested_surv_file <- file.path(input_dir, cohort_id, paste0("09.", cohort_lower, "_survival.csv"))
  flat_expr_file <- file.path(input_dir, paste0("07.", cohort_lower, "_expr.csv"))
  flat_surv_file <- file.path(input_dir, paste0("09.", cohort_lower, "_survival.csv"))

  list(
    expr_file = first_existing_file(c(nested_expr_file, flat_expr_file)),
    surv_file = first_existing_file(c(nested_surv_file, flat_surv_file)),
    expected_expr_file = nested_expr_file,
    expected_surv_file = nested_surv_file
  )
}

read_gene_list <- function(path, preferred_cols = c("gene", "gene_symbol", "symbol", "Gene", "genes")) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    return(character())
  }
  gene_col <- intersect(preferred_cols, colnames(df))[1] %||% colnames(df)[1]
  genes <- unique(trimws(as.character(df[[gene_col]])))
  genes[!is.na(genes) & genes != ""]
}

read_coef_file <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    stop(paste0("系数文件为空: ", path), call. = FALSE)
  }

  gene_col <- intersect(c("gene", "gene_symbol", "symbol", "Gene", "feature"), colnames(df))[1] %||% colnames(df)[1]
  coef_col <- intersect(c("coef", "coefficient", "Coefficient", "beta", "Coef"), colnames(df))[1]
  if (is.null(coef_col) || is.na(coef_col)) {
    numeric_cols <- colnames(df)[vapply(df, function(x) any(!is.na(suppressWarnings(as.numeric(x)))), logical(1))]
    numeric_cols <- setdiff(numeric_cols, gene_col)
    coef_col <- numeric_cols[1]
  }
  if (is.null(coef_col) || is.na(coef_col)) {
    stop(paste0("无法在系数文件中识别 coef 列: ", path), call. = FALSE)
  }

  out <- data.frame(
    gene = trimws(as.character(df[[gene_col]])),
    coef = suppressWarnings(as.numeric(df[[coef_col]])),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$gene) & out$gene != "" & !is.na(out$coef), , drop = FALSE]
  out <- out[out$coef != 0, , drop = FALSE]
  out <- out[!duplicated(out$gene), , drop = FALSE]
  if (nrow(out) == 0) {
    stop(paste0("系数文件未得到有效 gene + coef: ", path), call. = FALSE)
  }
  out
}

prepare_expr_matrix <- function(expr_file, gene_vec = NULL, dataset_name = "dataset") {
  if (!file.exists(expr_file)) {
    stop(paste0("表达矩阵文件不存在: ", expr_file), call. = FALSE)
  }
  expr_df <- read.csv(expr_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(expr_df) == 0) {
    stop(paste0("表达矩阵为空: ", expr_file), call. = FALSE)
  }

  expr_df$gene_symbol <- trimws(as.character(rownames(expr_df)))
  expr_df <- expr_df[!is.na(expr_df$gene_symbol) & expr_df$gene_symbol != "", , drop = FALSE]
  expr_df <- expr_df[!duplicated(expr_df$gene_symbol), , drop = FALSE]

  if (!is.null(gene_vec)) {
    gene_vec <- unique(trimws(as.character(gene_vec)))
    gene_vec <- gene_vec[!is.na(gene_vec) & gene_vec != ""]
    expr_df <- expr_df[expr_df$gene_symbol %in% gene_vec, , drop = FALSE]
    if (nrow(expr_df) == 0) {
      stop(paste0("表达矩阵中未匹配到任何目标基因: ", dataset_name), call. = FALSE)
    }
  }

  matched_genes <- expr_df$gene_symbol
  missing_genes <- if (is.null(gene_vec)) character() else setdiff(gene_vec, matched_genes)
  rownames(expr_df) <- expr_df$gene_symbol
  expr_df$gene_symbol <- NULL
  expr_mat <- as.data.frame(t(expr_df), check.names = FALSE)
  expr_mat$sample <- rownames(expr_mat)
  rownames(expr_mat) <- NULL
  expr_mat$sample <- trimws(as.character(expr_mat$sample))

  gene_cols <- setdiff(colnames(expr_mat), "sample")
  for (g in gene_cols) {
    expr_mat[[g]] <- suppressWarnings(as.numeric(expr_mat[[g]]))
  }

  list(
    expr = expr_mat,
    matched_genes = matched_genes,
    missing_genes = missing_genes,
    dataset_name = dataset_name
  )
}

prepare_surv_data <- function(surv_file) {
  if (!file.exists(surv_file)) {
    stop(paste0("生存信息文件不存在: ", surv_file), call. = FALSE)
  }
  surv_df <- read.csv(surv_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(surv_df) == 0) {
    stop(paste0("生存信息表为空: ", surv_file), call. = FALSE)
  }
  colnames(surv_df)[1] <- "sample"
  required_cols <- c("sample", "OS.time", "OS")
  missing_cols <- setdiff(required_cols, colnames(surv_df))
  if (length(missing_cols) > 0) {
    stop(paste0("生存信息缺少必需列: ", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  surv_df$sample <- trimws(as.character(surv_df$sample))
  surv_df$OS.time <- suppressWarnings(as.numeric(surv_df$OS.time))
  surv_df$OS <- suppressWarnings(as.numeric(surv_df$OS))
  surv_df <- surv_df[!is.na(surv_df$OS.time) & !is.na(surv_df$OS), , drop = FALSE]
  surv_df <- surv_df[surv_df$OS.time > 0 & surv_df$OS %in% c(0, 1), , drop = FALSE]
  surv_df <- surv_df[!duplicated(surv_df$sample), , drop = FALSE]
  if (nrow(surv_df) == 0) {
    stop(paste0("生存信息清洗后为空: ", surv_file), call. = FALSE)
  }
  surv_df
}

merge_expr_surv <- function(expr_list, surv_df) {
  merged <- merge(surv_df, expr_list$expr, by = "sample", all = FALSE)
  if (nrow(merged) == 0) {
    stop("表达矩阵与生存信息合并后为空，请检查 sample ID。", call. = FALSE)
  }
  merged
}

drop_invalid_gene_cols <- function(data_df, gene_cols) {
  for (g in gene_cols) {
    data_df[[g]] <- suppressWarnings(as.numeric(data_df[[g]]))
  }
  valid <- gene_cols[colSums(!is.na(data_df[, gene_cols, drop = FALSE])) > 0]
  sd_vec <- vapply(data_df[, valid, drop = FALSE], sd, numeric(1), na.rm = TRUE)
  valid <- valid[!is.na(sd_vec) & sd_vec > 0]
  list(data = data_df[, c("sample", "OS.time", "OS", valid), drop = FALSE], genes = valid)
}

calculate_risk_score <- function(data_df, coef_df) {
  gene_vec <- coef_df$gene
  coef_vec <- coef_df$coef
  names(coef_vec) <- gene_vec
  missing_genes <- setdiff(gene_vec, colnames(data_df))
  if (length(missing_genes) > 0) {
    stop(paste0("以下风险模型基因在数据中缺失: ", paste(missing_genes, collapse = ", ")), call. = FALSE)
  }
  as.numeric(as.matrix(data_df[, gene_vec, drop = FALSE]) %*% coef_vec)
}

make_risk_df <- function(data_df, coef_df, cutoff = NULL, cutoff_rule = "cohort_median") {
  risk_score <- calculate_risk_score(data_df, coef_df)
  cutoff_value <- cutoff %||% median(risk_score, na.rm = TRUE)
  risk_df <- data_df
  risk_df$riskScore <- risk_score
  risk_df$cutoff <- cutoff_value
  risk_df$cutoff_rule <- cutoff_rule
  risk_df$riskGroup <- ifelse(risk_df$riskScore > cutoff_value, "High risk", "Low risk")
  risk_df$riskGroup <- factor(risk_df$riskGroup, levels = c("Low risk", "High risk"))
  risk_df
}

days_to_years <- function(days) {
  yr <- days / 365
  ifelse(yr == floor(yr), paste0(floor(yr), "y"), paste0(round(yr, 1), "y"))
}

extract_auc <- function(roc_res, time_points = c(365, 1095, 1825)) {
  auc <- rep(NA_real_, length(time_points))
  if (!is.null(roc_res) && !is.null(roc_res$AUC)) {
    auc[seq_along(roc_res$AUC)] <- as.numeric(roc_res$AUC)
  }
  names(auc) <- paste0("auc_", days_to_years(time_points[seq_along(auc)]))
  auc
}

plot_timeROC <- function(risk_df, cohort_name, out_prefix, output_dir, time_points = c(365, 1095, 1825)) {
  ensure_dir(output_dir)
  roc_res <- timeROC::timeROC(
    T = risk_df$OS.time,
    delta = risk_df$OS,
    marker = risk_df$riskScore,
    cause = 1,
    weighting = "marginal",
    times = time_points,
    iid = TRUE
  )
  yr_labels <- paste0(days_to_years(time_points), "-year AUC = ", sprintf("%.3f", roc_res$AUC))
  colors <- c("#f4be7e", "#f47e84", "#6694e9")

  for (ext in c("pdf", "png")) {
    out_file <- file.path(output_dir, paste0(out_prefix, "_timeROC.", ext))
    if (ext == "pdf") {
      pdf(out_file, width = 6, height = 5)
    } else {
      png(out_file, width = 6, height = 5, units = "in", res = 600)
    }
    plot(roc_res, time = time_points[1], col = colors[1], lwd = 2, title = FALSE)
    if (length(time_points) >= 2) plot(roc_res, time = time_points[2], col = colors[2], lwd = 2, add = TRUE, title = FALSE)
    if (length(time_points) >= 3) plot(roc_res, time = time_points[3], col = colors[3], lwd = 2, add = TRUE, title = FALSE)
    abline(0, 1, lty = 2, col = "grey50")
    legend("bottomright", legend = yr_labels, col = colors[seq_along(time_points)], lwd = 2, bty = "n")
    title(main = cohort_name)
    dev.off()
  }
  roc_res
}

safe_timeROC <- function(risk_df, cohort_name, out_prefix, output_dir, log_dir, time_points = c(365, 1095, 1825)) {
  tryCatch(
    plot_timeROC(risk_df, cohort_name, out_prefix, output_dir, time_points),
    error = function(e) {
      write_note(c(paste0("Cohort: ", cohort_name), paste0("timeROC failed: ", e$message)),
                 file.path(log_dir, paste0(out_prefix, "_timeROC_note.txt")))
      NULL
    }
  )
}

plot_km <- function(risk_df, cohort_name, out_prefix, output_dir, log_dir) {
  ensure_dir(output_dir)
  if (length(unique(risk_df$riskGroup)) < 2) {
    write_note(c(
      paste0("Cohort: ", cohort_name),
      "KM analysis skipped because only one risk group is present."
    ), file.path(log_dir, paste0(out_prefix, "_KM_note.txt")))
    return(invisible(NULL))
  }
  fit <- survival::survfit(survival::Surv(OS.time, OS) ~ riskGroup, data = risk_df)
  legend_labs <- levels(droplevels(as.factor(risk_df$riskGroup)))
  palette_map <- c("Low risk" = "#6694e9", "High risk" = "#f47e84")
  p_km <- survminer::ggsurvplot(
    fit,
    data = risk_df,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.col = "strata",
    legend.title = "Risk group",
    legend.labs = legend_labs,
    palette = unname(palette_map[legend_labs]),
    ggtheme = ggplot2::theme_bw(),
    title = cohort_name,
    xlab = "Time (day)",
    ylab = "Survival probability"
  )
  pdf(file.path(output_dir, paste0(out_prefix, "_KM.pdf")), width = 6, height = 6)
  print(p_km)
  dev.off()
  png(file.path(output_dir, paste0(out_prefix, "_KM.png")), width = 6, height = 6, units = "in", res = 600)
  print(p_km)
  dev.off()
  invisible(p_km)
}

plot_risk_distribution <- function(risk_df, cohort_name, out_prefix, output_dir) {
  ensure_dir(output_dir)
  risk_df <- risk_df[order(risk_df$riskScore), , drop = FALSE]
  risk_df$order_id <- seq_len(nrow(risk_df))
  p1 <- ggplot2::ggplot(risk_df, ggplot2::aes(x = order_id, y = riskScore, color = riskGroup)) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_hline(yintercept = unique(risk_df$cutoff)[1], linetype = 2) +
    ggplot2::geom_vline(xintercept = sum(risk_df$riskGroup == "Low risk") + 0.5, linetype = 2) +
    ggplot2::scale_color_manual(values = c("Low risk" = "#6694e9", "High risk" = "#f47e84"), drop = FALSE) +
    ggplot2::labs(title = cohort_name, x = "Patients (increasing risk score)", y = "Risk Score") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.title = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  p2 <- ggplot2::ggplot(risk_df, ggplot2::aes(x = order_id, y = OS.time / 365, color = factor(OS, levels = c(0, 1), labels = c("Alive", "Dead")))) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_vline(xintercept = sum(risk_df$riskGroup == "Low risk") + 0.5, linetype = 2) +
    ggplot2::scale_color_manual(values = c("Alive" = "#6694e9", "Dead" = "#f47e84")) +
    ggplot2::labs(x = "Patients (increasing risk score)", y = "Survival time (Years)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  p_all <- gridExtra::grid.arrange(p1, p2, ncol = 1)
  ggplot2::ggsave(file.path(output_dir, paste0(out_prefix, "_risk_survival_distribution.pdf")), p_all, width = 6, height = 7)
  ggplot2::ggsave(file.path(output_dir, paste0(out_prefix, "_risk_survival_distribution.png")), p_all, width = 6, height = 7, dpi = 600)
  invisible(p_all)
}

plot_heatmap <- function(risk_df, gene_vec, cohort_name, out_prefix, output_dir) {
  ensure_dir(output_dir)
  missing_genes <- setdiff(gene_vec, colnames(risk_df))
  if (length(missing_genes) > 0) {
    stop(paste0("热图绘制失败，缺失基因: ", paste(missing_genes, collapse = ", ")), call. = FALSE)
  }
  expr_mat <- as.matrix(risk_df[, gene_vec, drop = FALSE])
  rownames(expr_mat) <- risk_df$sample
  expr_mat <- t(scale(t(expr_mat)))
  expr_mat[is.na(expr_mat)] <- 0
  ann_col <- data.frame(RiskGroup = factor(risk_df$riskGroup, levels = c("Low risk", "High risk")))
  rownames(ann_col) <- risk_df$sample
  ann_colors <- list(RiskGroup = c("Low risk" = "#6694e9", "High risk" = "#f47e84"))
  for (ext in c("pdf", "png")) {
    out_file <- file.path(output_dir, paste0(out_prefix, "_heatmap.", ext))
    if (ext == "pdf") {
      pdf(out_file, width = 8, height = 6)
    } else {
      png(out_file, width = 8, height = 6, units = "in", res = 600)
    }
    ok <- tryCatch({
      pheatmap::pheatmap(t(expr_mat), annotation_col = ann_col,
        annotation_colors = ann_colors,
        show_colnames = FALSE, cluster_cols = FALSE, main = cohort_name)
      TRUE
    }, error = function(e) { message("heatmap warning: ", e$message); FALSE })
    dev.off()
    if (!ok) try(file.remove(out_file), silent = TRUE)
  }
  invisible(NULL)
}

plot_risk_set <- function(risk_df, coef_df, cohort_name, out_prefix, output_dir, log_dir, time_points = c(365, 1095, 1825)) {
  ensure_dir(output_dir)
  write.csv(risk_df, file.path(output_dir, paste0(out_prefix, "_risk_score.csv")), row.names = FALSE, quote = FALSE)
  roc_res <- safe_timeROC(risk_df, cohort_name, out_prefix, output_dir, log_dir, time_points)
  plot_km(risk_df, cohort_name, out_prefix, output_dir, log_dir)
  plot_risk_distribution(risk_df, cohort_name, out_prefix, output_dir)
  plot_heatmap(risk_df, coef_df$gene, cohort_name, out_prefix, output_dir)
  roc_res
}

plot_cox_forest <- function(res_df, output_dir, out_prefix, title_text = "Cox model", ci_label = "95%", max_genes = 30) {
  ensure_dir(output_dir)
  hr_l_col <- if ("HR.L" %in% colnames(res_df)) "HR.L" else "HR.95L"
  hr_h_col <- if ("HR.H" %in% colnames(res_df)) "HR.H" else "HR.95H"
  forest_df <- res_df[!is.na(res_df$HR) & !is.na(res_df[[hr_l_col]]) & !is.na(res_df[[hr_h_col]]) & !is.na(res_df$pvalue), , drop = FALSE]
  forest_df <- forest_df[is.finite(forest_df$HR) & is.finite(forest_df[[hr_l_col]]) & is.finite(forest_df[[hr_h_col]]), , drop = FALSE]
  if (nrow(forest_df) == 0) return(invisible(NULL))
  forest_df <- forest_df[order(forest_df$HR), , drop = FALSE]
  forest_df <- utils::head(forest_df, max_genes)

  tabletext <- cbind(
    c("Gene", forest_df$gene),
    c("HR", format(round(forest_df$HR, 5), nsmall = 3)),
    c(paste0("lower ", ci_label, " CI"), format(round(forest_df[[hr_l_col]], 5), nsmall = 3)),
    c(paste0("upper ", ci_label, " CI"), format(round(forest_df[[hr_h_col]], 5), nsmall = 3)),
    c("pvalue", format(round(forest_df$pvalue, 5), nsmall = 3))
  )

  x_min <- floor(min(forest_df[[hr_l_col]], na.rm = TRUE) * 10) / 10
  x_max <- ceiling(max(forest_df[[hr_h_col]], na.rm = TRUE) * 10) / 10
  x_min <- min(x_min, 1)
  x_max <- max(x_max, 1)
  clip_min <- max(0, x_min - 0.2)
  clip_max <- x_max + 0.2
  xticks_vec <- pretty(c(clip_min, clip_max), n = 5)
  xticks_vec <- xticks_vec[xticks_vec >= 0]
  n_plot_rows <- nrow(forest_df) + 1
  lineheight_cm <- min(1.8, max(0.8, 10 / n_plot_rows))
  plot_height_in <- min(12, max(4.5, 0.9 + 0.55 * n_plot_rows))

  fp <- forestplot::forestplot(
    labeltext = tabletext,
    graph.pos = 4,
    col = forestplot::fpColors(box = "red", lines = "royalblue", zero = "gray50"),
    mean = c(NA, forest_df$HR),
    lower = c(NA, forest_df[[hr_l_col]]),
    upper = c(NA, forest_df[[hr_h_col]]),
    boxsize = 0.1,
    lwd.ci = 3,
    ci.vertices.height = 0.08,
    ci.vertices = TRUE,
    zero = 1,
    lwd.zero = 0.5,
    colgap = grid::unit(5, "mm"),
    xticks = xticks_vec,
    clip = c(clip_min, clip_max),
    lwd.xaxis = 2,
    lineheight = grid::unit(lineheight_cm, "cm"),
    graphwidth = grid::unit(0.6, "npc"),
    cex = 0.9,
    fn.ci_norm = forestplot::fpDrawCircleCI,
    hrzl_lines = list("2" = grid::gpar(col = "black", lty = 1, lwd = 2)),
    txt_gp = forestplot::fpTxtGp(
      label = grid::gpar(cex = 1),
      ticks = grid::gpar(cex = 0.8, fontface = "bold"),
      xlab = grid::gpar(cex = 1, fontface = "bold"),
      title = grid::gpar(cex = 1.25, fontface = "bold")
    ),
    xlab = "Hazard Ratio",
    grid = TRUE,
    title = title_text
  )

  pdf(file.path(output_dir, paste0(out_prefix, ".pdf")), width = 16, height = plot_height_in, family = "Times", onefile = FALSE)
  print(fp)
  dev.off()
  png(file.path(output_dir, paste0(out_prefix, ".png")), width = 16, height = plot_height_in, units = "in", res = 600, family = "Times")
  print(fp)
  dev.off()
  invisible(NULL)
}

plot_zph_panels <- function(zph_obj, output_dir, out_prefix, title_prefix = "Cox PH diagnostics") {
  ensure_dir(output_dir)
  term_names <- rownames(zph_obj$table)
  gene_terms <- term_names[term_names != "GLOBAL"]
  if (length(gene_terms) == 0) {
    return(invisible(NULL))
  }
  n_panel <- length(gene_terms)
  ncol_plot <- if (n_panel <= 4) 2 else 3
  nrow_plot <- ceiling(n_panel / ncol_plot)
  draw_panel <- function() {
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mfrow = c(nrow_plot, ncol_plot), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
    for (term in gene_terms) {
      plot(zph_obj, var = term, main = paste0(term, " Schoenfeld residuals"), xlab = "Time", ylab = "Scaled Schoenfeld residuals")
      abline(h = 0, lty = 2, col = "grey50")
    }
    mtext(title_prefix, outer = TRUE, font = 2, cex = 1.2)
  }
  pdf(file.path(output_dir, paste0(out_prefix, ".pdf")), width = 12, height = max(6, 3.5 * nrow_plot))
  draw_panel()
  dev.off()
  png(file.path(output_dir, paste0(out_prefix, ".png")), width = 12, height = max(6, 3.5 * nrow_plot), units = "in", res = 600)
  draw_panel()
  dev.off()
  invisible(NULL)
}

write_summary_row <- function(path, row_df) {
  ensure_dir(dirname(path))
  if (file.exists(path)) {
    old <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    out <- rbind(old, row_df)
  } else {
    out <- row_df
  }
  write.csv(out, path, row.names = FALSE, quote = FALSE)
  invisible(path)
}
