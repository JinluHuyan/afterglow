# Differential expression routines and the high-level afterglow workflow.
# limma is the default method because it is broadly applicable to microarray and
# other log-scale expression matrices. limma-voom, edgeR and DESeq2 are available
# for high-throughput sequencing count matrices when users explicitly request
# count-scale methods.

#' Differential expression on a matrix
#'
#' @param mat Numeric matrix, genes x samples. For `method = "limma"`, this is
#'   usually log2-scale microarray or normalized expression. For
#'   `method = "limma_voom"`, `"edgeR"`, or `"deseq2"`, this must represent
#'   non-negative sequencing counts or count-like values.
#' @param group Data frame with `sample` and `group`, or a named group vector.
#' @param control_group Control group label.
#' @param treatment_group Treatment group label.
#' @param method Differential expression method. `"limma"` is recommended for
#'   microarray or log2 normalized expression; `"limma_voom"` and `"edgeR"`
#'   support RNA-seq count matrices; `"deseq2"` supports integer count matrices;
#'   `"welch"` is a lightweight base-R fallback.
#' @param assay_type Optional assay description used for user guidance. Use
#'   `"microarray"` for chip/log-expression matrices and `"rna_seq"` for raw
#'   or estimated sequencing counts.
#' @param log_scale Whether `mat` is already on log scale. If `FALSE`, log2
#'   fold change is computed from group means with `pseudocount`.
#' @param pseudocount Small positive number used for linear-scale logFC and
#'   count rounding safeguards.
#' @param deseq_round_counts Logical. DESeq2 requires integer counts. If `TRUE`,
#'   non-integer turnover-corrected values are rounded after non-negative
#'   validation. This is mathematically convenient but should be reported as an
#'   approximation.
#' @param verbose Logical. If `TRUE`, print progress and method diagnostics.
#' @return Data frame ordered by adjusted p-value.
#' @export
differential_expression <- function(mat,
                                    group,
                                    control_group,
                                    treatment_group,
                                    method = c("limma", "limma_voom", "edgeR", "deseq2", "welch"),
                                    assay_type = c("auto", "microarray", "rna_seq"),
                                    log_scale = TRUE,
                                    pseudocount = 1e-8,
                                    deseq_round_counts = TRUE,
                                    verbose = TRUE) {
  method <- .match_arg(method[[1L]], c("limma", "limma_voom", "edgeR", "deseq2", "welch"), "method")
  assay_type <- .match_arg(assay_type[[1L]], c("auto", "microarray", "rna_seq"), "assay_type")

  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  .stop_if(is.null(rownames(mat)) || is.null(colnames(mat)), "mat must have row and column names.")

  if (is.vector(group) && !is.data.frame(group)) {
    .stop_if(is.null(names(group)), "Named group vector must use sample names as names(group).")
    group <- data.frame(sample = names(group), group = as.character(group), stringsAsFactors = FALSE)
  }
  .stop_if(!all(c("sample", "group") %in% colnames(group)), "group must contain sample and group columns.")
  missing_samples <- setdiff(colnames(mat), group$sample)
  .stop_if(length(missing_samples) > 0L,
           sprintf("Group table is missing samples: %s", paste(missing_samples, collapse = ", ")))
  group <- group[match(colnames(mat), group$sample), , drop = FALSE]
  control_idx <- which(group$group == control_group)
  treatment_idx <- which(group$group == treatment_group)
  .stop_if(length(control_idx) < 2L || length(treatment_idx) < 2L,
           "Each group needs at least two samples for differential testing.")

  .inform(verbose, "afterglow: differential expression method = '%s'.", method)
  if (method == "limma") {
    .inform(verbose, "afterglow: limma is the default and is appropriate for microarray/chip or log-scale normalized expression matrices.")
    .stop_if(assay_type == "rna_seq" && !isTRUE(log_scale),
             "limma with assay_type='rna_seq' expects log-scale RNA-seq input (for example log2 CPM/TPM or voom-transformed values). Use log_scale=TRUE after transformation, or choose method='deseq2' for count-scale input.")
  }
  if (method == "deseq2") {
    .inform(verbose, "afterglow: DESeq2 is intended for high-throughput sequencing count matrices. Provide raw or estimated counts, not log2 intensities.")
  }
  if (method == "limma_voom") {
    .inform(verbose, "afterglow: limma_voom is intended for RNA-seq count matrices and uses edgeR normalization followed by limma::voom.")
  }
  if (method == "edgeR") {
    .inform(verbose, "afterglow: edgeR uses calcNormFactors, estimateDisp and glmQLFit/glmQLFTest on count-scale data.")
  }

  control_mean <- rowMeans(mat[, control_idx, drop = FALSE], na.rm = TRUE)
  treatment_mean <- rowMeans(mat[, treatment_idx, drop = FALSE], na.rm = TRUE)
  log_fc <- if (isTRUE(log_scale)) {
    treatment_mean - control_mean
  } else {
    log2((treatment_mean + pseudocount) / (control_mean + pseudocount))
  }

  if (method == "limma") {
    .stop_if(!requireNamespace("limma", quietly = TRUE),
             "method='limma' requires the Bioconductor package 'limma'. Install with BiocManager::install('limma').")
    design_group <- factor(group$group, levels = c(control_group, treatment_group))
    design <- stats::model.matrix(~ 0 + design_group)
    # Use syntactically safe contrast labels even when user-supplied group names
    # contain spaces, punctuation, or non-ASCII characters.
    colnames(design) <- c("control", "treatment")
    fit <- limma::lmFit(mat, design)
    contrast <- limma::makeContrasts(contrasts = "treatment-control", levels = design)
    fit2 <- limma::contrasts.fit(fit, contrast)
    fit2 <- limma::eBayes(fit2)
    tab <- limma::topTable(fit2, number = Inf, sort.by = "none")
    out <- data.frame(
      gene_id = rownames(mat),
      log2_fc = tab$logFC,
      average_expression = tab$AveExpr,
      statistic = tab$t,
      p_value = tab$P.Value,
      adjusted_p_value = tab$adj.P.Val,
      method = "limma",
      stringsAsFactors = FALSE
    )
  } else if (method == "limma_voom") {
    .stop_if(isTRUE(log_scale),
             "limma_voom requires count-scale data. Use expression_scale='linear' and provide raw or estimated RNA-seq counts.")
    .stop_if(any(is.na(mat)), "limma_voom input contains NA values after turnover correction.")
    .stop_if(any(!is.finite(mat)), "limma_voom input contains non-finite values after turnover correction.")
    .stop_if(any(mat < 0), "limma_voom input contains negative values after turnover correction.")
    .stop_if(!requireNamespace("limma", quietly = TRUE),
             "method='limma_voom' requires the Bioconductor package 'limma'. Install with BiocManager::install('limma').")
    .stop_if(!requireNamespace("edgeR", quietly = TRUE),
             "method='limma_voom' requires the Bioconductor package 'edgeR'. Install with BiocManager::install('edgeR').")
    non_integer <- any(abs(mat - round(mat)) > sqrt(.Machine$double.eps))
    if (non_integer) {
      .inform(verbose, "afterglow: limma_voom received non-integer count-like values; interpret results as an approximation.")
    }
    condition <- factor(group$group, levels = c(control_group, treatment_group))
    design <- stats::model.matrix(~ 0 + condition)
    colnames(design) <- c("control", "treatment")
    dge <- edgeR::DGEList(counts = mat, group = condition)
    dge <- edgeR::calcNormFactors(dge)
    v <- limma::voom(dge, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    contrast <- limma::makeContrasts(contrasts = "treatment-control", levels = design)
    fit2 <- limma::contrasts.fit(fit, contrast)
    fit2 <- limma::eBayes(fit2)
    tab <- limma::topTable(fit2, number = Inf, sort.by = "none")
    out <- data.frame(
      gene_id = rownames(mat),
      log2_fc = tab$logFC,
      average_expression = tab$AveExpr,
      statistic = tab$t,
      p_value = tab$P.Value,
      adjusted_p_value = tab$adj.P.Val,
      method = if (non_integer) "limma_voom_count_like" else "limma_voom",
      stringsAsFactors = FALSE
    )
  } else if (method == "edgeR") {
    .stop_if(isTRUE(log_scale),
             "edgeR requires count-scale data. Use expression_scale='linear' and provide raw or estimated RNA-seq counts.")
    .stop_if(any(is.na(mat)), "edgeR input contains NA values after turnover correction.")
    .stop_if(any(!is.finite(mat)), "edgeR input contains non-finite values after turnover correction.")
    .stop_if(any(mat < 0), "edgeR input contains negative values after turnover correction.")
    .stop_if(!requireNamespace("edgeR", quietly = TRUE),
             "method='edgeR' requires the Bioconductor package 'edgeR'. Install with BiocManager::install('edgeR').")
    non_integer <- any(abs(mat - round(mat)) > sqrt(.Machine$double.eps))
    if (non_integer) {
      .inform(verbose, "afterglow: edgeR received non-integer count-like values; interpret results as an approximation.")
    }
    condition <- factor(group$group, levels = c(control_group, treatment_group))
    design <- stats::model.matrix(~ 0 + condition)
    colnames(design) <- c("control", "treatment")
    dge <- edgeR::DGEList(counts = mat, group = condition)
    dge <- edgeR::calcNormFactors(dge)
    robust_fit <- requireNamespace("statmod", quietly = TRUE)
    if (!robust_fit) {
      .inform(verbose, "afterglow: edgeR robust QL fitting requires 'statmod'; using robust=FALSE.")
    }
    dge <- edgeR::estimateDisp(dge, design, robust = robust_fit)
    fit <- edgeR::glmQLFit(dge, design, robust = robust_fit)
    qlf <- edgeR::glmQLFTest(fit, contrast = c(-1, 1))
    tab <- qlf$table
    out <- data.frame(
      gene_id = rownames(mat),
      log2_fc = tab$logFC,
      average_expression = tab$logCPM,
      statistic = tab$F,
      p_value = tab$PValue,
      adjusted_p_value = .bh(tab$PValue),
      method = if (non_integer) "edgeR_QLF_count_like" else "edgeR_QLF",
      stringsAsFactors = FALSE
    )
  } else if (method == "deseq2") {
    .stop_if(isTRUE(log_scale),
             "DESeq2 requires count-scale data. Use expression_scale='linear' and provide raw or estimated RNA-seq counts.")
    .stop_if(any(is.na(mat)), "DESeq2 input contains NA values after turnover correction.")
    .stop_if(any(!is.finite(mat)), "DESeq2 input contains non-finite values after turnover correction.")
    .stop_if(any(mat < 0), "DESeq2 input contains negative values after turnover correction.")
    .stop_if(any(mat > .Machine$integer.max),
             "DESeq2 input contains values larger than the maximum R integer after turnover correction. Use method='limma_voom' or 'edgeR', or filter high-amplification genes before DESeq2.")

    count_mat <- mat
    non_integer <- any(abs(count_mat - round(count_mat)) > sqrt(.Machine$double.eps))
    if (non_integer) {
      .stop_if(!isTRUE(deseq_round_counts),
               "DESeq2 requires integer counts. Set deseq_round_counts=TRUE or use method='limma'.")
      .inform(verbose, "afterglow: rounding non-integer turnover-corrected counts for DESeq2.")
      count_mat <- round(count_mat)
    }
    count_mat[count_mat < 0] <- 0
    storage.mode(count_mat) <- "integer"

    .stop_if(!requireNamespace("DESeq2", quietly = TRUE),
             "method='deseq2' requires the Bioconductor package 'DESeq2'. Install with BiocManager::install('DESeq2').")

    col_data <- data.frame(
      condition = factor(group$group, levels = c(control_group, treatment_group)),
      row.names = group$sample
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = count_mat,
      colData = col_data,
      design = ~ condition
    )
    dds <- DESeq2::DESeq(dds, quiet = !isTRUE(verbose))
    res <- DESeq2::results(dds, contrast = c("condition", treatment_group, control_group))
    res <- as.data.frame(res)
    out <- data.frame(
      gene_id = rownames(res),
      log2_fc = res$log2FoldChange,
      average_expression = res$baseMean,
      statistic = res$stat,
      p_value = res$pvalue,
      adjusted_p_value = res$padj,
      method = if (non_integer) "DESeq2_rounded_turnover_counts" else "DESeq2",
      stringsAsFactors = FALSE
    )
  } else {
    p_value <- rep(NA_real_, nrow(mat))
    statistic <- rep(NA_real_, nrow(mat))

    for (i in seq_len(nrow(mat))) {
      x <- mat[i, treatment_idx]
      y <- mat[i, control_idx]
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
      if (length(x) >= 2L && length(y) >= 2L && stats::sd(x) + stats::sd(y) > 0) {
        tt <- stats::t.test(x = x, y = y, var.equal = FALSE)
        p_value[[i]] <- tt$p.value
        statistic[[i]] <- unname(tt$statistic)
      }
    }

    out <- data.frame(
      gene_id = rownames(mat),
      log2_fc = log_fc,
      average_expression = rowMeans(mat, na.rm = TRUE),
      statistic = statistic,
      p_value = p_value,
      adjusted_p_value = .bh(p_value),
      method = "welch_t_test",
      stringsAsFactors = FALSE
    )
  }

  out <- out[order(out$adjusted_p_value, -abs(out$log2_fc), na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  .inform(verbose, "afterglow: differential expression completed for %d genes.", nrow(out))
  .inform(verbose, "afterglow: %d genes have adjusted_p_value < 0.05.", sum(out$adjusted_p_value < 0.05, na.rm = TRUE))
  out
}

#' Full half-life-aware differential expression workflow
#'
#' @param expr Numeric expression matrix, genes x samples.
#' @param group Group table or named vector.
#' @param half_life_ref Gene-level half-life reference table.
#' @param time_hours Time from stimulus to sampling.
#' @param control_group Control group label.
#' @param treatment_group Treatment group label.
#' @param expression_scale `"auto"`, `"log2"`, or `"linear"`.
#' @param model `"pulse"` or `"synthesis_rate"`.
#' @param de_method `"limma"`, `"limma_voom"`, `"edgeR"`, `"deseq2"`, or
#'   `"welch"`.
#' @param assay_type `"auto"`, `"microarray"`, or `"rna_seq"`. This is used for
#'   method guidance and documentation in the result object.
#' @param pseudocount Pseudocount for transformations.
#' @param negative_action Handling of non-positive inferred values.
#' @param deseq_round_counts Logical. Round non-integer corrected counts when
#'   `de_method = "deseq2"`.
#' @param write_corrected_matrix Logical. If `TRUE`, write the corrected
#'   expression matrix to `corrected_matrix_file`.
#' @param corrected_matrix_file CSV path for the corrected expression matrix.
#' @param write_results Logical. If `TRUE`, write the result table to
#'   `results_file`.
#' @param results_file CSV path for the result table.
#' @param verbose Logical. If `TRUE`, print progress and result summaries.
#' @return An `afterglow_de` list with results and intermediate matrices.
#' @export
run_afterglow_de <- function(expr,
                             group,
                             half_life_ref,
                             time_hours,
                             control_group,
                             treatment_group,
                             expression_scale = c("auto", "log2", "linear"),
                             model = c("pulse", "synthesis_rate"),
                             de_method = c("limma", "limma_voom", "edgeR", "deseq2", "welch"),
                             assay_type = c("auto", "microarray", "rna_seq"),
                             pseudocount = 1e-8,
                             negative_action = c("floor", "NA", "keep"),
                             deseq_round_counts = TRUE,
                             write_corrected_matrix = FALSE,
                             corrected_matrix_file = "afterglow_corrected_expression.csv",
                             write_results = FALSE,
                             results_file = "afterglow_results.csv",
                             verbose = TRUE) {
  expression_scale <- .match_arg(expression_scale[[1L]], c("auto", "log2", "linear"), "expression_scale")
  model <- .match_arg(model[[1L]], c("pulse", "synthesis_rate"), "model")
  de_method <- .match_arg(de_method[[1L]], c("limma", "limma_voom", "edgeR", "deseq2", "welch"), "de_method")
  assay_type <- .match_arg(assay_type[[1L]], c("auto", "microarray", "rna_seq"), "assay_type")
  negative_action <- .match_arg(negative_action[[1L]], c("floor", "NA", "keep"), "negative_action")

  .inform(verbose, "afterglow: starting half-life-aware differential expression workflow.")
  .inform(verbose, "afterglow: for microarray/chip or log2 normalized expression, use de_method='limma'.")
  .inform(verbose, "afterglow: for RNA-seq raw counts, use de_method='limma_voom', 'edgeR' or 'deseq2' with expression_scale='linear'.")

  inferred <- infer_turnover_expression(
    expr = expr,
    group = group,
    half_life_ref = half_life_ref,
    time_hours = time_hours,
    control_group = control_group,
    treatment_group = treatment_group,
    expression_scale = expression_scale,
    model = model,
    pseudocount = pseudocount,
    negative_action = negative_action,
    verbose = verbose
  )

  de_log_scale <- inferred$settings$expression_scale == "log2"

  results <- differential_expression(
    mat = inferred$corrected_matrix,
    group = inferred$group,
    control_group = control_group,
    treatment_group = treatment_group,
    method = de_method,
    assay_type = assay_type,
    log_scale = de_log_scale,
    pseudocount = pseudocount,
    deseq_round_counts = deseq_round_counts,
    verbose = verbose
  )

  results <- merge(results, inferred$diagnostics, by = "gene_id", all.x = TRUE, sort = FALSE)
  results <- results[order(results$adjusted_p_value, -abs(results$log2_fc), na.last = TRUE), , drop = FALSE]
  rownames(results) <- NULL

  if (isTRUE(write_corrected_matrix)) {
    corrected_out <- data.frame(gene_id = rownames(inferred$corrected_matrix), inferred$corrected_matrix, check.names = FALSE)
    utils::write.csv(corrected_out, corrected_matrix_file, row.names = FALSE)
    .inform(verbose, "afterglow: wrote corrected expression matrix to %s.", corrected_matrix_file)
  }
  if (isTRUE(write_results)) {
    utils::write.csv(results, results_file, row.names = FALSE)
    .inform(verbose, "afterglow: wrote differential expression results to %s.", results_file)
  }

  out <- c(inferred, list(results = results))
  out$settings$de_method <- de_method
  out$settings$assay_type <- assay_type
  out$settings$deseq_round_counts <- deseq_round_counts
  out$settings$corrected_matrix_file <- if (isTRUE(write_corrected_matrix)) corrected_matrix_file else NA_character_
  out$settings$results_file <- if (isTRUE(write_results)) results_file else NA_character_
  class(out) <- "afterglow_de"

  .inform(verbose, "afterglow: workflow finished.")
  out
}

#' Print method for afterglow_de objects
#'
#' @param x An `afterglow_de` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.afterglow_de <- function(x, ...) {
  cat("afterglow result\n")
  cat("  turnover model:", x$settings$model, "\n")
  cat("  differential expression method:", x$settings$de_method, "\n")
  cat("  assay type:", x$settings$assay_type, "\n")
  cat("  time_hours:", x$settings$time_hours, "\n")
  cat("  expression_scale:", x$settings$expression_scale, "\n")
  cat("  genes tested:", nrow(x$results), "\n")
  cat("  significant genes at adjusted_p_value < 0.05:", sum(x$results$adjusted_p_value < 0.05, na.rm = TRUE), "\n")
  cat("  non-positive inferred values before handling:", x$settings$n_non_positive_inferred, "\n")
  if (!is.na(x$settings$corrected_matrix_file)) {
    cat("  corrected expression matrix:", x$settings$corrected_matrix_file, "\n")
  }
  if (!is.na(x$settings$results_file)) {
    cat("  results table:", x$settings$results_file, "\n")
  }
  invisible(x)
}
