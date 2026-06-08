# Half-life perturbation sensitivity analysis.

#' Half-life perturbation sensitivity analysis
#'
#' Re-runs the turnover inversion after perturbing each gene's half-life estimate
#' by log-normal noise. This gives a reproducible, user-facing way to quantify
#' how stable corrected fold-changes are to TTDB or user-supplied half-life
#' uncertainty.
#'
#' @param expr Numeric matrix, genes x samples. Can be linear expression or
#'   log2-scale expression.
#' @param group Data frame with `sample` and `group`, or a named group vector.
#' @param half_life_ref Data frame with `gene_id` and `half_life_hours`.
#' @param time_hours Time from stimulus at t0 to sampling at t1, in hours.
#' @param control_group Control group label.
#' @param treatment_group Stimulated/treatment group label.
#' @param expression_scale `"linear"`, `"log2"`, or `"auto"`.
#' @param model `"pulse"` or `"synthesis_rate"`.
#' @param perturbation_sd Numeric vector of log-normal standard deviations.
#'   Values of 0.05, 0.10, 0.20 and 0.50 approximately represent 5%, 10%, 20%
#'   and 50% multiplicative uncertainty for small perturbations.
#' @param n_iter Number of random perturbation replicates per non-zero
#'   perturbation level.
#' @param seed Optional random seed for reproducibility.
#' @param truth_log2_fc Optional named numeric vector of true log2 fold-changes,
#'   used in simulation studies to report recovery metrics.
#' @param pseudocount Pseudocount used for fold-change calculations.
#' @param negative_action Handling of non-positive inferred expression.
#' @param verbose Logical. If `TRUE`, print progress messages.
#' @return An `afterglow_sensitivity` list with `summary`, `gene_level`, and
#'   `settings` components.
#' @export
half_life_sensitivity <- function(expr,
                                  group,
                                  half_life_ref,
                                  time_hours,
                                  control_group,
                                  treatment_group,
                                  expression_scale = c("auto", "log2", "linear"),
                                  model = c("pulse", "synthesis_rate"),
                                  perturbation_sd = c(0, 0.05, 0.10, 0.20),
                                  n_iter = 20L,
                                  seed = NULL,
                                  truth_log2_fc = NULL,
                                  pseudocount = 1e-8,
                                  negative_action = c("floor", "NA", "keep"),
                                  verbose = TRUE) {
  expression_scale <- .match_arg(expression_scale[[1L]], c("auto", "log2", "linear"), "expression_scale")
  model <- .match_arg(model[[1L]], c("pulse", "synthesis_rate"), "model")
  negative_action <- .match_arg(negative_action[[1L]], c("floor", "NA", "keep"), "negative_action")
  .stop_if(length(time_hours) != 1L || !is.finite(time_hours) || time_hours < 0,
           "time_hours must be one finite non-negative number.")
  .stop_if(length(n_iter) != 1L || is.na(n_iter) || n_iter < 1L,
           "n_iter must be a positive integer.")
  .stop_if(any(!is.finite(perturbation_sd)) || any(perturbation_sd < 0),
           "perturbation_sd must contain non-negative finite values.")
  n_iter <- as.integer(n_iter)

  if (!is.null(seed)) set.seed(seed)

  nominal <- infer_turnover_expression(
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

  lfc_from_ag <- function(ag) {
    ctrl_idx <- which(ag$group$group == ag$settings$control_group)
    trt_idx <- which(ag$group$group == ag$settings$treatment_group)
    ctrl_mean <- rowMeans(ag$corrected_linear[, ctrl_idx, drop = FALSE], na.rm = TRUE)
    trt_mean <- rowMeans(ag$corrected_linear[, trt_idx, drop = FALSE], na.rm = TRUE)
    out <- log2((trt_mean + pseudocount) / (ctrl_mean + pseudocount))
    names(out) <- rownames(ag$corrected_linear)
    out
  }

  nominal_lfc <- lfc_from_ag(nominal)
  base_ref <- nominal$half_life_ref
  base_ref$half_life_hours <- .safe_numeric(base_ref$half_life_hours)

  if (!is.null(truth_log2_fc)) {
    .stop_if(is.null(names(truth_log2_fc)), "truth_log2_fc must be a named numeric vector.")
    truth_log2_fc <- as.numeric(truth_log2_fc[base_ref$gene_id])
  }

  rows <- list()
  k <- 1L
  levels <- unique(perturbation_sd)
  for (sd_i in levels) {
    reps <- if (sd_i == 0) 1L else n_iter
    for (iter_i in seq_len(reps)) {
      ref_i <- base_ref
      if (sd_i > 0) {
        ref_i$half_life_hours <- pmax(ref_i$half_life_hours * exp(stats::rnorm(nrow(ref_i), 0, sd_i)), .Machine$double.eps)
      }
      ag_i <- infer_turnover_expression(
        expr = nominal$original_linear,
        group = nominal$group,
        half_life_ref = ref_i[, intersect(c("gene_id", "half_life_hours", "n_records"), names(ref_i)), drop = FALSE],
        time_hours = time_hours,
        control_group = control_group,
        treatment_group = treatment_group,
        expression_scale = "linear",
        model = model,
        pseudocount = pseudocount,
        negative_action = negative_action,
        verbose = FALSE
      )
      lfc_i <- lfc_from_ag(ag_i)
      gene_df <- data.frame(
        perturbation_sd = sd_i,
        iteration = iter_i,
        gene_id = base_ref$gene_id,
        half_life_hours = base_ref$half_life_hours,
        perturbed_half_life_hours = ref_i$half_life_hours,
        nominal_log2_fc = nominal_lfc[base_ref$gene_id],
        perturbed_log2_fc = lfc_i[base_ref$gene_id],
        model_amplification_factor = ag_i$diagnostics$model_amplification_factor,
        stringsAsFactors = FALSE
      )
      gene_df$delta_log2_fc <- gene_df$perturbed_log2_fc - gene_df$nominal_log2_fc
      if (!is.null(truth_log2_fc)) {
        gene_df$truth_log2_fc <- truth_log2_fc
        gene_df$error_vs_truth <- gene_df$perturbed_log2_fc - gene_df$truth_log2_fc
      }
      rows[[k]] <- gene_df
      k <- k + 1L
    }
  }

  gene_level <- do.call(rbind, rows)
  summary_rows <- lapply(split(gene_level, gene_level$perturbation_sd), function(x) {
    ok <- is.finite(x$perturbed_log2_fc) & is.finite(x$nominal_log2_fc)
    out <- data.frame(
      perturbation_sd = unique(x$perturbation_sd),
      n_iter = length(unique(x$iteration)),
      n_genes = length(unique(x$gene_id)),
      median_abs_delta_log2_fc = stats::median(abs(x$delta_log2_fc[ok]), na.rm = TRUE),
      p95_abs_delta_log2_fc = as.numeric(stats::quantile(abs(x$delta_log2_fc[ok]), 0.95, na.rm = TRUE)),
      pearson_vs_nominal = stats::cor(x$perturbed_log2_fc[ok], x$nominal_log2_fc[ok], method = "pearson"),
      spearman_vs_nominal = stats::cor(x$perturbed_log2_fc[ok], x$nominal_log2_fc[ok], method = "spearman"),
      stringsAsFactors = FALSE
    )
    if ("truth_log2_fc" %in% names(x)) {
      ok_truth <- ok & is.finite(x$truth_log2_fc)
      out$pearson_vs_truth <- stats::cor(x$perturbed_log2_fc[ok_truth], x$truth_log2_fc[ok_truth], method = "pearson")
      out$rmse_vs_truth <- sqrt(mean((x$perturbed_log2_fc[ok_truth] - x$truth_log2_fc[ok_truth])^2))
    }
    out
  })
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL

  out <- list(
    summary = summary,
    gene_level = gene_level,
    nominal = nominal,
    settings = list(
      time_hours = time_hours,
      control_group = control_group,
      treatment_group = treatment_group,
      expression_scale = nominal$settings$expression_scale,
      model = model,
      perturbation_sd = levels,
      n_iter = n_iter,
      seed = seed,
      negative_action = negative_action
    )
  )
  class(out) <- "afterglow_sensitivity"
  out
}

#' Print method for afterglow_sensitivity objects
#'
#' @param x An `afterglow_sensitivity` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.afterglow_sensitivity <- function(x, ...) {
  cat("afterglow half-life sensitivity result\n")
  cat("  turnover model:", x$settings$model, "\n")
  cat("  time_hours:", x$settings$time_hours, "\n")
  cat("  perturbation levels:", paste(x$settings$perturbation_sd, collapse = ", "), "\n")
  cat("  iterations per non-zero level:", x$settings$n_iter, "\n")
  cat("  genes tested:", length(unique(x$gene_level$gene_id)), "\n")
  print(x$summary)
  invisible(x)
}
