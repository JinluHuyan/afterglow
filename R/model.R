# Mathematical core of afterglow.
# The equations are kept explicit rather than algebraically compressed so users
# can audit the model assumptions and modify individual steps if needed.

#' Validate expression, group, and half-life inputs
#'
#' @param expr Numeric matrix, genes x samples.
#' @param group Data frame with sample and group columns, or named vector.
#' @param half_life_ref Data frame with `gene_id` and `half_life_hours`.
#' @param control_group Control group label.
#' @param treatment_group Treatment/stimulus group label.
#' @return Normalized list used internally by `run_afterglow_de()`.
#' @export
validate_afterglow_inputs <- function(expr,
                                      group,
                                      half_life_ref,
                                      control_group,
                                      treatment_group) {
  .stop_if(!is.matrix(expr) && !is.data.frame(expr), "expr must be a numeric matrix or data frame.")
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  .stop_if(is.null(rownames(expr)), "expr must have gene IDs as row names.")
  .stop_if(is.null(colnames(expr)), "expr must have sample IDs as column names.")
  .stop_if(anyDuplicated(rownames(expr)) > 0L, "expr row names must be unique gene IDs.")
  .stop_if(anyDuplicated(colnames(expr)) > 0L, "expr column names must be unique sample IDs.")

  if (is.vector(group) && !is.data.frame(group)) {
    .stop_if(is.null(names(group)), "Named group vector must use sample names as names(group).")
    group <- data.frame(sample = names(group), group = as.character(group), stringsAsFactors = FALSE)
  }
  .stop_if(!all(c("sample", "group") %in% colnames(group)), "group must contain sample and group columns.")

  missing_samples <- setdiff(colnames(expr), group$sample)
  .stop_if(length(missing_samples) > 0L,
           sprintf("Group table is missing samples: %s", paste(missing_samples, collapse = ", ")))

  group <- group[match(colnames(expr), group$sample), , drop = FALSE]
  .stop_if(!control_group %in% group$group, sprintf("control_group '%s' is not present.", control_group))
  .stop_if(!treatment_group %in% group$group, sprintf("treatment_group '%s' is not present.", treatment_group))

  .stop_if(!all(c("gene_id", "half_life_hours") %in% colnames(half_life_ref)),
           "half_life_ref must contain gene_id and half_life_hours columns.")

  half_life_ref$gene_id <- as.character(half_life_ref$gene_id)
  half_life_ref$half_life_hours <- .safe_numeric(half_life_ref$half_life_hours)
  half_life_ref <- half_life_ref[!is.na(half_life_ref$half_life_hours) & half_life_ref$half_life_hours > 0, , drop = FALSE]
  .stop_if(nrow(half_life_ref) == 0L, "No positive half-life values in half_life_ref.")
  .stop_if(anyDuplicated(half_life_ref$gene_id) > 0L, "half_life_ref must have one row per gene_id.")

  common <- intersect(rownames(expr), half_life_ref$gene_id)
  .stop_if(length(common) == 0L, "No overlap between expression gene IDs and half_life_ref$gene_id.")

  list(
    expr = expr[common, , drop = FALSE],
    group = group,
    half_life_ref = half_life_ref[match(common, half_life_ref$gene_id), , drop = FALSE],
    control_idx = which(group$group == control_group),
    treatment_idx = which(group$group == treatment_group),
    common_genes = common
  )
}

#' Infer stimulus-time expression using mRNA half-life
#'
#' @param expr Numeric matrix, genes x samples. Can be linear expression or
#'   log2-scale expression.
#' @param group Data frame with `sample` and `group`, or a named group vector.
#' @param half_life_ref Data frame with `gene_id` and `half_life_hours`.
#' @param time_hours Time from stimulus at t0 to sampling at t1, in hours.
#' @param control_group Control group label.
#' @param treatment_group Stimulated/treatment group label.
#' @param expression_scale `"linear"`, `"log2"`, or `"auto"`. Auto treats data
#'   as log2 when all values are non-negative and the 95th percentile is below
#'   50; for raw counts, use `"linear"` explicitly.
#' @param model `"pulse"` or `"synthesis_rate"`.
#' @param pseudocount Pseudocount used before log2 transformation and for
#'   denominator protection.
#' @param negative_action How to handle inferred non-positive expression:
#'   `"floor"` replaces with `pseudocount`, `"NA"` marks as missing, and
#'   `"keep"` keeps mathematically inferred values.
#' @param verbose Logical. If `TRUE`, print progress and model diagnostics.
#' @return List containing corrected matrix, linear-scale intermediate values,
#'   matched half-life table, and diagnostics.
#' @export
infer_turnover_expression <- function(expr,
                                      group,
                                      half_life_ref,
                                      time_hours,
                                      control_group,
                                      treatment_group,
                                      expression_scale = c("auto", "log2", "linear"),
                                      model = c("pulse", "synthesis_rate"),
                                      pseudocount = 1e-8,
                                      negative_action = c("floor", "NA", "keep"),
                                      verbose = TRUE) {
  expression_scale <- .match_arg(expression_scale[[1L]], c("auto", "log2", "linear"), "expression_scale")
  model <- .match_arg(model[[1L]], c("pulse", "synthesis_rate"), "model")
  negative_action <- .match_arg(negative_action[[1L]], c("floor", "NA", "keep"), "negative_action")

  .stop_if(length(time_hours) != 1L || !is.finite(time_hours) || time_hours < 0,
           "time_hours must be one finite non-negative number.")
  .stop_if(pseudocount <= 0, "pseudocount must be positive.")

  .inform(verbose, "afterglow: validating expression, group, and half-life inputs.")
  x <- validate_afterglow_inputs(expr, group, half_life_ref, control_group, treatment_group)
  expr <- x$expr
  group <- x$group
  ref <- x$half_life_ref
  control_idx <- x$control_idx
  treatment_idx <- x$treatment_idx

  if (expression_scale == "auto") {
    # Conservative scale detection. Values like 14-16 are typical of log2
    # microarray or normalized expression data. Raw RNA-seq counts should be
    # supplied with expression_scale = "linear".
    q95 <- as.numeric(stats::quantile(expr, probs = 0.95, na.rm = TRUE))
    expression_scale <- if (min(expr, na.rm = TRUE) >= 0 && q95 < 50) "log2" else "linear"
    .inform(verbose, "afterglow: expression_scale='auto' resolved to '%s'.", expression_scale)
  }

  expr_linear <- if (expression_scale == "log2") {
    # Do not subtract pseudocount during inverse transformation because common
    # log2 intensities are not guaranteed to be log2(x + pseudocount). Treat the
    # input as log2 relative abundance.
    2^expr
  } else {
    expr
  }

  .stop_if(any(expr_linear < 0, na.rm = TRUE),
           "Linear expression contains negative values. Use log2 scale or transform your data first.")

  lambda <- log(2) / ref$half_life_hours
  decay_fraction <- exp(-lambda * time_hours)
  growth_factor <- exp(lambda * time_hours)
  pulse_amplification_factor <- growth_factor
  synthesis_rate_amplification_factor <- if (time_hours == 0) {
    rep(Inf, length(decay_fraction))
  } else {
    1 / (1 - decay_fraction)
  }
  model_amplification_factor <- if (model == "pulse") {
    pulse_amplification_factor
  } else {
    synthesis_rate_amplification_factor
  }

  control_baseline <- rowMeans(expr_linear[, control_idx, drop = FALSE], na.rm = TRUE)

  corrected_linear <- expr_linear
  .inform(
    verbose,
    "afterglow: applying '%s' turnover model to %d matched genes over %.4g hour(s).",
    model,
    nrow(expr),
    time_hours
  )

  if (model == "pulse") {
    # Pulse model:
    #   B_g = control steady-state expression of gene g
    #   P_{g,s}(t0) = instantaneous net pulse caused by stimulus in treatment sample s
    #   E_{g,s}(t1) = B_g + P_{g,s}(t0) * exp(-lambda_g * delta_t)
    # Solving for the stimulus-time pulse gives:
    #   P_{g,s}(t0) = (E_{g,s}(t1) - B_g) * exp(lambda_g * delta_t)
    #   E_{g,s}^{inferred}(t0 after pulse) = B_g + P_{g,s}(t0)
    for (j in treatment_idx) {
      corrected_linear[, j] <- control_baseline + (expr_linear[, j] - control_baseline) * growth_factor
    }
  } else {
    # Synthesis-rate model:
    #   dM/dt = k_s - lambda * M
    #   M(t1) = k_s/lambda + (M(t0) - k_s/lambda) * exp(-lambda * delta_t)
    # Assuming pre-stimulus M(t0)=B_g, solve for the post-stimulus synthesis rate:
    #   k_s = lambda * (M(t1) - B_g * exp(-lambda * delta_t)) / (1 - exp(-lambda * delta_t))
    # To keep the output on an expression-like scale, report equivalent steady
    # expression k_s/lambda:
    #   M_equiv = (M(t1) - B_g * exp(-lambda * delta_t)) / (1 - exp(-lambda * delta_t))
    .stop_if(time_hours == 0, "synthesis_rate model is undefined when time_hours is 0.")
    denom <- 1 - decay_fraction
    for (j in treatment_idx) {
      corrected_linear[, j] <- (expr_linear[, j] - control_baseline * decay_fraction) / denom
    }
  }

  n_non_positive <- sum(corrected_linear <= 0, na.rm = TRUE)
  if (n_non_positive > 0L && negative_action == "floor") {
    corrected_linear[corrected_linear <= 0] <- pseudocount
  } else if (n_non_positive > 0L && negative_action == "NA") {
    corrected_linear[corrected_linear <= 0] <- NA_real_
  }
  .inform(verbose, "afterglow: detected %d non-positive inferred value(s) before negative_action handling.", n_non_positive)

  corrected_matrix <- if (expression_scale == "log2") {
    log2(corrected_linear + pseudocount)
  } else {
    corrected_linear
  }

  rownames(corrected_matrix) <- rownames(expr)
  colnames(corrected_matrix) <- colnames(expr)

  diagnostics <- data.frame(
    gene_id = rownames(expr),
    half_life_hours = ref$half_life_hours,
    decay_rate_per_hour = lambda,
    decay_fraction_t1_over_t0 = decay_fraction,
    correction_growth_factor = growth_factor,
    pulse_amplification_factor = pulse_amplification_factor,
    synthesis_rate_amplification_factor = synthesis_rate_amplification_factor,
    model_amplification_factor = model_amplification_factor,
    control_baseline_linear = control_baseline,
    stringsAsFactors = FALSE
  )

  list(
    corrected_matrix = corrected_matrix,
    corrected_linear = corrected_linear,
    original_matrix = expr,
    original_linear = expr_linear,
    group = group,
    half_life_ref = ref,
    diagnostics = diagnostics,
    settings = list(
      time_hours = time_hours,
      control_group = control_group,
      treatment_group = treatment_group,
      expression_scale = expression_scale,
      model = model,
      pseudocount = pseudocount,
      negative_action = negative_action,
      n_non_positive_inferred = n_non_positive
    )
  )
}

# Backward-compatible internal alias used by early local scripts. It is not exported.
validate_turnover_inputs <- validate_afterglow_inputs
