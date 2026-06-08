# Confidence scoring for afterglow results.
# Quantifies per-gene reliability of the inferred transcriptional fold-change.

#' Estimate measurement CV from replicate expression data
#'
#' Computes the median coefficient of variation across genes within each group,
#' giving users an empirical CV to feed into \code{\link{compute_confidence}}.
#'
#' @param expr Numeric matrix (genes x samples), linear scale.
#' @param group Data frame with \code{sample} and \code{group} columns, or a
#'   named character vector.
#' @return Named numeric vector of median CV per group, plus an overall
#'   \code{"pooled"} estimate (median across all groups).
#' @details
#' CV is defined as sd/mean for each gene within each group. The function
#' reports the median across genes (robust to outliers). Typical values:
#' \itemize{
#'   \item Technical replicates: 0.01-0.03
#'   \item Biological replicates (bulk RNA-seq, n=3): 0.05-0.15
#'   \item Single-cell pseudobulk: 0.10-0.30
#' }
#' @export
estimate_cv <- function(expr, group) {
  if (is.vector(group) && !is.data.frame(group)) {
    group <- data.frame(sample = names(group), group = as.character(group),
                        stringsAsFactors = FALSE)
  }
  expr <- as.matrix(expr)
  grps <- unique(group$group)
  cv_per_group <- vapply(grps, function(g) {
    idx <- which(group$group == g)
    if (length(idx) < 2L) return(NA_real_)
    means <- rowMeans(expr[, idx, drop = FALSE], na.rm = TRUE)
    sds <- apply(expr[, idx, drop = FALSE], 1, sd, na.rm = TRUE)
    # Only compute CV for genes with non-zero mean
    valid <- means > 0
    median(sds[valid] / means[valid], na.rm = TRUE)
  }, numeric(1L))
  names(cv_per_group) <- grps
  c(cv_per_group, pooled = median(cv_per_group, na.rm = TRUE))
}

#' Compute per-gene afterglow confidence scores
#'
#' Assigns a confidence score (0-1) and tier to each gene based on how much
#' the correction amplifies measurement noise. Genes with high correction
#' factors (short half-life relative to experiment duration) receive lower
#' confidence because small measurement errors are magnified.
#'
#' @param half_life_hours Numeric vector of mRNA half-lives in hours.
#' @param time_hours Scalar: time between stimulus and RNA measurement.
#' @param observed_fc Optional numeric vector of observed log2 fold-changes.
#'   When provided, signal-to-noise is factored into the score.
#' @param n_records Optional integer vector: number of TTDB records per gene.
#'   More records indicate a more reliable half-life estimate.
#' @param cv_measurement Assumed coefficient of variation of expression
#'   measurements (default 0.05).
#' @param model Kinetic model used to define the error amplification factor:
#'   `"pulse"` uses `exp(lambda * time_hours)`, while `"synthesis_rate"`
#'   uses `1 / (1 - exp(-lambda * time_hours))`.
#' @return Data frame with columns: amplification_factor, correction_factor
#'   (backward-compatible alias), model, confidence_tier, confidence_score, snr,
#'   n_records, score_cf, score_snr, score_evidence.
#' @details
#' The confidence score combines three components:
#' \itemize{
#'   \item \strong{Model-specific amplification factor} (weight 0.6):
#'     pulse uses CF = exp(ln2 * t / hl), whereas synthesis_rate uses
#'     CF = 1 / (1 - exp(-ln2 * t / hl)). This is the primary determinant:
#'     it controls how much measurement error is amplified during model inversion.
#'   \item \strong{Signal-to-noise} (weight 0.25): Whether the observed FC is
#'     large enough to exceed the amplified noise floor.
#'   \item \strong{Half-life evidence} (weight 0.15): Number of independent
#'     measurements supporting the half-life estimate.
#' }
#'
#' Confidence tiers (empirically calibrated from simulation):
#' \itemize{
#'   \item High (CF 1-2): error amplification <2x, recovery r > 0.99
#'   \item Good (CF 2-5): amplification 2-5x, recovery r > 0.97
#'   \item Moderate (CF 5-10): amplification 5-10x, recovery r > 0.90
#'   \item Low (CF 10-50): amplification 10-50x, interpret with caution
#'   \item Unreliable (CF > 50): should be filtered from results
#' }
#' @export
compute_confidence <- function(half_life_hours,
                               time_hours,
                               observed_fc = NULL,
                               n_records = NULL,
                               cv_measurement = 0.05,
                               model = c("pulse", "synthesis_rate")) {
  model <- .match_arg(model[[1L]], c("pulse", "synthesis_rate"), "model")
  n <- length(half_life_hours)
  .stop_if(length(time_hours) != 1L || !is.finite(time_hours) || time_hours <= 0,
           "time_hours must be a single positive finite number.")
  .stop_if(any(half_life_hours <= 0, na.rm = TRUE),
           "half_life_hours must be positive.")

  lambda <- log(2) / half_life_hours
  decay_fraction <- exp(-lambda * time_hours)
  CF <- if (model == "pulse") {
    exp(lambda * time_hours)
  } else {
    1 / (1 - decay_fraction)
  }

  # Confidence tier based on model-specific amplification factor

  tier <- ifelse(CF <= 2, "High",
           ifelse(CF <= 5, "Good",
            ifelse(CF <= 10, "Moderate",
             ifelse(CF <= 50, "Low", "Unreliable"))))

  # Numeric score from amplification factor: logistic decay
  # Calibrated: CF=1 -> ~1.0, CF=5 -> ~0.7, CF=10 -> 0.5, CF=50 -> ~0.1
  CF_half <- 10
  k <- 1.5
  score_cf <- 1 / (1 + (CF / CF_half)^k)


  # Signal-to-noise component
  if (!is.null(observed_fc)) {
    .stop_if(length(observed_fc) != n, "observed_fc must have same length as half_life_hours.")
    noise_log2 <- CF * cv_measurement / log(2)
    snr <- abs(observed_fc) / pmax(noise_log2, 0.01)
    score_snr <- 1 / (1 + exp(-2 * (snr - 2)))
  } else {
    score_snr <- rep(0.5, n)  # neutral when not provided
    snr <- rep(NA_real_, n)
  }

  # Half-life evidence component
  if (!is.null(n_records)) {
    .stop_if(length(n_records) != n, "n_records must have same length as half_life_hours.")
    n_records[is.na(n_records)] <- 1L
    score_evidence <- pmin(1, sqrt(n_records / 5))
  } else {
    score_evidence <- rep(0.5, n)  # neutral when not provided
  }

  # Combined score (weighted)
  confidence <- 0.6 * score_cf + 0.25 * score_snr + 0.15 * score_evidence
  confidence <- pmin(1, pmax(0, confidence))

  data.frame(
    model = model,
    amplification_factor = CF,
    correction_factor = CF,
    confidence_tier = factor(tier, levels = c("High", "Good", "Moderate", "Low", "Unreliable")),
    confidence_score = confidence,
    snr = snr,
    n_records = if (!is.null(n_records)) n_records else rep(NA_integer_, n),
    score_cf = score_cf,
    score_snr = score_snr,
    score_evidence = score_evidence,
    stringsAsFactors = FALSE
  )
}

#' Add confidence scores to afterglow results
#'
#' Convenience wrapper that extracts half-life and time information from an
#' afterglow result object and appends confidence scores to the diagnostics.
#'
#' @param ag_result List returned by \code{\link{infer_turnover_expression}}.
#' @param observed_fc Optional log2 FC vector used for the signal-to-noise
#'   component of the confidence score. If NULL, the raw observed treatment-vs-
#'   control log2FC is computed from the uncorrected matrix. The returned table
#'   always includes both `raw_observed_log2_fc` and `corrected_log2_fc` so users
#'   can see which signal entered the SNR calculation.
#' @param cv_measurement Assumed measurement CV (default 0.05).
#' @return The input ag_result with an added \code{confidence} data frame.
#' @export
add_confidence <- function(ag_result, observed_fc = NULL, cv_measurement = 0.05) {
  .stop_if(is.null(ag_result$diagnostics), "ag_result must contain diagnostics (from infer_turnover_expression).")
  .stop_if(is.null(ag_result$settings$time_hours), "ag_result must contain settings$time_hours.")

  hl <- ag_result$diagnostics$half_life_hours
  t_h <- ag_result$settings$time_hours
  model <- ag_result$settings$model
  n_rec <- ag_result$half_life_ref$n_records  # may be NULL

  ctrl_idx <- which(ag_result$group$group == ag_result$settings$control_group)
  trt_idx <- which(ag_result$group$group == ag_result$settings$treatment_group)

  raw_ctrl_mean <- rowMeans(ag_result$original_linear[, ctrl_idx, drop = FALSE], na.rm = TRUE)
  raw_trt_mean <- rowMeans(ag_result$original_linear[, trt_idx, drop = FALSE], na.rm = TRUE)
  raw_observed_fc <- log2((raw_trt_mean + 1e-8) / (raw_ctrl_mean + 1e-8))

  corrected_ctrl_mean <- rowMeans(ag_result$corrected_linear[, ctrl_idx, drop = FALSE], na.rm = TRUE)
  corrected_trt_mean <- rowMeans(ag_result$corrected_linear[, trt_idx, drop = FALSE], na.rm = TRUE)
  corrected_fc <- log2((corrected_trt_mean + 1e-8) / (corrected_ctrl_mean + 1e-8))

  confidence_input_fc <- if (is.null(observed_fc)) raw_observed_fc else observed_fc

  conf <- compute_confidence(
    half_life_hours = hl,
    time_hours = t_h,
    observed_fc = confidence_input_fc,
    n_records = n_rec,
    cv_measurement = cv_measurement,
    model = model
  )
  conf$gene_id <- ag_result$diagnostics$gene_id
  conf$raw_observed_log2_fc <- raw_observed_fc
  conf$corrected_log2_fc <- corrected_fc
  conf$confidence_input_log2_fc <- confidence_input_fc
  conf <- conf[, c("gene_id", "model", "amplification_factor", "correction_factor",
                   "confidence_tier", "confidence_score", "snr", "n_records",
                   "raw_observed_log2_fc", "corrected_log2_fc", "confidence_input_log2_fc",
                   "score_cf", "score_snr", "score_evidence")]

  ag_result$confidence <- conf
  ag_result
}
