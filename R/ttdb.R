# TTDB processing helpers.
# Raw TTDB stability tables can contain multiple records per gene across
# studies, samples, cell types, and conditions. This file converts those records
# into a gene-level half-life reference suitable for afterglow modeling.

#' Build a gene-level half-life reference table
#'
#' @param stability_data A TTDB stability data frame, or a list returned by
#'   `read_ttdb_downloads()`. If a list is supplied, the function prefers
#'   `species_stability_combined`, then `species_stability`.
#' @param species Optional species filter, e.g. `"Human"` or `"Mouse"`.
#' @param cell_type Optional cell type filter. Uses exact matching.
#' @param condition Optional condition filter. Uses exact matching.
#' @param technique Optional technique filter. Uses exact matching.
#' @param gene_col_preference Candidate gene identifier columns, in priority
#'   order.
#' @param half_life_col Column containing half-life values in hours.
#' @param decay_rate_col Column containing decay rates per hour.
#' @param aggregate One of `"median"`, `"mean"`, or `"geometric_mean"`. Median
#'   is the default because TTDB integrates heterogeneous studies and half-life
#'   estimates can be skewed.
#' @param min_r_squared Optional lower bound for `r_squared`, when available.
#' @param use_decay_rate_if_half_life_missing If `TRUE`, missing half-life is
#'   computed as `log(2) / decay_rate` when a positive decay rate exists.
#' @param collapse_case Either `"as_is"` or `"upper"` for gene IDs.
#' @param verbose Logical. If `TRUE`, print progress and filtering diagnostics.
#' @return Data frame with one row per gene and columns:
#'   `gene_id`, `half_life_hours`, `decay_rate_per_hour`, `n_records`.
#' @export
make_half_life_reference <- function(stability_data,
                                     species = NULL,
                                     cell_type = NULL,
                                     condition = NULL,
                                     technique = NULL,
                                     gene_col_preference = c("gene_name_x", "gene_name_y", "ensembl_id", "gene_id"),
                                     half_life_col = "half_life",
                                     decay_rate_col = "decay_rate",
                                     aggregate = c("median", "mean", "geometric_mean"),
                                     min_r_squared = NULL,
                                     use_decay_rate_if_half_life_missing = TRUE,
                                     collapse_case = c("as_is", "upper"),
                                     verbose = TRUE) {
  aggregate <- .match_arg(aggregate[[1L]], c("median", "mean", "geometric_mean"), "aggregate")
  collapse_case <- .match_arg(collapse_case[[1L]], c("as_is", "upper"), "collapse_case")

  if (is.list(stability_data) && !is.data.frame(stability_data)) {
    if ("species_stability_combined" %in% names(stability_data)) {
      stability_data <- stability_data$species_stability_combined
    } else if ("species_stability" %in% names(stability_data)) {
      stability_data <- stability_data$species_stability
    } else {
      stop("TTDB list must contain species_stability_combined or species_stability.", call. = FALSE)
    }
  }

  .stop_if(!is.data.frame(stability_data), "stability_data must be a data frame or TTDB list.")

  gene_col <- .first_existing_col(stability_data, gene_col_preference, required = TRUE, what = "gene identifier column")
  .stop_if(!half_life_col %in% colnames(stability_data) && !decay_rate_col %in% colnames(stability_data),
           "stability_data must contain half_life and/or decay_rate.")

  dat <- stability_data
  n_before_filtering <- nrow(dat)

  # Exact matching is used by default to avoid unintended matches across species,
  # cell types, or experimental conditions.
  if (!is.null(species) && "species_name" %in% colnames(dat)) {
    dat <- dat[dat$species_name %in% species, , drop = FALSE]
  }
  if (!is.null(cell_type) && "cell_type" %in% colnames(dat)) {
    dat <- dat[dat$cell_type %in% cell_type, , drop = FALSE]
  }
  if (!is.null(condition) && "condition" %in% colnames(dat)) {
    dat <- dat[dat$condition %in% condition, , drop = FALSE]
  }
  if (!is.null(technique) && "technique" %in% colnames(dat)) {
    dat <- dat[dat$technique %in% technique, , drop = FALSE]
  }
  if (!is.null(min_r_squared) && "r_squared" %in% colnames(dat)) {
    r2 <- .safe_numeric(dat$r_squared)
    dat <- dat[is.na(r2) | r2 >= min_r_squared, , drop = FALSE]
  }

  .stop_if(nrow(dat) == 0L, "No TTDB rows remain after filtering.")

  gene_id <- .collapse_gene_id(dat[[gene_col]], collapse_case)
  half_life <- if (half_life_col %in% colnames(dat)) .safe_numeric(dat[[half_life_col]]) else rep(NA_real_, nrow(dat))
  decay_rate <- if (decay_rate_col %in% colnames(dat)) .safe_numeric(dat[[decay_rate_col]]) else rep(NA_real_, nrow(dat))

  # TTDB defines t_1/2 = ln(2) / k. Missing half-life values can therefore be
  # inferred from a positive decay rate. Zero or negative decay rates are not
  # physically meaningful in this first-order decay model.
  if (isTRUE(use_decay_rate_if_half_life_missing)) {
    needs_fill <- is.na(half_life) & !is.na(decay_rate) & decay_rate > 0
    half_life[needs_fill] <- log(2) / decay_rate[needs_fill]
  }

  keep <- !is.na(gene_id) & gene_id != "" & !is.na(half_life) & half_life > 0
  dat2 <- data.frame(
    gene_id = gene_id[keep],
    half_life_hours = half_life[keep],
    stringsAsFactors = FALSE
  )

  .stop_if(nrow(dat2) == 0L, "No valid positive half-life values are available.")

  split_hl <- split(dat2$half_life_hours, dat2$gene_id)
  fun <- switch(
    aggregate,
    median = function(x) stats::median(x, na.rm = TRUE),
    mean = function(x) mean(x, na.rm = TRUE),
    geometric_mean = .geometric_mean
  )

  ref <- data.frame(
    gene_id = names(split_hl),
    half_life_hours = vapply(split_hl, fun, numeric(1L)),
    n_records = vapply(split_hl, length, integer(1L)),
    stringsAsFactors = FALSE
  )

  # Recalculate decay_rate_per_hour from the aggregated half-life so the two
  # columns remain mathematically consistent.
  ref$decay_rate_per_hour <- log(2) / ref$half_life_hours
  ref <- ref[order(ref$gene_id), c("gene_id", "half_life_hours", "decay_rate_per_hour", "n_records")]
  rownames(ref) <- NULL
  .inform(
    verbose,
    "afterglow: built half-life reference for %d genes from %d retained TTDB rows (%d rows before filtering).",
    nrow(ref),
    nrow(dat2),
    n_before_filtering
  )
  ref
}
