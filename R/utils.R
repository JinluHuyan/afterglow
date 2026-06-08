# Utility functions used internally by afterglow.
# The package intentionally avoids tidyverse dependencies so that installation
# remains lightweight and reproducible in minimal R environments.

#' Configure an external TTDB directory
#'
#' For release builds, TTDB CSV files are expected to live outside the package.
#' Use this helper once per R session, or set the environment variable
#' `AFTERGLOW_TTDB_DIR`, to point afterglow to a local TTDB download.
#'
#' @param path Directory containing TTDB CSV files. Use `NULL` to clear the
#'   session option.
#' @return The previous option value, invisibly.
#' @export
set_ttdb_dir <- function(path = NULL) {
  old <- getOption("afterglow.ttdb_dir", NULL)
  if (is.null(path)) {
    options(afterglow.ttdb_dir = NULL)
  } else {
    .stop_if(!dir.exists(path), sprintf("TTDB directory does not exist: %s", path))
    options(afterglow.ttdb_dir = normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  invisible(old)
}

#' Return the active TTDB directory
#'
#' Directory resolution order is: explicit `ttdb_dir` argument in
#' `read_ttdb_downloads()`, `options(afterglow.ttdb_dir = ...)`, environment
#' variable `AFTERGLOW_TTDB_DIR`, and finally package `extdata/ttdb` if present.
#'
#' @return Character scalar. The active TTDB directory path, or the package
#'   extdata placeholder path when no external directory is configured.
#' @export
get_ttdb_dir <- function() {
  opt <- getOption("afterglow.ttdb_dir", NULL)
  if (!is.null(opt) && nzchar(opt)) return(normalizePath(opt, winslash = "/", mustWork = FALSE))
  env <- Sys.getenv("AFTERGLOW_TTDB_DIR", unset = "")
  if (nzchar(env)) return(normalizePath(env, winslash = "/", mustWork = FALSE))
  ttdb_extdata_dir()
}

#' Return the installed TTDB extdata directory
#'
#' @return Character scalar. The directory where the package expects TTDB CSV
#'   files when they are shipped inside `inst/extdata/ttdb`.
#' @export
ttdb_extdata_dir <- function() {
  # Locate package-installed files under inst/extdata/ttdb.
  system.file("extdata", "ttdb", package = "afterglow", mustWork = FALSE)
}

#' List bundled TTDB-like CSV files
#'
#' @param path Directory to inspect. Defaults to the package TTDB extdata folder.
#' @return A named character vector of existing CSV file paths.
#' @export
available_ttdb_files <- function(path = get_ttdb_dir()) {
  expected <- c(
    computed_feature = "computed_feature.csv",
    species_stability = "species_stability_0927.csv",
    species_stability_combined = "species_stability_combined_0927_2_no_threshold.csv",
    species_top_100_gene_correlations = "species_top_100_gene_correlations.csv",
    study_info = "stuty_info_0913.csv"
  )

  files <- file.path(path, expected)
  files <- files[file.exists(files)]
  names(files) <- names(expected)[file.exists(file.path(path, expected))]
  files
}

.stop_if <- function(condition, message) {
  # Small wrapper to keep validation code readable and consistent.
  if (isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

.match_arg <- function(x, choices, name = deparse(substitute(x))) {
  # Lightweight alternative to match.arg() with explicit error text.
  if (length(x) != 1L || !x %in% choices) {
    stop(sprintf("%s must be one of: %s", name, paste(choices, collapse = ", ")), call. = FALSE)
  }
  x
}

.first_existing_col <- function(data, candidates, required = TRUE, what = "column") {
  # Return the first candidate column that exists in the data frame.
  found <- candidates[candidates %in% colnames(data)]
  if (length(found) > 0L) {
    return(found[[1L]])
  }
  if (isTRUE(required)) {
    stop(
      sprintf("Cannot find %s. Tried: %s", what, paste(candidates, collapse = ", ")),
      call. = FALSE
    )
  }
  NA_character_
}

.safe_numeric <- function(x) {
  # Convert character/factor-like vectors to numeric. Non-numeric values become
  # NA and are handled by downstream validation.
  suppressWarnings(as.numeric(as.character(x)))
}

.collapse_gene_id <- function(x, collapse_case = c("as_is", "upper")) {
  # Gene symbol case can be biologically meaningful; keep it unchanged by
  # default and only uppercase when explicitly requested.
  collapse_case <- .match_arg(collapse_case[[1L]], c("as_is", "upper"), "collapse_case")
  x <- trimws(as.character(x))
  if (collapse_case == "upper") {
    x <- toupper(x)
  }
  x
}

.geometric_mean <- function(x, na.rm = TRUE) {
  # Half-life estimates can be right-skewed across studies. Geometric mean is
  # provided as an optional aggregation strategy, while median remains default.
  x <- x[x > 0]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  exp(mean(log(x), na.rm = na.rm))
}

.row_var <- function(mat) {
  # Base-R row variance helper. For very large matrices, this can be replaced
  # by matrixStats::rowVars() in a future performance-focused release.
  n <- ncol(mat)
  if (n < 2L) {
    return(rep(NA_real_, nrow(mat)))
  }
  row_means <- rowMeans(mat, na.rm = TRUE)
  rowSums((mat - row_means)^2, na.rm = TRUE) / (n - 1L)
}

.bh <- function(p) {
  # Centralized Benjamini-Hochberg/FDR adjustment helper.
  stats::p.adjust(p, method = "BH")
}

.inform <- function(verbose, ...) {
  # User-facing progress messages are controlled by a single verbose flag.
  if (isTRUE(verbose)) {
    message(sprintf(...))
  }
}
