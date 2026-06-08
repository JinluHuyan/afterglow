# Input/output helpers for user-provided expression, group, and TTDB files.
# File parsing is intentionally separated from statistical modeling so users can
# inspect imported data before running the full analysis.

#' Read expression matrix from CSV
#'
#' The expected format is genes in rows and samples in columns. The first column
#' can be an unnamed gene column, `gene`, `Gene`, `gene_id`, or similar.
#'
#' @param file Path to expression CSV.
#' @param gene_col Optional gene column name or column index. If `NULL`, the
#'   first column is used when it is non-numeric or unnamed.
#' @param check_names Passed to `utils::read.csv`; default keeps sample names
#'   unchanged.
#' @param verbose Logical. If `TRUE`, print progress and import diagnostics.
#' @return Numeric matrix with gene IDs as row names and samples as columns.
#' @export
read_expression_csv <- function(file, gene_col = NULL, check_names = FALSE, verbose = TRUE) {
  dat <- utils::read.csv(file, check.names = check_names, stringsAsFactors = FALSE)
  .stop_if(ncol(dat) < 2L, "Expression CSV must contain one gene column and at least one sample column.")

  # Some expression CSV files have a blank or whitespace-only top-left cell.
  # read.csv() converts such a header to names like "" or "X". We explicitly
  # treat the first column as the gene ID column in that situation instead of
  # attempting to coerce it into numeric sample values.
  colnames_trimmed <- trimws(colnames(dat))
  gene_col_index <- NA_integer_

  if (is.null(gene_col)) {
    possible <- c("", "X", "gene", "Gene", "gene_id", "GeneID", "symbol", "Symbol")
    gene_col <- .first_existing_col(dat, possible, required = FALSE, what = "gene column")
    if (is.na(gene_col)) {
      gene_name_candidates <- tolower(colnames_trimmed) %in% c("", "x", "gene", "genes", "geneid", "gene_id", "symbol", "id")
      if (any(gene_name_candidates)) {
        gene_col_index <- which(gene_name_candidates)[[1L]]
      } else {
        first_column_numeric <- !any(is.na(.safe_numeric(dat[[1L]])))
        if (isFALSE(first_column_numeric)) {
          gene_col_index <- 1L
        } else {
          stop(
            "Could not automatically identify the gene ID column. Please set gene_col explicitly.",
            call. = FALSE
          )
        }
      }
    }

    if (is.na(gene_col) && is.na(gene_col_index)) {
      gene_col_index <- 1L
    }
  }

  if (is.numeric(gene_col) && !is.null(gene_col)) {
    gene_col_index <- as.integer(gene_col[[1L]])
  } else if (!is.na(gene_col_index)) {
    # Keep the index identified above.
    gene_col_index <- as.integer(gene_col_index)
  } else {
    exact_match <- which(colnames(dat) == gene_col)
    trimmed_match <- which(colnames_trimmed == trimws(gene_col))
    if (length(exact_match) > 0L) {
      gene_col_index <- exact_match[[1L]]
    } else if (length(trimmed_match) > 0L) {
      gene_col_index <- trimmed_match[[1L]]
    } else {
      stop(sprintf("gene_col '%s' not found in expression CSV.", gene_col), call. = FALSE)
    }
  }

  .stop_if(is.na(gene_col_index) || gene_col_index < 1L || gene_col_index > ncol(dat),
           "gene_col does not identify a valid column.")
  gene_ids <- dat[[gene_col_index]]
  dat <- dat[, -gene_col_index, drop = FALSE]

  mat <- as.matrix(data.frame(lapply(dat, .safe_numeric), check.names = FALSE))
  rownames(mat) <- trimws(as.character(gene_ids))

  .stop_if(anyNA(rownames(mat)) || any(rownames(mat) == ""), "Expression matrix contains empty gene IDs.")
  .stop_if(anyDuplicated(rownames(mat)) > 0L, "Expression matrix contains duplicated gene IDs.")
  .stop_if(any(is.infinite(mat), na.rm = TRUE), "Expression matrix contains Inf or -Inf.")
  .inform(verbose, "afterglow: imported expression matrix with %d genes and %d samples.", nrow(mat), ncol(mat))
  .inform(verbose, "afterglow: expression import detected %d missing numeric values.", sum(is.na(mat)))
  mat
}

#' Read sample group table from CSV
#'
#' @param file Path to group CSV.
#' @param sample_col Name of sample/accession column.
#' @param group_col Name of group column.
#' @param verbose Logical. If `TRUE`, print progress and import diagnostics.
#' @return Data frame with `sample` and `group` columns.
#' @export
read_group_csv <- function(file, sample_col = "Accession", group_col = "Grouptype", verbose = TRUE) {
  dat <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  .stop_if(!sample_col %in% colnames(dat), sprintf("sample_col '%s' not found.", sample_col))
  .stop_if(!group_col %in% colnames(dat), sprintf("group_col '%s' not found.", group_col))

  out <- data.frame(
    sample = as.character(dat[[sample_col]]),
    group = as.character(dat[[group_col]]),
    stringsAsFactors = FALSE
  )

  .stop_if(anyDuplicated(out$sample) > 0L, "Group table contains duplicated sample IDs.")
  .stop_if(anyNA(out$sample) || any(out$sample == ""), "Group table contains empty sample IDs.")
  .stop_if(anyNA(out$group) || any(out$group == ""), "Group table contains empty group labels.")
  .inform(verbose, "afterglow: imported group table with %d samples and %d groups.", nrow(out), length(unique(out$group)))
  out
}

#' Read TTDB downloaded CSV files
#'
#' This function is intentionally permissive: it reads every known TTDB CSV that
#' exists in `ttdb_dir`, and silently skips files that are not present. This lets
#' users point afterglow to either bundled example tables or a full local TTDB
#' download without changing the downstream workflow.
#'
#' @param ttdb_dir Directory containing TTDB CSV files. Defaults to
#'   `get_ttdb_dir()`, which checks `options(afterglow.ttdb_dir)`, the
#'   `AFTERGLOW_TTDB_DIR` environment variable, and package extdata.
#' @param verbose Logical. If `TRUE`, print progress and import diagnostics.
#' @return Named list of data frames.
#' @export
read_ttdb_downloads <- function(ttdb_dir = get_ttdb_dir(), verbose = TRUE) {
  files <- available_ttdb_files(ttdb_dir)
  .stop_if(length(files) == 0L, sprintf(
    "No TTDB CSV files found in: %s. Download TTDB from https://sysbio.gzzoc.com/ttdb/download.html and then call set_ttdb_dir('path/to/ttdb') or read_ttdb_downloads('path/to/ttdb').",
    ttdb_dir
  ))

  out <- lapply(files, function(path) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  })
  names(out) <- names(files)
  .inform(verbose, "afterglow: imported %d TTDB file(s) from %s.", length(out), ttdb_dir)
  out
}
