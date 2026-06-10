\name{make_half_life_reference}
\alias{make_half_life_reference}
\title{Build a gene-level half-life reference table}
\usage{
make_half_life_reference(
  stability_data,
  species = NULL,
  cell_type = NULL,
  tissue = NULL,
  condition = NULL,
  technique = NULL,
  study_id = NULL,
  sample_id = NULL,
  gene_col_preference = c("gene_name_x", "gene_name_y", "ensembl_id", "gene_id"),
  half_life_col = "half_life",
  decay_rate_col = "decay_rate",
  aggregate = c("median", "mean", "geometric_mean"),
  min_r_squared = NULL,
  use_decay_rate_if_half_life_missing = TRUE,
  collapse_case = c("as_is", "upper"),
  verbose = TRUE
)
}
\arguments{
\item{stability_data}{A TTDB stability data frame, or a list returned by
\code{\link{read_ttdb_downloads}}.}
\item{species}{Optional species filter, for example \code{"Human"}.}
\item{cell_type}{Optional exact cell type filter.}
\item{tissue}{Optional tissue or cell-type filter. TTDB stores this metadata
in the \code{cell_type} column, so \code{tissue} is a user-facing alias for
\code{cell_type}.}
\item{condition}{Optional exact condition filter.}
\item{technique}{Optional exact technique filter.}
\item{study_id}{Optional exact TTDB study identifier filter.}
\item{sample_id}{Optional exact TTDB sample identifier filter.}
\item{gene_col_preference}{Candidate gene identifier columns, in priority order.}
\item{half_life_col}{Column containing half-life values in hours.}
\item{decay_rate_col}{Column containing decay rates per hour.}
\item{aggregate}{One of \code{"median"}, \code{"mean"}, or
\code{"geometric_mean"}.}
\item{min_r_squared}{Optional lower bound for \code{r_squared}, when available.}
\item{use_decay_rate_if_half_life_missing}{If \code{TRUE}, missing half-life is
computed as \code{log(2) / decay_rate} when decay rate is positive.}
\item{collapse_case}{Either \code{"as_is"} or \code{"upper"} for gene IDs.}
\item{verbose}{Logical. If \code{TRUE}, print progress and filtering diagnostics.}
}
\value{
Data frame with one row per gene and columns \code{gene_id},
\code{half_life_hours}, \code{decay_rate_per_hour}, and \code{n_records}.
}
\description{
Convert TTDB stability records into a one-row-per-gene half-life reference table.
}
\details{
TTDB can contain multiple records per gene across studies, samples, cell types,
and conditions. This function filters the rows and aggregates half-life values.
When a matching cell type or tissue-like entry is available, prefer filtering by
\code{cell_type} or its alias \code{tissue} before falling back to a broader
species-level reference.
The decay rate is recalculated as \code{log(2) / half_life_hours} to keep the
two quantities mathematically consistent.
}
