\name{make_half_life_reference}
\alias{make_half_life_reference}
\title{Build a gene-level half-life reference table}
\usage{
make_half_life_reference(
  stability_data,
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
  verbose = TRUE
)
}
\arguments{
\item{stability_data}{A TTDB stability data frame, or a list returned by
\code{\link{read_ttdb_downloads}}.}
\item{species}{Optional species filter, for example \code{"Human"}.}
\item{cell_type}{Optional exact cell type filter.}
\item{condition}{Optional exact condition filter.}
\item{technique}{Optional exact technique filter.}
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
The decay rate is recalculated as \code{log(2) / half_life_hours} to keep the
two quantities mathematically consistent.
}
