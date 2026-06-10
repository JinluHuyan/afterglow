\name{ttdb_filter_options}
\alias{ttdb_filter_options}
\title{List TTDB metadata values available for filtering}
\usage{
ttdb_filter_options(
  stability_data,
  species = NULL,
  cell_type = NULL,
  condition = NULL,
  technique = NULL,
  study_id = NULL,
  sample_id = NULL,
  fields = c("species_name", "cell_type", "condition", "technique", "study_id", "sample_id"),
  gene_col_preference = c("gene_name_x", "gene_name_y", "ensembl_id", "gene_id")
)
}
\arguments{
\item{stability_data}{A TTDB stability data frame, or a list returned by
\code{\link{read_ttdb_downloads}}.}
\item{species, cell_type, condition, technique, study_id, sample_id}{Optional exact filters applied before listing options.}
\item{fields}{Metadata columns to summarize.}
\item{gene_col_preference}{Candidate gene identifier columns, in priority order for counting distinct genes per metadata value.}
}
\value{
Data frame with columns \code{field}, \code{value}, \code{n_records} and
\code{n_genes}.
}
\description{
List the species, cell-type, condition, technique, study and sample values
available in a TTDB stability table after optional exact filtering.
}
\details{
Use this helper before \code{\link{make_half_life_reference}} to identify
cell-type or study-specific TTDB subsets. TTDB stores cell type and sample
source labels in the \code{cell_type} column, and the available labels can also
be browsed on the TTDB website.
}
\examples{
ttdb <- data.frame(
  gene_name_x = c("A", "B"),
  species_name = c("Human", "Human"),
  cell_type = c("HepG2", "K562"),
  condition = c("WT", "WT"),
  half_life = c(1, 2)
)
ttdb_filter_options(ttdb, species = "Human", fields = "cell_type")
}
