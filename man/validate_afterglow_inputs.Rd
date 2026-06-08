\name{validate_afterglow_inputs}
\alias{validate_afterglow_inputs}
\title{Validate expression, group, and half-life inputs}
\usage{
validate_afterglow_inputs(
  expr,
  group,
  half_life_ref,
  control_group,
  treatment_group
)
}
\arguments{
\item{expr}{Numeric matrix, genes by samples. Rows must be genes and columns
must be samples. Row names are matched to \code{half_life_ref$gene_id}; column
names are matched to \code{group$sample}.}
\item{group}{Data frame with \code{sample} and \code{group} columns, or a named
vector whose names are sample IDs and whose values are group labels.}
\item{half_life_ref}{Data frame with \code{gene_id} and
\code{half_life_hours}. Each gene must appear once and half-life values must be
positive.}
\item{control_group}{The label in \code{group$group} identifying untreated or
baseline samples.}
\item{treatment_group}{The label in \code{group$group} identifying stimulated or
treated samples.}
}
\value{
A normalized list used internally by \code{\link{run_afterglow_de}}.
}
\description{
Check required structure, sample matching, group labels, half-life positivity,
and gene overlap before afterglow modeling.
}
