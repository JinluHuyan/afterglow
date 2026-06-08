\name{read_group_csv}
\alias{read_group_csv}
\title{Read sample group table from CSV}
\usage{
read_group_csv(file, sample_col = "Accession", group_col = "Grouptype", verbose = TRUE)
}
\arguments{
\item{file}{Path to group CSV.}
\item{sample_col}{Name of sample/accession column.}
\item{group_col}{Name of group column.}
\item{verbose}{Logical. If \code{TRUE}, print import diagnostics.}
}
\value{
Data frame with columns \code{sample} and \code{group}.
}
\description{
Read a sample group table and standardize it to the columns expected by
\code{\link{run_afterglow_de}}.
}
