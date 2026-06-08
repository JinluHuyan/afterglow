\name{read_expression_csv}
\alias{read_expression_csv}
\title{Read expression matrix from CSV}
\usage{
read_expression_csv(file, gene_col = NULL, check_names = FALSE, verbose = TRUE)
}
\arguments{
\item{file}{Path to expression CSV.}
\item{gene_col}{Optional gene column name or column index. If \code{NULL}, the
first likely gene column is used.}
\item{check_names}{Passed to \code{utils::read.csv}; default keeps sample names
unchanged.}
\item{verbose}{Logical. If \code{TRUE}, print import diagnostics.}
}
\value{
Numeric matrix with gene IDs as row names and samples as columns.
}
\description{
Read a gene-by-sample expression matrix from CSV. The first column can be an
unnamed gene column or a named gene identifier column.
}
\details{
If the top-left CSV header cell is blank or contains only whitespace, the first
column is treated as the gene ID column. This prevents gene IDs from being
coerced into numeric expression values.
}
