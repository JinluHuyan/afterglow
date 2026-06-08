\name{differential_expression}
\alias{differential_expression}
\title{Differential expression on an expression matrix}
\usage{
differential_expression(
  mat,
  group,
  control_group,
  treatment_group,
  method = c("limma", "limma_voom", "edgeR", "deseq2", "welch"),
  assay_type = c("auto", "microarray", "rna_seq"),
  log_scale = TRUE,
  pseudocount = 1e-08,
  deseq_round_counts = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{mat}{Numeric matrix with genes as rows and samples as columns. For
\code{method = "limma"}, this is usually log2-scale microarray or normalized
expression. For \code{"limma_voom"}, \code{"edgeR"}, or \code{"deseq2"}, this
must represent non-negative sequencing counts or count-like values.}
\item{group}{Data frame with \code{sample} and \code{group} columns, or a named
group vector.}
\item{control_group}{Control group label.}
\item{treatment_group}{Treatment group label.}
\item{method}{Differential expression method. Supported values are
\code{"limma"}, \code{"limma_voom"}, \code{"edgeR"}, \code{"deseq2"}, and
\code{"welch"}.}
\item{assay_type}{Optional assay description used for user guidance:
\code{"auto"}, \code{"microarray"}, or \code{"rna_seq"}.}
\item{log_scale}{Whether \code{mat} is already on log scale. Count-scale
methods require \code{log_scale = FALSE}.}
\item{pseudocount}{Small positive number used for linear-scale log fold change
and count safeguards.}
\item{deseq_round_counts}{Logical. DESeq2 requires integer counts. If
\code{TRUE}, non-integer turnover-corrected count-like values are rounded after
non-negative validation.}
\item{verbose}{Logical. If \code{TRUE}, print progress and method diagnostics.}
}
\value{
Data frame ordered by adjusted p-value with columns including
\code{gene_id}, \code{log2_fc}, \code{average_expression}, \code{statistic},
\code{p_value}, \code{adjusted_p_value}, and \code{method}.
}
\description{
Run two-group differential expression testing on an expression matrix.
}
\details{
\code{method = "limma"} is recommended for microarray/chip data or log2
normalized expression matrices. \code{method = "limma_voom"} uses
\pkg{edgeR} normalization followed by \pkg{limma} voom and is appropriate for
RNA-seq count matrices. \code{method = "edgeR"} uses edgeR quasi-likelihood
testing. \code{method = "deseq2"} uses DESeq2 and requires integer counts;
rounded turnover-corrected counts should be reported as an approximation.
\code{method = "welch"} is a lightweight base-R fallback for quick checks.
}
