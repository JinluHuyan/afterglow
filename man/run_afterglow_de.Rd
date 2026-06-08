\name{run_afterglow_de}
\alias{run_afterglow_de}
\title{Full half-life-aware differential expression workflow}
\usage{
run_afterglow_de(
  expr,
  group,
  half_life_ref,
  time_hours,
  control_group,
  treatment_group,
  expression_scale = c("auto", "log2", "linear"),
  model = c("pulse", "synthesis_rate"),
  de_method = c("limma", "limma_voom", "edgeR", "deseq2", "welch"),
  assay_type = c("auto", "microarray", "rna_seq"),
  pseudocount = 1e-08,
  negative_action = c("floor", "NA", "keep"),
  deseq_round_counts = TRUE,
  write_corrected_matrix = FALSE,
  corrected_matrix_file = "afterglow_corrected_expression.csv",
  write_results = FALSE,
  results_file = "afterglow_results.csv",
  verbose = TRUE
)
}
\arguments{
\item{expr}{Numeric expression matrix with genes as rows and samples as
columns. Row names must match \code{half_life_ref$gene_id}; column names must
match \code{group$sample}.}
\item{group}{Data frame with \code{sample} and \code{group} columns, or a named
group vector.}
\item{half_life_ref}{Gene-level half-life reference table with \code{gene_id}
and \code{half_life_hours}.}
\item{time_hours}{Time from stimulation at t0 to sampling at t1, in hours.}
\item{control_group}{Group label for untreated or baseline samples.}
\item{treatment_group}{Group label for stimulated or treated samples.}
\item{expression_scale}{\code{"log2"} for log2 expression,
\code{"linear"} for counts or linear expression, or \code{"auto"} for
heuristic detection. Use \code{"linear"} explicitly for raw counts.}
\item{model}{\code{"pulse"} for an instantaneous or transient net response, or
\code{"synthesis_rate"} for a sustained altered RNA synthesis rate.}
\item{de_method}{Differential expression method: \code{"limma"},
\code{"limma_voom"}, \code{"edgeR"}, \code{"deseq2"}, or \code{"welch"}.}
\item{assay_type}{Assay description: \code{"microarray"}, \code{"rna_seq"}, or
\code{"auto"}.}
\item{pseudocount}{Small positive value used for log transformations and
non-positive value handling.}
\item{negative_action}{How to handle non-positive inferred expression values:
\code{"floor"}, \code{"NA"}, or \code{"keep"}.}
\item{deseq_round_counts}{Logical. Round non-integer corrected counts when
\code{de_method = "deseq2"}.}
\item{write_corrected_matrix}{Logical. If \code{TRUE}, write the corrected
expression matrix to CSV.}
\item{corrected_matrix_file}{CSV path for the corrected expression matrix.}
\item{write_results}{Logical. If \code{TRUE}, write the result table to CSV.}
\item{results_file}{CSV path for the result table.}
\item{verbose}{Logical. If \code{TRUE}, print progress and result summaries.}
}
\value{
An \code{afterglow_de} list with corrected matrices, diagnostics, settings, and
a \code{results} data frame.
}
\description{
Run the complete workflow: validate inputs, infer half-life-corrected
stimulus-time expression, perform differential expression testing, and
optionally write output CSV files.
}
\examples{
expr <- matrix(
  c(10, 10.1, 11, 11.2),
  nrow = 1,
  dimnames = list("GENE1", c("C1", "C2", "T1", "T2"))
)
group <- data.frame(
  sample = c("C1", "C2", "T1", "T2"),
  group = c("C", "C", "T", "T")
)
half_life_ref <- data.frame(gene_id = "GENE1", half_life_hours = 1)
fit <- run_afterglow_de(
  expr = expr,
  group = group,
  half_life_ref = half_life_ref,
  time_hours = 1,
  control_group = "C",
  treatment_group = "T",
  expression_scale = "linear",
  model = "pulse",
  de_method = "welch",
  verbose = FALSE
)
head(fit$results)
}
