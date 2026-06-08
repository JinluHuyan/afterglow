\name{half_life_sensitivity}
\alias{half_life_sensitivity}
\title{Half-life perturbation sensitivity analysis}
\usage{
half_life_sensitivity(
  expr,
  group,
  half_life_ref,
  time_hours,
  control_group,
  treatment_group,
  expression_scale = c("auto", "log2", "linear"),
  model = c("pulse", "synthesis_rate"),
  perturbation_sd = c(0, 0.05, 0.1, 0.2),
  n_iter = 20L,
  seed = NULL,
  truth_log2_fc = NULL,
  pseudocount = 1e-08,
  negative_action = c("floor", "NA", "keep"),
  verbose = TRUE
)
}
\arguments{
\item{expr}{Numeric matrix with genes as rows and samples as columns.}
\item{group}{Data frame with \code{sample} and \code{group} columns, or a named
group vector.}
\item{half_life_ref}{Data frame with \code{gene_id} and
\code{half_life_hours}.}
\item{time_hours}{Time from stimulus at t0 to sampling at t1, in hours.}
\item{control_group}{Control group label.}
\item{treatment_group}{Stimulated or treatment group label.}
\item{expression_scale}{\code{"linear"}, \code{"log2"}, or \code{"auto"}.}
\item{model}{\code{"pulse"} or \code{"synthesis_rate"}.}
\item{perturbation_sd}{Numeric vector of log-normal standard deviations used to
perturb half-life values. Values 0.05, 0.10, 0.20, and 0.50 approximately
represent 5, 10, 20, and 50 percent multiplicative uncertainty for small
perturbations.}
\item{n_iter}{Number of random perturbation replicates per non-zero
perturbation level.}
\item{seed}{Optional random seed for reproducibility.}
\item{truth_log2_fc}{Optional named numeric vector of true log2 fold-changes,
used in simulation studies to report recovery metrics.}
\item{pseudocount}{Small positive value used for fold-change calculations.}
\item{negative_action}{How to handle inferred non-positive expression values:
\code{"floor"}, \code{"NA"}, or \code{"keep"}.}
\item{verbose}{Logical. If \code{TRUE}, print progress messages.}
}
\value{
An \code{afterglow_sensitivity} list with \code{summary}, \code{gene_level},
\code{nominal}, and \code{settings} components.
}
\description{
Re-runs the selected afterglow kinetic inversion after perturbing each gene's
half-life estimate by reproducible log-normal noise. The summary table reports
fold-change stability against the nominal afterglow result and, when
\code{truth_log2_fc} is supplied, against known truth.
}
\examples{
expr <- matrix(
  c(10, 10, 12, 12, 20, 20, 21, 21),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("GENE1", "GENE2"), c("C1", "C2", "T1", "T2"))
)
group <- data.frame(sample = colnames(expr), group = c("C", "C", "T", "T"))
half_life_ref <- data.frame(gene_id = rownames(expr), half_life_hours = c(1, 6))
sens <- half_life_sensitivity(
  expr = expr,
  group = group,
  half_life_ref = half_life_ref,
  time_hours = 1,
  control_group = "C",
  treatment_group = "T",
  expression_scale = "linear",
  model = "pulse",
  perturbation_sd = c(0, 0.1),
  n_iter = 2,
  seed = 1,
  verbose = FALSE
)
sens$summary
}
