\name{infer_turnover_expression}
\alias{infer_turnover_expression}
\title{Infer stimulus-time expression using mRNA half-life}
\usage{
infer_turnover_expression(
  expr,
  group,
  half_life_ref,
  time_hours,
  control_group,
  treatment_group,
  expression_scale = c("auto", "log2", "linear"),
  model = c("pulse", "synthesis_rate"),
  pseudocount = 1e-08,
  negative_action = c("floor", "NA", "keep"),
  verbose = TRUE
)
}
\arguments{
\item{expr}{Numeric matrix, genes by samples.}
\item{group}{Data frame with \code{sample} and \code{group}, or a named group
vector.}
\item{half_life_ref}{Data frame with \code{gene_id} and
\code{half_life_hours}.}
\item{time_hours}{Time from stimulus at t0 to sampling at t1, in hours.}
\item{control_group}{Control group label.}
\item{treatment_group}{Stimulated/treatment group label.}
\item{expression_scale}{\code{"linear"}, \code{"log2"}, or \code{"auto"}.}
\item{model}{\code{"pulse"} or \code{"synthesis_rate"}.}
\item{pseudocount}{Pseudocount used before log2 transformation and denominator
protection.}
\item{negative_action}{How to handle inferred non-positive expression:
\code{"floor"}, \code{"NA"}, or \code{"keep"}.}
\item{verbose}{Logical. If \code{TRUE}, print progress and model diagnostics.}
}
\value{
List containing corrected matrix, linear-scale intermediate values, matched
half-life table, and diagnostics.
}
\description{
Infer the expression matrix at stimulus time using gene-specific mRNA half-life.
}
\details{
For the pulse model, the function uses
\deqn{E^{inferred}_{g,s}(t0)=B_g+(E_{g,s}(t1)-B_g)e^{\lambda_g \Delta t}}
where \eqn{B_g} is the control baseline and
\eqn{\lambda_g=\log(2)/t_{1/2,g}}.

For the synthesis-rate model, the function uses
\deqn{M^{equiv}_g=\frac{M_g(t1)-B_g e^{-\lambda_g\Delta t}}{1-e^{-\lambda_g\Delta t}}}
under the assumption that stimulation changes synthesis rate but not decay rate.
}
