\name{compute_confidence}
\alias{compute_confidence}
\title{Compute model-aware afterglow confidence scores}
\usage{
compute_confidence(
  half_life_hours,
  time_hours,
  observed_fc = NULL,
  n_records = NULL,
  cv_measurement = 0.05,
  model = c("pulse", "synthesis_rate")
)
}
\arguments{
\item{half_life_hours}{Numeric vector of mRNA half-lives in hours.}
\item{time_hours}{Scalar time between stimulus and RNA measurement, in hours.}
\item{observed_fc}{Optional numeric vector of observed log2 fold-changes used for the signal-to-noise component.}
\item{n_records}{Optional integer vector: number of independent half-life measurements per gene.}
\item{cv_measurement}{Assumed coefficient of variation of expression measurements.}
\item{model}{Kinetic model. \code{"pulse"} uses \eqn{\exp(\lambda t)}; \code{"synthesis_rate"} uses \eqn{1/(1-\exp(-\lambda t))}.}
}
\value{
A data frame with model, amplification_factor, correction_factor, confidence tier, confidence score, SNR, evidence score, and related components.
}
\description{
Computes per-gene confidence scores using model-specific error amplification factors. The pulse model is unstable for short-half-life genes sampled long after a transient pulse; the synthesis-rate model is unstable for long-half-life genes sampled early after a sustained rate change.
}
