\name{add_confidence}
\alias{add_confidence}
\title{Add model-aware confidence scores to an afterglow result}
\usage{
add_confidence(ag_result, observed_fc = NULL, cv_measurement = 0.05)
}
\arguments{
\item{ag_result}{A list returned by \code{infer_turnover_expression()} or \code{run_afterglow_de()}.}
\item{observed_fc}{Optional log2 fold-change vector used for the signal-to-noise part of confidence scoring. If \code{NULL}, raw observed treatment-vs-control log2 fold-change is used.}
\item{cv_measurement}{Assumed coefficient of variation of expression measurements.}
}
\value{
The input afterglow result with an added \code{confidence} data frame. The table includes model, amplification_factor, raw_observed_log2_fc, corrected_log2_fc, and confidence_input_log2_fc so users can distinguish the measured signal from the corrected signal.
}
\description{
Convenience wrapper around \code{compute_confidence()} that extracts the model, time, half-lives, and raw/corrected fold-changes from an afterglow result object.
}
