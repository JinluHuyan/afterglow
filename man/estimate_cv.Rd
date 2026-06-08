\name{estimate_cv}
\alias{estimate_cv}
\title{Estimate expression measurement coefficient of variation}
\usage{
estimate_cv(expr, group)
}
\arguments{
\item{expr}{Numeric expression matrix with genes as rows and samples as columns, on a linear scale.}
\item{group}{Data frame with \code{sample} and \code{group} columns, or a named group vector.}
}
\value{
Named numeric vector of median CV per group plus a pooled estimate.
}
\description{
Estimates replicate-level measurement variability for use in afterglow confidence scoring.
}
