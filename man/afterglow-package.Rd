\name{afterglow-package}
\alias{afterglow-package}
\alias{afterglow}
\title{mRNA half-life-aware differential expression analysis}
\description{
The \pkg{afterglow} package provides a workflow for differential expression
analysis that accounts for gene-specific mRNA half-life. It reads TTDB
downloaded tables or user-provided half-life tables, infers stimulus-time
expression under a pulse model or synthesis-rate model, and runs differential
testing on the inferred matrix.
}
\details{
The default pulse model assumes that the control group is steady-state and that
the treatment group receives an instantaneous net pulse at stimulus time. The
main user-facing function is \code{\link{run_afterglow_de}}.
}
