\name{set_ttdb_dir}
\alias{set_ttdb_dir}
\title{Configure an external TTDB directory}
\usage{
set_ttdb_dir(path = NULL)
}
\arguments{
\item{path}{Directory containing TTDB CSV files. Use \code{NULL} to clear the session option.}
}
\value{
The previous option value, invisibly.
}
\description{
Configure the directory where afterglow should look for local TTDB CSV files. For release builds, TTDB is expected to be downloaded separately and kept outside the package source. This helper sets \code{options(afterglow.ttdb_dir = ...)} for the current R session.
}
\examples{
\dontrun{
set_ttdb_dir("~/data/ttdb")
ttdb <- read_ttdb_downloads()
}
}
