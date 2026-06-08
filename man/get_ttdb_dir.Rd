\name{get_ttdb_dir}
\alias{get_ttdb_dir}
\title{Return the active TTDB directory}
\usage{
get_ttdb_dir()
}
\value{
Character scalar giving the active TTDB directory. Resolution order is \code{options(afterglow.ttdb_dir)}, \code{AFTERGLOW_TTDB_DIR}, then package \file{extdata/ttdb}.
}
\description{
Return the directory afterglow will use when \code{read_ttdb_downloads()} is called without an explicit path.
}
