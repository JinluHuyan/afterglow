\name{read_ttdb_downloads}
\alias{read_ttdb_downloads}
\title{Read TTDB downloaded CSV files}
\usage{
read_ttdb_downloads(ttdb_dir = get_ttdb_dir(), verbose = TRUE)
}
\arguments{
\item{ttdb_dir}{Directory containing TTDB CSV files. Defaults to \code{get_ttdb_dir()}, which checks \code{options(afterglow.ttdb_dir)}, \code{AFTERGLOW_TTDB_DIR}, then package extdata.}
\item{verbose}{Logical. If \code{TRUE}, print import diagnostics.}
}
\value{
Named list of data frames.
}
\description{
Read every known TTDB CSV file that exists in a directory. Public afterglow releases are expected to keep TTDB outside the package source because the full CSV download is large. Download TTDB from \url{https://sysbio.gzzoc.com/ttdb/download.html}, then call \code{set_ttdb_dir("path/to/ttdb")} or pass the directory explicitly.
}
