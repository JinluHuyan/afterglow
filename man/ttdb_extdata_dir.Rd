\name{ttdb_extdata_dir}
\alias{ttdb_extdata_dir}
\title{Return the installed TTDB extdata directory}
\usage{
ttdb_extdata_dir()
}
\value{
Character scalar giving the installed \file{extdata/ttdb} directory.
}
\description{
Return the package extdata TTDB directory. In public release builds this directory usually contains only a README placeholder because TTDB CSV files are large and are expected to be downloaded separately. For analysis, prefer \code{set_ttdb_dir()}, \code{AFTERGLOW_TTDB_DIR}, or an explicit path passed to \code{read_ttdb_downloads()}.
}
