\name{available_ttdb_files}
\alias{available_ttdb_files}
\title{List bundled or external TTDB-like CSV files}
\usage{
available_ttdb_files(path = get_ttdb_dir())
}
\arguments{
\item{path}{Directory to inspect. Defaults to \code{get_ttdb_dir()}, which can point to an external TTDB download.}
}
\value{
A named character vector of existing TTDB CSV file paths.
}
\description{
List known TTDB CSV files present in a directory. Missing files are skipped.
}
