% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGEO.R
\name{getGEO}
\alias{getGEO}
\title{Get GEO metadata}
\usage{
getGEO(searchterm, minnum = 0, maxnum = 1000)
}
\arguments{
\item{searchterm}{input searching terms as GEO database keywords, multiple terms are seperated by blanks}

\item{minnum}{The minimum return records, default is 0}

\item{maxnum}{The maximum return records, default is 1000}
}
\value{
A `data.frame`.
}
\description{
Get GEO metadata
}
\examples{
\donttest{
meta <- getGEO('ezh2 knockout',maxnum = 5)
}
}
