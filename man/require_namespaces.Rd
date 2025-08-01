% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swissArmyKnife.R
\name{require_namespaces}
\alias{require_namespaces}
\title{Vectorised version of requireNamespace}
\usage{
require_namespaces(pkgs)
}
\arguments{
\item{pkgs}{Vector of character strings with required package name(s)}
}
\description{
Carbon copy of the helper function described by Raphael Sonabend-Friend in the package
\code{survivalmodels}. Vectorises the \code{requireNamespace} function and returns \code{TRUE} if all
packages, \code{x}, are available and \code{FALSE} otherwise. For internal use only.
https://github.com/RaphaelS1/survivalmodels/blob/main/R/helpers.R
}
\keyword{internal}
