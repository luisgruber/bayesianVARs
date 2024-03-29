% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{print.bayesianVARs_bvar}
\alias{print.bayesianVARs_bvar}
\title{Pretty printing of a bvar object}
\usage{
\method{print}{bayesianVARs_bvar}(x, ...)
}
\arguments{
\item{x}{Object of class \code{bayesianVARs_bvar}, usually resulting from a call
of \code{\link[=bvar]{bvar()}}.}

\item{...}{Ignored.}
}
\value{
Returns \code{x} invisibly.
}
\description{
Pretty printing of a bvar object
}
\examples{
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Print model
mod

}
