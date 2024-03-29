% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{summary.bayesianVARs_bvar}
\alias{summary.bayesianVARs_bvar}
\title{Summary method for bayesianVARs_bvar objects}
\usage{
\method{summary}{bayesianVARs_bvar}(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)
}
\arguments{
\item{object}{A \code{bayesianVARs_bvar} object obtained via \code{\link[=bvar]{bvar()}}.}

\item{quantiles}{numeric vector which quantiles to compute.}

\item{...}{Currently ignored!}
}
\value{
An object of type \code{summary.bayesianVARs_bvar}.
}
\description{
Summary method for \code{bayesianVARs_bvar} objects.
}
\examples{
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate model
mod <- bvar(data, quiet = TRUE)

# Summary
sum <- summary(mod)
}
