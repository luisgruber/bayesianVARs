% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bvar_wrapper.R
\name{[.bayesianVARs_draws}
\alias{[.bayesianVARs_draws}
\title{Extract or Replace Parts of a bayesianVARs_draws object}
\usage{
\method{[}{bayesianVARs_draws}(x, i, j, ...)
}
\arguments{
\item{x}{An object of type \code{bayesianVARs_draws}.}

\item{i}{indices}

\item{j}{indices}

\item{...}{further indices}
}
\value{
An object of type \code{bayesianVARs_draws}.
}
\description{
Extract or replace parts of a \code{bayesianVARs_draws} object.
}
\examples{
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Extract coefficients, which are of class bayesianVARs_draws
phi <- coef(mod)
phi[1,1,1]
}
