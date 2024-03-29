% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bvar_wrapper.R
\name{stable_bvar}
\alias{stable_bvar}
\title{Stable posterior draws}
\usage{
stable_bvar(object, quiet = FALSE)
}
\arguments{
\item{object}{A \code{bayesianVARs_bvar} object obtained via \code{\link[=bvar]{bvar()}}.}

\item{quiet}{logical indicating whether informative output should be omitted.}
}
\value{
An object of type \code{bayesianVARs_bvar}.
}
\description{
\code{stable_bvar()} detects and discards all posterior draws of an
\code{bayesianVARs_bvar} object that do not fulfill the stability condition:
A VAR(p) model is considered as stable only if the eigenvalues of the
companion form matrix lie inside the unit circle.
}
\examples{
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Discard "unstable" draws
stable_mod <- stable_bvar(mod)

}
