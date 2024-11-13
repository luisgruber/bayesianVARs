#' Data from the US-economy
#'
#' 21 selected quarterly time-series from 1953:Q1 to 2021:Q2. From FRED-QD data
#' base (McCracken and Ng, 2021). Release date 2021-07. Data is transformed to be interpreted as
#' growth-rates (first log-differences with the exception of interest rates,
#' which are already growth rates).
#'
#' @format A matrix with 247 rows and 21 columns.
#'
#' @source Raw (untransformed) data available at
#'   \url{https://www.stlouisfed.org/research/economists/mccracken/fred-databases},
#'   \url{https://files.stlouisfed.org/files/htdocs/fred-md/quarterly/2021-07.csv}.
#' @references McCracken, M. W. and Ng, S. (2021). FRED-QD: A Quarterly
#'   Database for Macroeconomic Research, \emph{Review, Federal Reserve Bank of St.
#'   Louis}, \bold{103}(1), 1--44, \doi{10.20955/r.103.1-44}.
"usmacro_growth"


