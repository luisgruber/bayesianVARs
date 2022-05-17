#' Data from the US-economy
#'
#' 21 selected quarterly time-series from 1953:Q1 to 2021:Q2. From FRED-QD data
#' base. Release date 2021-07.
#'
#' @format A 'xts' object with 250 rows and 21 variables.
#'
#' @source \url{https://research.stlouisfed.org/econ/mccracken/fred-databases/}
"data_raw"

#' Data from the US-economy
#'
#' 21 selected quarterly time-series from 1953:Q1 to 2021:Q2. From FRED-QD data
#' base. Release date 2021-07. Data is transformed to be interpreted as
#' growth-rates (first log-differences with the exception of interest rates,
#' which are already growth rates).
#'
#' @format A 'xts' object with 247 rows and 21 variables.
#'
#' @source \url{https://research.stlouisfed.org/econ/mccracken/fred-databases/}
"dat_growth"

#' Data from the US-economy
#'
#' 21 selected quarterly time-series from 1953:Q1 to 2021:Q2. From FRED-QD data
#' base. Release date 2021-07. Data is transformed to be approximately
#' stationary (using the transformations proposed in
#' \url{https://doi.org/10.20955/r.103.1-44}.
#'
#' @format A 'xts' object with 247 rows and 21 variables.
#'
#' @source \url{https://research.stlouisfed.org/econ/mccracken/fred-databases/}
"dat"
