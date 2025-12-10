# Data from the US-economy

21 selected quarterly time-series from 1953:Q1 to 2021:Q2. From FRED-QD
data base (McCracken and Ng, 2021). Release date 2021-07. Data is
transformed to be interpreted as growth-rates (first log-differences
with the exception of interest rates, which are already growth rates).

## Usage

``` r
usmacro_growth
```

## Format

A matrix with 247 rows and 21 columns.

## Source

Raw (untransformed) data available at
<https://www.stlouisfed.org/research/economists/mccracken/fred-databases>,
<https://www.stlouisfed.org/-/media/project/frbstl/stlouisfed/research/fred-md/historical-vintages-of-fred-qd-2018-05-to-2024-12.zip>.

## References

McCracken, M. W. and Ng, S. (2021). FRED-QD: A Quarterly Database for
Macroeconomic Research, *Review, Federal Reserve Bank of St. Louis*,
**103**(1), 1–44,
[doi:10.20955/r.103.1-44](https://doi.org/10.20955/r.103.1-44) .
