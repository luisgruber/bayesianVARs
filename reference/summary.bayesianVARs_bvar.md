# Summary method for bayesianVARs_bvar objects

Summary method for `bayesianVARs_bvar` objects.

## Usage

``` r
# S3 method for class 'bayesianVARs_bvar'
summary(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), digits = 3L, ...)
```

## Arguments

- object:

  A `bayesianVARs_bvar` object obtained via
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- quantiles:

  numeric vector which quantiles to compute.

- digits:

  Single integer indicating the number of decimal places to be used for
  rounding the summary statistics. Negative values are not allowed.

- ...:

  Currently ignored!

## Value

An object of type `summary.bayesianVARs_bvar`.

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate model
mod <- bvar(data, quiet = TRUE)

# Summary
sum <- summary(mod)
```
