# Extract VAR coefficients

Extracts posterior draws of the VAR coefficients from a VAR model
estimated with
[`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

## Usage

``` r
# S3 method for class 'bayesianVARs_bvar'
coef(object, ...)
```

## Arguments

- object:

  A `bayesianVARs_bvar` object obtained from
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- ...:

  Currently ignored.

## Value

Returns a numeric array of dimension \\M \times K \times draws\\, where
M is the number of time-series, K is the number of covariates per
equation (including the intercept) and draws is the number of stored
posterior draws.

## See also

[`summary.bayesianVARs_draws()`](https://luisgruber.github.io/bayesianVARs/reference/summary.bayesianVARs_draws.md),
[`vcov.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/vcov.bayesianVARs_bvar.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Extract posterior draws of VAR coefficients
bvar_coefs <- coef(mod)
```
