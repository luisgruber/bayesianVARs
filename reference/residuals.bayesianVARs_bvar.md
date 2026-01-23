# Extract Model Residuals

Extract model residuals, defined as the difference between the observed
time-series and the in-sample predictions of the VAR model. Because
in-sample prediction is subject to uncertainty of the VAR parameter
estimates, this uncertainty carries over to the model residuals.

## Usage

``` r
# S3 method for class 'bayesianVARs_bvar'
residuals(object, ...)
```

## Arguments

- object:

  A `bayesianVARs_bvar` object estimated via
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- ...:

  Passed to
  [`fitted.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/fitted.bayesianVARs_bvar.md).

## Value

An object of class `bayesianVARs_residuals`.

## See also

[`fitted.bayesianVARs_bvar`](https://luisgruber.github.io/bayesianVARs/reference/fitted.bayesianVARs_bvar.md)

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

mod.resids <- residuals(mod)
plot(mod.resids)
```
