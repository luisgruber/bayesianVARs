# Extract posterior draws of the (time-varying) variance-covariance matrix for a VAR model

Returns the posterior draws of the possibly time-varying
variance-covariance matrix of a VAR estimated via
[`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).
Returns the full paths if `sv_keep="all"` when calling
[`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).
Otherwise, the draws of the variance-covariance matrix for the last
observation are returned, only.

## Usage

``` r
# S3 method for class 'bayesianVARs_bvar'
vcov(object, t = seq_len(nrow(object$logvar)), ...)
```

## Arguments

- object:

  An object of class `bayesianVARs_bvar` obtained via
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- t:

  Vector indicating which points in time should be extracted, defaults
  to all.

- ...:

  Currently ignored.

## Value

An array of class `bayesianVARs_draws` of dimension \\T \times M \times
M \times draws\\, where \\T\\ is the number of observations, \\M\\ the
number of time-series and \\draws\\ the number of stored posterior
draws.

## See also

[`summary.bayesianVARs_draws`](https://luisgruber.github.io/bayesianVARs/reference/summary.bayesianVARs_draws.md),
[`coef.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/coef.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Extract posterior draws of the variance-covariance matrix
bvar_vcov <- vcov(mod)
```
