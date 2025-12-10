# Summary statistics for bayesianVARs posterior draws.

Summary statistics for bayesianVARs posterior draws.

## Usage

``` r
# S3 method for class 'bayesianVARs_draws'
summary(object, quantiles = c(0.25, 0.5, 0.75), ...)
```

## Arguments

- object:

  An object of class `bayesianVARs_draws` usually obtained through
  extractors like
  [`coef.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/coef.md)
  and
  [`vcov.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/vcov.bayesianVARs_bvar.md).

- quantiles:

  A vector of quantiles to evaluate.

- ...:

  Currently ignored.

## Value

A list object of class `bayesianVARs_draws_summary` holding

- `mean`: Vector or matrix containing the posterior mean.

- `sd`: Vector or matrix containing the posterior standard deviation .

- `quantiles`: Array containing the posterior quantiles.

## See also

Available extractors:
[`coef.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/coef.md),
[`vcov.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/vcov.bayesianVARs_bvar.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Extract posterior draws of VAR coefficients
bvar_coefs <- coef(mod)

# Compute summary statistics
summary_stats <- summary(bvar_coefs)

# Compute summary statistics of VAR coefficients without using coef()
summary_stats <- summary(mod$PHI)

# Test which list elements of 'mod' are of class 'bayesianVARs_draws'.
names(mod)[sapply(names(mod), function(x) inherits(mod[[x]], "bayesianVARs_draws"))]
#> [1] "PHI"     "U"       "logvar"  "sv_para" "facload" "fac"    
```
