# Print method for summary.bayesianVARs_bvar objects

Print method for `summary.bayesianVARs_bvar` objects.

## Usage

``` r
# S3 method for class 'summary.bayesianVARs_bvar'
print(x, ...)
```

## Arguments

- x:

  A `summary.bayesianVARs_bvar` object obtained via
  [`summary.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/summary.bayesianVARs_bvar.md).

- ...:

  Currently ignored!

## Value

Returns `x` invisibly!

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate model
mod <- bvar(data, quiet = TRUE)

# Print summary
summary(mod)
#> 
#> Posterior median of reduced-form coefficients:
#>              GDPC1 CPIAUCSL FEDFUNDS
#> GDPC1.l1     0.230    0.003    0.043
#> CPIAUCSL.l1 -0.072    0.636   -0.003
#> FEDFUNDS.l1  0.013    0.037    0.998
#> intercept    0.006    0.001    0.000
#> 
#> Posterior interquartile range of of reduced-form coefficients:
#>             GDPC1 CPIAUCSL FEDFUNDS
#> GDPC1.l1    0.076    0.024    0.023
#> CPIAUCSL.l1 0.120    0.086    0.016
#> FEDFUNDS.l1 0.017    0.015    0.007
#> intercept   0.001    0.001    0.000
#> 
#> Posterior median of factor loadings:
#>          factor1
#> GDPC1      0.002
#> CPIAUCSL   0.000
#> FEDFUNDS   0.000
#> 
#> Posterior interquartile range of factor loadings:
#>          factor1
#> GDPC1      0.001
#> CPIAUCSL   0.000
#> FEDFUNDS   0.000
```
