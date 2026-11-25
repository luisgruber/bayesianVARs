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
#> GDPC1.l1     0.230    0.005    0.021
#> CPIAUCSL.l1 -0.045    0.611   -0.004
#> FEDFUNDS.l1  0.008    0.037    1.001
#> intercept    0.006    0.001    0.000
#> 
#> Posterior interquartile range of of reduced-form coefficients:
#>             GDPC1 CPIAUCSL FEDFUNDS
#> GDPC1.l1    0.092    0.026    0.022
#> CPIAUCSL.l1 0.110    0.093    0.015
#> FEDFUNDS.l1 0.016    0.014    0.006
#> intercept   0.001    0.001    0.000
#> 
#> Posterior median of factor loadings:
#>          factor1
#> GDPC1      0.000
#> CPIAUCSL   0.001
#> FEDFUNDS   0.000
#> 
#> Posterior interquartile range of factor loadings:
#>          factor1
#> GDPC1      0.000
#> CPIAUCSL   0.002
#> FEDFUNDS   0.000
```
