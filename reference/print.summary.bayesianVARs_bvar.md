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
#> GDPC1.l1     0.222    0.017    0.013
#> CPIAUCSL.l1 -0.053    0.634   -0.004
#> FEDFUNDS.l1  0.012    0.038    1.001
#> intercept    0.006    0.001    0.000
#> 
#> Posterior interquartile range of of reduced-form coefficients:
#>             GDPC1 CPIAUCSL FEDFUNDS
#> GDPC1.l1    0.101    0.037    0.021
#> CPIAUCSL.l1 0.111    0.096    0.014
#> FEDFUNDS.l1 0.019    0.015    0.006
#> intercept   0.001    0.001    0.000
#> 
#> Posterior median of factor loadings:
#>          factor1
#> GDPC1      0.000
#> CPIAUCSL   0.002
#> FEDFUNDS   0.000
#> 
#> Posterior interquartile range of factor loadings:
#>          factor1
#> GDPC1      0.000
#> CPIAUCSL   0.002
#> FEDFUNDS   0.000
```
