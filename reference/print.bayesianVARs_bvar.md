# Pretty printing of a bvar object

Pretty printing of a bvar object

## Usage

``` r
# S3 method for class 'bayesianVARs_bvar'
print(x, ...)
```

## Arguments

- x:

  Object of class `bayesianVARs_bvar`, usually resulting from a call of
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- ...:

  Ignored.

## Value

Returns `x` invisibly.

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Print model
mod
#> 
#> Fitted bayesianVARs_bvar object with
#>   -       3 series
#>   -       1 lag(s)
#>   -     246 used observations
#>   -     247 total observations
#>   -    1000 MCMC draws
#>   -       1 thinning
#>   -    1000 burn-in
#> 
```
