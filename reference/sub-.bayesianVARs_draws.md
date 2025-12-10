# Extract or Replace Parts of a bayesianVARs_draws object

Extract or replace parts of a `bayesianVARs_draws` object.

## Usage

``` r
# S3 method for class 'bayesianVARs_draws'
x[i, j, ...]
```

## Arguments

- x:

  An object of type `bayesianVARs_draws`.

- i:

  indices

- j:

  indices

- ...:

  further indices

## Value

An object of type `bayesianVARs_draws`.

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Extract coefficients, which are of class bayesianVARs_draws
phi <- coef(mod)
phi[1,1,1]
#> [1] 0.2225462
#> attr(,"class")
#> [1] "bayesianVARs_coef"  "bayesianVARs_draws"
```
