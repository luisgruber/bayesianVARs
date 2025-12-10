# Stable posterior draws

`stable_bvar()` detects and discards all posterior draws of an
`bayesianVARs_bvar` object that do not fulfill the stability condition:
A VAR(p) model is considered as stable only if the eigenvalues of the
companion form matrix lie inside the unit circle.

## Usage

``` r
stable_bvar(object, quiet = FALSE)
```

## Arguments

- object:

  A `bayesianVARs_bvar` object obtained via
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- quiet:

  logical indicating whether informative output should be omitted.

## Value

An object of type `bayesianVARs_bvar`.

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Discard "unstable" draws
stable_mod <- stable_bvar(mod)
#> 
#> Original 'bayesianVARs_bvar' object consists of 1000 posterior draws.
#> 
#> Detected 552 unstable draws.
#> 
#> Remaining draws: 448 !
```
