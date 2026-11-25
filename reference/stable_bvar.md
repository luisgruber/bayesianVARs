# Stable posterior draws

`stable_bvar()` detects and discards all posterior draws of an
`bayesianVARs_bvar` object that do not fulfill the stability condition:
A VAR(p) model is considered as stable only if all eigenvalues of the
companion form matrix lie inside the unit circle, i.e. have modulus less
than 1, see Lütkepohl (2005).

## Usage

``` r
stable_bvar(object, lim = 1, quiet = FALSE)
```

## Arguments

- object:

  A `bayesianVARs_bvar` object obtained via
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- lim:

  positive numerical less-than or equal to 1. Posterior draws for which
  all eigenvalues of the companion form matrix have modulus less than
  `lim` are kept and the remaining ones are discarded.

- quiet:

  logical indicating whether informative output should be omitted.

## Value

An object of type `bayesianVARs_bvar`.

## References

Lütkepohl, H. (2005). Stable Vector Autoregressive Processes. In: New
Introduction to Multiple Time Series Analysis. Springer, Berlin,
Heidelberg.
[doi:doi.org/10.1007/978-3-540-27752-1_2](https://doi.org/doi.org/10.1007/978-3-540-27752-1_2)
.

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Discard "unstable" draws
stable_mod <- stable_bvar(mod)
#> Original 'bayesianVARs_bvar' object consists of 1000 posterior draws.
#> Detected 618 unstable draws.
#> Remaining draws: 382!
```
