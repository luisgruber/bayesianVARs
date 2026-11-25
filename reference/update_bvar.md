# Update MCMC-chain obtained via bvar()

`update_bvar()` simulates from the joint posterior distribution of the
parameters and latent variables using the latest state of an existing
`bayesianVARs_bvar` object obtained via
[`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md)
as starting values and returns the posterior draws.

## Usage

``` r
update_bvar(object, draws = 1000L, quiet = FALSE)
```

## Arguments

- object:

  `bayesianVARs_bvar` object obtained via
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- draws:

  single integer indicating the number of draws.

- quiet:

  logical value indicating whether information about the progress during
  sampling should be displayed during sampling.

## Value

An object of type 'bayesianVARs_bvar', see
[`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, draws = 10L, burnin = 0L, quiet = TRUE)

# Update model
mod_new <- update_bvar(mod, draws = 10L, quiet = TRUE)
```
