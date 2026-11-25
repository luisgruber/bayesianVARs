# Specify starting values

`get_startvals()` can be used to obtain starting values for
[`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).
If you want to specify your own starting values, it is recommended to
first run `get_startvals()`, and then adapt the resulting
`bayesianVARs_startvals` object accordingly.

## Usage

``` r
get_startvals(
  prior_phi,
  prior_sigma,
  intercept,
  data = NULL,
  nobs = nrow(data)
)
```

## Arguments

- prior_phi:

  `bayesianVARs_prior_phi` object. Best use constructor
  [`specify_prior_phi()`](https://luisgruber.github.io/bayesianVARs/reference/specify_prior_phi.md).

- prior_sigma:

  `bayesianVARs_prior_sigma` object. Best use constructor
  [`specify_prior_sigma()`](https://luisgruber.github.io/bayesianVARs/reference/specify_prior_sigma.md).

- intercept:

  single logical or numeric indicating whether a constant term
  (intercept) will be included. In case of numeric input, a value of 0
  is interpreted as `FALSE` and a positive non-zero value as `TRUE`.

- data:

  Optional. Data matrix.

- nobs:

  single integer, indicating the number of observations of `data`.

## Value

Object of class `bayesianVARs_startvals.`

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Specify priors
prior_phi <- specify_prior_phi(data, lags = 1L)
prior_sigma <- specify_prior_sigma(data, quiet = TRUE)

# Get starting values
startvals <- get_startvals(prior_phi, prior_sigma, intercept = 10, data = data)
# Estimate a model
mod <- bvar(data, lags = 1L, prior_intercept = 10, draws = 10L, burnin = 0L, 
            startvals = startvals, quiet = TRUE)
```
