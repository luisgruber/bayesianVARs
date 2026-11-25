# Predict method for Bayesian VARs

Simulates from (out-of-sample) predictive density for Bayesian VARs
estimated via
[`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md)
and computes log predictive likelhoods if ex-post observed data is
supplied.

## Usage

``` r
# S3 method for class 'bayesianVARs_bvar'
predict(
  object,
  ahead = 1L,
  each = 1L,
  stable = TRUE,
  stable_lim = 1,
  simulate_predictive = TRUE,
  LPL = FALSE,
  Y_obs = NA,
  LPL_VoI = NA,
  quiet = FALSE,
  ...
)
```

## Arguments

- object:

  A `bayesianVARs_bvar` object, obtained from
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- ahead:

  Integer vector (or coercible to such), indicating the number of steps
  ahead at which to predict.

- each:

  Single integer (or coercible to such) indicating how often should be
  drawn from the posterior predictive distribution for each draw that
  has been stored during MCMC sampling.

- stable:

  logical indicating whether to consider only those draws from the
  posterior that fulfill the 'stable' criterion. Default is `TRUE`.

- stable_lim:

  postive numerical less-than or equal to 1. Only those posterior draws
  are used for prediction for which all eigenvalues of the companion
  form matrix have modulus less than `stable_lim`, see
  [`stable_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/stable_bvar.md).
  Only used if `stable=TRUE`.

- simulate_predictive:

  logical, indicating whether the posterior predictive distribution
  should be simulated.

- LPL:

  logical indicating whether `ahead`-step-ahead log predictive
  likelihoods should be computed. If `LPL=TRUE`, `Y_obs` has to be
  specified.

- Y_obs:

  Data matrix of observed values for computation of log predictive
  likelihood. Each of `ncol(object$Yraw)` columns is assumed to contain
  a single time-series of length `length(ahead)`.

- LPL_VoI:

  either integer vector or character vector of column-names indicating
  for which subgroup of time-series in `object$Yraw` a joint log
  predictive likelihood shall be computed.

- quiet:

  logical indicating whether informative output should be omitted.

- ...:

  Currently ignored!

## Value

Object of class `bayesianVARs_predict`, a list that may contain the
following elements:

- `predictions` array of dimensions
  `c(length(ahead), ncol(object$Yraw), each * dim(object$PHI)[3])`
  containing the simulations from the predictive density (if
  `simulate_predictive=TRUE`).

- `LPL` vector of length `length(ahead)` containing the
  log-predictive-likelihoods (taking into account the joint distribution
  of all variables) (if `LPL=TRUE`).

- `LPL_univariate` matrix of dimension
  `c(length(ahead), ncol(object$Yraw)` containing the marginalized
  univariate log-predictive-likelihoods of each series (if `LPL=TRUE`).

- `LPL_VoI` vector of length `length(ahead)` containing the
  log-predictive-likelihoods for a subset of variables (if `LPL=TRUE`
  and `LPL_VoI != NA`).

- `Yraw` matrix containing the data used for the estimation of the VAR.

- `LPL_draws` matrix containing the simulations of the
  log-predictive-likelihood (if `LPL=TRUE`).

- `PL_univariate_draws` array containing the simulations of the
  univariate predictive-likelihoods (if `LPL=TRUE`).

- `LPL_sub_draws` matrix containing the simulations of the
  log-predictive-likelihood for a subset of variables (if `LPL=TRUE` and
  `LPL_VoI != NA`).

## See also

[`stable_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/stable_bvar.md),
[`plot.bayesianVARs_predict()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_predict.md),
[`pairs.bayesianVARs_predict()`](https://luisgruber.github.io/bayesianVARs/reference/pairs_predict.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Split data in train and test
train <- data[1:(nrow(data)-4),]
test <- data[-c(1:(nrow(data)-4)),]

# Estimate model using train data only
mod <- bvar(train, quiet = TRUE)

# Simulate from 1-step to 4-steps ahead posterior predictive and compute
# log-predictive-likelihoods
predictions <- predict(mod, ahead = 1:4, LPL = TRUE, Y_obs = test)
#> 'stable=TRUE': Calling 'stable_bvar()' to discard those posterior
#>           draws that do not fulfill the stable criterion.
#> Using the remaining 542 stable posterior draws for prediction!

# Summary
summary(predictions)
#> 
#> LPL:
#>   t+1   t+2   t+3   t+4 
#> 2.964 8.977 8.829 6.499 
#> 
#> Marginal univariate LPLs:
#>       GDPC1 CPIAUCSL FEDFUNDS
#> t+1 -0.9551   0.4695    3.702
#> t+2  2.7591   2.7089    3.396
#> t+3  2.9269   2.5829    3.229
#> t+4  3.0538   0.4387    3.081
#> 
#> Prediction quantiles:
#> , , GDPC1
#> 
#>          t+1       t+2       t+3       t+4
#> 5%  -0.06967 -0.048658 -0.036069 -0.025065
#> 50% -0.02097  0.001588  0.007057  0.009054
#> 95%  0.03566  0.039089  0.046252  0.041513
#> 
#> , , CPIAUCSL
#> 
#>           t+1       t+2       t+3       t+4
#> 5%  -0.019291 -0.017333 -0.016687 -0.013098
#> 50% -0.008315 -0.005376 -0.003321 -0.001078
#> 95%  0.002702  0.007175  0.010066  0.012650
#> 
#> , , FEDFUNDS
#> 
#>           t+1       t+2       t+3       t+4
#> 5%  -0.020487 -0.027833 -0.035777 -0.040580
#> 50% -0.003742 -0.004332 -0.005792 -0.004706
#> 95%  0.014443  0.022497  0.023405  0.036294
#> 

# Visualize via fan-charts
plot(predictions)


# \donttest{
# In order to evaluate the joint predictive density of a subset of the
# variables (variables of interest), consider specifying 'LPL_VoI':
predictions <- predict(mod, ahead = 1:4, LPL = TRUE, Y_obs = test, LPL_VoI = c("GDPC1","FEDFUNDS"))
#> 'stable=TRUE': Calling 'stable_bvar()' to discard those posterior
#>           draws that do not fulfill the stable criterion.
#> Using the remaining 542 stable posterior draws for prediction!
predictions$LPL_VoI
#>      t+1      t+2      t+3      t+4 
#> 2.612975 6.185825 6.148418 6.127385 
# }
```
