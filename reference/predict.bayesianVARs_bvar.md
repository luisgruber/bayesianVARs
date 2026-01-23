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
  simulate_predictive = TRUE,
  LPL = FALSE,
  Y_obs = NA,
  LPL_VoI = NA,
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
#> 
#> 492/1000 stable posterior draws remaining for prediction!

# Summary
summary(predictions)
#> 
#> LPL:
#>   t+1   t+2   t+3   t+4 
#> 3.234 8.839 8.629 6.442 
#> 
#> Marginal univariate LPLs:
#>      GDPC1 CPIAUCSL FEDFUNDS
#> t+1 -1.046   0.4946    3.587
#> t+2  2.658   2.7642    3.284
#> t+3  2.780   2.6174    3.123
#> t+4  2.893   0.5619    2.996
#> 
#> Prediction quantiles:
#> , , GDPC1
#> 
#>          t+1       t+2       t+3       t+4
#> 5%  -0.08087 -0.036587 -0.023680 -0.020865
#> 50% -0.02243 -0.001922  0.004844  0.006433
#> 95%  0.03125  0.037331  0.031066  0.032412
#> 
#> , , CPIAUCSL
#> 
#>           t+1       t+2       t+3       t+4
#> 5%  -0.017982 -0.019123 -0.018902 -0.013825
#> 50% -0.007775 -0.005189 -0.002953 -0.001493
#> 95%  0.003257  0.006441  0.010182  0.012025
#> 
#> , , FEDFUNDS
#> 
#>           t+1       t+2       t+3       t+4
#> 5%  -0.024906 -0.033167 -0.037943 -0.046417
#> 50% -0.004319 -0.004624 -0.006759 -0.004937
#> 95%  0.015106  0.019991  0.030742  0.042475
#> 

# Visualize via fan-charts
plot(predictions)


# \donttest{
# In order to evaluate the joint predictive density of a subset of the
# variables (variables of interest), consider specifying 'LPL_VoI':
predictions <- predict(mod, ahead = 1:4, LPL = TRUE, Y_obs = test, LPL_VoI = c("GDPC1","FEDFUNDS"))
#> 'stable=TRUE': Calling 'stable_bvar()' to discard those posterior
#>           draws that do not fulfill the stable criterion.
#> 
#> 492/1000 stable posterior draws remaining for prediction!
predictions$LPL_VoI
#>      t+1      t+2      t+3      t+4 
#> 2.485835 6.018957 5.981567 5.954947 
# }
```
