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
#> 661/1000 stable posterior draws remaining for prediction!

# Summary
summary(predictions)
#> 
#> LPL:
#>   t+1   t+2   t+3   t+4 
#> 3.179 8.841 8.710 6.504 
#> 
#> Marginal univariate LPLs:
#>       GDPC1 CPIAUCSL FEDFUNDS
#> t+1 -0.7522   0.5097    3.586
#> t+2  2.6422   2.7711    3.296
#> t+3  2.8231   2.6009    3.139
#> t+4  2.9553   0.5632    3.009
#> 
#> Prediction quantiles:
#> , , GDPC1
#> 
#>          t+1       t+2       t+3       t+4
#> 5%  -0.08795 -0.044676 -0.031418 -0.025550
#> 50% -0.02285 -0.002203  0.004421  0.007084
#> 95%  0.03916  0.046564  0.037960  0.044622
#> 
#> , , CPIAUCSL
#> 
#>           t+1       t+2       t+3      t+4
#> 5%  -0.018304 -0.019707 -0.018021 -0.01523
#> 50% -0.007872 -0.005365 -0.002934 -0.00120
#> 95%  0.003491  0.008279  0.009535  0.01203
#> 
#> , , FEDFUNDS
#> 
#>           t+1       t+2       t+3       t+4
#> 5%  -0.022097 -0.030451 -0.040058 -0.048608
#> 50% -0.004467 -0.004472 -0.006268 -0.005041
#> 95%  0.012831  0.020932  0.032039  0.045410
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
#> 661/1000 stable posterior draws remaining for prediction!
predictions$LPL_VoI
#>      t+1      t+2      t+3      t+4 
#> 2.702619 6.078233 6.054422 6.015822 
# }
```
