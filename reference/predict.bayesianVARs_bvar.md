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
#> 552/1000 stable posterior draws remaining for prediction!

# Summary
summary(predictions)
#> 
#> LPL:
#>   t+1   t+2   t+3   t+4 
#> 4.234 9.586 9.057 6.262 
#> 
#> Marginal univariate LPLs:
#>       GDPC1 CPIAUCSL FEDFUNDS
#> t+1 -0.6468   0.5340    3.592
#> t+2  2.7997   2.8042    3.416
#> t+3  2.8024   2.5535    3.322
#> t+4  2.8997   0.6063    3.277
#> 
#> Prediction quantiles:
#> , , GDPC1
#> 
#>          t+1       t+2       t+3       t+4
#> 5%  -0.07360 -0.044123 -0.025555 -0.024188
#> 50% -0.01864  0.001661  0.005585  0.007187
#> 95%  0.05851  0.056716  0.040457  0.041119
#> 
#> , , CPIAUCSL
#> 
#>           t+1      t+2       t+3        t+4
#> 5%  -0.018125 -0.01744 -0.015966 -0.0170495
#> 50% -0.006756 -0.00406 -0.001766 -0.0003074
#> 95%  0.003689  0.00796  0.011154  0.0115700
#> 
#> , , FEDFUNDS
#> 
#>           t+1       t+2       t+3       t+4
#> 5%  -0.016168 -0.024425 -0.029679 -0.036915
#> 50% -0.004118 -0.004735 -0.005671 -0.004132
#> 95%  0.009494  0.013771  0.020476  0.024615
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
#> 552/1000 stable posterior draws remaining for prediction!
predictions$LPL_VoI
#>      t+1      t+2      t+3      t+4 
#> 2.752971 6.625609 6.435607 6.326910 
# }
```
