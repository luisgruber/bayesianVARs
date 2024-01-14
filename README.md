
# bayesianVARs <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/luisgruber/bayesianVARs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/luisgruber/bayesianVARs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Estimation of Bayesian vectorautoregressions with/without stochastic
volatility.

Implements several modern hierarchical shrinkage priors, amongst them
Dirichlet-Laplace prior (DL), hierarchical Minnesota prior (HM),
Horseshoe prior (HS), normal-gamma prior (NG),
$R^2$-induced-Dirichlet-decomposition prior (R2D2) and stochastic search
variable selection prior (SSVS).

Concerning the error-term, the user can either specify an
order-invariant factor structure or an order-variant cholesky structure.

# Installation

Install directly from GitHub.

``` r
devtools::install_github("luisgruber/bayesianVARs")
```

# Usage

The main workhorse to conduct Bayesian inference for
vectorautoregression models in this package is the function `bvar()`.

Some features:

- Prediction, plotting, extraction of model parameters and extraction of
  fitted values with the usual generic functions `predict()`, `plot()`,
  `coef()`, `vcov()` and `fitted()`.
- Configure prior distributions with helper functions
  `specify_prior_phi()` and `specify_prior_sigma()`.

# Demonstration

``` r
set.seed(537)
# load package
library(bayesianVARs)

# Load data
train_data <-100 * usmacro_growth[1:237,c("GDPC1", "PCECC96", "GPDIC1", "AWHMAN", "GDPCTPI", "CES2000000008x", "FEDFUNDS", "GS10", "EXUSUKx", "S&P 500")]
test_data <-100 * usmacro_growth[238:241,c("GDPC1", "PCECC96", "GPDIC1", "AWHMAN", "GDPCTPI", "CES2000000008x", "FEDFUNDS", "GS10", "EXUSUKx", "S&P 500")]
                                   
# Estimate model using default prior settings
mod <- bvar(train_data, lags = 2L, draws = 2000, burnin = 1000, sv_keep = "all")

# Out of sample prediction and log-predictive-likelihood evaluation
pred <- predict(mod, ahead = 1:4, LPL = TRUE, Y_obs = test_data)

# Visualize in-sample fit plus out-of-sample prediction intervals
plot(mod, predictions = pred)
```

# Documentation

``` r
browseVignettes("bayesianVARs")
```
