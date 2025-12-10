# Markov Chain Monte Carlo Sampling for Bayesian Vectorautoregressions

`bvar` simulates from the joint posterior distribution of the parameters
and latent variables and returns the posterior draws.

## Usage

``` r
bvar(
  data,
  lags = 1L,
  draws = 1000L,
  burnin = 1000L,
  thin = 1L,
  prior_intercept = 10,
  prior_phi = specify_prior_phi(data = data, lags = lags, prior = "HS"),
  prior_sigma = specify_prior_sigma(data = data, type = "factor", quiet = TRUE),
  sv_keep = "last",
  quiet = FALSE,
  startvals = list(),
  expert = list()
)
```

## Arguments

- data:

  Data matrix (can be a time series object). Each of \\M\\ columns is
  assumed to contain a single time-series of length \\T\\.

- lags:

  Integer indicating the order of the VAR, i.e. the number of lags of
  the dependent variables included as predictors.

- draws:

  single integer indicating the number of draws after the burnin

- burnin:

  single integer indicating the number of draws discarded as burnin

- thin:

  single integer. Every \\thin\\th draw will be stored. Default is
  `thin=1L`.

- prior_intercept:

  Either `prior_intercept=FALSE` and no constant term (intercept) will
  be included. Or a numeric vector of length \\M\\ indicating the
  (fixed) prior standard deviations on the constant term. A single
  number will be recycled accordingly. Default is `prior_intercept=10`.

- prior_phi:

  `bayesianVARs_prior_phi` object specifying prior for the reduced form
  VAR coefficients. Best use constructor
  [`specify_prior_phi`](https://luisgruber.github.io/bayesianVARs/reference/specify_prior_phi.md).

- prior_sigma:

  `bayesianVARs_prior_sigma` object specifying prior for the
  variance-covariance matrix of the VAR. Best use constructor
  [`specify_prior_sigma`](https://luisgruber.github.io/bayesianVARs/reference/specify_prior_sigma.md).

- sv_keep:

  String equal to `"all"` or `"last"`. In case of `sv_keep = "last"`,
  the default, only draws for the very last log-variance \\h_T\\ are
  stored.

- quiet:

  logical value indicating whether information about the progress during
  sampling should be displayed during sampling (default is `TRUE`).

- startvals:

  optional list with starting values.

- expert:

  optional list with expert settings.

## Value

An object of type `bayesianVARs_bvar`, a list containing the following
objects:

- `PHI`: A `bayesianVARs_coef` object, an array, containing the
  posterior draws of the VAR coefficients (including the intercept).

- `U`: A `bayesianVARs_draws` object, a matrix, containing the posterior
  draws of the contemporaneous coefficients (if cholesky decomposition
  for sigma is specified).

- `logvar`: A `bayesianVARs_draws` object containing the log-variance
  draws.

- `sv_para`: A `baysesianVARs_draws` object containing the posterior
  draws of the stochastic volatility related parameters.

- `phi_hyperparameter`: A matrix containing the posterior draws of the
  hyperparameters of the conditional normal prior on the VAR
  coefficients.

- `u_hyperparameter`: A matrix containing the posterior draws of the
  hyperparameters of the conditional normal prior on U (if cholesky
  decomposition for sigma is specified).

- `bench`: `proc_time` object containing the run time of the sampler.

- `V_prior`: An array containing the posterior draws of the variances of
  the conditional normal prior on the VAR coefficients.

- `facload`: A `bayesianVARs_draws` object, an array, containing draws
  from the posterior distribution of the factor loadings matrix (if
  factor decomposition for sigma is specified).

- `fac`: A `bayesianVARs_draws` object, an array, containing factor
  draws from the posterior distribution (if factor decomposition for
  sigma is specified).

- `Y`: Matrix containing the dependent variables used for estimation.

- `X` matrix containing the lagged values of the dependent variables,
  i.e. the covariates.

- `lags`: Integer indicating the lag order of the VAR.

- `intercept`: Logical indicating whether a constant term is included.

- `heteroscedastic` logical indicating whether heteroscedasticity is
  assumed.

- `Yraw`: Matrix containing the dependent variables, including the
  initial 'lags' observations.

- `Traw`: Integer indicating the total number of observations.

- `sigma_type`: Character specifying the decomposition of the
  variance-covariance matrix.

- `datamat`: Matrix containing both 'Y' and 'X'.

- `config`: List containing information on configuration parameters.

## Details

The VAR(p) model is of the following form: \\ \boldsymbol{y}^\prime_t =
\boldsymbol{\iota}^\prime + \boldsymbol{x}^\prime_t\boldsymbol{\Phi} +
\boldsymbol{\epsilon}^\prime_t\\, where \\\boldsymbol{y}\_t\\ is a
\\M\\-dimensional vector of dependent variables and
\\\boldsymbol{\epsilon}\_t\\ is the error term of the same dimension.
\\\boldsymbol{x}\_t\\ is a \\K=pM\\-dimensional vector containing
lagged/past values of the dependent variables \\\boldsymbol{y}\_{t-l}\\
for \\l=1,\dots,p\\ and \\\boldsymbol{\iota}\\ is a constant term
(intercept) of dimension \\M\times 1\\. The reduced-form coefficient
matrix \\\boldsymbol{\Phi}\\ is of dimension \\K \times M\\.

`bvar` offers two different specifications for the errors: The user can
choose between a factor stochastic volatility structure or a cholesky
stochastic volatility structure. In both cases the disturbances
\\\boldsymbol{\epsilon}\_t\\ are assumed to follow a \\M\\-dimensional
multivariate normal distribution with zero mean and variance-covariance
matrix \\\boldsymbol{\Sigma}\_t\\. In case of the cholesky specification
\\\boldsymbol{\Sigma}\_t = \boldsymbol{U}^{\prime -1} \boldsymbol{D}\_t
\boldsymbol{U}^{-1}\\, where \\\boldsymbol{U}^{-1}\\ is upper
unitriangular (with ones on the diagonal). The diagonal matrix
\\\boldsymbol{D}\_t\\ depends upon latent log-variances, i.e.
\\\boldsymbol{D}\_t=diag(exp(h\_{1t}),\dots, exp(h\_{Mt})\\. The
log-variances follow a priori independent autoregressive processes
\\h\_{it}\sim N(\mu_i + \phi_i(h\_{i,t-1}-\mu_i),\sigma_i^2)\\ for
\\i=1,\dots,M\\. In case of the factor structure,
\\\boldsymbol{\Sigma}\_t = \boldsymbol{\Lambda} \boldsymbol{V}\_t
\boldsymbol{\Lambda}^\prime + \boldsymbol{G}\_t\\. The diagonal matrices
\\\boldsymbol{V}\_t\\ and \\\boldsymbol{G}\_t\\ depend upon latent
log-variances, i.e. \\\boldsymbol{G}\_t=diag(exp(h\_{1t}),\dots,
exp(h\_{Mt})\\ and \\\boldsymbol{V}\_t=diag(exp(h\_{M+1,t}),\dots,
exp(h\_{M+r,t})\\. The log-variances follow a priori independent
autoregressive processes \\h\_{it}\sim N(\mu_i +
\phi_i(h\_{i,t-1}-\mu_i),\sigma_i^2)\\ for \\i=1,\dots,M\\ and
\\h\_{M+j,t}\sim N(\phi_ih\_{M+j,t-1},\sigma\_{M+j}^2)\\ for
\\j=1,\dots,r\\.

## MCMC algorithm

To sample efficiently the reduced-form VAR coefficients assuming a
**factor structure for the errors**, the equation per equation algorithm
in Kastner & Huber (2020) is implemented. All parameters and latent
variables associated with the factor-structure are sampled using package
[`factorstochvol-package`](https://rdrr.io/pkg/factorstochvol/man/factorstochvol-package.html)'s
function `update_fsv` callable on the C-level only.

To sample efficiently the reduced-form VAR coefficients, assuming a
**cholesky-structure for the errors**, the corrected triangular
algorithm in Carriero et al. (2021) is implemented. The SV parameters
and latent variables are sampled using package
[`stochvol`](https://gregorkastner.github.io/stochvol/reference/stochvol-package.html)'s
[`update_fast_sv`](https://gregorkastner.github.io/stochvol/reference/update_fast_sv.html)
function. The precision parameters, i.e. the free off-diagonal elements
in \\\boldsymbol{U}\\, are sampled as in Cogley and Sargent (2005).

## References

Gruber, L. and Kastner, G. (2025). Forecasting macroeconomic data with
Bayesian VARs: Sparse or dense? It depends! *International Journal of
Forecasting*.
[doi:10.1016/j.ijforecast.2025.02.001](https://doi.org/10.1016/j.ijforecast.2025.02.001)
.

Kastner, G. and Huber, F. Sparse (2020). Bayesian vector autoregressions
in huge dimensions. *Journal of Forecasting*. **39**, 1142–1165,
[doi:10.1002/for.2680](https://doi.org/10.1002/for.2680) .

Kastner, G. (2019). Sparse Bayesian Time-Varying Covariance Estimation
in Many Dimensions *Journal of Econometrics*, **210**(1), 98–115,
[doi:10.1016/j.jeconom.2018.11.007](https://doi.org/10.1016/j.jeconom.2018.11.007)
.

Carriero, A. and Chan, J. and Clark, T. E. and Marcellino, M. (2021).
Corrigendum to “Large Bayesian vector autoregressions with stochastic
volatility and non-conjugate priors” \[J. Econometrics 212 (1) (2019)
137–154\]. *Journal of Econometrics*,
[doi:10.1016/j.jeconom.2021.11.010](https://doi.org/10.1016/j.jeconom.2021.11.010)
.

Cogley, S. and Sargent, T. (2005). Drifts and volatilities: monetary
policies and outcomes in the post WWII US. *Review of Economic
Dynamics*, **8**, 262–302,
[doi:10.1016/j.red.2004.10.009](https://doi.org/10.1016/j.red.2004.10.009)
.

Hosszejni, D. and Kastner, G. (2021). Modeling Univariate and
Multivariate Stochastic Volatility in R with stochvol and
factorstochvol. *Journal of Statistical Software*, *100*, 1–-34.
[doi:10.18637/jss.v100.i12](https://doi.org/10.18637/jss.v100.i12) .

## See also

- Helpers for prior configuration:
  [`specify_prior_phi()`](https://luisgruber.github.io/bayesianVARs/reference/specify_prior_phi.md),
  [`specify_prior_sigma()`](https://luisgruber.github.io/bayesianVARs/reference/specify_prior_sigma.md).

- Plotting:
  [`plot.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_bvar.md).

- Extractors:
  [`coef.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/coef.md),
  [`vcov.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/vcov.bayesianVARs_bvar.md).

- 'stable' bvar:
  [`stable_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/stable_bvar.md).

- summary method:
  [`summary.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/summary.bayesianVARs_bvar.md).

- predict method:
  [`predict.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/predict.bayesianVARs_bvar.md).

- fitted method:
  [`fitted.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/fitted.bayesianVARs_bvar.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Plot
plot(mod)


# Summary
summary(mod)
#> 
#> Posterior median of reduced-form coefficients:
#>              GDPC1 CPIAUCSL FEDFUNDS
#> GDPC1.l1     0.153   -0.003    0.020
#> CPIAUCSL.l1 -0.098    0.607   -0.001
#> FEDFUNDS.l1  0.014    0.043    1.001
#> intercept    0.006    0.001    0.000
#> 
#> Posterior interquartile range of of reduced-form coefficients:
#>             GDPC1 CPIAUCSL FEDFUNDS
#> GDPC1.l1    0.111    0.028    0.018
#> CPIAUCSL.l1 0.140    0.086    0.017
#> FEDFUNDS.l1 0.020    0.014    0.008
#> intercept   0.001    0.001    0.000
#> 
#> Posterior median of factor loadings:
#>          factor1
#> GDPC1      0.000
#> CPIAUCSL   0.000
#> FEDFUNDS   0.001
#> 
#> Posterior interquartile range of factor loadings:
#>          factor1
#> GDPC1      0.000
#> CPIAUCSL   0.000
#> FEDFUNDS   0.001
```
