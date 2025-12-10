# Specify prior on Sigma

Configures prior on the variance-covariance of the VAR.

## Usage

``` r
specify_prior_sigma(
  data = NULL,
  M = ncol(data),
  type = c("factor", "cholesky"),
  factor_factors = 1L,
  factor_restrict = c("none", "upper"),
  factor_priorfacloadtype = c("rowwiseng", "colwiseng", "normal"),
  factor_priorfacload = 0.1,
  factor_facloadtol = 1e-18,
  factor_priorng = c(1, 1),
  factor_priormu = c(0, 10),
  factor_priorphiidi = c(10, 3),
  factor_priorphifac = c(10, 3),
  factor_priorsigmaidi = 1,
  factor_priorsigmafac = 1,
  factor_priorh0idi = "stationary",
  factor_priorh0fac = "stationary",
  factor_heteroskedastic = TRUE,
  factor_priorhomoskedastic = NA,
  factor_interweaving = 4,
  cholesky_U_prior = c("HS", "DL", "R2D2", "NG", "SSVS", "normal", "HMP"),
  cholesky_U_tol = 1e-18,
  cholesky_heteroscedastic = TRUE,
  cholesky_priormu = c(0, 100),
  cholesky_priorphi = c(20, 1.5),
  cholesky_priorsigma2 = c(0.5, 0.5),
  cholesky_priorh0 = "stationary",
  cholesky_priorhomoscedastic = as.numeric(NA),
  cholesky_DL_a = "1/n",
  cholesky_DL_tol = 0,
  cholesky_R2D2_a = 0.4,
  cholesky_R2D2_b = 0.5,
  cholesky_R2D2_tol = 0,
  cholesky_NG_a = 0.5,
  cholesky_NG_b = 0.5,
  cholesky_NG_c = 0.5,
  cholesky_NG_tol = 0,
  cholesky_SSVS_c0 = 0.001,
  cholesky_SSVS_c1 = 1,
  cholesky_SSVS_p = 0.5,
  cholesky_HMP_lambda3 = c(0.01, 0.01),
  cholesky_normal_sds = 10,
  expert_sv_offset = 0,
  quiet = FALSE,
  ...
)
```

## Arguments

- data:

  Optional. Data matrix (can be a time series object). Each of \\M\\
  columns is assumed to contain a single time-series of length \\T\\.

- M:

  positive integer indicating the number of time-series of the VAR.

- type:

  character, one of `"factor"` (the default) or `"cholesky"`, indicating
  which decomposition to be applied to the covariance-matrix.

- factor_factors:

  Number of latent factors to be estimated. Only required if
  `type="factor"`.

- factor_restrict:

  Either "upper" or "none", indicating whether the factor loadings
  matrix should be restricted to have zeros above the diagonal ("upper")
  or whether all elements should be estimated from the data ("none").
  Setting `restrict` to "upper" often stabilizes MCMC estimation and can
  be important for identifying the factor loadings matrix, however, it
  generally is a strong prior assumption. Setting `restrict` to "none"
  is usually the preferred option if identification of the factor
  loadings matrix is of less concern but covariance estimation or
  prediction is the goal. Only required if `type="factor"`.

- factor_priorfacloadtype:

  Can be `"normal"`, `"rowwiseng"`, `"colwiseng"`. Only required if
  `type="factor"`.

  `"normal"`:

  :   Normal prior. The value of `priorfacload` is interpreted as the
      standard deviations of the Gaussian prior distributions for the
      factor loadings.

  `"rowwiseng"`:

  :   Row-wise Normal-Gamma prior. The value of `priorfacload` is
      interpreted as the shrinkage parameter `a`.

  `"colwiseng"`:

  :   Column-wise Normal-Gamma prior. The value of `priorfacload` is
      interpreted as the shrinkage parameter `a`.

  For details please see Kastner (2019).

- factor_priorfacload:

  Either a matrix of dimensions `M` times `factor_factors` with positive
  elements or a single number (which will be recycled accordingly). Only
  required if `type="factor"`. The meaning of `factor_priorfacload`
  depends on the setting of `factor_priorfacloadtype` and is explained
  there.

- factor_facloadtol:

  Minimum number that the absolute value of a factor loadings draw can
  take. Prevents numerical issues that can appear when strong shrinkage
  is enforced if chosen to be greater than zero. Only required if
  `type="factor"`.

- factor_priorng:

  Two-element vector with positive entries indicating the Normal-Gamma
  prior's hyperhyperparameters `c` and `d` (cf. Kastner (2019)). Only
  required if `type="factor"`.

- factor_priormu:

  Vector of length 2 denoting prior mean and standard deviation for
  unconditional levels of the idiosyncratic log variance processes. Only
  required if `type="factor"`.

- factor_priorphiidi:

  Vector of length 2, indicating the shape parameters for the Beta prior
  distributions of the transformed parameters `(phi+1)/2`, where `phi`
  denotes the persistence of the idiosyncratic log variances. Only
  required if `type="factor"`.

- factor_priorphifac:

  Vector of length 2, indicating the shape parameters for the Beta prior
  distributions of the transformed parameters `(phi+1)/2`, where `phi`
  denotes the persistence of the factor log variances. Only required if
  `type="factor"`.

- factor_priorsigmaidi:

  Vector of length `M` containing the prior volatilities of log
  variances. If `factor_priorsigmaidi` has exactly one element, it will
  be recycled for all idiosyncratic log variances. Only required if
  `type="factor"`.

- factor_priorsigmafac:

  Vector of length `factor_factors` containing the prior volatilities of
  log variances. If `factor_priorsigmafac` has exactly one element, it
  will be recycled for all factor log variances. Only required if
  `type="factor"`.

- factor_priorh0idi:

  Vector of length 1 or `M`, containing information about the Gaussian
  prior for the initial idiosyncratic log variances. Only required if
  `type="factor"`. If an element of `factor_priorh0idi` is a nonnegative
  number, the conditional prior of the corresponding initial log
  variance h0 is assumed to be Gaussian with mean 0 and standard
  deviation `factor_priorh0idi` times \\sigma\\. If an element of
  `factor_priorh0idi` is the string 'stationary', the prior of the
  corresponding initial log volatility is taken to be from the
  stationary distribution, i.e. h0 is assumed to be Gaussian with mean 0
  and variance \\sigma^2/(1-phi^2)\\.

- factor_priorh0fac:

  Vector of length 1 or `factor_factors`, containing information about
  the Gaussian prior for the initial factor log variances. Only required
  if `type="factor"`. If an element of `factor_priorh0fac` is a
  nonnegative number, the conditional prior of the corresponding initial
  log variance h0 is assumed to be Gaussian with mean 0 and standard
  deviation `factor_priorh0fac` times \\sigma\\. If an element of
  `factor_priorh0fac` is the string 'stationary', the prior of the
  corresponding initial log volatility is taken to be from the
  stationary distribution, i.e. h0 is assumed to be Gaussian with mean 0
  and variance \\sigma^2/(1-phi^2)\\.

- factor_heteroskedastic:

  Vector of length 1, 2, or `M + factor_factors`, containing logical
  values indicating whether time-varying
  (`factor_heteroskedastic = TRUE`) or constant
  (`factor_heteroskedastic = FALSE`) variance should be estimated. If
  `factor_heteroskedastic` is of length 2 it will be recycled
  accordingly, whereby the first element is used for all idiosyncratic
  variances and the second element is used for all factor variances.
  Only required if `type="factor"`.

- factor_priorhomoskedastic:

  Only used if at least one element of `factor_heteroskedastic` is set
  to `FALSE`. In that case, `factor_priorhomoskedastic` must be a matrix
  with positive entries and dimension c(M, 2). Values in column 1 will
  be interpreted as the shape and values in column 2 will be interpreted
  as the rate parameter of the corresponding inverse gamma prior
  distribution of the idiosyncratic variances. Only required if
  `type="factor"`.

- factor_interweaving:

  The following values for interweaving the factor loadings are accepted
  (Only required if `type="factor"`):

  0:

  :   No interweaving.

  1:

  :   Shallow interweaving through the diagonal entries.

  2:

  :   Deep interweaving through the diagonal entries.

  3:

  :   Shallow interweaving through the largest absolute entries in each
      column.

  4:

  :   Deep interweaving through the largest absolute entries in each
      column.

  For details please see Kastner et al. (2017). A value of 4 is the
  highly recommended default.

- cholesky_U_prior:

  character, one of `"HS"`, `"R2D2"`, `"NG"`, `"DL"`, `"SSVS"`, `"HMP"`
  or `"normal"`. Only required if `type="cholesky"`.

- cholesky_U_tol:

  Minimum number that the absolute value of an free off-diagonal element
  of an \\U\\-draw can take. Prevents numerical issues that can appear
  when strong shrinkage is enforced if chosen to be greater than zero.
  Only required if `type="cholesky"`.

- cholesky_heteroscedastic:

  single logical indicating whether time-varying
  (`cholesky_heteroscedastic = TRUE`) or constant
  (`cholesky_heteroscedastic = FALSE`) variance should be estimated.
  Only required if `type="cholesky"`.

- cholesky_priormu:

  Vector of length 2 denoting prior mean and standard deviation for
  unconditional levels of the log variance processes. Only required if
  `type="cholesky"`.

- cholesky_priorphi:

  Vector of length 2, indicating the shape parameters for the Beta prior
  distributions of the transformed parameters `(phi+1)/2`, where `phi`
  denotes the persistence of the log variances. Only required if
  `type="cholesky"`.

- cholesky_priorsigma2:

  Vector of length 2, indicating the shape and the rate for the Gamma
  prior distributions on the variance of the log variance processes.
  (Currently only one global setting for all \\M\\ processes is
  supported). Only required if `type="cholesky"`.

- cholesky_priorh0:

  Vector of length 1 or `M`, containing information about the Gaussian
  prior for the initial idiosyncratic log variances. Only required if
  `type="cholesky"`. If an element of `cholesky_priorh0` is a
  nonnegative number, the conditional prior of the corresponding initial
  log variance h0 is assumed to be Gaussian with mean 0 and standard
  deviation `cholesky_priorh0` times \\sigma\\. If an element of
  `cholesky_priorh0` is the string 'stationary', the prior of the
  corresponding initial log volatility is taken to be from the
  stationary distribution, i.e. h0 is assumed to be Gaussian with mean 0
  and variance \\sigma^2/(1-phi^2)\\.

- cholesky_priorhomoscedastic:

  Only used if `cholesky_heteroscedastic=FALSE`. In that case,
  `cholesky_priorhomoscedastic` must be a matrix with positive entries
  and dimension c(M, 2). Values in column 1 will be interpreted as the
  shape and values in column 2 will be interpreted as the scale
  parameter of the corresponding inverse gamma prior distribution of the
  variances. Only required if `type="cholesky"`.

- cholesky_DL_a:

  (Single) positive real number. The value is interpreted as the
  concentration parameter for the local scales. Smaller values enforce
  heavier shrinkage. A matrix of dimension `c(s,2)` specifies a discrete
  hyperprior, where the first column contains s support points and the
  second column contains the associated prior probabilities.
  `cholesky_DL_a` has only to be specified if `cholesky_U_prior="DL"`.

- cholesky_DL_tol:

  Minimum number that a parameter draw of one of the shrinking
  parameters of the Dirichlet Laplace prior can take. Prevents numerical
  issues that can appear when strong shrinkage is enforced if chosen to
  be greater than zero. `DL_tol` has only to be specified if
  `cholesky_U_prior="DL"`.

- cholesky_R2D2_a:

  (Single) positive real number. The value is interpreted as the
  concentration parameter for the local scales. Smaller values enforce
  heavier shrinkage. A matrix of dimension `c(s,2)` specifies a discrete
  hyperprior, where the first column contains s support points and the
  second column contains the associated prior probabilities.
  cholesky_R2D2_a has only to be specified if `cholesky_U_prior="R2D2"`.

- cholesky_R2D2_b:

  single positive number, where greater values indicate heavier
  regularization. `cholesky_R2D2_b` has only to be specified if
  `cholesky_U_prior="R2D2"`.

- cholesky_R2D2_tol:

  Minimum number that a parameter draw of one of the shrinking
  parameters of the R2D2 prior can take. Prevents numerical issues that
  can appear when strong shrinkage is enforced if chosen to be greater
  than zero. `cholesky_R2D2_tol` has only to be specified if
  `cholesky_U_prior="R2D2"`.

- cholesky_NG_a:

  (Single) positive real number. The value is interpreted as the
  concentration parameter for the local scales. Smaller values enforce
  heavier shrinkage. A matrix of dimension `c(s,2)` specifies a discrete
  hyperprior, where the first column contains s support points and the
  second column contains the associated prior probabilities.
  `cholesky_NG_a` has only to be specified if `cholesky_U_prior="NG"`.

- cholesky_NG_b:

  (Single) positive real number. The value indicates the shape parameter
  of the inverse gamma prior on the global scales. `cholesky_NG_b` has
  only to be specified if `cholesky_U_prior="NG"`.

- cholesky_NG_c:

  (Single) positive real number. The value indicates the scale parameter
  of the inverse gamma prior on the global scales. Expert option would
  be to set the scale parameter proportional to NG_a. E.g. in the case
  where a discrete hyperprior for NG_a is chosen, a desired proportion
  of let's say 0.2 is achieved by setting NG_c="0.2a" (character
  input!). `cholesky_NG_c` has only to be specified if
  `cholesky_U_prior="NG"`.

- cholesky_NG_tol:

  Minimum number that a parameter draw of one of the shrinking
  parameters of the normal-gamma prior can take. Prevents numerical
  issues that can appear when strong shrinkage is enforced if chosen to
  be greater than zero. `cholesky_NG_tol` has only to be specified if
  `cholesky_U_prior="NG"`.

- cholesky_SSVS_c0:

  single positive number indicating the (unscaled) standard deviation of
  the spike component. `cholesky_SSVS_c0` has only to be specified if
  `choleksy_U_prior="SSVS"`. It should be that \\SSVS\_{c0} \ll
  SSVS\_{c1}\\!

- cholesky_SSVS_c1:

  single positive number indicating the (unscaled) standard deviation of
  the slab component. `cholesky_SSVS_c1` has only to be specified if
  `choleksy_U_prior="SSVS"`. It should be that \\SSVS\_{c0} \ll
  SSVS\_{c1}\\!

- cholesky_SSVS_p:

  Either a single positive number in the range `(0,1)` indicating the
  (fixed) prior inclusion probability of each coefficient. Or numeric
  vector of length 2 with positive entries indicating the shape
  parameters of the Beta distribution. In that case a Beta hyperprior is
  placed on the prior inclusion probability. `cholesky_SSVS_p` has only
  to be specified if `choleksy_U_prior="SSVS"`.

- cholesky_HMP_lambda3:

  numeric vector of length 2. Both entries must be positive. The first
  indicates the shape and the second the rate of the Gamma hyperprior on
  the contemporaneous coefficients. `cholesky_HMP_lambda3` has only to
  be specified if `choleksy_U_prior="HMP"`.

- cholesky_normal_sds:

  numeric vector of length \\\frac{M^2-M}{2}\\, indicating the prior
  variances for the free off-diagonal elements in \\U\\. A single number
  will be recycled accordingly! Must be positive. `cholesky_normal_sds`
  has only to be specified if `choleksy_U_prior="normal"`.

- expert_sv_offset:

  ... Do not use!

- quiet:

  logical indicating whether informative output should be omitted.

- ...:

  Do not use!

## Value

Object of class `bayesianVARs_prior_sigma`.

## Details

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

## References

Kastner, G. (2019). Sparse Bayesian Time-Varying Covariance Estimation
in Many Dimensions *Journal of Econometrics*, **210**(1), 98–115,
[doi:10.1016/j.jeconom.2018.11.007](https://doi.org/10.1016/j.jeconom.2018.11.007)

Kastner, G., Frühwirth-Schnatter, S., and Lopes, H.F. (2017). Efficient
Bayesian Inference for Multivariate Factor Stochastic Volatility Models.
*Journal of Computational and Graphical Statistics*, **26**(4), 905–917,
[doi:10.1080/10618600.2017.1322091](https://doi.org/10.1080/10618600.2017.1322091)
.

## See also

[`specify_prior_phi()`](https://luisgruber.github.io/bayesianVARs/reference/specify_prior_phi.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# examples with stochastic volatility (heteroscedasticity) -----------------
# factor-decomposition with 2 factors and colwise normal-gamma prior on the loadings
sigma_factor_cng_sv <- specify_prior_sigma(data = data, type = "factor",
factor_factors = 2L, factor_priorfacloadtype = "colwiseng", factor_heteroskedastic = TRUE)
#> 
#> Since argument 'type' is specified with 'factor', all arguments starting with 'cholesky_' are being ignored.

# cholesky-decomposition with Dirichlet-Laplace prior on U
sigma_cholesky_dl_sv <- specify_prior_sigma(data = data, type = "cholesky",
cholesky_U_prior = "DL", cholesky_DL_a = 0.5, cholesky_heteroscedastic = TRUE)
#> 
#> Since argument 'type' is specified with 'cholesky', all arguments starting with 'factor_' are being ignored.

# examples without stochastic volatility (homoscedasticity) ----------------
# factor-decomposition with 2 factors and colwise normal-gamma prior on the loadings
sigma_factor_cng <- specify_prior_sigma(data = data, type = "factor",
factor_factors = 2L, factor_priorfacloadtype = "colwiseng",
factor_heteroskedastic = FALSE, factor_priorhomoskedastic = matrix(c(0.5,0.5),
ncol(data), 2))
#> 
#> Since argument 'type' is specified with 'factor', all arguments starting with 'cholesky_' are being ignored.
#> 
#> Cannot do deep factor_interweaving if (some) factor_factors are homoskedastic. Setting 'factor_interweaving' to 3.

# cholesky-decomposition with Horseshoe prior on U
sigma_cholesky_dl <- specify_prior_sigma(data = data, type = "cholesky",
cholesky_U_prior = "HS", cholesky_heteroscedastic = FALSE)
#> 
#> Since argument 'type' is specified with 'cholesky', all arguments starting with 'factor_' are being ignored.
#> 
#> Argument 'cholesky_priorhomoscedastic' not specified. Setting both shape and rate of inverse gamma prior equal to 0.01.

# \donttest{
# Estimate model with your prior configuration of choice
mod <- bvar(data, prior_sigma = sigma_factor_cng_sv, quiet = TRUE)
# }
```
