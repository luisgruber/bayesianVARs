# Specify prior on PHI

Configures prior on PHI, the matrix of reduced-form VAR coefficients.

## Usage

``` r
specify_prior_phi(
  data = NULL,
  M = ncol(data),
  lags = 1L,
  prior = "HS",
  priormean = 0,
  PHI_tol = 1e-18,
  DL_a = "1/K",
  DL_tol = 0,
  R2D2_a = 0.1,
  R2D2_b = 0.5,
  R2D2_tol = 0,
  NG_a = 0.1,
  NG_b = 1,
  NG_c = 1,
  NG_tol = 0,
  SSVS_c0 = 0.01,
  SSVS_c1 = 100,
  SSVS_semiautomatic = TRUE,
  SSVS_p = 0.5,
  HMP_lambda1 = c(0.01, 0.01),
  HMP_lambda2 = c(0.01, 0.01),
  normal_sds = 10,
  global_grouping = "global",
  ...
)
```

## Arguments

- data:

  Optional. Data matrix (can be a time series object). Each of \\M\\
  columns is assumed to contain a single time-series of length \\T\\.

- M:

  positive integer indicating the number of time-series of the VAR.

- lags:

  positive integer indicating the order of the VAR, i.e. the number of
  lags of the dependent variables included as predictors.

- prior:

  character, one of `"HS"`, `"R2D2"`, `"NG"`, `"DL"`, `"SSVS"`, `"HMP"`
  or `"normal"`.

- priormean:

  real numbers indicating the prior means of the VAR coefficients. One
  single number means that the prior mean of all own-lag coefficients
  w.r.t. the first lag equals `priormean` and `0` else. A vector of
  length M means that the prior mean of the own-lag coefficients w.r.t.
  the first lag equals `priormean` and `0` else. If `priormean` is a
  matrix of dimension `c(lags*M,M)`, then each of the \\M\\ columns is
  assumed to contain `lags*M` prior means for the VAR coefficients of
  the respective VAR equations.

- PHI_tol:

  Minimum number that the absolute value of a VAR coefficient draw can
  take. Prevents numerical issues that can appear when strong shrinkage
  is enforced if chosen to be greater than zero.

- DL_a:

  (Single) positive real number. The value is interpreted as the
  concentration parameter for the local scales. Smaller values enforce
  heavier shrinkage. If the argument `global_grouping` specifies e.g.
  `k` groups, then `DL_a` can be a numeric vector of length `k` and the
  elements indicate the shrinkage in each group. A matrix of dimension
  `c(s,2)` specifies a discrete hyperprior, where the first column
  contains `s` support points and the second column contains the
  associated prior probabilities. `DL_a` has only to be specified if
  `prior="DL"`.

- DL_tol:

  Minimum number that a parameter draw of one of the shrinking
  parameters of the Dirichlet Laplace prior can take. Prevents numerical
  issues that can appear when strong shrinkage is enforced if chosen to
  be greater than zero. `DL_tol` has only to be specified if
  `prior="DL"`.

- R2D2_a:

  (Single) positive real number. The value is interpreted as the
  concentration parameter for the local scales. Smaller values enforce
  heavier shrinkage. If the argument `global_grouping` specifies e.g.
  `k` groups, then `R2D2_a` can be a numeric vector of length `k` and
  the elements indicate the shrinkage in each group. A matrix of
  dimension `c(s,2)` specifies a discrete hyperprior, where the first
  column contains `s` support points and the second column contains the
  associated prior probabilities. `R2D2_a` has only to be specified if
  `prior="R2D2"`.

- R2D2_b:

  (Single) positive real number. The value indicates the shape parameter
  of the inverse gamma prior on the (semi-)global scales. If the
  argument `global_grouping` specifies e.g. `k` groups, then `NG_b` can
  be a numeric vector of length `k` and the elements determine the shape
  parameter in each group. `R2D2_b` has only to be specified if
  `prior="R2D2"`.

- R2D2_tol:

  Minimum number that a parameter draw of one of the shrinking
  parameters of the R2D2 prior can take. Prevents numerical issues that
  can appear when strong shrinkage is enforced if chosen to be greater
  than zero. `R2D2_tol` has only to be specified if `prior="R2D2"`.

- NG_a:

  (Single) positive real number. The value is interpreted as the
  concentration parameter for the local scales. Smaller values enforce
  heavier shrinkage. If the argument `global_grouping` specifies e.g.
  `k` groups, then `NG_a` can be a numeric vector of length `k` and the
  elements indicate the shrinkage in each group. A matrix of dimension
  `c(s,2)` specifies a discrete hyperprior, where the first column
  contains `s` support points and the second column contains the
  associated prior probabilities. `NG_a` has only to be specified if
  `prior="NG"`.

- NG_b:

  (Single) positive real number. The value indicates the shape parameter
  of the inverse gamma prior on the (semi-)global scales. If the
  argument `global_grouping` specifies e.g. `k` groups, then `NG_b` can
  be a numeric vector of length `k` and the elements determine the shape
  parameter in each group. `NG_b` has only to be specified if
  `prior="NG"`.

- NG_c:

  (Single) positive real number. The value indicates the scale parameter
  of the inverse gamma prior on the (semi-)global scales. If the
  argument `global_grouping` specifies e.g. `k` groups, then `NG_c` can
  be a numeric vector of length `k` and the elements determine the scale
  parameter in each group. Expert option would be to set the scale
  parameter proportional to `NG_a`. E.g. in the case where a discrete
  hyperprior for `NG_a` is chosen, a desired proportion of let's say
  `0.2` is achieved by setting `NG_c="0.2*a"` (character input!). `NG_c`
  has only to be specified if `prior="NG"`.

- NG_tol:

  Minimum number that a parameter draw of one of the shrinking
  parameters of the normal-gamma prior can take. Prevents numerical
  issues that can appear when strong shrinkage is enforced if chosen to
  be greater than zero. `NG_tol` has only to be specified if
  `prior="NG"`.

- SSVS_c0:

  single positive number indicating the (unscaled) standard deviation of
  the spike component. `SSVS_c0` has only to be specified if
  `prior="SSVS"`. It should be that \\SSVS\_{c0} \ll SSVS\_{c1}\\!
  `SSVS_c0` has only to be specified if `prior="SSVS"`.

- SSVS_c1:

  single positive number indicating the (unscaled) standard deviation of
  the slab component. `SSVS_c0` has only to be specified if
  `prior="SSVS"`. It should be that \\SSVS\_{c0} \ll SSVS\_{c1}\\!

- SSVS_semiautomatic:

  logical. If `SSVS_semiautomatic=TRUE` both `SSVS_c0` and `SSVS_c1`
  will be scaled by the variances of the posterior of PHI under a FLAT
  conjugate (dependent Normal-Wishart prior). `SSVS_semiautomatic` has
  only to be specified if `prior="SSVS"`.

- SSVS_p:

  Either a single positive number in the range `(0,1)` indicating the
  (fixed) prior inclusion probability of each coefficient. Or numeric
  vector of length 2 with positive entries indicating the shape
  parameters of the Beta distribution. In that case a Beta hyperprior is
  placed on the prior inclusion probability. `SSVS_p` has only to be
  specified if `prior="SSVS"`.

- HMP_lambda1:

  numeric vector of length 2. Both entries must be positive. The first
  indicates the shape and the second the rate of the Gamma hyperprior on
  own-lag coefficients. `HMP_lambda1` has only to be specified if
  `prior="HMP"`.

- HMP_lambda2:

  numeric vector of length 2. Both entries must be positive. The first
  indicates the shape and the second the rate of the Gamma hyperprior on
  cross-lag coefficients. `HMP_lambda2` has only to be specified if
  `prior="HMP"`.

- normal_sds:

  numeric vector of length \\n\\, where \\n = lags M^2\\ is the number
  of all VAR coefficients (excluding the intercept), indicating the
  prior variances. A single number will be recycled accordingly! Must be
  positive. `normal_sds` has only to be specified if `prior="normal"`.

- global_grouping:

  One of `"global"`, `"equation-wise"`, `"covariate-wise"`,
  `"olcl-lagwise"` `"fol"` indicating the sub-groups of the
  semi-global(-local) modifications to HS, R2D2, NG, DL and SSVS prior.
  Works also with user-specified indicator matrix of dimension
  `c(lags*M,M)`. Only relevant if `prior="HS"`, `prior="DL"`,
  `prior="R2D2"`, `prior="NG"` or `prior="SSVS"`.

- ...:

  Do not use!

## Value

A `baysianVARs_prior_phi`-object.

## Details

For details concerning prior-elicitation for VARs please see Gruber &
Kastner (2025).

Currently one can choose between six hierarchical shrinkage priors and a
normal prior: `prior="HS"` stands for the Horseshoe-prior, `prior="R2D2`
for the R\\^2\\-induced-Dirichlet-decompostion-prior, `prior="NG"` for
the normal-gamma-prior, `prior="DL"` for the Dirichlet-Laplace-prior,
`prior="SSVS"` for the stochastic-search-variable-selection-prior,
`prior="HMP"` for the semi-hierarchical Minnesota prior and
`prior=normal` for the normal-prior.

Semi-global shrinkage, i.e. group-specific shrinkage for pre-specified
subgroups of the coefficients, can be achieved through the argument
`global_grouping`.

## References

Gruber, L. and Kastner, G. (2025). Forecasting macroeconomic data with
Bayesian VARs: Sparse or dense? It depends! *International Journal of
Forecasting*.
[doi:10.1016/j.ijforecast.2025.02.001](https://doi.org/10.1016/j.ijforecast.2025.02.001)
.

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Horseshoe prior for a VAR(2)
phi_hs <- specify_prior_phi(data = data, lags = 2L ,prior = "HS")

# Semi-global-local Horseshoe prior for a VAR(2) with semi-global shrinkage parameters for
# cross-lag and own-lag coefficients in each lag
phi_hs_sg <- specify_prior_phi(data = data, lags = 2L, prior = "HS",
global_grouping = "olcl-lagwise")

# Semi-global-local Horseshoe prior for a VAR(2) with equation-wise shrinkage
# construct indicator matrix for equation-wise shrinkage
semi_global_mat <- matrix(1:ncol(data), 2*ncol(data), ncol(data),
byrow = TRUE)
phi_hs_ew <- specify_prior_phi(data = data, lags = 2L, prior = "HS",
global_grouping = semi_global_mat)
# (for equation-wise shrinkage one can also use 'global_grouping = "equation-wise"')

# \donttest{
# Estimate model with your prior configuration of choice
mod <- bvar(data, lags = 2L, prior_phi = phi_hs_sg, quiet = TRUE)
# }
```
