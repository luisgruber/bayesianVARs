bayesianVARs: Hierarchical shrinkage priors
================
Luis Gruber
16 5 2022

Estimation of Bayesian vectorautoregressions. Implements several modern
hierarchical shrinkage priors, amongst them
![R^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R%5E2 "R^2")-induced-Dirichlet-decomposition
(R2D2) prior, Dirichlet-Laplace (DL) prior, Stochastic Search Variable
Selection (SSVS) and the Hierarchical Minnesota prior.

# Installation

Install directly from GitHub.

``` r
devtools::install_github("luisgruber/bayesianVARs")
```

# Getting started

``` r
# load package
library(bayesianVARs)

# Some data (from FRED-MD database)
data <- dat_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Data for estimation
Y_est <- data[1:100,]

# Specify prior for reduced-form VAR coefficients (with default settings)
priorPHI <- specify_priorPHI(prior = "R2D2")

# Specify prior for L (Decomposition of variance-covariance matrix in the form of t(L^(-1))%*%D_t%*%L^(-1), where L is upper triangular)
priorL <- specify_priorL(prior = "DL")

# Estimate VAR(2) with stochastic volatitlity
mod <- bvar(Yraw = Y_est, p = 2, draws = 5000, burnin = 1000,
            priorPHI = priorPHI, priorL = priorL, SV = TRUE, progressbar = TRUE)

# Posterior summary of PHI
summary(mod)
```

    ## 
    ## Posterior median of PHI:
    ##                  GDPC1   CPIAUCSL   FEDFUNDS
    ## GDPC1.l1     1.638e-11  8.360e-18  6.678e-09
    ## CPIAUCSL.l1 -1.109e-10  6.907e-01  1.478e-07
    ## FEDFUNDS.l1 -2.713e-16  4.066e-02  1.287e+00
    ## GDPC1.l2     4.869e-16 -1.247e-13  4.574e-13
    ## CPIAUCSL.l2 -5.187e-24  7.750e-02  3.442e-12
    ## FEDFUNDS.l2 -1.011e-01 -1.954e-05 -3.341e-01
    ## intercept    1.651e-02  7.356e-04  1.500e-03
    ## 
    ## Posterior interquartile range of PHI:
    ##                 GDPC1  CPIAUCSL  FEDFUNDS
    ## GDPC1.l1    4.550e-04 3.274e-04 3.012e-02
    ## CPIAUCSL.l1 2.778e-03 2.501e-01 2.314e-01
    ## FEDFUNDS.l1 2.536e-03 1.228e-01 1.754e-01
    ## GDPC1.l2    6.069e-05 7.844e-06 7.049e-05
    ## CPIAUCSL.l2 2.126e-07 2.925e-01 1.311e-03
    ## FEDFUNDS.l2 1.155e-01 9.568e-02 1.730e-01
    ## intercept   3.159e-03 1.319e-03 2.088e-03

``` r
# Traceplot of global shrinkage parameter
ts.plot(mod$phi_hyperparameter$zeta1)
```

![](read_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Simulate from predictive density and compare to ex-post realized value by means
# log predictive likelihood
Y_obs <- data[101:104,]
pred <- predict(mod, nsteps = 4, LPL = TRUE, Y_obs = Y_obs ,LPL_VoI = c("CPIAUCSL", "FEDFUNDS"))

# Histograms of predictive densities
par(mfrow=c(4,3))
for (i in paste0("t+",1:4)) {
  for(j in c("GDPC1" ,"CPIAUCSL", "FEDFUNDS")){
  hist(pred$predictions[,i,j], main = paste0(i,": ",j), xlab = "")
}
}
```

![](read_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
# Summary of predictive evaluation
summary(pred)
```

    ## 
    ## LPL:
    ##   t+1   t+2   t+3   t+4 
    ## 9.205 9.187 9.279 8.829 
    ## 
    ## Marginal joint LPL of CPIAUCSL & FEDFUNDS:
    ##   t+1   t+2   t+3   t+4 
    ## 5.637 5.726 5.702 5.668 
    ## 
    ## Marginal univariate LPLs:
    ##     GDPC1 CPIAUCSL FEDFUNDS
    ## t+1 3.613    4.316    1.260
    ## t+2 3.496    4.095    1.507
    ## t+3 3.534    3.895    1.635
    ## t+4 3.036    3.535    1.835
    ## 
    ## Prediction quantiles:
    ## , , GDPC1
    ## 
    ##           t+1       t+2       t+3       t+4
    ## 5%  -0.011935 -0.013955 -0.014017 -0.014216
    ## 50%  0.005112  0.004069  0.003913  0.004034
    ## 95%  0.022006  0.021375  0.021346  0.022128
    ## 
    ## , , CPIAUCSL
    ## 
    ##          t+1       t+2        t+3      t+4
    ## 5%  0.002013 0.0004063 -0.0009667 -0.00150
    ## 50% 0.010601 0.0115145  0.0122505  0.01274
    ## 95% 0.019727 0.0232591  0.0261385  0.02839
    ## 
    ## , , FEDFUNDS
    ## 
    ##         t+1     t+2     t+3     t+4
    ## 5%  0.09688 0.08249 0.06949 0.05796
    ## 50% 0.11483 0.11369 0.11152 0.10963
    ## 95% 0.12995 0.13943 0.14718 0.15575
