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

Data used in the following example is from Michael W. McCracken and
Serena Ng, “FRED-QD: A Quarterly Database for Macroeconomic Research,”
Federal Reserve Bank of St. Louis Review, First Quarter 2021, pp. 1-44.
<https://doi.org/10.20955/r.103.1-44>.

``` r
set.seed(537)
# load package
library(bayesianVARs)

# Some data (from FRED-MD database)
data <- dat_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Data for estimation
Y_est <- data[1:100,]

# Specify prior for reduced-form VAR coefficients (with default settings)
prior <- "R2D2" # or "DL", "SSVS", "HMP", "normal"
priorPHI <- specify_priorPHI(prior = "R2D2")

# Specify prior for L (Decomposition of variance-covariance matrix in the form of t(L^(-1))%*%D_t%*%L^(-1), where L is upper triangular)
priorL <- specify_priorL(prior = "DL")

# Estimate VAR(2) with stochastic volatility
mod <- bvar(Yraw = Y_est, p = 2, draws = 5000, burnin = 1000,
            priorPHI = priorPHI, priorL = priorL, SV = TRUE, progressbar = TRUE)

# Posterior summary of PHI
summary(mod)
```

    ## 
    ## Posterior median of PHI:
    ##                  GDPC1   CPIAUCSL   FEDFUNDS
    ## GDPC1.l1     2.600e-22  5.664e-14  1.146e-04
    ## CPIAUCSL.l1 -3.823e-12  7.143e-01  7.155e-03
    ## FEDFUNDS.l1 -2.986e-07  2.230e-02  1.265e+00
    ## GDPC1.l2     2.325e-10 -1.249e-16  1.223e-29
    ## CPIAUCSL.l2 -1.625e-09  6.716e-02  3.589e-11
    ## FEDFUNDS.l2 -7.502e-02 -2.066e-07 -3.139e-01
    ## intercept    1.646e-02  7.708e-04  1.305e-03
    ## 
    ## Posterior interquartile range of PHI:
    ##                 GDPC1  CPIAUCSL  FEDFUNDS
    ## GDPC1.l1    3.123e-07 0.0001682 0.0720803
    ## CPIAUCSL.l1 8.327e-03 0.2364082 0.2927033
    ## FEDFUNDS.l1 3.478e-02 0.0957405 0.1961262
    ## GDPC1.l2    5.686e-04 0.0001618 0.0001116
    ## CPIAUCSL.l2 3.252e-02 0.2697496 0.0002177
    ## FEDFUNDS.l2 1.209e-01 0.0614795 0.1823603
    ## intercept   3.168e-03 0.0012899 0.0022168

``` r
# Traceplot of global shrinkage parameter
ts.plot(mod$phi_hyperparameter$zeta1)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Simulate from predictive density and compare to ex-post realized value by 
# means of log predictive likelihood

# get ex-post observed data
Y_obs <- data[101:104,]
# predict
pred <- predict(mod, nsteps = 4, LPL = TRUE, Y_obs = Y_obs ,LPL_VoI = c("CPIAUCSL", "FEDFUNDS"))

# Histograms of predictive densities
par(mfrow=c(4,3))
for (i in paste0("t+",1:4)) {
  for(j in c("GDPC1" ,"CPIAUCSL", "FEDFUNDS")){
  hist(pred$predictions[,i,j], main = paste0(i,": ",j), xlab = "")
}
}
```

![](README_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
# Summary of predictive evaluation
summary(pred)
```

    ## 
    ## LPL:
    ##   t+1   t+2   t+3   t+4 
    ## 9.307 9.306 9.336 8.940 
    ## 
    ## Marginal joint LPL of CPIAUCSL & FEDFUNDS:
    ##   t+1   t+2   t+3   t+4 
    ## 5.690 5.784 5.728 5.708 
    ## 
    ## Marginal univariate LPLs:
    ##     GDPC1 CPIAUCSL FEDFUNDS
    ## t+1 3.636    4.334    1.291
    ## t+2 3.537    4.109    1.554
    ## t+3 3.561    3.919    1.650
    ## t+4 3.091    3.570    1.880
    ## 
    ## Prediction quantiles:
    ## , , GDPC1
    ## 
    ##           t+1       t+2       t+3       t+4
    ## 5%  -0.012521 -0.013541 -0.013760 -0.013461
    ## 50%  0.005351  0.004844  0.004739  0.004448
    ## 95%  0.023620  0.022066  0.022527  0.022722
    ## 
    ## , , CPIAUCSL
    ## 
    ##          t+1        t+2       t+3       t+4
    ## 5%  0.001726 -0.0001635 -0.001052 -0.002234
    ## 50% 0.010477  0.0111110  0.011822  0.012336
    ## 95% 0.019874  0.0232649  0.025888  0.028110
    ## 
    ## , , FEDFUNDS
    ## 
    ##         t+1     t+2     t+3     t+4
    ## 5%  0.09626 0.08008 0.06766 0.05618
    ## 50% 0.11501 0.11409 0.11192 0.11050
    ## 95% 0.12911 0.13768 0.14495 0.15230
