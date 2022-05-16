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

# Estimate VAR(2) with stochastic volatitlity
mod <- bvar(Yraw = Y_est, p = 2, draws = 5000, burnin = 1000,
            priorPHI = priorPHI, priorL = priorL, SV = TRUE, progressbar = TRUE)

# Posterior summary of PHI
summary(mod)
```

    ## 
    ## Posterior median of PHI:
    ##                  GDPC1   CPIAUCSL   FEDFUNDS
<<<<<<< HEAD
    ## GDPC1.l1     2.600e-22  5.664e-14  1.146e-04
    ## CPIAUCSL.l1 -3.823e-12  7.143e-01  7.155e-03
    ## FEDFUNDS.l1 -2.986e-07  2.230e-02  1.265e+00
    ## GDPC1.l2     2.325e-10 -1.249e-16  1.223e-29
    ## CPIAUCSL.l2 -1.625e-09  6.716e-02  3.589e-11
    ## FEDFUNDS.l2 -7.502e-02 -2.066e-07 -3.139e-01
    ## intercept    1.646e-02  7.708e-04  1.305e-03
||||||| d488290
    ## GDPC1.l1     1.638e-11  8.360e-18  6.678e-09
    ## CPIAUCSL.l1 -1.109e-10  6.907e-01  1.478e-07
    ## FEDFUNDS.l1 -2.713e-16  4.066e-02  1.287e+00
    ## GDPC1.l2     4.869e-16 -1.247e-13  4.574e-13
    ## CPIAUCSL.l2 -5.187e-24  7.750e-02  3.442e-12
    ## FEDFUNDS.l2 -1.011e-01 -1.954e-05 -3.341e-01
    ## intercept    1.651e-02  7.356e-04  1.500e-03
=======
    ## GDPC1.l1     4.705e-21  2.092e-12  5.331e-11
    ## CPIAUCSL.l1 -3.209e-07  7.077e-01  1.142e-01
    ## FEDFUNDS.l1 -1.659e-06  2.936e-02  1.270e+00
    ## GDPC1.l2     1.105e-11 -6.636e-22  2.031e-13
    ## CPIAUCSL.l2 -1.006e-08  9.085e-02  5.309e-11
    ## FEDFUNDS.l2 -2.987e-02 -6.341e-07 -3.263e-01
    ## intercept    1.611e-02  8.160e-04  1.675e-03
>>>>>>> 7f7e9003da7aa2d47be741346ed6d184dd924ee3
    ## 
    ## Posterior interquartile range of PHI:
    ##                 GDPC1  CPIAUCSL  FEDFUNDS
<<<<<<< HEAD
    ## GDPC1.l1    3.123e-07 0.0001682 0.0720803
    ## CPIAUCSL.l1 8.327e-03 0.2364082 0.2927033
    ## FEDFUNDS.l1 3.478e-02 0.0957405 0.1961262
    ## GDPC1.l2    5.686e-04 0.0001618 0.0001116
    ## CPIAUCSL.l2 3.252e-02 0.2697496 0.0002177
    ## FEDFUNDS.l2 1.209e-01 0.0614795 0.1823603
    ## intercept   3.168e-03 0.0012899 0.0022168
||||||| d488290
    ## GDPC1.l1    4.550e-04 3.274e-04 3.012e-02
    ## CPIAUCSL.l1 2.778e-03 2.501e-01 2.314e-01
    ## FEDFUNDS.l1 2.536e-03 1.228e-01 1.754e-01
    ## GDPC1.l2    6.069e-05 7.844e-06 7.049e-05
    ## CPIAUCSL.l2 2.126e-07 2.925e-01 1.311e-03
    ## FEDFUNDS.l2 1.155e-01 9.568e-02 1.730e-01
    ## intercept   3.159e-03 1.319e-03 2.088e-03
=======
    ## GDPC1.l1    7.218e-05 5.237e-04 0.0179714
    ## CPIAUCSL.l1 2.067e-01 2.431e-01 0.3371642
    ## FEDFUNDS.l1 2.860e-02 1.103e-01 0.1907421
    ## GDPC1.l2    3.455e-04 3.063e-05 0.0005711
    ## CPIAUCSL.l2 2.947e-02 2.785e-01 0.0011570
    ## FEDFUNDS.l2 1.118e-01 8.357e-02 0.1812153
    ## intercept   3.283e-03 1.242e-03 0.0021808
>>>>>>> 7f7e9003da7aa2d47be741346ed6d184dd924ee3

``` r
# Traceplot of global shrinkage parameter
ts.plot(mod$phi_hyperparameter$zeta1)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

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

![](README_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
# Summary of predictive evaluation
summary(pred)
```

    ## 
    ## LPL:
    ##   t+1   t+2   t+3   t+4 
<<<<<<< HEAD
    ## 9.307 9.306 9.336 8.940 
||||||| d488290
    ## 9.205 9.187 9.279 8.829 
=======
    ## 9.437 9.385 9.428 9.110 
>>>>>>> 7f7e9003da7aa2d47be741346ed6d184dd924ee3
    ## 
    ## Marginal joint LPL of CPIAUCSL & FEDFUNDS:
    ##   t+1   t+2   t+3   t+4 
<<<<<<< HEAD
    ## 5.690 5.784 5.728 5.708 
||||||| d488290
    ## 5.637 5.726 5.702 5.668 
=======
    ## 5.833 5.856 5.816 5.812 
>>>>>>> 7f7e9003da7aa2d47be741346ed6d184dd924ee3
    ## 
    ## Marginal univariate LPLs:
    ##     GDPC1 CPIAUCSL FEDFUNDS
<<<<<<< HEAD
    ## t+1 3.636    4.334    1.291
    ## t+2 3.537    4.109    1.554
    ## t+3 3.561    3.919    1.650
    ## t+4 3.091    3.570    1.880
||||||| d488290
    ## t+1 3.613    4.316    1.260
    ## t+2 3.496    4.095    1.507
    ## t+3 3.534    3.895    1.635
    ## t+4 3.036    3.535    1.835
=======
    ## t+1 3.647    4.344    1.457
    ## t+2 3.570    4.119    1.683
    ## t+3 3.588    3.936    1.778
    ## t+4 3.181    3.580    1.952
>>>>>>> 7f7e9003da7aa2d47be741346ed6d184dd924ee3
    ## 
    ## Prediction quantiles:
    ## , , GDPC1
    ## 
<<<<<<< HEAD
    ##           t+1       t+2       t+3       t+4
    ## 5%  -0.012521 -0.013541 -0.013760 -0.013461
    ## 50%  0.005351  0.004844  0.004739  0.004448
    ## 95%  0.023620  0.022066  0.022527  0.022722
||||||| d488290
    ##           t+1       t+2       t+3       t+4
    ## 5%  -0.011935 -0.013955 -0.014017 -0.014216
    ## 50%  0.005112  0.004069  0.003913  0.004034
    ## 95%  0.022006  0.021375  0.021346  0.022128
=======
    ##           t+1      t+2       t+3       t+4
    ## 5%  -0.011507 -0.01241 -0.012545 -0.012727
    ## 50%  0.006661  0.00553  0.005558  0.005395
    ## 95%  0.023956  0.02330  0.023269  0.023926
>>>>>>> 7f7e9003da7aa2d47be741346ed6d184dd924ee3
    ## 
    ## , , CPIAUCSL
    ## 
<<<<<<< HEAD
    ##          t+1        t+2       t+3       t+4
    ## 5%  0.001726 -0.0001635 -0.001052 -0.002234
    ## 50% 0.010477  0.0111110  0.011822  0.012336
    ## 95% 0.019874  0.0232649  0.025888  0.028110
||||||| d488290
    ##          t+1       t+2        t+3      t+4
    ## 5%  0.002013 0.0004063 -0.0009667 -0.00150
    ## 50% 0.010601 0.0115145  0.0122505  0.01274
    ## 95% 0.019727 0.0232591  0.0261385  0.02839
=======
    ##          t+1       t+2       t+3       t+4
    ## 5%  0.001736 0.0004717 -0.001208 -0.002108
    ## 50% 0.010499 0.0111807  0.011831  0.012239
    ## 95% 0.019816 0.0229956  0.025614  0.027502
>>>>>>> 7f7e9003da7aa2d47be741346ed6d184dd924ee3
    ## 
    ## , , FEDFUNDS
    ## 
<<<<<<< HEAD
    ##         t+1     t+2     t+3     t+4
    ## 5%  0.09626 0.08008 0.06766 0.05618
    ## 50% 0.11501 0.11409 0.11192 0.11050
    ## 95% 0.12911 0.13768 0.14495 0.15230
||||||| d488290
    ##         t+1     t+2     t+3     t+4
    ## 5%  0.09688 0.08249 0.06949 0.05796
    ## 50% 0.11483 0.11369 0.11152 0.10963
    ## 95% 0.12995 0.13943 0.14718 0.15575
=======
    ##         t+1    t+2    t+3     t+4
    ## 5%  0.09502 0.0777 0.0622 0.04942
    ## 50% 0.11431 0.1125 0.1102 0.10784
    ## 95% 0.13000 0.1398 0.1482 0.15611
>>>>>>> 7f7e9003da7aa2d47be741346ed6d184dd924ee3
