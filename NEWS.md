# bayesianVARs 0.1.5

* Bug fix which makes estimation possible of a VAR with factor stochastic volatility structure on the errors with the restriction that the loadings matrix has zeros above the diagonal. Thanks to Stefan Haan for reporting the bug.

# bayesianVARs 0.1.4

* For consistency with other functions, from now on prior_intercept in bvar() specifies standard deviations instead of variances.
* Bug fix concerning 'additional check' valgrind. Seems to pass the check now. Thanks to Brian Ripley for reporting the bug.

# bayesianVARs 0.1.3

* bugfix concerning VAR with factor structure on errors with homoscedastic factors.

# bayesianVARs 0.1.2

* Added minimum version to factorstochvol in the Imports field of the DESCRIPTION file in order to avoid unnecessary building errors. Thanks to Sergey Fedorov for pointing this out.
* vcov.bayesianVARs_bvar method now can be specified for specific time-points.
* bugfix in cpp function which constructs variance-covariance matrices. If a Cholesky structure for the errors had been specified, exported functions such as vcov, predict and fitted were affected.

# bayesianVARs 0.1.1

* Fixed clang-UBSAN issue.
* Fixed undefined figure references in vignette.

# bayesianVARs 0.1.0

* Initial CRAN submission.
