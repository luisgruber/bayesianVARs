# bayesianVARs 0.1.2

* Added minimum version to factorstochvol in the Imports field of the DESCRIPTION file in order to avoid unnecessary building errors. Thanks to Sergey Fedorov for pointing this out.
* vcov.bayesianVARs_bvar method now can be specified for specific time-points.
* bugfix in cpp function which constructs variance-covariance matrices. If a Cholesky structure for the errors had been specified, exported functions such as vcov, predict and fitted were affected.

# bayesianVARs 0.1.1

* Fixed clang-UBSAN issue.
* Fixed undefined figure references in vignette.

# bayesianVARs 0.1.0

* Initial CRAN submission.
