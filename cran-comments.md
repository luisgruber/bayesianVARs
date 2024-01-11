## Resubmission
This is a resubmission. As pointed out by Benjamin Altmann in this version I have:

* either removed \dontrun{} or replaced \dontrun{} with \donttest{} in examples.
* checked that in examples, vignettes and demos users's options() are always reset, e.g.: -> inst/doc/bayesianVARs-vignette.R
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

## Resubmission
This is a resubmission. In this version I have:

* Removed $(SHLIB_OPENMP_CXXFLAGS) in the makevars files.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* IMHO all possibly misspelled words in DESCRIPTION are correctly spelled.
