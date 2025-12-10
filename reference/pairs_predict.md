# Pairwise visualization of out-of-sample posterior predictive densities.

Pairwise visualization of out-of-sample posterior predictive densities.

## Usage

``` r
# S3 method for class 'bayesianVARs_predict'
pairs(x, vars, ahead, ...)
```

## Arguments

- x:

  An object of class `bayesianVARs_predict` obtained via
  [`predict.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/predict.bayesianVARs_bvar.md).

- vars:

  Integer vector (or coercible to such) indicating which variables to
  plot.

- ahead:

  Integer vector (or coercible to such) indicating which step ahead to
  plot. `max(ahead)` must be smaller equal to `dim(x$predictions)[1]`.

- ...:

  Currently ignored!

## Value

Returns `x` invisibly.

## Note

Note that that `bayesianVARs_predict` can also be used withing
[`plot.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_bvar.md).

## See also

Other plotting
[`plot.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_bvar.md),
[`plot.bayesianVARs_fitted()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_fitted.md),
[`plot.bayesianVARs_predict()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_predict.md)
[`posterior_heatmap()`](https://luisgruber.github.io/bayesianVARs/reference/posterior_heatmap.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Simulate from posterior predictive
predictions <- predict(mod, ahead = 1:3)
#> 'stable=TRUE': Calling 'stable_bvar()' to discard those posterior
#>           draws that do not fulfill the stable criterion.
#> 
#> 522/1000 stable posterior draws remaining for prediction!

# Visualize
pairs(predictions, vars = 1:3, ahead = 1:3)


```
