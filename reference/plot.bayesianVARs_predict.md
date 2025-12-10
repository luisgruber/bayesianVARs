# Fan chart

Visualization of (out-of-sample) predictive distribution.

## Usage

``` r
# S3 method for class 'bayesianVARs_predict'
plot(
  x,
  dates = NULL,
  vars = "all",
  ahead = NULL,
  quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
  n_col = 1L,
  first_obs = 1L,
  ...
)
```

## Arguments

- x:

  An object of type `bayesianVARs_predict` obtained via
  [`predict.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/predict.bayesianVARs_bvar.md).

- dates:

  optional vector of dates for labeling the x-axis. The default values
  is `NULL`; in this case, the axis will be labeled with numbers.

- vars:

  character vector containing the names of the variables to be
  visualized. The default is `"all"` indicating that all variables are
  visualized.

- ahead:

  Integer vector (or coercible to such) indicating which step ahead to
  plot. `max(ahead)` must be smaller equal to `dim(x$predictions)[1]`.

- quantiles:

  numeric vector indicating which quantiles to plot.

- n_col:

  integer indicating the number of columns to use for plotting.

- first_obs:

  integer indicating the first observation to be used for plotting.

- ...:

  Currently ignored!

## Value

Returns `x` invisibly!

## See also

Other plotting
[`plot.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_bvar.md),
[`plot.bayesianVARs_fitted()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_fitted.md),
[`pairs.bayesianVARs_predict()`](https://luisgruber.github.io/bayesianVARs/reference/pairs_predict.md)
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
#> 477/1000 stable posterior draws remaining for prediction!

# Visualize
plot(predictions, vars = 1:3, ahead = 1:3)
```
