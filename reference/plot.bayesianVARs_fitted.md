# Visualization of in-sample fit of an estimated VAR.

Visualization of in-sample fit of an estimated VAR.

## Usage

``` r
# S3 method for class 'bayesianVARs_fitted'
plot(
  x,
  dates = NULL,
  vars = "all",
  quantiles = c(0.05, 0.5, 0.95),
  n_col = 1L,
  ...
)
```

## Arguments

- x:

  A `bayesianVARs_fitted` object.

- dates:

  optional vector of dates for labelling the x-axis. The default values
  is `NULL`; in this case, the axis will be labeled with numbers.

- vars:

  character vector containing the names of the variables to be
  visualized. The default is `"all"` indicating that the fit of all
  variables is visualized.

- quantiles:

  numeric vector indicating which quantiles to plot.

- n_col:

  integer indicating the number of columns to use for plotting.

- ...:

  Currently ignored.

## Value

returns `x` invisibly

## See also

- fitted method for class 'bayesianVARs_bvar':
  [`fitted.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/fitted.bayesianVARs_bvar.md).

- Other plotting
  [`plot.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_bvar.md),
  `plot.bayesianVARs_fitted()`,
  [`plot.bayesianVARs_predict()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_predict.md),
  [`pairs.bayesianVARs_predict()`](https://luisgruber.github.io/bayesianVARs/reference/pairs_predict.md),
  [`posterior_heatmap()`](https://luisgruber.github.io/bayesianVARs/reference/posterior_heatmap.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Simulate predicted historical values including the error term.
pred <- fitted(mod, error_term = TRUE)

# Visualize
plot(pred)
```
