# Visualization of the residuals of an estimated VAR.

Visualization of the residuals of an estimated VAR.

## Usage

``` r
# S3 method for class 'bayesianVARs_residuals'
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

  A `bayesianVARs_residuals` object.

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

- residuals method for class 'bayesianVARs_bvar':
  [`residuals.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/residuals.bayesianVARs_bvar.md).

- Other plotting
  [`plot.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_bvar.md),
  [`plot.bayesianVARs_fitted()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_fitted.md),
  [`plot.bayesianVARs_predict()`](https://luisgruber.github.io/bayesianVARs/reference/plot.bayesianVARs_predict.md),
  [`pairs.bayesianVARs_predict()`](https://luisgruber.github.io/bayesianVARs/reference/pairs_predict.md),
  [`posterior_heatmap()`](https://luisgruber.github.io/bayesianVARs/reference/posterior_heatmap.md).

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

mod.resids <- residuals(mod)

# Visualize
plot(mod.resids)
```
