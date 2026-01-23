# Impulse Responses Plot

Visualization of the impulse responses. Responses are plotted on a grid,
where rows correspond to variables and columns correspond to shocks.

## Usage

``` r
# S3 method for class 'bayesianVARs_irf'
plot(
  x,
  vars = "all",
  quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
  default_hair_color = "#FF000003",
  true_irf = NULL,
  ...
)
```

## Arguments

- x:

  An object of type `bayesianVARs_irf` obtained via
  [`irf`](https://luisgruber.github.io/bayesianVARs/reference/irf.md).

- vars:

  character vector containing the names of the variables to be
  visualized. The default is `"all"` indicating that all variables are
  visualized.

- quantiles:

  numeric vector indicating which quantiles to plot. If `hairy=TRUE` was
  specified when calling
  [`irf`](https://luisgruber.github.io/bayesianVARs/reference/irf.md), a
  proportion of `max(quantiles)` IRFs will be plotted. Specify `0` to
  plot a point-estimate only. If `hairy=FALSE` was specified (the
  default), point-wise quantiles will be plotted. Note that the curve of
  point-wise medians is not necessarily in the set of IRFs (see Inoue
  2022).

- default_hair_color:

  the color of the IRF samples, if `hairy=TRUE` was specified.

- true_irf:

  If the true IRFs are known (because the data was simulated) they can
  be plotted alongside the estimates, such that the quality of the
  estimates may be judged. `true_irf` should be a numeric array with
  dimensions variables, shocks and time, in that order.

- ...:

  Currently ignored!

## References

Inoue, A. and Kilian, L. (2022). Joint Bayesian inference about impulse
responses in VAR models. *Journal of Econometrics*,
[doi:10.1016/j.jeconom.2021.05.010](https://doi.org/10.1016/j.jeconom.2021.05.010)
.

## See also

[`irf`](https://luisgruber.github.io/bayesianVARs/reference/irf.md)

## Author

Stefan Haan <sthaan@edu.aau.at>
