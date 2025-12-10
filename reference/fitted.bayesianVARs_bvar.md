# Simulate fitted/predicted historical values for an estimated VAR model

Simulates the fitted/predicted (in-sample) values for an estimated VAR
model.

## Usage

``` r
# S3 method for class 'bayesianVARs_bvar'
fitted(object, error_term = TRUE, ...)
```

## Arguments

- object:

  A `bayesianVARs_bvar` object estimated via
  [`bvar()`](https://luisgruber.github.io/bayesianVARs/reference/bvar.md).

- error_term:

  logical indicating whether to include the error term or not.

- ...:

  Currently ignored.

## Value

An object of class `bayesianVARs_fitted`.

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Estimate a model
mod <- bvar(data, sv_keep = "all", quiet = TRUE)

# Simulate predicted historical values including the error term.
pred <- fitted(mod, error_term = TRUE)

# Simulate fitted historical values not including the error term.
fit <- fitted(mod, error_term = FALSE)

# Visualize
plot(pred)

plot(fit)
```
