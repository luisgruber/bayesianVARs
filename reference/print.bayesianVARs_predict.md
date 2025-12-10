# Print method for bayesianVARs_predict objects

Print method for bayesianVARs_predict objects.

## Usage

``` r
# S3 method for class 'bayesianVARs_predict'
print(x, ...)
```

## Arguments

- x:

  A `bayesianVARs_predict` object obtained via
  [`predict.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/predict.bayesianVARs_bvar.md).

- ...:

  Currently ignored!

## Value

Returns `x` invisibly.

## Examples

``` r
# Access a subset of the usmacro_growth dataset
data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

# Split data in train and test
train <- data[1:(nrow(data)-4),]
test <- data[-c(1:(nrow(data)-4)),]

# Estimate model using train data only
mod <- bvar(train, quiet = TRUE)

# Simulate from 1-step ahead posterior predictive
predictions <- predict(mod, ahead = 1L)
#> 'stable=TRUE': Calling 'stable_bvar()' to discard those posterior
#>           draws that do not fulfill the stable criterion.
#> 
#> 643/1000 stable posterior draws remaining for prediction!
print(predictions)
#> 
#> Generic functions for bayesianVARs_predict objects:
#>   - summary.bayesianVARs_predict(),
#>   - pairs.bayesianVARs_predict(),
#>   - plot.bayesianVARs_predict() (alias for pairs.bayesianVARs_predict()).
```
