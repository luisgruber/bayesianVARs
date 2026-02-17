# Print method for summary.bayesianVARs_predict objects

Print method for `summary.bayesianVARs_predict` objects.

## Usage

``` r
# S3 method for class 'summary.bayesianVARs_predict'
print(x, ...)
```

## Arguments

- x:

  A `summary.bayesianVARs_predict` object obtained via
  [`summary.bayesianVARs_predict()`](https://luisgruber.github.io/bayesianVARs/reference/summary.bayesianVARs_predict.md).

- ...:

  Currently ignored.

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
#> 599/1000 stable posterior draws remaining for prediction!
sum <- summary(predictions)
print(sum)
#> 
#> Prediction quantiles:
#> , , GDPC1
#> 
#>          t+1
#> 5%  -0.09297
#> 50% -0.01958
#> 95%  0.05450
#> 
#> , , CPIAUCSL
#> 
#>           t+1
#> 5%  -0.019524
#> 50% -0.006549
#> 95%  0.004611
#> 
#> , , FEDFUNDS
#> 
#>           t+1
#> 5%  -0.020377
#> 50% -0.003513
#> 95%  0.013922
#> 
```
