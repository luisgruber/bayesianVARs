# Summary method for bayesianVARs_predict objects

Summary method for `bayesianVARs_predict` objects.

## Usage

``` r
# S3 method for class 'bayesianVARs_predict'
summary(object, ...)
```

## Arguments

- object:

  A `bayesianVARs_predict` object obtained via
  [`predict.bayesianVARs_bvar()`](https://luisgruber.github.io/bayesianVARs/reference/predict.bayesianVARs_bvar.md).

- ...:

  Currently ignored!

## Value

A `summary.bayesianVARs_predict` object.

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
#> Using the remaining 627 stable posterior draws for prediction!
summary(predictions)
#> 
#> Prediction quantiles:
#> , , GDPC1
#> 
#>          t+1
#> 5%  -0.07407
#> 50% -0.01820
#> 95%  0.02693
#> 
#> , , CPIAUCSL
#> 
#>           t+1
#> 5%  -0.018565
#> 50% -0.007542
#> 95%  0.003227
#> 
#> , , FEDFUNDS
#> 
#>           t+1
#> 5%  -0.021216
#> 50% -0.002924
#> 95%  0.016862
#> 
```
