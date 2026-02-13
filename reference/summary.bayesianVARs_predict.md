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
#> 
#> 556/1000 stable posterior draws remaining for prediction!
summary(predictions)
#> 
#> Prediction quantiles:
#> , , GDPC1
#> 
#>          t+1
#> 5%  -0.06684
#> 50% -0.02114
#> 95%  0.02705
#> 
#> , , CPIAUCSL
#> 
#>           t+1
#> 5%  -0.020639
#> 50% -0.008359
#> 95%  0.002530
#> 
#> , , FEDFUNDS
#> 
#>           t+1
#> 5%  -0.019689
#> 50% -0.002872
#> 95%  0.019201
#> 
```
