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
#> 558/1000 stable posterior draws remaining for prediction!
sum <- summary(predictions)
print(sum)
#> 
#> Prediction quantiles:
#> , , GDPC1
#> 
#>          t+1
#> 5%  -0.06715
#> 50% -0.02025
#> 95%  0.03612
#> 
#> , , CPIAUCSL
#> 
#>           t+1
#> 5%  -0.019310
#> 50% -0.007503
#> 95%  0.003487
#> 
#> , , FEDFUNDS
#> 
#>           t+1
#> 5%  -0.021193
#> 50% -0.002992
#> 95%  0.018641
#> 
```
