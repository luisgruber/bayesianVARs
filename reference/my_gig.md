# Draw from generalized inverse Gaussian

Vectorized version of [`rgig`](https://rdrr.io/pkg/GIGrvg/man/rgig.html)

## Usage

``` r
my_gig(n, lambda, chi, psi)
```

## Arguments

- n:

  A single integer indicating the number of draws to generate.

- lambda:

  vector of shape parameters.

- chi:

  vector of shape/scale parameters. Must be nonnegative for positive
  lambdas and positive else.

- psi:

  vector of shape/scale parameters. Must be nonnegative for negative
  lambdas and positive else.

## Value

Matrix of dimension `c(n,m)`, where `m` is the maximum length of
`lambda`, `psi` and `chi`.

## Examples

``` r
gigsamples <- my_gig(2, c(1,1), c(1,1), c(1,1))
```
