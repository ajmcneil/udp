# udp

<!-- badges: start -->
[![R-CMD-check](https://github.com/ajmcneil/udp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ajmcneil/udp/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`udp` provides tools for constructing and applying **uniform distribution
preserving** transformations: maps of the unit interval (and unit hypercube)
to itself that leave the uniform distribution invariant. It also includes
utilities for checking the distribution-preserving property.

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("ajmcneil/udp")
```

## Example

A uniform distribution preserving (udp) transformation folds the unit interval
onto itself, so `T(U)` stays uniform when `U` is. The map is many-to-one:
`udpsi()` inverts it stochastically, picking one pre-image at random in the
proportion that keeps the uniform property running backwards.

``` r
library(udp)

set.seed(1)
U <- runif(1000)

T4 <- udplegendre(5)      # a degree-5 shifted-Legendre udp function
plot(T4)

V <- udptrans(T4, U)      # still uniform
U2 <- udpsi(T4, V)        # a random pre-image of V, also uniform

pcoincide(T4)             # P(udpsi lands back on the original U)
#> [1] 0.3390592
```

Shuffles are piecewise-linear udp bijections built from a permutation and a
vector of signs. `aceshuffle()` searches for a shuffle of each margin of a
bivariate sample that maximises the linear correlation of the transformed
pair, which can expose dependence that is invisible to the ordinary
correlation:

``` r
GC <- copula::gumbelCopula(3)
W <- copula::rCopula(1000, GC)
X <- cbind(udpsi(T4, W[, 1]), udpsi(udpcosine(6), W[, 2]))

cor(X)[1, 2]              # near zero
#> [1] -0.01131886

fit <- aceshuffle(X[, 1], X[, 2], m = 50)
fit$correlation          # recovered
#> [1] 0.7974555
```

See `vignette("udp")` for the full tour.

## References

- McNeil, A. J. (2021). Modelling volatile time series with v-transforms and
  copulas. *Risks*, **9**(1), 14.
  <https://doi.org/10.3390/risks9010014>
- McNeil, A. J., Nešlehová, J. G. and Smith, A. D. (2025). Measures and models
  of non-monotonic dependence. arXiv:2512.10828.
  <https://arxiv.org/abs/2512.10828>

## License

MIT © Alexander J. McNeil
