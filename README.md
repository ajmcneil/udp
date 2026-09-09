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

``` r
library(udp)

# TODO: add a short worked example
```

## References

- McNeil, A. J. (2021). Modelling volatile time series with v-transforms and
  copulas. *Risks*, **9**(1), 14.
  <https://doi.org/10.3390/risks9010014>
- McNeil, A. J., Nešlehová, J. G. and Smith, A. D. (2025). Measures and models
  of non-monotonic dependence. arXiv:2512.10828.
  <https://arxiv.org/abs/2512.10828>

## License

MIT © Alexander J. McNeil
