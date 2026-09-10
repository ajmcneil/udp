# udp 0.1.1

* `aceshuffle()` now returns `V`, a two-column matrix of the shuffled data
  (`udptrans(shuffle1, U1)`, `udptrans(shuffle2, U2)`), in place of the echoed
  input.
* `aceshuffle()` resolves the joint reflection
  `(shuffle1, shuffle2) -> (1 - shuffle1, 1 - shuffle2)` toward the
  upward-trending orientation, measured by
  `E[shuffle(U) | U > 1/2] - E[shuffle(U) | U < 1/2]` summed over the pair.
* The transformation curve in every udp `plot()` method is drawn at line width
  1.5 (was 2).
* Documentation switched to Oxford (`-ize`) spelling.
* Added a worked example to the README and a package vignette,
  `vignette("udp")`.
* Added a pkgdown site at <https://ajmcneil.github.io/udp/>.

# udp 0.1.0

* First development release.
* `udp` virtual class with the `udptrans()`, `udpsi()`, `udpinverse()` and
  `pcoincide()` generics.
* udp transformation families: `shuffle()` (permutation-and-sign bijections of
  the unit interval), `udpcosine()` (triangle waves) and `udplegendre()`
  (shifted-Legendre compositions).
* v-transforms: `vlinear()`, `v2p()`, `v3p()`, `v2b()`, `v3b()`, with
  `vtransform()` methods for evaluation, gradient and numeric inversion.
* `aceshuffle()`: alternating conditional expectation search for the pair of
  shuffles that maximises the linear correlation of a bivariate sample.
* `plot()` methods for every class.
