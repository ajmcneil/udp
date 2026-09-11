# udp 0.1.1

* Added an internal (unexported) `udpbreakcdf()`: for a set of query points
  (typically `udpbreaks(x)`), the probability that `udpsi()` selects a
  pre-image at or below each one -- the running sum of `udpinverse()`'s
  selection probabilities over the pre-images not exceeding it.
* Added an internal (unexported) `udpbreaks()` generic returning the
  partition of `[0, 1]`, including `0` and `1`, on which a udp
  transformation is piecewise continuously differentiable. For
  `udplegendre` this is the turning points of the underlying
  shifted-Legendre polynomial together with the transversal pre-images
  of their critical values -- matches the known piece counts for degree
  2 to 6 (2, 5, 6, 13, 12).
* Added `udpderiv()`, the derivative `T'(u)` of a udp transformation, as a
  generic over every class. At the finitely many points where `T` is not
  differentiable it returns the left derivative (the right derivative at
  `u = 0`). For `udplegendre`, `T'` is computed from the same exact polynomial
  root-finding as `udpinverse()`; at a turning point of the underlying
  shifted-Legendre polynomial the one-sided slope is `2` or `-2` times the
  number of turning points sharing that critical value (mirror pairs under
  the symmetry of even degree contribute together).
* `aceshuffle()` now returns `V`, a two-column matrix of the shuffled data
  (`udptrans(shuffle1, U1)`, `udptrans(shuffle2, U2)`), in place of the echoed
  input.
* `aceshuffle()` resolves the joint reflection
  `(shuffle1, shuffle2) -> (1 - shuffle1, 1 - shuffle2)` toward the
  upward-trending orientation, measured by
  `E[shuffle(U) | U > 1/2] - E[shuffle(U) | U < 1/2]` summed over the pair.
* The transformation curve in every udp `plot()` method is drawn at line width
  1.5 (was 2).
* `plot()` gridlines for the `shuffle`, `udpcosine` and `udplegendre` classes
  are drawn with `segments()` rather than `abline()`, so they no longer extend
  past the panel under `par()` multi-panel layouts.
* `plot()` methods gain an `embellish` argument controlling the gridlines and
  the v-transform inadmissible-zone shading: `"colour"` (default, red), `"bw"`
  (grey) or `"none"`. This replaces the `shading` argument of the
  `vtransform` method.
* The `vtransform` `plot()` method gains `...`, forwarding further graphical
  parameters to `graphics::plot()` as the other `plot()` methods already do.
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
