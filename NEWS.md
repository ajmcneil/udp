# udp 0.1.1

* `udpinverse()` now records a pre-image at a turning point of `T` once per
  branch meeting there, for every class: each column is always one branch,
  and at a turning value the selection probabilities are the one-sided
  limits. Previously v-transforms did this, but `udpzigzag` and `udpcosine`
  merged a shared peak or trough into one root (for `udpcosine` with the
  wrong probability, since the merged root counted as a single branch), and
  the polynomial classes merged a turning point into one root whose
  probability could fall back to an equal split over all of the row's
  pre-images when `g'` rounded to exactly `0`. `finalise_prob()` now shares
  the probability among infinite weights, and the polynomial classes mark
  every pre-image at a turning point of `g` as infinitely weighted. The
  internal `udpzigzaginverse()`, which existed only to merge roots, is gone.
* Fixed `udpderiv()` near turning points for `udplegendre`,
  `udplegendrebex` and `udpcosinebex`: within `1e-5` of a turning point it
  returned the left-hand slope on both sides, so the sign was wrong just to
  the right; for `udpcosinebex` with more than one turning point, each was
  also paired with another's `h''`, so the sign could be wrong on either
  side. The slope is now the left one at and left of the turning point and
  the right one to its right.
* Added `fitbsicopula()`, maximum likelihood estimation of a `bsicopula`
  from copula data -- pseudo-observations, or probability integral
  transforms from fitted margins (`pseudo = FALSE`) -- returning a new
  `fitbsicopula` class (with `show`, `coef`, `logLik` and `vcov` methods).
  The base copula must be an `rvinecopulib` `bicop_dist`; `randomizermod`
  may be `NULL` or a `randsdvine`, whose parametric copulas are all
  estimated (by default only `copZ1Z2_V1V2`, the tree-2 copulas being
  parameter-free independence copulas), from a two-stage fit by default.
  Unless `udpfix = TRUE`, the udp parameters are estimated too: those of
  v-transforms, the interior breakpoints of `udpzigzag` objects (their
  number of pieces and `up` held fixed), and the weights of
  `udplegendrebex` and `udpcosinebex` objects, normalized to unit length
  (their degree held fixed); parameter-free udps (shuffles, `vsymmetric()`,
  fixed-degree `udpcosine`/`udplegendre`) contribute none. Because the
  likelihood has a kink or cusp wherever a udp breakpoint crosses an
  observation, udp parameters are fitted by continuation, first with the
  carrier values coarsely clamped, which stops the optimizer stalling short
  of the main peak. Standard errors are optional: a parametric bootstrap
  of the whole procedure (`se = "bootstrap"`), or a quicker Hessian
  (`se = "hessian"`) that is unreliable for udp breakpoints. For base
  copulas whose density is unbounded at some corner of the unit square,
  carrier values are clamped into `[1 / (2n), 1 - 1 / (2n)]` in the base
  density during fitting (`vfloor = "auto"`), which stops the optimizer
  being trapped where both udps send one observation to that corner.
* `udpinverse()` for `udpzigzag` objects is vectorized, and several hundred
  times faster.
* `udplegendre()`, `udplegendrebex()` and `udpcosinebex()` build `F` from a
  vectorized sublevel-set measure (bisection plus safeguarded Newton on each
  monotone piece of `g`) instead of one polynomial root-finding per grid
  point, and `udpinverse()` for the two expansion classes finds all
  pre-images the same way: construction is about 6-8 times faster and
  `udpinverse()` about 10 times faster, with results unchanged to within
  the classes' interpolation accuracy.
* Fixed two failures of `dbsicopula()` for `randsdvine` models: pre-image
  matching was stricter (`1e-6`) than the spline-based classes' own
  accuracy, and cumulative selection probabilities could exceed `1` by
  rounding, which `rvinecopulib` rejected.
* License changed from MIT to GPL-3, matching `basiscor`, which depends on
  `udp`.
* Every `plot()` method whose axes are both `[0, 1]` (`shuffle`, `udpcosine`,
  `udpcosinebex`, `udplegendre`, `udplegendrebex`, and the `vtransform`
  `"transform"`/`"inverse"`/`"pdown"` panels) now draws with `asp = 1`, so
  the unit square renders as an actual square rather than being stretched
  to the plotting device's own aspect ratio.
* Added `udpcdf()` and `udpquantile()`, exposing `F` and `F^{-1}` directly
  for udp classes built as `T(u) = F(g(u))` (currently `udplegendre`; `F` was
  previously only reachable composed with `g`, via `udptrans()`). Supports
  computing quantities of `g(U)` itself, such as the maximal/minimal
  correlation attainable between two shifted-Legendre polynomials of a
  comonotonic/countermonotonic pair.
* Added `boundaryadjust()`, a general-purpose utility that nudges values
  exactly at `0` or `1` strictly inside `(0, 1)` -- useful for the output of
  `udptrans()` or of any other fitted transformation that should lie in the
  open unit interval but can attain the closed boundary exactly. The default
  tolerance is `1 / (2 * n)`, where `n` is the number of rows for a matrix
  (one tolerance per observation, however many columns) and `length(u)`
  otherwise.
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
