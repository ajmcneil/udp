# The base-copula interface of a bsicopula.
#
# bsicopula() needs six things from its base copula (and from the copulas
# of a randmixture), which the internal generics below supply. Each takes the
# copula as its first argument. Methods exist for rvinecopulib's bicop_dist
# and, through the "ANY" method, for the copula package's parCopula; both
# packages are Suggests, so neither class is referred to formally. Another
# family plugs in by defining an S4 class and methods for it:
#
#   basecopula_supported(cop)             TRUE
#   basecopula_sample(cop, n)             n x 2 matrix of draws
#   basecopula_density(cop, v1, v2)       copula density
#   basecopula_cdf(cop, v1, v2)           copula CDF
#   basecopula_h(cop, given, x, given_var)
#       P(X <= x | G = given), X the coordinate other than 'given_var' (1 or
#       2) and G the one with that index, at interior points x in (0, 1)
#   basecopula_hinv(cop, given, q, given_var)
#       its inverse in x, for levels q in (0, 1)
#   basecopula_fast_hinv(cop)             TRUE if basecopula_hinv() is cheap
#                                         (a closed form or a library call),
#                                         FALSE if it iterates
#
# The wrappers basecopula_cond_cdf() and basecopula_cond_quantile() (R/hbsicopula.R)
# add the exact values at the ends of [0, 1] and the clipping to [0, 1] that
# the callers rely on, so a method need only handle the interior.

# TRUE if x is a parCopula object from the copula package -- FALSE, not an
# error, when copula isn't installed.
is_parCopula <- function(x) {
  requireNamespace("copula", quietly = TRUE) && methods::is(x, "parCopula")
}

# TRUE if x is a bicop_dist object from rvinecopulib. inherits() only looks
# at the class attribute, so this needs no namespace load.
is_bicop_dist <- function(x) inherits(x, "bicop_dist")

setOldClass("bicop_dist")

setGeneric("basecopula_supported", function(basecopula) standardGeneric("basecopula_supported"))
setGeneric("basecopula_sample", function(basecopula, n) standardGeneric("basecopula_sample"))
setGeneric("basecopula_density", function(basecopula, v1, v2) standardGeneric("basecopula_density"))
setGeneric("basecopula_cdf", function(basecopula, v1, v2) standardGeneric("basecopula_cdf"))
setGeneric("basecopula_h", function(basecopula, given, x, given_var) standardGeneric("basecopula_h"))
setGeneric("basecopula_hinv", function(basecopula, given, q, given_var) standardGeneric("basecopula_hinv"))
setGeneric("basecopula_fast_hinv", function(basecopula) standardGeneric("basecopula_fast_hinv"))

# Anything that is neither a bicop_dist nor a parCopula is not a base copula.
unsupported_basecopula <- function(basecopula) {
  stop(sprintf(
    "'%s' is not a supported base copula: use a bicop_dist (rvinecopulib), a parCopula (copula) or another family with the basecopula_*() methods.",
    class(basecopula)[1L]
  ), call. = FALSE)
}

# ---- bicop_dist (rvinecopulib) ---------------------------------------------
# hbicop() conventions: with cond_var = 1 the copula's first variable is
# conditioned on, with 2 the second; arguments are never transposed, so
# non-exchangeable copulas (90 and 270 degree rotations) stay correct.

setMethod("basecopula_supported", "bicop_dist", function(basecopula) TRUE)
setMethod("basecopula_sample", "bicop_dist", function(basecopula, n) rvinecopulib::rbicop(n, basecopula))
setMethod("basecopula_density", "bicop_dist", function(basecopula, v1, v2) {
  rvinecopulib::dbicop(cbind(v1, v2), basecopula)
})
setMethod("basecopula_cdf", "bicop_dist", function(basecopula, v1, v2) {
  rvinecopulib::pbicop(cbind(v1, v2), basecopula)
})
setMethod("basecopula_h", "bicop_dist", function(basecopula, given, x, given_var) {
  if (given_var == 1L) {
    rvinecopulib::hbicop(cbind(given, x), cond_var = 1, family = basecopula)
  } else {
    rvinecopulib::hbicop(cbind(x, given), cond_var = 2, family = basecopula)
  }
})
setMethod("basecopula_hinv", "bicop_dist", function(basecopula, given, q, given_var) {
  if (given_var == 1L) {
    rvinecopulib::hbicop(cbind(given, q), cond_var = 1, family = basecopula, inverse = TRUE)
  } else {
    rvinecopulib::hbicop(cbind(q, given), cond_var = 2, family = basecopula, inverse = TRUE)
  }
})
setMethod("basecopula_fast_hinv", "bicop_dist", function(basecopula) TRUE)

# ---- parCopula (copula) -----------------------------------------------------
# Reached through the "ANY" method, since copula is only Suggested. The
# h-function is a central difference of pCopula() in the conditioning
# argument: copula::cCopula() is not used, since it can condition only on the
# first coordinate, has no inverse for rotated copulas, and for a rotation
# that flips the conditioned coordinate disagrees with the derivative of
# pCopula() (copula 1.1.7). The inverse is found by bisection, the
# conditional CDF being increasing.

setMethod("basecopula_supported", "ANY", function(basecopula) is_parCopula(basecopula))
setMethod("basecopula_sample", "ANY", function(basecopula, n) {
  if (!is_parCopula(basecopula)) unsupported_basecopula(basecopula)
  copula::rCopula(n, basecopula)
})
setMethod("basecopula_density", "ANY", function(basecopula, v1, v2) {
  if (!is_parCopula(basecopula)) unsupported_basecopula(basecopula)
  copula::dCopula(cbind(v1, v2), basecopula)
})
setMethod("basecopula_cdf", "ANY", function(basecopula, v1, v2) {
  if (!is_parCopula(basecopula)) unsupported_basecopula(basecopula)
  copula::pCopula(cbind(v1, v2), basecopula)
})
setMethod("basecopula_h", "ANY", function(basecopula, given, x, given_var) {
  if (!is_parCopula(basecopula)) unsupported_basecopula(basecopula)
  lo <- pmax(given - 1e-5, 0)
  hi <- pmin(given + 1e-5, 1)
  cdf <- function(gg) {
    copula::pCopula(if (given_var == 1L) cbind(gg, x) else cbind(x, gg), basecopula)
  }
  (cdf(hi) - cdf(lo)) / (hi - lo)
})
setMethod("basecopula_hinv", "ANY", function(basecopula, given, q, given_var) {
  if (!is_parCopula(basecopula)) unsupported_basecopula(basecopula)
  lo <- rep(0, length(q))
  hi <- rep(1, length(q))
  for (i in seq_len(50L)) {
    mid <- (lo + hi) / 2
    below <- basecopula_cond_cdf(basecopula, given, mid, given_var) < q
    lo[below] <- mid[below]
    hi[!below] <- mid[!below]
  }
  (lo + hi) / 2
})
setMethod("basecopula_fast_hinv", "ANY", function(basecopula) FALSE)
