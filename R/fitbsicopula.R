# Maximum likelihood estimation of bsicopula objects from pseudo-observations.
#
# Every free parameter is optimised on an unconstrained scale (a scaled logit
# for parameters bounded on both sides, a shifted log for those bounded only
# below), but that scale never leaves this file: estimates, standard errors
# and the fitted object are all on the natural scale (delta, kappa, xi, rho,
# ...). Parameters are gathered into "blocks", one per estimable component
# (the base copula, each udp, each estimated randsdvine copula), each
# carrying its natural values, bounds and a setter that writes new values
# back into a bsicopula.

# Parameter names and bounds of the rvinecopulib parametric families, as
# enforced by rvinecopulib::bicop_dist() (version 0.7.x). The bounds are the
# same for every rotation.
bicop_par_info <- function(family) {
  info <- switch(family,
    indep    = list(names = character(0), lower = numeric(0), upper = numeric(0)),
    gaussian = list(names = "rho", lower = -1, upper = 1),
    t        = list(names = c("rho", "nu"), lower = c(-1, 2), upper = c(1, 50)),
    clayton  = list(names = "theta", lower = 1e-10, upper = 28),
    gumbel   = list(names = "theta", lower = 1, upper = 50),
    frank    = list(names = "theta", lower = -35, upper = 35),
    joe      = list(names = "theta", lower = 1, upper = 30),
    bb1      = list(names = c("theta", "delta"), lower = c(0, 1), upper = c(7, 7)),
    bb6      = list(names = c("theta", "delta"), lower = c(1, 1), upper = c(6, 8)),
    bb7      = list(names = c("theta", "delta"), lower = c(1, 0.01), upper = c(6, 25)),
    bb8      = list(names = c("theta", "delta"), lower = c(1, 1e-4), upper = c(8, 1)),
    NULL
  )
  if (is.null(info)) {
    stop(sprintf("estimation is not supported for the '%s' copula family.", family),
      call. = FALSE
    )
  }
  info
}

# TRUE when a bicop_dist base copula's density is unbounded at any corner of
# the unit square. udp breakpoints map interior points to V = 0 (a
# v-transform's fulcrum) or to V = 0 and 1 alternately (a zigzag's
# breakpoints), so the likelihood spikes to infinity when both udps send the
# same observation to an unbounded corner -- and the spike is log-singular,
# so even a slow blow-up (Gaussian, plain Gumbel, BB6, all like eps^-0.6)
# draws Nelder-Mead straight into it. fitbsicopula() therefore clamps the
# carrier values for every such copula. Only independence, Frank and BB8
# are bounded at all four corners. (Every other family is unbounded at some
# corner for every rotation and parameter value -- the Gaussian at (0, 0)
# and (1, 1) for rho > 0 and at the other two for rho < 0 -- so the rule
# needs neither, and the objective stays one function when rho changes sign
# mid-fit.)
setGeneric("basecopula_needs_vfloor", function(cop) standardGeneric("basecopula_needs_vfloor"))

setMethod("basecopula_needs_vfloor", "bicop_dist", function(cop) {
  !(cop$family %in% c("indep", "frank", "bb8"))
})

# The absolute spherical t copula has an asymptote at (1, 1) (upper tail
# dependence), and a zigzag's breakpoints map interior points to 1.
setMethod("basecopula_needs_vfloor", "astcopula", function(cop) TRUE)

# Any other base copula: be safe.
setMethod("basecopula_needs_vfloor", "ANY", function(cop) TRUE)

# A parameter block for a bicop_dist held somewhere in a bsicopula: get()
# extracts it, set(obj, cop) stores a replacement.
bicop_block <- function(label, cop, set) {
  info <- bicop_par_info(cop$family)
  list(
    label = label, names = info$names, value = as.numeric(cop$parameters),
    maps = bounded_maps(info$lower, info$upper),
    set = function(obj, value) {
      set(obj, rvinecopulib::bicop_dist(cop$family, cop$rotation, value))
    }
  )
}

# The parameter block of a base copula held in a bsicopula, a generic so that
# a new base-copula family can be estimated by defining a method. set(obj, cop)
# stores a replacement copula in obj.
setGeneric("basecopula_fit_block", function(cop, label, set) standardGeneric("basecopula_fit_block"))

setMethod("basecopula_fit_block", "bicop_dist", function(cop, label, set) bicop_block(label, cop, set))

# The absolute spherical t copula has the single parameter nu > 0, estimated
# on the log scale.
setMethod("basecopula_fit_block", "astcopula", function(cop, label, set) {
  list(
    label = label, names = "nu", value = cop@nu,
    maps = bounded_maps(0, Inf),
    set = function(obj, value) set(obj, astcopula(value))
  )
})

setMethod("basecopula_fit_block", "ANY", function(cop, label, set) {
  stop(sprintf(
    "estimation is not supported for a base copula of class '%s'; use a bicop_dist (rvinecopulib) or an astcopula.",
    class(cop)[1L]
  ),
  call. = FALSE
  )
})

# A one-line description of a base copula, for show().
basecopula_label <- function(cop) {
  if (is_bicop_dist(cop)) {
    sprintf("%s (rotation %s)", cop$family, cop$rotation)
  } else if (methods::is(cop, "astcopula")) {
    "absolute spherical t (ast)"
  } else {
    class(cop)[1L]
  }
}

# udp_fitpars() gives either box bounds (lower/upper) or, for parameters
# with a joint constraint such as a zigzag's ordered breakpoints, its own
# 'maps'.
udp_block <- function(label, x, set) {
  fp <- udp_fitpars(x)
  list(
    label = label, names = names(fp$value), value = unname(fp$value),
    maps = if (is.null(fp$maps)) bounded_maps(fp$lower, fp$upper) else fp$maps,
    set = function(obj, value) set(obj, udp_setfitpars(x, value))
  )
}

# The parameter blocks to estimate, dropping any with no free parameters --
# so a randsdvine's tree-2 copulas, independence by default, are estimated
# only when the user has made them parametric.
fit_blocks <- function(object, udpfix) {
  blocks <- list(basecopula_fit_block(object@basecopula, "basecopula", function(obj, cop) {
    obj@basecopula <- cop
    obj
  }))
  if (!udpfix) {
    blocks <- c(blocks, list(
      udp_block("udp1", object@udp1, function(obj, x) {
        obj@udp1 <- x
        obj
      }),
      udp_block("udp2", object@udp2, function(obj, x) {
        obj@udp2 <- x
        obj
      })
    ))
  }
  rm <- object@randomizermod
  if (!is.null(rm)) {
    blocks <- c(blocks, list(
      bicop_block("copZ1Z2_V1V2", rm@copZ1Z2_V1V2, function(obj, cop) {
        obj@randomizermod@copZ1Z2_V1V2 <- cop
        obj
      }),
      bicop_block("copZ1V2_V1", rm@copZ1V2_V1, function(obj, cop) {
        obj@randomizermod@copZ1V2_V1 <- cop
        obj
      }),
      bicop_block("copV1Z2_V2", rm@copV1Z2_V2, function(obj, cop) {
        obj@randomizermod@copV1Z2_V2 <- cop
        obj
      })
    ))
  }
  Filter(function(b) length(b$value) > 0L, blocks)
}

# Natural <-> unconstrained scale, elementwise. The map back is squeezed
# strictly inside the bounds, since rvinecopulib rejects values at or beyond
# them and a v-transform's delta must stay inside (0, 1).
to_free <- function(x, lower, upper) {
  both <- is.finite(lower) & is.finite(upper)
  below <- is.finite(lower) & !is.finite(upper)
  y <- x
  y[both] <- stats::qlogis((x[both] - lower[both]) / (upper[both] - lower[both]))
  y[below] <- log(x[below] - lower[below])
  y
}

from_free <- function(y, lower, upper) {
  both <- is.finite(lower) & is.finite(upper)
  below <- is.finite(lower) & !is.finite(upper)
  x <- y
  p <- pmin(pmax(stats::plogis(y[both]), 1e-10), 1 - 1e-10)
  x[both] <- lower[both] + (upper[both] - lower[both]) * p
  x[below] <- lower[below] + pmax(exp(y[below]), 1e-10)
  x
}

# Move starting values off (or away from) their bounds: a start sitting on a
# bound maps to an infinite unconstrained value, and one very close to it
# leaves Nelder-Mead's initial simplex almost no room to move.
nudge_start <- function(x, lower, upper) {
  both <- is.finite(lower) & is.finite(upper)
  below <- is.finite(lower) & !is.finite(upper)
  margin <- 0.01 * (upper[both] - lower[both])
  x[both] <- pmin(pmax(x[both], lower[both] + margin), upper[both] - margin)
  x[below] <- pmax(x[below], lower[below] + 0.01)
  x
}

# A block's maps between the natural and unconstrained scales, and its
# starting-value nudge, for elementwise box bounds.
bounded_maps <- function(lower, upper) {
  list(
    to_free = function(x) to_free(x, lower, upper),
    from_free = function(y) from_free(y, lower, upper),
    nudge = function(x) nudge_start(x, lower, upper)
  )
}

# Negative log-likelihood of a bsicopula at the data. Any failure along the
# way -- a zero or non-finite density, or an error from udpinverse()'s
# pre-image matching at an awkward parameter value -- returns a large finite
# penalty rather than stopping the optimiser.
bsicopula_negll <- function(U, object, vfloor) {
  d <- tryCatch(
    suppressWarnings(dbsicopula_eval(U[, 1], U[, 2], object, vfloor)),
    error = function(e) NA_real_
  )
  if (anyNA(d) || any(!is.finite(d)) || any(d <= 0)) {
    return(1e10)
  }
  -sum(log(d))
}

# One maximum likelihood fit of 'object' over the parameters in 'blocks',
# with Hessian-based standard errors when 'hessian' is TRUE. A block's
# natural values can outnumber its free parameters -- a udplegendrebex's
# unit-length weights have one fewer degree of freedom than weights -- so
# natural and free vectors are split into blocks separately, and 'npar'
# counts the free parameters.
fit_stage <- function(U, object, blocks, vfloor, hessian, method, control) {
  nat_sizes <- vapply(blocks, function(b) length(b$value), integer(1))
  free_sizes <- vapply(blocks, function(b) length(b$maps$to_free(b$value)), integer(1))
  idx <- rep(seq_along(blocks), nat_sizes)
  free_idx <- rep(seq_along(blocks), free_sizes)
  labels <- unlist(lapply(blocks, function(b) paste(b$label, b$names, sep = ".")))
  by_block <- function(v, i) split(v, factor(i, levels = seq_along(blocks)))
  all_to_free <- function(nat) {
    unlist(Map(function(b, v) b$maps$to_free(v), blocks, by_block(nat, idx)), use.names = FALSE)
  }
  all_from_free <- function(th) {
    unlist(Map(function(b, v) b$maps$from_free(v), blocks, by_block(th, free_idx)),
      use.names = FALSE
    )
  }
  start <- unlist(lapply(blocks, function(b) b$maps$nudge(b$value)), use.names = FALSE)
  npar <- sum(free_sizes)
  nnat <- length(start)

  build <- function(nat) {
    obj <- object
    for (i in seq_along(blocks)) obj <- blocks[[i]]$set(obj, nat[idx == i])
    obj
  }
  negll_nat <- function(nat) {
    obj <- tryCatch(build(nat), error = function(e) NULL)
    if (is.null(obj)) 1e10 else bsicopula_negll(U, obj, vfloor)
  }

  if (npar == 0L) {
    return(list(
      object = object, estimate = numeric(0), se = numeric(0), boot = NULL,
      vcov = matrix(numeric(0), 0, 0), loglik = -bsicopula_negll(U, object, vfloor),
      npar = 0L, convergence = 0L, message = "no free parameters",
      counts = c(`function` = NA_integer_, gradient = NA_integer_), method = NA_character_
    ))
  }

  if (is.null(method)) method <- if (npar == 1L) "BFGS" else "Nelder-Mead"
  if (method == "Nelder-Mead" && is.null(control$maxit)) control$maxit <- 2000
  fn <- function(th) negll_nat(all_from_free(th))
  opt <- stats::optim(all_to_free(start), fn, method = method, control = control)
  # The likelihood has kinks (at each observation, as a udp breakpoint
  # passes it), on which Nelder-Mead's simplex can collapse (code 10) short
  # of the optimum. Restarting from where it stopped, with a fresh simplex,
  # usually gets it moving again.
  restarts <- 0L
  while (method == "Nelder-Mead" && opt$convergence == 10L && restarts < 5L) {
    counts <- opt$counts
    opt <- stats::optim(opt$par, fn, method = method, control = control)
    opt$counts <- opt$counts + counts
    restarts <- restarts + 1L
  }
  est <- stats::setNames(all_from_free(opt$par), labels)

  V <- matrix(NA_real_, nnat, nnat, dimnames = list(labels, labels))
  sds <- stats::setNames(rep(NA_real_, nnat), labels)
  if (hessian) {
    # Inverse Hessian on the free scale, carried to the natural scale by the
    # delta method (J the Jacobian of the free-to-natural map).
    H <- stats::optimHess(opt$par, fn)
    Vfree <- tryCatch(solve(H), error = function(e) NULL)
    if (is.null(Vfree) || any(diag(Vfree) <= 0)) {
      warning("Hessian is not positive definite; standard errors are NA.", call. = FALSE)
    } else {
      h <- 1e-6
      J <- vapply(seq_len(npar), function(k) {
        e <- replace(numeric(npar), k, h)
        (all_from_free(opt$par + e) - all_from_free(opt$par - e)) / (2 * h)
      }, numeric(nnat))
      V[] <- J %*% Vfree %*% t(J)
      sds[] <- sqrt(pmax(diag(V), 0))
    }
  }

  if (opt$convergence != 0L) {
    warning(sprintf("optim() did not converge (code %d).", opt$convergence), call. = FALSE)
  }
  list(
    object = build(unname(est)), estimate = est, se = sds, vcov = V,
    loglik = -opt$value, npar = npar, convergence = as.integer(opt$convergence),
    message = if (is.null(opt$message)) "" else opt$message,
    counts = opt$counts, method = method
  )
}

#' Class of fitted bivariate stochastic inversion copulas
#'
#' The result of [fitbsicopula()].
#'
#' @slot bsicopula the fitted \linkS4class{bsicopula}.
#' @slot estimate named vector of the estimated parameters, on their natural
#'   scale. Names are `component.parameter`, e.g. `basecopula.rho`,
#'   `udp1.delta`, `copZ1Z2_V1V2.rho`.
#' @slot se standard errors matching `estimate`; all `NA` when the fit was
#'   run with `se = FALSE`.
#' @slot vcov the matching covariance matrix -- the inverse observed
#'   information (`se = "hessian"`) or the covariance of the bootstrap
#'   estimates (`se = "bootstrap"`); all `NA` when `se = FALSE`.
#' @slot se_method `"none"`, `"hessian"` or `"bootstrap"`.
#' @slot boot for `se = "bootstrap"`, the `B`-row matrix of bootstrap
#'   estimates, one column per parameter (a row of `NA`s for a replicate
#'   whose refit failed); otherwise a matrix with no rows. Useful for
#'   percentile intervals.
#' @slot loglik the maximized log-likelihood.
#' @slot nobs number of observations.
#' @slot npar number of free parameters, as used by AIC and BIC. It is one
#'   fewer than the number of weights reported for each
#'   \linkS4class{udplegendrebex} or \linkS4class{udpcosinebex}, whose
#'   weights are constrained to unit length.
#' @slot convergence the [stats::optim()] convergence code (`0` for success).
#' @slot message the [stats::optim()] message, if any.
#' @slot counts the [stats::optim()] function/gradient evaluation counts.
#' @slot method the [stats::optim()] method used.
#' @slot vfloor the clamp `f` applied to the carrier values (into `[f, 1 -
#'   f]`) before evaluating the base copula density, or `NA` if none was
#'   applied.
#' @slot stage1 for a two-stage fit, the first-stage fit (with independent
#'   randomizers) as a `fitbsicopula` object; otherwise `NULL`.
#'
#' @seealso [fitbsicopula()].
#' @export
setClass("fitbsicopula", slots = list(
  bsicopula = "bsicopula", estimate = "numeric", se = "numeric", vcov = "matrix",
  se_method = "character", boot = "matrix", loglik = "numeric", nobs = "integer", npar = "integer", convergence = "integer",
  message = "character", counts = "integer", method = "character",
  vfloor = "numeric", stage1 = "ANY"
))

new_fitbsicopula <- function(res, nobs, vfloor, se_method = "none", stage1 = NULL) {
  boot <- if (is.null(res$boot)) {
    matrix(numeric(0), 0, length(res$estimate), dimnames = list(NULL, names(res$estimate)))
  } else {
    res$boot
  }
  new("fitbsicopula",
    bsicopula = res$object, estimate = res$estimate, se = res$se, vcov = res$vcov,
    se_method = se_method, boot = boot, loglik = res$loglik, nobs = nobs, npar = res$npar, convergence = res$convergence,
    message = res$message, counts = res$counts, method = res$method,
    vfloor = if (is.null(vfloor)) NA_real_ else vfloor, stage1 = stage1
  )
}

#' Fit a bivariate stochastic inversion copula by maximum likelihood
#'
#' Estimates the parameters of a \linkS4class{bsicopula} from bivariate
#' copula data -- pseudo-observations (ranks), or probability integral
#' transforms from fitted marginal models (`pseudo = FALSE`) -- by maximizing
#' the likelihood built from [dbsicopula()]. `object` supplies both the model structure -- the
#' copula families and rotations, the udp classes, the randomizer model --
#' and the starting values for the optimization.
#'
#' The base copula must be a parametric `bicop_dist` object from
#' \pkg{rvinecopulib} or an \linkS4class{astcopula} (whose one parameter,
#' `nu`, is estimated on the log scale), and every copula of a
#' \linkS4class{randsdvine} a parametric `bicop_dist`; the family and
#' rotation are held fixed and only the parameters are estimated.
#' `randomizermod` must be `NULL` or a \linkS4class{randsdvine}.
#'
#' **udp parameters.** With `udpfix = FALSE`, the continuous parameters of
#' `udp1` and `udp2` are estimated alongside the copula parameters. At
#' present that is supported for v-transforms (`delta`, `kappa`, `xi`) and
#' for \linkS4class{udpzigzag} objects, whose interior breakpoints are
#' estimated (reported as `break1`, `break2`, ...) while the number of
#' pieces and the direction of the first piece, `up`, are held at their
#' given values -- set them from a priori knowledge, [aceshuffle()] output
#' or a plot -- and for \linkS4class{udplegendrebex} and
#' \linkS4class{udpcosinebex} objects, whose weights are estimated
#' (reported as `coef1`, `coef2`, ...) with the degree held fixed. Only the
#' direction of the weight vector matters (`T` is unchanged by a positive
#' rescaling), so the weights are normalized to unit length and have one
#' fewer free parameter than weights, which is what `npar`, AIC and BIC
#' count; starting weights from `basisexpand()` in the \pkg{basiscor}
#' package are already normalized. The fit keeps the sign of the
#' largest-magnitude starting weight, so it stays in the starting
#' orientation: negating every weight gives the reflection `1 - T`, the
#' counterpart of flipping a zigzag's `up`. The fitted object is rebuilt with
#' the constructor's default `ngrid`. udps with no continuous parameters -- [vsymmetric()], shuffles such as
#' [udpid()] and [udpflip()], and \linkS4class{udpcosine} or
#' \linkS4class{udplegendre} objects, whose degree is regarded as fixed --
#' simply contribute nothing. For other classes use `udpfix = TRUE`, which
#' holds both udps at their given values: the natural choice when they have
#' been determined externally, e.g. by [aceshuffle()].
#'
#' **Randomizer.** When `randomizermod` is a \linkS4class{randsdvine}, the
#' parameters of all three of its copulas are estimated. The two tree-2
#' copulas default to the independence copula, which has no parameters, so
#' by default only the tree-3 copula `copZ1Z2_V1V2` is estimated; to
#' estimate a tree-2 copula too, give it a parametric family (a Gaussian
#' with a small correlation, say) as its starting value. With
#' `twostage = TRUE`, the model is first fitted with independent
#' randomizers, and that fit's base copula and udps are used as starting
#' values for the full fit. Since the independent-randomizer model is
#' misspecified for such data, those values can lie in the basin of a poor
#' local maximum (with an absolute spherical t base copula the likelihood was
#' 30 to 50 units lower in simulations), so the full fit is also run from the
#' starting values in `object` and the better of the two is kept; this
#' doubles the work, and `twostage = FALSE` gives the single fit from `object`.
#'
#' **Boundary clamp.** udp breakpoints map interior points to `0` or `1`: a
#' v-transform's fulcrum `delta` to `0`, a zigzag's breakpoints to `0` and
#' `1` alternately. When the base copula's density is unbounded at a corner
#' of the unit square, the likelihood becomes unbounded as breakpoints of
#' both udps approach the same observation and send it to that corner, and
#' even a slow blow-up is enough to trap the optimizer there.
#' `vfloor = "auto"` guards against this by clamping the carrier values
#' into `[1 / (2n), 1 - 1 / (2n)]` before evaluating the base copula
#' density, for every base copula unbounded at some corner -- all but the
#' independence, Frank and BB8 copulas. Observations with a `u` outside
#' `[1 / (2n), 1 - 1 / (2n)]` (possible with `pseudo = FALSE`) are exempt:
#' every udp sends `u = 0` and `u = 1` to fixed values, so their extreme
#' carrier values are data, not the product of a moving breakpoint.
#' `vfloor = TRUE` or `FALSE` forces the clamp on or off. The clamp enters
#' only the fitting; [dbsicopula()] is exact.
#'
#' **Optimization and standard errors.** Parameters are optimized with
#' [stats::optim()] on an internal unconstrained scale; every result is
#' reported on the natural scale. The likelihood is only piecewise smooth
#' in a udp's breakpoint parameters (a v-transform's `delta`): it has a
#' kink or cusp wherever the breakpoint crosses an observation, and for many
#' base copulas the profile likelihood in `delta` (or a zigzag breakpoint) is
#' a fine sawtooth. The optimizer may stop on a neighbouring tooth, which
#' moves `delta` by about `1 / n` and the log-likelihood by a unit or two.
#' Worse, from a neutral start the optimizer can stall on one of these
#' features well short of the main peak. So when udp parameters are
#' estimated, the model is fitted by continuation: first with the carrier
#' values coarsely clamped into `[0.05, 0.95]`, which smooths the kinks
#' away, then into `[0.01, 0.99]`, then as requested by `vfloor`, each fit
#' starting the next. Nelder-Mead is also restarted automatically when its
#' simplex collapses on a kink.
#'
#' `se = "bootstrap"` (or `TRUE`) gives parametric bootstrap standard
#' errors: `B` samples of size `n` are drawn from the fitted model with
#' [rbsicopula()], converted to pseudo-observations when `pseudo = TRUE`,
#' and put through the whole estimation procedure again, from the same
#' starting values as the original fit. This accounts for the rank
#' transformation and is the recommended choice. With `pseudo = FALSE` the
#' simulated uniforms are refitted as they are; the standard errors then
#' ignore the uncertainty from estimating the marginal models, as in
#' inference functions for margins. `se = "hessian"` is much quicker, using the
#' numerically differentiated Hessian of the log-likelihood on the natural
#' scale, but it ignores the rank transformation. It is also unreliable for
#' `delta`: when the base copula density vanishes along the edge `V = 0`
#' (Gaussian with positive correlation, Clayton, Gumbel, ...), the
#' likelihood has a cusp at every observation and the Hessian can understate
#' the standard error of `delta` many times over; when the edge density is
#' finite and nonzero (Joe, survival Clayton), it still understates it,
#' by a smaller factor.
#'
#' @param U a two-column numeric matrix (or data frame) of copula data, with
#'   every value strictly inside `(0, 1)`; [boundaryadjust()] moves values
#'   of exactly `0` or `1` inside.
#' @param object a \linkS4class{bsicopula} giving the model structure and
#'   starting values.
#' @param pseudo logical; are `U` pseudo-observations (ranks scaled into
#'   `(0, 1)`)? Set `FALSE` for probability integral transforms from fitted
#'   marginal models. Affects only the bootstrap.
#' @param udpfix logical; hold the parameters of `udp1` and `udp2` fixed at
#'   their given values?
#' @param twostage logical; when `randomizermod` is a
#'   \linkS4class{randsdvine}, obtain starting values from a first fit with
#'   independent randomizers (the full fit is also run from `object` and the
#'   better kept; see Details)? Ignored otherwise.
#' @param se `FALSE` (the default) for no standard errors, `"bootstrap"` (or
#'   `TRUE`) for parametric bootstrap standard errors, or `"hessian"` for
#'   Hessian-based ones; see Details.
#' @param B number of bootstrap replicates for `se = "bootstrap"`.
#' @param vfloor `"auto"`, `TRUE` or `FALSE`; see Details.
#' @param method the [stats::optim()] method. Defaults to `"Nelder-Mead"`, or
#'   `"BFGS"` when there is a single free parameter.
#' @param control a list of control parameters passed to [stats::optim()].
#'   For Nelder-Mead, `maxit` defaults to `2000`.
#'
#' @return An object of class \linkS4class{fitbsicopula}.
#' @references
#' McNeil, A. J. and Nešlehová, J. G. (2026). Stochastic inversion of
#' multivariate uniform-distribution-preserving transformations.
#' \href{https://arxiv.org/abs/2607.07174}{arXiv:2607.07174}
#'
#' Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and
#' D-vines for stochastic volatility. *Journal of Multivariate Analysis*,
#' **217**, 105695. \doi{10.1016/j.jmva.2026.105695}
#' @include bsicopula.R astcopula.R
#' @export
#'
#' @examples
#' if (requireNamespace("rvinecopulib", quietly = TRUE)) {
#'   set.seed(1)
#'   truth <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.6), vlinear(0.4), vlinear(0.6))
#'   U <- rvinecopulib::pseudo_obs(rbsicopula(500, truth))
#'   start <- bsicopula(rvinecopulib::bicop_dist("gaussian", 0, 0.3), vlinear(0.5), vlinear(0.5))
#'   fit <- fitbsicopula(U, start)
#'   fit
#'   \donttest{
#'   fitbsicopula(U, start, se = "bootstrap", B = 100)
#'   }
#' }
fitbsicopula <- function(U, object, pseudo = TRUE, udpfix = FALSE, twostage = TRUE,
                         se = FALSE, B = 200, vfloor = "auto", method = NULL,
                         control = list()) {
  if (!methods::is(object, "bsicopula")) {
    stop("'object' must be an object of class 'bsicopula'.", call. = FALSE)
  }
  U <- as.matrix(U)
  if (!is.numeric(U) || ncol(U) != 2L || nrow(U) < 2L || anyNA(U)) {
    stop("'U' must be a two-column numeric matrix without missing values.", call. = FALSE)
  }
  if (any(U <= 0 | U >= 1)) {
    stop(
      "every value of 'U' must be strictly inside (0, 1); boundaryadjust() moves values of exactly 0 or 1 inside.",
      call. = FALSE
    )
  }
  for (arg in c("pseudo", "udpfix", "twostage")) {
    val <- get(arg)
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop(sprintf("'%s' must be TRUE or FALSE.", arg), call. = FALSE)
    }
  }
  se_method <- if (isFALSE(se)) {
    "none"
  } else if (isTRUE(se)) {
    "bootstrap"
  } else if (identical(se, "hessian") || identical(se, "bootstrap")) {
    se
  } else {
    stop("'se' must be FALSE, TRUE, \"hessian\" or \"bootstrap\".", call. = FALSE)
  }
  if (se_method == "bootstrap" &&
    (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 2 || B != round(B))) {
    stop("'B' must be a single integer of at least 2.", call. = FALSE)
  }
  rm <- object@randomizermod
  if (!is.null(rm) && !methods::is(rm, "randsdvine")) {
    stop("fitbsicopula() supports only randomizermod = NULL or a 'randsdvine' object.",
      call. = FALSE
    )
  }

  n <- nrow(U)
  use_floor <- if (identical(vfloor, "auto")) {
    basecopula_needs_vfloor(object@basecopula)
  } else if (isTRUE(vfloor) || isFALSE(vfloor)) {
    vfloor
  } else {
    stop("'vfloor' must be \"auto\", TRUE or FALSE.", call. = FALSE)
  }
  floor_value <- if (use_floor) 1 / (2 * n) else NULL

  fit <- fit_procedure(U, object, udpfix, twostage, floor_value, se_method == "hessian",
    method, control
  )
  res <- fit$res
  stage1 <- if (is.null(fit$res1)) NULL else new_fitbsicopula(fit$res1, n, floor_value)
  if (se_method == "bootstrap" && res$npar > 0L) {
    res <- bootstrap_se(res, object, n, B, pseudo, udpfix, twostage, floor_value,
      method, control
    )
  }
  new_fitbsicopula(res, n, floor_value, se_method, stage1)
}

# The whole estimation procedure from the starting object: for a randsdvine
# with twostage = TRUE, a first fit with independent randomizers whose base
# copula and udps then start the full fit. The independent-randomizer model
# is misspecified for the data, so its estimates can lie in the basin of a
# poor local maximum of the full likelihood (with an absolute spherical t base
# copula, 30 to 50 log-likelihood units below the best in simulations), and
# the full fit is therefore also run from the user's own starting values, the
# better of the two being kept. Returns the full fit's result ('res') and the
# first stage's ('res1', NULL for a single-stage fit).
fit_procedure <- function(U, object, udpfix, twostage, vfloor, hessian, method, control) {
  res1 <- NULL
  if (!is.null(object@randomizermod) && twostage) {
    object1 <- object
    object1@randomizermod <- NULL
    res1 <- fit_continued(U, object1, udpfix, vfloor, FALSE, method, control)
    object2 <- object
    object2@basecopula <- res1$object@basecopula
    object2@udp1 <- res1$object@udp1
    object2@udp2 <- res1$object@udp2
    res <- fit_continued(U, object2, udpfix, vfloor, hessian, method, control)
    direct <- tryCatch(
      fit_continued(U, object, udpfix, vfloor, hessian, method, control),
      error = function(e) NULL
    )
    if (!is.null(direct) && direct$loglik > res$loglik) res <- direct
    return(list(res = res, res1 = res1))
  }
  list(res = fit_continued(U, object, udpfix, vfloor, hessian, method, control), res1 = res1)
}

# Carrier-value clamps for the preliminary fits of fit_continued().
continuation_levels <- c(0.05, 0.01)

# One fit by continuation. The likelihood has a kink or cusp wherever a udp
# breakpoint crosses an observation, all produced by carrier values near 0
# or 1; from a neutral start Nelder-Mead often stalls on one well short of
# the main peak (in simulations, 10 to 20 log-likelihood units short in
# most fits). Clamping the carrier values coarsely smooths those features
# away, so when udp parameters are being estimated the model is fitted
# first with each clamp in continuation_levels, each fit starting the next,
# and only then with the requested clamp 'vfloor'. With the udps fixed the
# surface is smooth and one fit suffices.
fit_continued <- function(U, object, udpfix, vfloor, hessian, method, control) {
  blocks <- fit_blocks(object, udpfix)
  if (any(vapply(blocks, function(b) b$label %in% c("udp1", "udp2"), logical(1)))) {
    for (f in continuation_levels) {
      object <- suppressWarnings(
        fit_stage(U, object, fit_blocks(object, udpfix), f, FALSE, method, control)
      )$object
    }
  }
  fit_stage(U, object, fit_blocks(object, udpfix), vfloor, hessian, method, control)
}

# Parametric bootstrap: B samples of size n from the fitted model, each
# ranked to pseudo-observations if the real data were, and put through the
# whole estimation procedure again from the user's starting object -- not
# from the fitted values: the likelihood is rugged in udp breakpoints, and
# refits started at the fitted values tend to stay on the nearby tooth,
# understating the spread. Returns 'res' with se, vcov and boot filled in.
# A refit that errors leaves a row of NAs and is left out of se/vcov.
bootstrap_se <- function(res, start, n, B, pseudo, udpfix, twostage, vfloor, method,
                         control) {
  fitted <- res$object
  est <- matrix(NA_real_, B, length(res$estimate), dimnames = list(NULL, names(res$estimate)))
  for (b in seq_len(B)) {
    Ub <- rbsicopula(n, fitted)
    if (pseudo) Ub <- apply(Ub, 2, rank) / (n + 1)
    rb <- tryCatch(
      suppressWarnings(
        fit_procedure(Ub, start, udpfix, twostage, vfloor, FALSE, method, control)$res
      ),
      error = function(e) NULL
    )
    if (!is.null(rb)) est[b, ] <- rb$estimate
  }
  ok <- stats::complete.cases(est)
  if (sum(ok) < 2L) {
    stop("fewer than two bootstrap refits succeeded.", call. = FALSE)
  }
  if (!all(ok)) {
    warning(sprintf("%d of %d bootstrap refits failed and were dropped.", sum(!ok), B),
      call. = FALSE
    )
  }
  res$vcov <- stats::cov(est[ok, , drop = FALSE])
  res$se <- sqrt(diag(res$vcov))
  res$boot <- est
  res
}

#' @describeIn fitbsicopula-class Show method for fitbsicopula objects.
#' @param object an object of class \linkS4class{fitbsicopula}.
#' @export
setMethod("show", "fitbsicopula", function(object) {
  bc <- object@bsicopula
  udpname <- function(x) if (methods::is(x, "vtransform")) x@name else class(x)
  rm <- bc@randomizermod
  cat("Fitted bsicopula (maximum likelihood, n = ", object@nobs, ")\n", sep = "")
  cat("base copula: ", basecopula_label(bc@basecopula), "\n", sep = "")
  cat("udp1: ", udpname(bc@udp1), ", udp2: ", udpname(bc@udp2), "\n", sep = "")
  cat("randomizer: ", if (is.null(rm)) "independent" else "randsdvine", "\n", sep = "")
  if (object@npar > 0L) {
    tab <- cbind(estimate = object@estimate, se = object@se)
    if (all(is.na(object@se))) tab <- tab[, "estimate", drop = FALSE]
    cat("\n")
    print(signif(tab, 5))
    if (object@se_method == "bootstrap") {
      nok <- sum(stats::complete.cases(object@boot))
      cat("(parametric bootstrap standard errors, ", nok, " of ", nrow(object@boot),
        " refits)\n",
        sep = ""
      )
    } else if (object@se_method == "hessian") {
      cat("(Hessian standard errors; unreliable for udp breakpoints such as delta)\n")
    }
  }
  cat("\nlog-likelihood: ", format(object@loglik, digits = 6),
    "  AIC: ", format(-2 * object@loglik + 2 * object@npar, digits = 6),
    "  BIC: ", format(-2 * object@loglik + log(object@nobs) * object@npar, digits = 6),
    "\n",
    sep = ""
  )
  if (object@convergence != 0L) {
    cat("optim() did not converge (code ", object@convergence, ")\n", sep = "")
  }
  if (!is.na(object@vfloor)) {
    cat("carrier values clamped into [", format(object@vfloor, digits = 4), ", 1 - ",
      format(object@vfloor, digits = 4), "] in the base copula density\n",
      sep = ""
    )
  }
})

#' @describeIn fitbsicopula-class The named vector of estimates.
#' @export
setMethod("coef", "fitbsicopula", function(object) object@estimate)

#' @describeIn fitbsicopula-class The log-likelihood, as a `logLik` object,
#'   so [stats::AIC()] and [stats::BIC()] work.
#' @param ... unused.
#' @exportS3Method stats::logLik
logLik.fitbsicopula <- function(object, ...) {
  structure(object@loglik, df = object@npar, nobs = object@nobs, class = "logLik")
}

#' @describeIn fitbsicopula-class The covariance matrix of the estimates
#'   (requires the fit to have been run with standard errors).
#' @exportS3Method stats::vcov
vcov.fitbsicopula <- function(object, ...) {
  if (all(is.na(object@vcov))) {
    stop("no covariance matrix available; refit with se = \"bootstrap\" or \"hessian\".",
      call. = FALSE
    )
  }
  object@vcov
}
