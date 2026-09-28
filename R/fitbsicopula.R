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

# TRUE when a bicop_dist base copula's density is unbounded at the (0, 0)
# corner. Such copulas make the likelihood of a udp whose pre-image of 0 is
# an interior point (a v-transform's fulcrum, say) spike to infinity when
# both udps' pre-images of 0 sit on the same observation -- and the spike is
# log-singular, so even a slow blow-up (Gaussian with rho > 0, plain Gumbel,
# BB6, all like eps^-0.6) draws Nelder-Mead straight into it. fitbsicopula()
# therefore floors the carrier values for all of them. The Gaussian is
# floored whatever the sign of rho, since rho can change sign during a fit
# and the objective must stay one function; for rho < 0 the floor only moves
# the few points with V < 1 / (2n). Left alone: independence, Frank, Joe and
# BB8 (rotation 0), survival Clayton, and every 90/270 rotation, whose
# densities are bounded at (0, 0).
basecopula_needs_vfloor <- function(cop) {
  fam <- cop$family
  rot <- cop$rotation
  fam %in% c("gaussian", "t") ||
    (rot == 0 && fam %in% c("clayton", "gumbel", "bb1", "bb6", "bb7")) ||
    (rot == 180 && fam %in% c("gumbel", "joe", "bb1", "bb6", "bb7"))
}

# A parameter block for a bicop_dist held somewhere in a bsicopula: get()
# extracts it, set(obj, cop) stores a replacement.
bicop_block <- function(label, cop, set) {
  info <- bicop_par_info(cop$family)
  list(
    label = label, names = info$names, value = as.numeric(cop$parameters),
    lower = info$lower, upper = info$upper,
    set = function(obj, value) {
      set(obj, rvinecopulib::bicop_dist(cop$family, cop$rotation, value))
    }
  )
}

udp_block <- function(label, x, set) {
  fp <- udp_fitpars(x)
  list(
    label = label, names = names(fp$value), value = unname(fp$value),
    lower = fp$lower, upper = fp$upper,
    set = function(obj, value) set(obj, udp_setfitpars(x, value))
  )
}

# The parameter blocks to estimate, dropping any with no free parameters.
fit_blocks <- function(object, udpfix, sdvinesimple) {
  blocks <- list(bicop_block("basecopula", object@basecopula, function(obj, cop) {
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
    blocks <- c(blocks, list(bicop_block("copZ1Z2_V1V2", rm@copZ1Z2_V1V2, function(obj, cop) {
      obj@randomizermod@copZ1Z2_V1V2 <- cop
      obj
    })))
    if (!sdvinesimple) {
      blocks <- c(blocks, list(
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
# with Hessian-based standard errors when 'hessian' is TRUE.
fit_stage <- function(U, object, blocks, vfloor, hessian, method, control) {
  sizes <- vapply(blocks, function(b) length(b$value), integer(1))
  idx <- rep(seq_along(blocks), sizes)
  labels <- unlist(lapply(blocks, function(b) paste(b$label, b$names, sep = ".")))
  lower <- unlist(lapply(blocks, `[[`, "lower"))
  upper <- unlist(lapply(blocks, `[[`, "upper"))
  start <- nudge_start(unlist(lapply(blocks, `[[`, "value")), lower, upper)
  npar <- length(start)

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
      object = object, estimate = numeric(0), se = numeric(0),
      vcov = matrix(numeric(0), 0, 0), loglik = -bsicopula_negll(U, object, vfloor),
      npar = 0L, convergence = 0L, message = "no free parameters",
      counts = c(`function` = NA_integer_, gradient = NA_integer_), method = NA_character_
    ))
  }

  if (is.null(method)) method <- if (npar == 1L) "BFGS" else "Nelder-Mead"
  if (method == "Nelder-Mead" && is.null(control$maxit)) control$maxit <- 2000
  fn <- function(th) negll_nat(from_free(th, lower, upper))
  opt <- stats::optim(to_free(start, lower, upper), fn, method = method, control = control)
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
  est <- stats::setNames(from_free(opt$par, lower, upper), labels)

  V <- matrix(NA_real_, npar, npar, dimnames = list(labels, labels))
  sds <- stats::setNames(rep(NA_real_, npar), labels)
  if (hessian) {
    H <- stats::optimHess(unname(est), negll_nat)
    Vinv <- tryCatch(solve(H), error = function(e) NULL)
    if (is.null(Vinv) || any(diag(Vinv) <= 0)) {
      warning("Hessian is not positive definite; standard errors are NA.", call. = FALSE)
    } else {
      V[] <- Vinv
      sds[] <- sqrt(diag(Vinv))
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
#' @slot npar number of estimated parameters.
#' @slot convergence the [stats::optim()] convergence code (`0` for success).
#' @slot message the [stats::optim()] message, if any.
#' @slot counts the [stats::optim()] function/gradient evaluation counts.
#' @slot method the [stats::optim()] method used.
#' @slot vfloor the floor applied to the carrier values before evaluating the
#'   base copula density, or `NA` if none was applied.
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
    matrix(numeric(0), 0, res$npar, dimnames = list(NULL, names(res$estimate)))
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
#' pseudo-observations by maximizing the likelihood built from
#' [dbsicopula()]. `object` supplies both the model structure -- the
#' copula families and rotations, the udp classes, the randomizer model --
#' and the starting values for the optimization.
#'
#' The base copula (and every copula of a \linkS4class{randsdvine}) must be
#' a parametric `bicop_dist` object from \pkg{rvinecopulib}; its family and
#' rotation are held fixed and only its parameters are estimated.
#' `randomizermod` must be `NULL` or a \linkS4class{randsdvine}.
#'
#' **udp parameters.** With `udpfix = FALSE`, the continuous parameters of
#' `udp1` and `udp2` are estimated alongside the copula parameters. At
#' present that is supported for v-transforms (`delta`, `kappa`, `xi`).
#' udps with no continuous parameters -- [vsymmetric()], shuffles such as
#' [udpid()] and [udpflip()], and \linkS4class{udpcosine} or
#' \linkS4class{udplegendre} objects, whose degree is regarded as fixed --
#' simply contribute nothing. For other classes use `udpfix = TRUE`, which
#' holds both udps at their given values: the natural choice when they have
#' been determined externally, e.g. by [aceshuffle()].
#'
#' **Randomizer.** When `randomizermod` is a \linkS4class{randsdvine},
#' `sdvinesimple = TRUE` estimates only its tree-3 copula `copZ1Z2_V1V2`,
#' leaving the two tree-2 copulas at their given values (by default the
#' independence copula); `sdvinesimple = FALSE` estimates all three. With
#' `twostage = TRUE`, the model is first fitted with independent
#' randomizers, and that fit's base copula and udps are used as starting
#' values for the full fit.
#'
#' **Boundary floor.** For many udps an interior point (a v-transform's
#' fulcrum `delta`, for instance) maps to `0`. When the base copula's
#' density is unbounded at `(0, 0)`, the likelihood becomes unbounded as
#' such points approach an observation in both margins at once, and even a
#' slow blow-up is enough to trap the optimizer there. `vfloor = "auto"`
#' guards against this by flooring the carrier values at `1 / (2n)` before
#' evaluating the base copula density, but only for base copulas that need
#' it: the Gaussian (for either sign of the correlation) and t copulas;
#' Clayton, Gumbel, BB1, BB6 and BB7 (rotation 0); and survival Gumbel,
#' Joe, BB1, BB6 and BB7 (rotation 180). `vfloor = TRUE` or `FALSE` forces
#' the floor on or off. The floor enters only the fitting; [dbsicopula()]
#' is exact.
#'
#' **Optimization and standard errors.** Parameters are optimized with
#' [stats::optim()] on an internal unconstrained scale; every result is
#' reported on the natural scale. The likelihood is only piecewise smooth
#' in a udp's breakpoint parameters (a v-transform's `delta`): it has a
#' kink or cusp wherever the breakpoint crosses an observation, and for many
#' base copulas the profile likelihood in `delta` is a fine sawtooth. The
#' optimizer may stop on a neighbouring tooth, which moves `delta` by about
#' `1 / n` and the log-likelihood by a unit or two; Nelder-Mead is restarted
#' automatically when its simplex collapses on a kink.
#'
#' `se = "bootstrap"` (or `TRUE`) gives parametric bootstrap standard
#' errors: `B` samples of size `n` are drawn from the fitted model with
#' [rbsicopula()], converted to pseudo-observations, and refitted, starting
#' from the fitted values. This accounts for the rank transformation and is
#' the recommended choice. `se = "hessian"` is much quicker, using the
#' numerically differentiated Hessian of the log-likelihood on the natural
#' scale, but it ignores the rank transformation. It is also unreliable for
#' `delta`: when the base copula density vanishes along the edge `V = 0`
#' (Gaussian with positive correlation, Clayton, Gumbel, ...), the
#' likelihood has a cusp at every observation and the Hessian can understate
#' the standard error of `delta` many times over; when the edge density is
#' finite and nonzero (Joe, survival Clayton), it still understates it,
#' by a smaller factor.
#'
#' @param U a two-column numeric matrix (or data frame) of pseudo-observations,
#'   with every value strictly inside `(0, 1)`.
#' @param object a \linkS4class{bsicopula} giving the model structure and
#'   starting values.
#' @param udpfix logical; hold the parameters of `udp1` and `udp2` fixed at
#'   their given values?
#' @param sdvinesimple logical; when `randomizermod` is a
#'   \linkS4class{randsdvine}, estimate only `copZ1Z2_V1V2`?
#' @param twostage logical; when `randomizermod` is a
#'   \linkS4class{randsdvine}, obtain starting values from a first fit with
#'   independent randomizers? Ignored otherwise.
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
#' @include bsicopula.R
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
fitbsicopula <- function(U, object, udpfix = FALSE, sdvinesimple = TRUE, twostage = TRUE,
                         se = FALSE, B = 200, vfloor = "auto", method = NULL,
                         control = list()) {
  if (!methods::is(object, "bsicopula")) {
    stop("'object' must be an object of class 'bsicopula'.", call. = FALSE)
  }
  U <- as.matrix(U)
  if (!is.numeric(U) || ncol(U) != 2L || nrow(U) < 2L || anyNA(U) || any(U <= 0 | U >= 1)) {
    stop(
      "'U' must be a two-column numeric matrix of pseudo-observations strictly inside (0, 1).",
      call. = FALSE
    )
  }
  for (arg in c("udpfix", "sdvinesimple", "twostage")) {
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
  if (!is_bicop_dist(object@basecopula)) {
    stop("fitbsicopula() requires the base copula to be a bicop_dist object (rvinecopulib).",
      call. = FALSE
    )
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

  stage1 <- NULL
  if (!is.null(rm) && twostage) {
    object1 <- object
    object1@randomizermod <- NULL
    res1 <- fit_stage(U, object1, fit_blocks(object1, udpfix, sdvinesimple),
      floor_value, FALSE, method, control
    )
    stage1 <- new_fitbsicopula(res1, n, floor_value)
    object@basecopula <- res1$object@basecopula
    object@udp1 <- res1$object@udp1
    object@udp2 <- res1$object@udp2
  }

  res <- fit_stage(U, object, fit_blocks(object, udpfix, sdvinesimple),
    floor_value, se_method == "hessian", method, control
  )
  if (se_method == "bootstrap" && res$npar > 0L) {
    res <- bootstrap_se(res, n, B, udpfix, sdvinesimple, floor_value, method, control)
  }
  new_fitbsicopula(res, n, floor_value, se_method, stage1)
}

# Parametric bootstrap: B samples of size n from the fitted model, each
# ranked to pseudo-observations (as the real data were) and refitted from
# the fitted values -- which are close enough that a randsdvine model needs
# no first stage. Returns 'res' with se, vcov and boot filled in. A refit
# that errors leaves a row of NAs and is left out of se/vcov.
bootstrap_se <- function(res, n, B, udpfix, sdvinesimple, vfloor, method, control) {
  fitted <- res$object
  blocks <- fit_blocks(fitted, udpfix, sdvinesimple)
  est <- matrix(NA_real_, B, res$npar, dimnames = list(NULL, names(res$estimate)))
  for (b in seq_len(B)) {
    Ub <- apply(rbsicopula(n, fitted), 2, rank) / (n + 1)
    rb <- tryCatch(
      suppressWarnings(fit_stage(Ub, fitted, blocks, vfloor, FALSE, method, control)),
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
  cat("base copula: ", bc@basecopula$family, " (rotation ", bc@basecopula$rotation, ")\n",
    sep = ""
  )
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
    cat("carrier values floored at ", format(object@vfloor, digits = 4),
      " in the base copula density\n",
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
