test_that("udpzigzag() validates its arguments", {
  expect_error(udpzigzag(), "exactly one of")
  expect_error(udpzigzag(breaks = 0.5, widths = c(1, 1)), "exactly one of")

  expect_error(udpzigzag(breaks = c(0.5, 0.5)), "strictly increasing")
  expect_error(udpzigzag(breaks = c(0.6, 0.4)), "strictly increasing")
  expect_error(udpzigzag(breaks = 0), "strictly inside")
  expect_error(udpzigzag(breaks = 1), "strictly inside")
  expect_error(udpzigzag(breaks = NA_real_), "strictly increasing")

  expect_error(udpzigzag(widths = numeric(0)), "positive")
  expect_error(udpzigzag(widths = c(1, -1)), "positive")
  expect_error(udpzigzag(widths = c(1, 0)), "positive")
  expect_error(udpzigzag(widths = c(1, NA)), "positive")
  expect_error(udpzigzag(widths = c(1, Inf)), "positive")

  expect_error(udpzigzag(breaks = 0.5, up = NA), "single logical")
  expect_error(udpzigzag(breaks = 0.5, up = c(TRUE, FALSE)), "single logical")
  expect_error(udpzigzag(breaks = 0.5, up = "yes"), "single logical")
})

test_that("udpzigzag() stores the full breakpoint vector from either input", {
  x <- udpzigzag(breaks = c(0.3, 0.5))
  expect_s4_class(x, "udpzigzag")
  expect_equal(x@breaks, c(0, 0.3, 0.5, 1))
  expect_true(x@up)

  y <- udpzigzag(widths = c(1, 2, 1), up = FALSE)
  expect_equal(y@breaks, c(0, 0.25, 0.75, 1))
  expect_false(y@up)

  # widths need not already sum to 1
  z <- udpzigzag(widths = c(2, 2, 4))
  expect_equal(z@breaks, c(0, 0.25, 0.5, 1))

  # a single piece
  single <- udpzigzag(breaks = numeric(0))
  expect_equal(single@breaks, c(0, 1))
})

test_that("udptrans() matches known closed forms", {
  u <- seq(0, 1, length.out = 101)

  # a single piece reduces exactly to udpid() / udpflip()
  expect_equal(udptrans(udpzigzag(breaks = numeric(0), up = TRUE), u), udptrans(udpid(), u))
  expect_equal(udptrans(udpzigzag(breaks = numeric(0), up = FALSE), u), udptrans(udpflip(), u))

  # equal widths reduce to udpcosine(n), matching its own starting direction
  for (n in 1:6) {
    up <- (n %% 2L) == 1L
    expect_equal(udptrans(udpzigzag(widths = rep(1, n), up = up), u),
      udptrans(udpcosine(n), u),
      tolerance = 1e-10
    )
  }

  # two pieces reduce exactly to vlinear(delta)
  for (delta in c(0.2, 0.4, 0.5, 0.75)) {
    x <- udpzigzag(widths = c(delta, 1 - delta), up = FALSE)
    expect_equal(udptrans(x, u), udptrans(vlinear(delta), u), tolerance = 1e-10)
  }
})

test_that("udptrans() maps [0, 1] into [0, 1]", {
  u <- seq(0, 1, length.out = 2001)
  for (w in list(1, c(1, 1), c(1, 2, 1), c(3, 1, 1, 4))) {
    for (up in c(TRUE, FALSE)) {
      y <- udptrans(udpzigzag(widths = w, up = up), u)
      expect_false(anyNA(y))
      expect_true(all(y >= 0 & y <= 1))
    }
  }
})

test_that("udptrans() is piecewise linear with slope of magnitude 1 / width", {
  x <- udpzigzag(widths = c(1, 2, 1), up = TRUE)
  b <- x@breaks
  w <- diff(b)
  u1 <- b[-length(b)] + 0.25 * w
  u2 <- b[-length(b)] + 0.75 * w
  slope <- (udptrans(x, u2) - udptrans(x, u1)) / (u2 - u1)
  expect_equal(abs(slope), 1 / w, tolerance = 1e-8)
})

test_that("udptrans() preserves the shape and attributes of u", {
  u <- ts(seq(0.05, 0.95, length.out = 24), frequency = 4)
  y <- udptrans(udpzigzag(widths = c(1, 2, 1)), u)
  expect_s3_class(y, "ts")
  expect_identical(tsp(y), tsp(u))
})

test_that("zigzag udp transformations preserve the uniform distribution", {
  skip_on_cran()
  set.seed(20240908)
  u <- runif(1e5)
  for (w in list(1, c(1, 1), c(1, 2, 1), c(3, 1, 1, 4))) {
    x <- udpzigzag(widths = w)
    p <- suppressWarnings(stats::ks.test(udptrans(x, u), "punif")$p.value)
    expect_gte(p, 0.01)
  }
})

test_that("udpzigzaginverse() returns one root per piece, weighted by width", {
  x <- udpzigzag(widths = c(1, 2, 1), up = TRUE) # breaks c(0, .25, .75, 1)
  res <- udpzigzaginverse(x, c(0, 0.4, 1))

  expect_equal(res[[1]]$u, c(0, 0.75))
  expect_equal(res[[1]]$w, c(0.25, 0.75))
  expect_equal(res[[2]]$u, c(0.1, 0.55, 0.85))
  expect_equal(res[[2]]$w, c(0.25, 0.5, 0.25))
  expect_equal(res[[3]]$u, c(0.25, 1))
  expect_equal(res[[3]]$w, c(0.75, 0.25))

  # every reported root really is a pre-image
  for (j in seq_along(res)) {
    r <- res[[j]]$u
    expect_equal(udptrans(x, r), rep(c(0, 0.4, 1)[j], length(r)), tolerance = 1e-9)
  }
  # each row's summed weights always add to 1
  expect_equal(vapply(res, function(r) sum(r$w), numeric(1)), c(1, 1, 1))
})

test_that("udpzigzaginverse() gives one root per piece for v in (0, 1)", {
  set.seed(1)
  for (w in list(1, c(1, 1), c(1, 2, 1), c(3, 1, 1, 4))) {
    x <- udpzigzag(widths = w)
    n <- length(w)
    v <- runif(30)
    res <- udpzigzaginverse(x, v)
    for (j in seq_along(v)) {
      expect_length(res[[j]]$u, n)
      expect_false(is.unsorted(res[[j]]$u))
      expect_true(all(res[[j]]$u >= 0 & res[[j]]$u <= 1))
      expect_equal(udptrans(x, res[[j]]$u), rep(v[j], n), tolerance = 1e-9)
    }
  }
})

test_that("udpzigzaginverse() validates v", {
  x <- udpzigzag(widths = c(1, 1))
  expect_error(udpzigzaginverse(x, 1.5), "in \\[0, 1\\]")
  expect_error(udpzigzaginverse(x, -0.01), "in \\[0, 1\\]")
  expect_error(udpzigzaginverse(x, c(0.5, NA)), "in \\[0, 1\\]")
})

test_that("udpinverse() packs roots with a prob attribute summing to 1", {
  x <- udpzigzag(widths = c(1, 2, 1), up = TRUE)
  v <- c(0, 0.2, 0.5, 0.8, 1)
  M <- udpinverse(x, v, prob = TRUE)

  expect_equal(dim(M), c(length(v), 3L))
  present <- !is.na(M)
  P <- attr(M, "prob")
  expect_equal(dim(P), dim(M))
  expect_equal(is.na(P), !present)
  expect_equal(rowSums(P, na.rm = TRUE), rep(1, length(v)))

  # every present entry really is a pre-image of its row's v
  expect_equal(
    udptrans(x, M[present]),
    rep(v, ncol(M))[as.vector(present)],
    tolerance = 1e-9
  )

  # without prob = TRUE (the default) there is no "prob" attribute
  M2 <- udpinverse(x, v)
  expect_null(attr(M2, "prob"))
})

test_that("udpinverse() without prob = TRUE matches udpzigzaginverse()'s roots", {
  x <- udpzigzag(widths = c(3, 1, 1, 4))
  v <- c(0, 0.15, 0.6, 1)
  M <- udpinverse(x, v)
  res <- udpzigzaginverse(x, v)
  for (j in seq_along(v)) {
    row <- M[j, ]
    expect_equal(row[!is.na(row)], res[[j]]$u)
  }
})

test_that("udpderiv() matches finite differences within a piece", {
  x <- udpzigzag(widths = c(1, 2, 1), up = TRUE)
  u <- c(0.1, 0.4, 0.6, 0.9)
  eps <- 1e-6
  fd <- (udptrans(x, u + eps) - udptrans(x, u - eps)) / (2 * eps)
  expect_equal(udpderiv(x, u), fd, tolerance = 1e-5)
})

test_that("udpderiv() takes the left piece at a kink and the right piece at u = 0", {
  x <- udpzigzag(widths = c(1, 2, 1), up = TRUE) # breaks c(0, .25, .75, 1)
  # at u = 0.25, left piece (1, increasing) has slope 1 / 0.25 = 4
  expect_equal(udpderiv(x, 0.25), 4)
  # at u = 0.75, left piece (2, decreasing) has slope -1 / 0.5 = -2
  expect_equal(udpderiv(x, 0.75), -2)
  # at u = 0, no left piece: the right piece (1, increasing) is used
  expect_equal(udpderiv(x, 0), 4)
  # at u = 1, the left piece (3, increasing) has slope 1 / 0.25 = 4
  expect_equal(udpderiv(x, 1), 4)
})

test_that("udpderiv() preserves the shape and attributes of u", {
  u <- ts(seq(0.05, 0.95, length.out = 24), frequency = 4)
  y <- udpderiv(udpzigzag(widths = c(1, 2, 1)), u)
  expect_s3_class(y, "ts")
  expect_identical(tsp(y), tsp(u))
})

test_that("udpbreaks() returns the full breakpoint vector", {
  x <- udpzigzag(breaks = c(0.3, 0.5))
  expect_identical(udpbreaks(x), x@breaks)
  expect_identical(udpbreaks(udpzigzag(breaks = numeric(0))), c(0, 1))
})

test_that("pcoincide() is sum(widths^2) in closed form", {
  for (w in list(1, c(1, 1), c(1, 2, 1), c(3, 1, 1, 4), c(0.2, 0.8))) {
    x <- udpzigzag(widths = w)
    expect_equal(pcoincide(x), sum((w / sum(w))^2))
  }
  # agrees with udpcosine's closed form at equal widths
  for (n in 1:6) {
    expect_equal(pcoincide(udpzigzag(widths = rep(1, n))), pcoincide(udpcosine(n)))
  }
  # agrees with vlinear's closed form at two pieces
  for (delta in c(0.2, 0.4, 0.5, 0.75)) {
    x <- udpzigzag(widths = c(delta, 1 - delta), up = FALSE)
    expect_equal(pcoincide(x), pcoincide(vlinear(delta)))
  }
})

test_that("the default udp method for pcoincide() integrates to the closed form", {
  general <- selectMethod("pcoincide", "udp")
  for (w in list(c(1, 1), c(1, 2, 1), c(3, 1, 1, 4))) {
    x <- udpzigzag(widths = w)
    expect_equal(general(x), pcoincide(x), tolerance = 1e-6)
  }
})

test_that("pcoincide() matches a direct Monte Carlo estimate", {
  skip_on_cran()
  set.seed(2024)
  x <- udpzigzag(widths = c(3, 1, 1, 4), up = TRUE)
  u <- runif(2e5)
  v <- udptrans(x, u)
  u_back <- udpsi(x, v)
  mc <- mean(abs(u_back - u) < 1e-9)
  expect_equal(pcoincide(x), mc, tolerance = 0.01)
})

test_that("udpsi() returns one pre-image per v in (0, 1)", {
  set.seed(1)
  x <- udpzigzag(widths = c(1, 2, 1))
  v <- runif(1000)
  Z <- runif(1000)
  u <- udpsi(x, v, Z)

  expect_length(u, 1000L)
  expect_true(all(u >= 0 & u <= 1))
  expect_equal(udptrans(x, u), v, tolerance = 1e-9)
})

test_that("udpsi() at v = 0 and v = 1 samples among the shared roots", {
  set.seed(1)
  x <- udpzigzag(widths = c(1, 2, 1), up = TRUE)
  v <- rep(c(0, 1), each = 40)
  u <- udpsi(x, v, runif(length(v)))
  expect_true(all(u >= 0 & u <= 1))
  expect_equal(udptrans(x, u), v, tolerance = 1e-9)
  expect_true(all(u[v == 0] %in% udpzigzaginverse(x, 0)[[1]]$u))
  expect_true(all(u[v == 1] %in% udpzigzaginverse(x, 1)[[1]]$u))
})

test_that("udpsi() validates its arguments and edge inputs", {
  x <- udpzigzag(widths = c(1, 2, 1))
  expect_error(udpsi(x, c(0.3, 0.4), Z = 0.5), "same length")
  expect_error(udpsi(x, 1.2, 0.5), "in \\[0, 1\\]")
  expect_error(udpsi(x, c(0.5, NA), c(0.1, 0.1)), "in \\[0, 1\\]")
  expect_length(udpsi(x, numeric(0), numeric(0)), 0L)
})

test_that("udpsi() is reproducible and preserves attributes", {
  x <- udpzigzag(widths = c(1, 2, 1))
  set.seed(7)
  a <- udpsi(x, seq(0.1, 0.9, by = 0.1))
  set.seed(7)
  b <- udpsi(x, seq(0.1, 0.9, by = 0.1))
  expect_identical(a, b)

  v <- ts(seq(0.05, 0.95, length.out = 24), frequency = 4)
  out <- udpsi(x, v)
  expect_s3_class(out, "ts")
  expect_identical(tsp(out), tsp(v))
})

test_that("the stochastic inverse recovers a uniform input", {
  skip_on_cran()
  set.seed(20240908)
  x <- udpzigzag(widths = c(1, 2, 1))
  u <- runif(1e5)
  back <- udpsi(x, udptrans(x, u), runif(1e5))
  expect_gte(suppressWarnings(stats::ks.test(back, "punif")$p.value), 0.01)
})

test_that("plot() runs for a udpzigzag", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  for (w in list(1, c(1, 1), c(1, 2, 1), c(3, 1, 1, 4))) {
    expect_no_error(plot(udpzigzag(widths = w)))
  }
  expect_no_error(plot(udpzigzag(widths = c(1, 2, 1)), main = "zigzag", ylab = "y"))
})
