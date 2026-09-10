test_that("plot() runs for every panel type", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })

  for (x in list(vsymmetric(), v2p(delta = 0.4, kappa = 1.2))) {
    for (type in c("transform", "inverse", "gradient", "pdown")) {
      expect_no_error(plot(x, type = type))
    }
  }
  for (e in c("colour", "bw", "none")) {
    expect_no_error(plot(v3p(delta = 0.45, kappa = 0.8, xi = 1.2), embellish = e))
  }
  expect_error(plot(vsymmetric(), embellish = "grey"), "should be one of")
})

test_that("gridlines stay within the unit square", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  seg <- list()
  local_mocked_bindings(
    segments = function(x0, y0, x1, y1, ...) {
      seg[[length(seg) + 1L]] <<- c(x0, x1, y0, y1)
      invisible(NULL)
    },
    .package = "udp"
  )
  on.exit({
    dev.off()
    unlink(f)
  })

  # multi-panel layout with clipping disabled: the case where abline() spills
  op <- par(mfrow = c(2, 2), xpd = NA)
  on.exit(par(op), add = TRUE, after = FALSE)

  for (x in list(shuffle(c(3, 1, 2), signs = c(1, -1, 1)), udpcosine(5),
    udplegendre(6))) {
    seg <- list()
    plot(x)
    coords <- unlist(seg)
    expect_true(all(is.finite(coords)))
    expect_gte(min(coords), 0)
    expect_lte(max(coords), 1)

    # "none" draws no gridlines (the curve's own segments, if any, still count
    # for the shuffle, so only assert the count drops)
    seg <- list()
    plot(x, embellish = "none")
    n_none <- length(seg)
    seg <- list()
    plot(x, embellish = "colour")
    expect_gt(length(seg), n_none)
  }
})

test_that("plot() rejects an unknown type", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  expect_error(plot(vsymmetric(), type = "nonsense"), "Not a plot method")
})
