test_that("plot() runs for every panel type", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })

  for (x in list(Vsymmetric(), V2p(delta = 0.4, kappa = 1.2))) {
    for (type in c("transform", "inverse", "gradient", "pdown")) {
      expect_no_error(plot(x, type = type))
    }
  }
  expect_no_error(plot(V3p(delta = 0.45, kappa = 0.8, xi = 1.2), shading = FALSE))
})

test_that("plot() rejects an unknown type", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  expect_error(plot(Vsymmetric(), type = "nonsense"), "Not a plot method")
})
