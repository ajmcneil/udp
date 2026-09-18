# Every udp class's plot() method draws a [0, 1] x [0, 1] frame with
# asp = 1. On a device or panel that isn't exactly square, asp = 1 alone
# stretches one axis's displayed range past its data limits to satisfy the
# aspect ratio; pty = "s" (a square plotting region) is needed to prevent
# that. Use a deliberately wide, non-square device to catch a regression.

wide_png <- function() {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 700, height = 300)
  f
}

in_unit_interval <- function(usr) all(usr > -1e-6 & usr < 1 + 1e-6)

test_that("plot() keeps both axes within [0, 1] on a non-square device", {
  cases <- list(
    vsymmetric(),
    shuffle(c(3, 1, 2), signs = c(1, -1, 1)),
    udpcosine(4),
    udplegendre(3),
    udpcosinebex(c(1, 0.4)),
    udplegendrebex(c(1, 0.4)),
    udpzigzag(widths = c(1, 2, 1))
  )
  for (x in cases) {
    f <- wide_png()
    on.exit(unlink(f), add = TRUE)
    plot(x)
    expect_true(in_unit_interval(par("usr")), info = class(x))
    grDevices::dev.off()
  }
})

test_that("plot() of a vtransform's 'transform' type stays within [0, 1]", {
  f <- wide_png()
  on.exit(unlink(f), add = TRUE)
  plot(v2p(delta = 0.4, kappa = 0.8))
  expect_true(in_unit_interval(par("usr")))
  grDevices::dev.off()
})

test_that("plot() of a vtransform's other types keeps the x-axis within [0, 1]", {
  # 'inverse', 'pdown' and 'gradient' plot data with a natural range outside
  # [0, 1], so only their shared x-axis (v or u, always in [0, 1]) is checked.
  for (type in c("inverse", "pdown", "gradient")) {
    f <- wide_png()
    on.exit(unlink(f), add = TRUE)
    plot(v2p(delta = 0.4, kappa = 0.8), type = type)
    expect_true(in_unit_interval(par("usr")[1:2]), info = type)
    grDevices::dev.off()
  }
})
