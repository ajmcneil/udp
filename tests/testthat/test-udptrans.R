test_that("udp is a virtual class that cannot be instantiated", {
  expect_true(isVirtualClass("udp"))
  expect_error(new("udp"))
})

test_that("every udp transformation inherits from the udp class", {
  objects <- list(
    vsymmetric(),
    vlinear(delta = 0.4),
    v2p(delta = 0.4, kappa = 1.2),
    v3b(delta = 0.3, kappa = 1.4, xi = 1.1),
    shuffle(c(2, 1, 3)),
    udpcosine(3)
  )
  for (x in objects) {
    expect_s4_class(x, "udp")
  }
})

test_that("udptrans() dispatches to each construction", {
  u <- seq(0, 1, by = 0.25)

  # v-transform
  expect_equal(udptrans(vsymmetric(), u), abs(2 * u - 1))

  # shuffle: the identity shuffle leaves u unchanged
  expect_equal(udptrans(shuffle(c(1, 2)), u), u)

  # cosine udp transformation: degree 2 is the tent map
  expect_equal(udptrans(udpcosine(2), u), abs(2 * u - 1))
})

test_that("udptrans() errors for an object with no method", {
  expect_error(udptrans(42, 0.5))
})
