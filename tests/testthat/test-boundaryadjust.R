# boundaryadjust(): nudges exact 0/1 values strictly inside (0, 1).

test_that("boundaryadjust() leaves interior values untouched", {
  x <- c(0.1, 0.5, 0.9)
  expect_equal(boundaryadjust(x), x)
})

test_that("boundaryadjust() replaces exact 0/1 with the default tolerance", {
  x <- c(0, 0.3, 1)
  out <- boundaryadjust(x)
  expect_equal(out[2], 0.3)
  expect_equal(out[1], 1 / (2 * length(x)))
  expect_equal(out[3], 1 - 1 / (2 * length(x)))
})

test_that("boundaryadjust() honours a user-supplied btol", {
  x <- c(0, 0.3, 1)
  out <- boundaryadjust(x, btol = 0.01)
  expect_equal(out, c(0.01, 0.3, 0.99))
})

test_that("boundaryadjust() default tolerance for a matrix uses nrow(), not length()", {
  m <- matrix(c(0, 0.4, 0.6, 1), nrow = 2, ncol = 2)
  out <- boundaryadjust(m)
  expect_equal(out[1, 1], 1 / (2 * nrow(m)))
  expect_equal(out[2, 2], 1 - 1 / (2 * nrow(m)))
  expect_equal(out[2, 1], 0.4)
  expect_equal(out[1, 2], 0.6)
  expect_true(dim(out)[1] == nrow(m) && dim(out)[2] == ncol(m))
})
