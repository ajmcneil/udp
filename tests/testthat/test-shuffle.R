test_that("shuffle() validates perm and signs", {
  expect_s4_class(shuffle(c(2, 1, 3)), "shuffle")
  expect_error(shuffle(c(1, 2, 2)), "permutation")
  expect_error(shuffle(c(1, 2, 4)), "permutation")
  expect_error(shuffle(c(2, 1, 3), signs = c(1, -1)), "same length")
  expect_error(shuffle(c(2, 1, 3), signs = c(1, 0, 1)), "1 and -1")
})

test_that("shuffle() defaults to all-positive signs and stores an integer perm", {
  s <- shuffle(c(2, 1, 3))
  expect_identical(s@signs, c(1, 1, 1))
  expect_type(s@perm, "integer")
})

test_that("shtrans() matches the piecewise definition", {
  s <- shuffle(c(2, 1, 3)) # swap the first two thirds, third fixed
  expect_equal(shtrans(s, c(0.1, 0.2)), c(0.1, 0.2) + 1 / 3)
  expect_equal(shtrans(s, c(0.4, 0.6)), c(0.4, 0.6) - 1 / 3)
  expect_equal(shtrans(s, c(0.7, 0.9)), c(0.7, 0.9))

  sf <- shuffle(c(1, 2, 3), signs = c(-1, 1, 1)) # flip the first third
  expect_equal(shtrans(sf, c(0.1, 0.25)), 1 / 3 - c(0.1, 0.25))

  expect_equal(shtrans(s, c(0, 1)), c(1 / 3, 1))
})

test_that("shtrans() maps [0, 1] onto [0, 1] and preserves attributes", {
  set.seed(1)
  for (rep in seq_len(20)) {
    m <- sample.int(8, 1)
    s <- shuffle(sample.int(m), signs = sample(c(-1, 1), m, replace = TRUE))
    y <- shtrans(s, seq(0, 1, length.out = 101))
    expect_true(all(y >= -1e-12 & y <= 1 + 1e-12))
  }

  u <- ts(seq(0.05, 0.95, length.out = 12), frequency = 4)
  expect_s3_class(shtrans(shuffle(c(2, 1)), u), "ts")
})

test_that("shinverse() inverts shtrans() in both directions", {
  set.seed(2)
  u <- sort(runif(300))
  for (rep in seq_len(25)) {
    m <- sample(2:10, 1)
    s <- shuffle(sample.int(m), signs = sample(c(-1, 1), m, replace = TRUE))
    expect_equal(shinverse(s, shtrans(s, u)), u, tolerance = 1e-9)
    expect_equal(shtrans(s, shinverse(s, u)), u, tolerance = 1e-9)
  }
})

test_that("the inverse of a shuffle is the shuffle of the inverse permutation", {
  s <- shuffle(c(3, 1, 2), signs = c(1, -1, 1))
  o <- order(s@perm)
  inv <- shuffle(o, signs = s@signs[o])
  v <- sort(runif(101))
  expect_equal(shinverse(s, v), shtrans(inv, v))
})

test_that("plot() runs for a shuffle", {
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit({
    dev.off()
    unlink(f)
  })
  expect_no_error(plot(shuffle(c(3, 1, 2), signs = c(1, -1, 1))))
  expect_no_error(plot(shuffle(1:5)))
  expect_no_error(plot(shuffle(c(2, 1)), main = "swap"))
})

test_that("shuffles preserve the uniform distribution", {
  skip_on_cran()
  set.seed(20240908)
  u <- runif(1e5)
  for (perm in list(c(2, 1, 3), c(4, 2, 1, 3), c(1, 2, 3, 4, 5))) {
    m <- length(perm)
    s <- shuffle(perm, signs = rep(c(1, -1), length.out = m))
    p <- suppressWarnings(stats::ks.test(shtrans(s, u), "punif")$p.value)
    expect_gte(p, 0.01)
  }
})
