test_that("rZOIB_4p returns array with correct dimensions and finite values", {

  set.seed(123)
  n_mc <- 100
  alpha_point <- runif(1, min = 2, max = 10)
  beta_point <- runif(1, min = 2, max = 10)
  zoi_point <- runif(1, 0.1, 0.3)
  coi_point <- runif(1, 0.2, 0.8)
  weights <- runif(1, 0.5, 2)

  Y <- rZOIB_4p(n_mc = n_mc,
                alpha_point = alpha_point,
                beta_point = beta_point,
                zoi_point = zoi_point,
                coi_point = coi_point,
                weight = weights)

  expect_type(Y, "double")
  expect_true(is.vector(Y))
  expect_equal(length(Y), n_mc)
  expect_false(anyNA(Y))
  expect_true(all(is.finite(Y)))
  expect_true(all(Y >= 0))
  expect_true(all(Y <= max(weights)))
})

test_that("rZOIB_4p produces only zero/one inflation when zoi = 1", {

  set.seed(123)
  n_mc <- 100
  alpha_point <- runif(1, min = 2, max = 10)
  beta_point <- runif(1, min = 2, max = 10)
  zoi_point <- 1
  coi_point <- runif(1, 0.2, 0.8)
  weights <- runif(1, 0.5, 2)

  Y <- rZOIB_4p(n_mc = n_mc,
                alpha_point = alpha_point,
                beta_point = beta_point,
                zoi_point = zoi_point,
                coi_point = coi_point,
                weight = weights)

  # values must be exactly 0 or weights
  allowed <- c(0, weights)

  expect_true(all(Y %in% allowed))
})

test_that("rZOIB_4p reduces to Beta sampler when zoi = 0", {

  set.seed(123)
  n_mc <- 100
  alpha_point <- runif(1, min = 2, max = 10)
  beta_point <- runif(1, min = 2, max = 10)
  zoi_point <- 0
  coi_point <- runif(1, 0.2, 0.8)
  weights <- runif(1, 0.5, 2)

  Y <- rZOIB_4p(n_mc = n_mc,
                alpha_point = alpha_point,
                beta_point = beta_point,
                zoi_point = zoi_point,
                coi_point = coi_point,
                weight = weights)

  # no point masses at 0 or 1 expected
  expect_false(any(Y == 0))
  expect_false(any(Y == 1))
})

