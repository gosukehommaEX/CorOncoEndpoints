# Test correlation-related functions

# Test CorBoundResponseTTE ----
test_that("CorBoundResponseTTE returns valid bounds", {
  bounds <- CorBoundResponseTTE(p = 0.4)

  # Should return a vector of length 2
  expect_length(bounds, 2)

  # Lower bound should be less than upper bound
  expect_true(bounds[1] < bounds[2])

  # Bounds should be within [-1, 1]
  expect_true(bounds[1] >= -1 && bounds[1] <= 1)
  expect_true(bounds[2] >= -1 && bounds[2] <= 1)
})

test_that("CorBoundResponseTTE validates input", {
  # Invalid p (outside (0, 1))
  expect_error(
    CorBoundResponseTTE(p = 0),
    "'p' must be a single numeric value between 0 and 1 \\(exclusive\\)"
  )

  expect_error(
    CorBoundResponseTTE(p = 1),
    "'p' must be a single numeric value between 0 and 1 \\(exclusive\\)"
  )

  expect_error(
    CorBoundResponseTTE(p = 1.5),
    "'p' must be a single numeric value between 0 and 1 \\(exclusive\\)"
  )
})

test_that("CorBoundResponseTTE is symmetric around p = 0.5", {
  bounds1 <- CorBoundResponseTTE(p = 0.3)
  bounds2 <- CorBoundResponseTTE(p = 0.7)

  # Due to symmetry, bounds should have similar absolute values
  expect_equal(abs(bounds1[1]), abs(bounds2[2]), tolerance = 1e-10)
  expect_equal(abs(bounds1[2]), abs(bounds2[1]), tolerance = 1e-10)
})

# Test CorBoundResponsePFS ----
test_that("CorBoundResponsePFS returns valid bounds", {
  bounds <- CorBoundResponsePFS(
    p = 0.4,
    hazard_OS = 0.05,
    hazard_PFS = 0.08
  )

  # Should return a vector of length 2
  expect_length(bounds, 2)

  # Lower bound should be less than upper bound
  expect_true(bounds[1] < bounds[2])

  # Bounds should be within [-1, 1]
  expect_true(bounds[1] >= -1 && bounds[1] <= 1)
  expect_true(bounds[2] >= -1 && bounds[2] <= 1)
})

test_that("CorBoundResponsePFS validates hazard constraint", {
  # hazard_PFS must be > hazard_OS
  expect_error(
    CorBoundResponsePFS(p = 0.4, hazard_OS = 0.08, hazard_PFS = 0.05),
    "hazard_PFS.*must be strictly greater than.*hazard_OS"
  )
})

test_that("CorBoundResponsePFS bounds are narrower than general bounds", {
  p <- 0.4
  hazard_OS <- 0.05
  hazard_PFS <- 0.08

  general_bounds <- CorBoundResponseTTE(p = p)
  pfs_bounds <- CorBoundResponsePFS(p = p, hazard_OS = hazard_OS,
                                    hazard_PFS = hazard_PFS)

  # PFS-specific bounds should generally be narrower
  # (though this may not always be strictly true in all cases)
  expect_true(pfs_bounds[1] >= general_bounds[1] - 1e-10)
  expect_true(pfs_bounds[2] <= general_bounds[2] + 1e-10)
})

# Test CorResponsePFS ----
test_that("CorResponsePFS returns valid correlation", {
  rho_pfs <- CorResponsePFS(
    p = 0.4,
    hazard_OS = 0.05,
    hazard_PFS = 0.08,
    rho_OS_Response = 0.3,
    copula = "Clayton"
  )

  # Should return a single numeric value
  expect_length(rho_pfs, 1)
  expect_type(rho_pfs, "double")

  # Should be within [-1, 1]
  expect_true(rho_pfs >= -1 && rho_pfs <= 1)
})

test_that("CorResponsePFS validates inputs", {
  # Invalid rho_OS_Response (outside feasible range)
  bounds <- CorBoundResponseTTE(p = 0.4)

  expect_error(
    CorResponsePFS(
      p = 0.4,
      hazard_OS = 0.05,
      hazard_PFS = 0.08,
      rho_OS_Response = bounds[1] - 0.01,
      copula = "Clayton"
    ),
    "rho_OS_Response.*must be within the feasible range"
  )

  # Clayton copula cannot handle negative dependence
  expect_error(
    CorResponsePFS(
      p = 0.4,
      hazard_OS = 0.05,
      hazard_PFS = 0.08,
      rho_OS_Response = -0.2,
      copula = "Clayton"
    ),
    "Clayton copula cannot model negative dependence"
  )
})

test_that("CorResponsePFS is within PFS bounds", {
  p <- 0.4
  hazard_OS <- 0.05
  hazard_PFS <- 0.08
  rho_OS <- 0.3

  rho_pfs <- CorResponsePFS(
    p = p,
    hazard_OS = hazard_OS,
    hazard_PFS = hazard_PFS,
    rho_OS_Response = rho_OS,
    copula = "Clayton"
  )

  pfs_bounds <- CorBoundResponsePFS(
    p = p,
    hazard_OS = hazard_OS,
    hazard_PFS = hazard_PFS
  )

  # Calculated correlation should be within bounds
  expect_true(rho_pfs >= pfs_bounds[1])
  expect_true(rho_pfs <= pfs_bounds[2])
})

# Test CopulaParamResponseTTE ----
test_that("CopulaParamResponseTTE returns valid theta", {
  theta <- CopulaParamResponseTTE(
    p = 0.4,
    rho = 0.3,
    copula = "Clayton"
  )

  # Should return a single numeric value
  expect_length(theta, 1)
  expect_type(theta, "double")

  # For Clayton, theta should be positive
  expect_true(theta > 0)
})

test_that("CopulaParamResponseTTE validates rho range", {
  bounds <- CorBoundResponseTTE(p = 0.4)

  # rho outside feasible range
  expect_error(
    CopulaParamResponseTTE(p = 0.4, rho = bounds[1] - 0.01, copula = "Clayton"),
    "The value of rho must be in"
  )

  expect_error(
    CopulaParamResponseTTE(p = 0.4, rho = bounds[2] + 0.01, copula = "Clayton"),
    "The value of rho must be in"
  )
})

test_that("CopulaParamResponseTTE works with Frank copula", {
  theta_frank <- CopulaParamResponseTTE(
    p = 0.4,
    rho = 0.3,
    copula = "Frank"
  )

  expect_length(theta_frank, 1)
  expect_type(theta_frank, "double")

  # Frank copula can have positive theta for positive correlation
  expect_true(theta_frank > 0)
})

test_that("CopulaParamResponseTTE handles negative correlation with Frank", {
  theta_frank_neg <- CopulaParamResponseTTE(
    p = 0.4,
    rho = -0.2,
    copula = "Frank"
  )

  # Frank copula with negative correlation should have negative theta
  expect_true(theta_frank_neg < 0)
})
