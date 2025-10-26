# Test rOncoEndpoints function

test_that("rOncoEndpoints generates correct number of observations", {
  set.seed(123)
  data <- rOncoEndpoints(
    nsim = 10,
    group = c("Treatment", "Control"),
    n = c(50, 50),
    p = c(0.4, 0.3),
    hazard_OS = c(0.05, 0.07),
    rho_tte_resp = c(0.3, 0.2),
    copula = "Clayton"
  )

  # Check total number of rows
  expect_equal(nrow(data), 10 * (50 + 50))

  # Check number of simulations
  expect_equal(length(unique(data$simID)), 10)

  # Check groups
  expect_equal(unique(data$Group), c("Treatment", "Control"))
})

test_that("rOncoEndpoints enforces PFS <= OS constraint", {
  set.seed(456)
  data <- rOncoEndpoints(
    nsim = 100,
    group = "Treatment",
    n = 100,
    hazard_OS = 0.05,
    hazard_PFS = 0.08
  )

  # All PFS values should be <= OS values
  expect_true(all(data$PFS <= data$OS))
})

test_that("rOncoEndpoints generates only specified endpoints", {
  # OS only
  set.seed(789)
  data_os <- rOncoEndpoints(
    nsim = 10,
    n = 50,
    hazard_OS = 0.05
  )
  expect_true("OS" %in% names(data_os))
  expect_false("PFS" %in% names(data_os))
  expect_false("Response" %in% names(data_os))

  # PFS only
  set.seed(789)
  data_pfs <- rOncoEndpoints(
    nsim = 10,
    n = 50,
    hazard_PFS = 0.08
  )
  expect_false("OS" %in% names(data_pfs))
  expect_true("PFS" %in% names(data_pfs))
  expect_false("Response" %in% names(data_pfs))

  # Response only
  set.seed(789)
  data_resp <- rOncoEndpoints(
    nsim = 10,
    n = 50,
    p = 0.4
  )
  expect_false("OS" %in% names(data_resp))
  expect_false("PFS" %in% names(data_resp))
  expect_true("Response" %in% names(data_resp))
})

test_that("rOncoEndpoints validates input parameters", {
  # Invalid nsim
  expect_error(
    rOncoEndpoints(nsim = -1, n = 50, hazard_OS = 0.05),
    "'nsim' must be a single positive integer"
  )

  # Invalid n
  expect_error(
    rOncoEndpoints(nsim = 10, n = -50, hazard_OS = 0.05),
    "'n' must be a vector of positive integers"
  )

  # Invalid hazard_OS
  expect_error(
    rOncoEndpoints(nsim = 10, n = 50, hazard_OS = -0.05),
    "'hazard_OS' must be a numeric vector.*with all positive elements"
  )

  # Invalid p
  expect_error(
    rOncoEndpoints(nsim = 10, n = 50, p = 1.5),
    "'p' must be a numeric vector.*with all elements in \\(0, 1\\)"
  )

  # hazard_PFS <= hazard_OS
  expect_error(
    rOncoEndpoints(nsim = 10, n = 50, hazard_OS = 0.08, hazard_PFS = 0.05),
    "'hazard_PFS' must be strictly greater than 'hazard_OS'"
  )
})

test_that("rOncoEndpoints works with both Clayton and Frank copulas", {
  set.seed(111)
  data_clayton <- rOncoEndpoints(
    nsim = 50,
    n = 100,
    p = 0.4,
    hazard_OS = 0.05,
    rho_tte_resp = 0.3,
    copula = "Clayton"
  )

  set.seed(111)
  data_frank <- rOncoEndpoints(
    nsim = 50,
    n = 100,
    p = 0.4,
    hazard_OS = 0.05,
    rho_tte_resp = 0.3,
    copula = "Frank"
  )

  expect_true(all(c("OS", "Response") %in% names(data_clayton)))
  expect_true(all(c("OS", "Response") %in% names(data_frank)))
})

test_that("rOncoEndpoints generates Response as binary", {
  set.seed(222)
  data <- rOncoEndpoints(
    nsim = 100,
    n = 100,
    p = 0.4,
    hazard_OS = 0.05,
    rho_tte_resp = 0.3,
    copula = "Clayton"
  )

  # Response should only contain 0 and 1
  expect_true(all(data$Response %in% c(0, 1)))
})

test_that("rOncoEndpoints maintains group order", {
  set.seed(333)
  groups <- c("Arm A", "Arm B", "Arm C")
  data <- rOncoEndpoints(
    nsim = 10,
    group = groups,
    n = c(50, 60, 70),
    p = c(0.3, 0.4, 0.5),
    hazard_OS = c(0.05, 0.06, 0.07),
    rho_tte_resp = c(0.2, 0.3, 0.4),
    copula = "Clayton"
  )

  # Check that groups appear in correct order
  expect_equal(unique(data$Group), groups)
})
