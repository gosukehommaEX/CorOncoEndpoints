# Test CheckSimResults function

test_that("CheckSimResults returns correct structure", {
  set.seed(123)
  sim_data <- rOncoEndpoints(
    nsim = 100,
    group = c("Treatment", "Control"),
    n = c(50, 50),
    p = c(0.4, 0.3),
    hazard_OS = c(0.05, 0.07),
    rho_tte_resp = c(0.3, 0.2),
    copula = "Clayton"
  )

  result <- CheckSimResults(
    dataset = sim_data,
    p = c(Treatment = 0.4, Control = 0.3),
    hazard_OS = c(Treatment = 0.05, Control = 0.07),
    rho_tte_resp = c(Treatment = 0.3, Control = 0.2),
    copula = "Clayton"
  )

  # Check that result is a tibble
  expect_s3_class(result, "tbl_df")

  # Check required columns exist
  required_cols <- c("Group", "Endpoint", "Empirical", "Theoretical",
                     "Bias", "Relative_Bias", "SE", "MSE", "RMSE")
  expect_true(all(required_cols %in% names(result)))
})

test_that("CheckSimResults validates input parameters", {
  set.seed(456)
  sim_data <- rOncoEndpoints(
    nsim = 100,
    group = "Treatment",
    n = 50,
    p = 0.4,
    hazard_OS = 0.05,
    rho_tte_resp = 0.3,
    copula = "Clayton"
  )

  # Missing required parameter p
  expect_error(
    CheckSimResults(
      dataset = sim_data,
      hazard_OS = c(Treatment = 0.05),
      rho_tte_resp = c(Treatment = 0.3),
      copula = "Clayton"
    ),
    "'p' must be provided when Response endpoint is present"
  )

  # Mismatched group names
  expect_error(
    CheckSimResults(
      dataset = sim_data,
      p = c(Control = 0.4),
      hazard_OS = c(Treatment = 0.05),
      rho_tte_resp = c(Treatment = 0.3),
      copula = "Clayton"
    ),
    "is missing values for groups"
  )
})

test_that("CheckSimResults handles all three endpoints", {
  set.seed(789)
  sim_data <- rOncoEndpoints(
    nsim = 500,
    group = "Experimental",
    n = 100,
    p = 0.5,
    hazard_OS = 0.04,
    hazard_PFS = 0.08,
    rho_tte_resp = 0.4,
    copula = "Frank"
  )

  result <- CheckSimResults(
    dataset = sim_data,
    p = c(Experimental = 0.5),
    hazard_OS = c(Experimental = 0.04),
    hazard_PFS = c(Experimental = 0.08),
    rho_tte_resp = c(Experimental = 0.4),
    copula = "Frank"
  )

  # Should have endpoints for OS, PFS, Response, and correlations
  expect_true("Response" %in% result$Endpoint)
  expect_true("OS_Median" %in% result$Endpoint)
  expect_true("PFS_Median" %in% result$Endpoint)
  expect_true("Cor_OS_PFS" %in% result$Endpoint)
  expect_true("Cor_OS_Response" %in% result$Endpoint)
  expect_true("Cor_PFS_Response" %in% result$Endpoint)
})

test_that("CheckSimResults calculates bias correctly", {
  set.seed(111)
  sim_data <- rOncoEndpoints(
    nsim = 1000,
    group = "Treatment",
    n = 200,
    hazard_OS = 0.05
  )

  result <- CheckSimResults(
    dataset = sim_data,
    hazard_OS = c(Treatment = 0.05)
  )

  # Bias should be Empirical - Theoretical
  for (i in seq_len(nrow(result))) {
    expected_bias <- result$Empirical[i] - result$Theoretical[i]
    expect_equal(result$Bias[i], expected_bias, tolerance = 1e-10)
  }
})
