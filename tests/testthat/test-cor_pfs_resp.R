test_that("PFS and response correlation matches independent reference values", {
  # reference: inst/validation/python/reference.py, corr_pfs_r_std()
  expect_equal(.cor_pfs_resp(0.5, 0.3, -Inf)$cor, 0.370124648, tolerance = 1e-8)
  z_tau <- qnorm(1 - exp(-log(2) / 6 * 1.5))
  expect_equal(.cor_pfs_resp(0.3, 0.3, z_tau)$cor, 0.2772537111, tolerance = 1e-8)
})

test_that("small copula correlations give small non-zero correlations", {
  # the split point c_resp / theta lies far in the tail for these values
  expect_equal(.cor_pfs_resp(-0.001, 0.3, -Inf)$cor, -6.851620937e-4, tolerance = 1e-7)
  expect_equal(.cor_pfs_resp(0.0063, 0.6, -Inf)$cor, 4.484983549e-3, tolerance = 1e-7)
  expect_equal(.cor_pfs_resp(-0.0131, 0.3, -Inf)$cor, -8.956799692e-3, tolerance = 1e-7)
  expect_equal(.cor_pfs_resp(0, 0.3, -Inf)$cor, 0, tolerance = 1e-12)
})
