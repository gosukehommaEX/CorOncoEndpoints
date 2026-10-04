test_that("two-point rule has nodes at -1/sqrt(3) and 1/sqrt(3)", {
  gl <- .gauss_legendre(2)
  expect_equal(gl$x, c(-1, 1) / sqrt(3), tolerance = 1e-12)
  expect_equal(gl$w, c(1, 1), tolerance = 1e-12)
})

test_that("n-point rule integrates polynomials of degree 2n - 1 exactly", {
  gl <- .gauss_legendre(10)
  expect_equal(sum(gl$w), 2, tolerance = 1e-12)
  expect_equal(sum(gl$w * gl$x ^ 18), 2 / 19, tolerance = 1e-12)
  expect_equal(sum(gl$w * gl$x ^ 19), 0, tolerance = 1e-12)
  expect_equal(sum(gl$w * exp(gl$x)), exp(1) - exp(-1), tolerance = 1e-12)
})
