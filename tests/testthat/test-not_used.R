test_that("arguments that are not used give an error naming the first one", {
  expect_null(.not_used(a = 1, cond = FALSE, why = "here"))
  expect_null(.not_used(a = NULL, b = NULL, cond = TRUE, why = "here"))
  expect_error(.not_used(a = NULL, b = 2, c = 3, cond = TRUE, why = "without X"),
               "'b' is not used without X")
  # an argument named like the start of 'cond' is checked, not taken for 'cond'
  expect_error(.not_used(co = 1, cond = TRUE, why = "here"), "'co' is not used here")
})
