library(testthat)
source(testthat::test_path("..", "..", "MI_helper_functions_Semi.R"))

test_that("multipliers modify imputed values", {
  data <- data.frame(
    y = c(1, 2, 3),
    z = c(1, 2, 3),
    ps1 = c(NA, NA, NA),
    mult_ps1 = c(10, 20, 30)
  )
  imp <- mice::mice(data, m = 2, maxit = 0, printFlag = FALSE)
  exposures <- list(ps1 = list(dependent = "y", independent = "z"))
  imp <- ps_on_mice_imputation(imp, exposures, multipliers = "mult_ps1")
  expect_equal(imp$imp$ps1[, 1], c(10, 40, 90))
  expect_equal(imp$imp$ps1[, 2], c(10, 40, 90))
})
