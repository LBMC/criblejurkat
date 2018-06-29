library(CribleJurkat)
context("Run CribleJurkat analysis")

test_that("loading data", {
  expect_identical(
    paste0("gene", 1:100),
    paste0("gene", 1:100))
  expect_equal(
    10:6,
    10:6)
})
