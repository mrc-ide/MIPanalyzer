context("test-main")

test_that("dummy1 works", {
  expect_equal(dummy1(-5:5), (-5:5)^2)
  expect_equal(dummy1(), dummy1(1:5))
  expect_error(dummy1(NULL))
})
