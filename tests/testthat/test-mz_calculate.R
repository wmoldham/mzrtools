test_that("polarity is defined", {
  expect_warning(mz_calculate("C2", "fred"))
})

test_that("polarity works", {
  expect_equal(mz_calculate("C3"), 36)
  expect_equal(mz_calculate("C3", "positive"), 37.007276)
  expect_equal(mz_calculate("C3", "negative"), 34.992724)
})

test_that("math is correct", {
  expect_equal(mz_calculate("C3H4N2Si3"), 151.9682, tolerance = 1e-6)
})

