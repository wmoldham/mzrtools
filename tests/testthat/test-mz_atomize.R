test_that("molecular formula is correct", {
  expect_error(mz_atomize("Z1"))
  expect_error(mz_atomize("C1C1C1"))
})

test_that("output is correct", {
  expect_equal(mz_atomize("C14"), list(C = 14))
  expect_equal(mz_atomize("C"), list(C = 1))
})
