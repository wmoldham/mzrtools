test_that("molecular formula is correct", {
  expect_error(mz_atomize("Z1"))
  expect_error(mz_atomize("C1C1C1"))
})

test_that("output is correct", {
  expect_vector(mz_atomize("C14"))
  expect_named(mz_atomize("C14"), "C")
  expect_equal(mz_atomize("C"), 1 %>% rlang::set_names("C"))
})
