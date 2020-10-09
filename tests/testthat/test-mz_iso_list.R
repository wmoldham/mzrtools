test_that("output format correct", {
  expect_equal(ncol(mz_iso_list("C3H3N3O3")), 9)
  expect_equal(nrow(mz_iso_list("C3H3N3O3")), 2160)
})

test_that("non-isotope molecules return a value", {
  expect_equal(mz_iso_list("Cu4"), NA)
})

