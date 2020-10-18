test_that("polarity argument error works", {
  expect_error(mz_iso_annotate("C3", polarity = "X"))
})

test_that("output is correct", {
  expect_named(mz_iso_annotate("C3"), expected = c("elements", "isotopes", "polarity", "iso_list"))
})
