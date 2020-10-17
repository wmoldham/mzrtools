test_that("warning works", {
  expect_warning(mz_iso_annotate("C3", polarity = "x"))
})

test_that("output is correct", {
  expect_named(mz_iso_annotate("C3"), expected = c("elements", "isotopes", "polarity", "iso_list"))
})
