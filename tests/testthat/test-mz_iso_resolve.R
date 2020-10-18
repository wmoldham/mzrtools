test_that("format is correct", {
  expect_named(mz_iso_resolve("C3", tracer = "C"),
               expected = c("elements", "isotopes", "polarity", "iso_list", "tracer", "targets", "resolved"))
})
