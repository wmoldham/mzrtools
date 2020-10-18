test_that("format is correct", {
  expect_named(mz_iso_quant("C3", tracer = "C"),
               expected = c("elements", "isotopes", "polarity", "iso_list", "tracer", "targets", "resolved", "prob_matrix"))
})
