test_that("tracer argument errors work", {
  expect_error(mz_iso_target("C3", tracer = c("X", "C")))
  expect_error(mz_iso_target("C3", tracer = c("N")))
})

test_that("format is correct", {
  expect_named(mz_iso_target("C3", tracer = "C"),
               expected = c("elements", "isotopes", "polarity", "iso_list", "tracer", "targets"))
})


