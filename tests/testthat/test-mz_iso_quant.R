test_case <- mz_iso_quant("C3H3O3", tracer = c("C", "H"), purity = list(C = 0.99, H = 0.99))

test_that("format is correct", {
  expect_named(test_case,
               expected = c("elements", "isotopes", "polarity", "iso_list", "tracer", "targets", "resolved", "prob_matrix"))
})

test_that("math is correct", {
  expect_equal(sum(colSums(test_case[["prob_matrix"]])), ncol(test_case[["prob_matrix"]]))
})
