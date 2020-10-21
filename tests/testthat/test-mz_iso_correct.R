prob_matrix <- mz_iso_quant("C2")$prob_matrix

test_that("errors work", {
  expect_error(mz_iso_correct(prob_matrix, c(10, 1)))
  expect_error(mz_iso_correct(prob_matrix, c(1:5)))
})
