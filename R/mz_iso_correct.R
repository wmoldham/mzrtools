#' Correct intensity vector for natural isotope abundance
#'
#' Uses the `pracma::lsqlincon` function to correct the measured mass isotope
#' distribution for natural abundance with a lower bound of 0 for all corrected
#' values.
#'
#' @param prob_matrix Correction matrix output by `mz_iso_quant`
#' @param v Mass distribution vector of isotopologue intensities
#'
#' @return Non-negative mass distribution vector corrected for natural isotope abundance
#'
#' @export
#'
#' @examples
#' intensity_vector <- c(5710407, 847435, 66688649, 3426642, 8222080, 99981, 0)
#' v <- intensity_vector/sum(intensity_vector)
#' prob_matrix <- mz_iso_quant("C6H14O12P2")$prob_matrix
#' mz_iso_correct(prob_matrix, v)
#'
mz_iso_correct <- function(prob_matrix, v) {

  if (isFALSE(all.equal(sum(v), 1))) {
    stop("Sum of intensity vector is not 1")
  }

  if (length(v) != ncol(prob_matrix)) {
    stop("Incompatible dimensions")
  }

  pracma::lsqlincon(C = prob_matrix, d = v, lb = 0, ub = 1)

}
