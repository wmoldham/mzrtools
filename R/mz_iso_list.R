#' Generate a matrix of stable isotopes
#'
#' \code{mz_iso_list} generates a matrix enumerating all possible stable isotopes
#'     of a molecule.
#'
#' @inheritParams mz_atomize
#'
#' @return A matrix of isotope combinations
#'
#' @export
#'
#' @examples
#' mz_iso_list("C5H8O5")
#'
mz_iso_list <- function(molecule) {

  elements <- mz_atomize(molecule)

  isotopes <- elements[names(elements) %in% names(iso_info)]

  if (length(isotopes) == 0) return(NA)

  iso_names <- lapply(iso_info[names(isotopes)], "[[", "isotope")

  suppressMessages(
    matrices <-
      purrr::pmap(.l = list(lapply(iso_info[names(isotopes)], "[[", "number"),
                            isotopes,
                            iso_names),
                  .f = gtools::combinations,
                  repeats.allowed = TRUE) %>%
      purrr::map2(iso_names, summarize_combos)
  )


  if (length(matrices) == 1) {
    CHONS <- matrices[[1]]
  } else {
    CHONS <- matrices[[1]]
    for (i in 2:length(matrices)) {
      CHONS <- assemble_matrices(CHONS, matrices[[i]])
    }
  }

  colnames(CHONS) <- unlist(iso_names)

  CHONS

}

summarize_combos <- function(combo_matrix, isotope_names) {
  purrr::map_dfc(isotope_names, ~rowSums(combo_matrix == .x)) %>%
    stats::setNames(isotope_names) %>%
    as.matrix()
}

assemble_matrices <- function(matrix_a, matrix_b) {
  a2 <- rep(1, nrow(matrix_b)) %x% matrix_a
  b2 <- matrix_b %x% rep(1, nrow(matrix_a))
  cbind(a2, b2)
}
