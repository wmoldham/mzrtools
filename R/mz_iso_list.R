#' Generate a matrix of stable isotopes
#'
#' \code{mz_iso_list} generates a matrix enumerating all possible isotopes
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
mz_iso_list <- function(mol) {

  elements <- mz_atomize(mol)

  isotopes <- elements[names(elements) %in% names(iso_info)]

  if (length(isotopes) == 0) return(NA)

  iso_names <- lapply(iso_info[names(isotopes)], '[[', "isotope")

  matrices <-
    purrr::pmap(.l = list(isotopes,
                          lapply(iso_info[names(isotopes)], '[[', "number")),
                .f = gtools::combinations,
                repeats.allowed = TRUE) %>%
    purrr::map2(., iso_names, magrittr::set_colnames)


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

assemble_matrices <- function(matrix_a, matrix_b) {
  a2 <- rep(1, nrow(matrix_b)) %x% matrix_a
  b2 <- matrix_b %x% rep(1, nrow(matrix_a))
  cbind(a2, b2)
}
