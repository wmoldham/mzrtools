#' Annotate isotope table
#'
#' \code{mz_iso_annotate} generates a tibble enumerating all possible isotopes of
#'     a molecule and annotates this with the mass or m/\emph{z} of each isotope
#'     based on the polarity (pol) setting and the nominal mass shift of the isotope.
#'
#' @inheritParams mz_atomize
#' @inheritParams mz_calculate
#'
#' @return A list containing a tibble of annotated isotope information, the
#'     parsed molecular formula, and the isotopes present in the molecule.
#'
#' @export
#'
#' @examples
#' mz_iso_annotate("C5H8O5")
#' mz_iso_annotate("C5H8O5", pol = "negative")
#'
mz_iso_annotate <- function(mol, pol = "negative") {

  if (!(pol %in% c("neutral", "positive", "negative"))) {
    warning("Unrecognized polarity, neutral masses returned.")
  }

  elements <- mz_atomize(mol)
  isotopes <- elements[names(elements) %in% names(iso_info)]

  # get unlabeled mass
  no_label <- elements[names(elements) %nin% names(iso_info)]
  if (length(no_label) == 0) {
    unlabeled_mass <- 0
  } else {
    no_label <- stringr::str_c(names(no_label), no_label, collapse = "")
    unlabeled_mass <- mz_calculate(no_label, pol = pol)
  }

  # get isotope information

  # extractor function for iso_info
  annot <- function(sublist, iso_list) {
    out <- unlist(lapply(iso_info[names(isotopes)], "[[", sublist))
    names(out) <- unlist(lapply(iso_info[names(isotopes)], "[[", "isotope"))
    iso_list %*% out
  }

  iso_list <- mz_iso_list(mol)
  iso_mass <- c(annot("mass", iso_list) + unlabeled_mass)
  iso_shift <- c(annot("shift", iso_list))

  iso_list <-
    tibble::as_tibble(iso_list) %>%
    dplyr::bind_cols(mass = iso_mass, shift = iso_shift) %>%
    dplyr::arrange(mass)

  list(iso_list = iso_list, elements = elements, isotopes = isotopes)

}
