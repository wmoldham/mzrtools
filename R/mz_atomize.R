#' Parse a molecular formula into element and count
#'
#' \code{mz_atomize} separates a molecular formula into element and count,
#'     returning a named vector of element counts.
#'
#' @param molecule A string containing a molecular formula (\emph{e.g.}, "C2H7NO3S").
#'     Structural formulas containing parentheses are not acceptable. Charges
#'     may be included, but the charge count should follow the sign (\emph{e.g.},
#'     "C10H16N5O13P3-3").
#'
#' @return A vector of counts named by element.
#'
#' @export
#'
#' @examples
#' mz_atomize("C5H8O5-")
#' mz_atomize("C2H7NO3S")
#' mz_atomize("C10H16N5O13P3-3")
#'
mz_atomize <- function(molecule) {

  # check arguments
  correct_format <- "^(([A-Z]{1}[a-z]?)[0-9]*)((([A-Z]{1}[a-z]?)|(\\+|\\-))[0-9]*)*$"
  if (!stringr::str_detect(molecule, pattern = correct_format)) {
    stop("Incorrect format, provide molecular formula (e.g., \"C5H6O5\")")
  }

  atom_count <- "(([A-Z]{1}[a-z]?)|(\\+|\\-))[0-9]*"
  atoms <- stringr::str_extract_all(molecule, pattern = atom_count, simplify = TRUE)
  elements <- stringr::str_extract(atoms, "\\D*")

  # verify elements
  if(sum(elements %nin% names(atomic_mass)) > 0) {
    errs <- elements[elements %nin% names(atomic_mass)]
    stop(stringr::str_c("Unknown element(s) (",
                        stringr::str_c(errs, collapse = ", "),
                        ") used in supplied formula"))
  }

  # check for duplicates
  if (sum(duplicated(elements)) > 0) {
    errs <- unique(elements[duplicated(elements)])
    stop(stringr::str_c("Duplicated element(s) (",
                        stringr::str_c(z, collapse = ", "),
                        ") used in supplied formula"))
  }

  counts <-
    stringr::str_extract(atoms, "\\d+") %>%
    replace(is.na(.), "1") %>%
    as.integer()

  names(counts) <- elements
  counts

}
