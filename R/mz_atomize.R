mz_atomize <- function(mol) {

  # check arguments
  correct_format <- "^(([A-Z]{1}[a-z]?)[0-9]*)((([A-Z]{1}[a-z]?)|(\\+|\\-))[0-9]*)*$"
  if (!stringr::str_detect(mol, pattern = correct_format)) {
    stop("Incorrect format, provide molecular formula (e.g., \"C5H6O5\").")
  }

  atom_count <- "(([A-Z]{1}[a-z]?)|(\\+|\\-))[0-9]*"
  parsed <-
    stringr::str_extract_all(mol,
                             pattern = atom_count,
                             simplify = TRUE) %>%
    stringr::str_split_fixed(pattern = "(?<=\\D)(?=\\d)", n = 2) %>%
    magrittr::set_colnames(c("element", "count")) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(count = replace(.data$count, .data$count == "", 1),
                  count = as.double(.data$count))

  if (sum(!(parsed$element %in% c(names(atomic_mass), "+", "-"))) > 0) {
    z <- parsed$element[!(parsed$element %in% c(names(atomic_mass)))]
    stop(stringr::str_c("Unknown element(s) (",
                        stringr::str_c(z, collapse = ", "),
                        ") used in supplied formula."))
  }

  if (sum(duplicated(parsed$element)) > 0) {
    z <- parsed$element[duplicated(parsed$element)]
    stop(stringr::str_c("Duplicated element(s) (",
                        stringr::str_c(z, collapse = ", "),
                        ") used in supplied formula."))
  }

  if (sum(parsed$count == 0) > 0) {
    z <- parsed$element[parsed$count == 0]
    stop(stringr::str_c("Element(s) (",
                        stringr::str_c(z, collapse = ", "),
                        ") provided count of 0."))
  }

  parsed

}
