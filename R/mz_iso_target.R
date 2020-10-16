#' Identify LC-MS target isotopes
#'
#' \code{mz_iso_target} takes a molecular formula and tracer elements to
#'     determine which molecular isotopes should be measured by LC-MS.
#'
#' @inheritParams mz_iso_annotate
#' @param tracer A vector of elements labeled by the tracer.
#' @param ... Passes the arguments nominal_resolution (default 70,000) and
#'     nominal_mz (default 200) used to determine the delta m/\emph{z} below
#'     which peaks are considered unresolved.
#'
#' @return A list containing a tibble of annotated isotope targets and the output
#'     of `mz_iso_annotate`.
#' @export
#'
#' @examples
#' mz_iso_target("C3H4O3")
#'
mz_iso_target <- function(mol, tracer, pol = "negative", ...) {

  l <- mz_iso_annotate(mol, pol)

  iso <- l$iso_list

  labels <- unlist(lapply(iso_info[tracer], "[[", "label"))

  target_isos <- which(rowSums(iso[, labels, drop = FALSE]) == iso[, "shift"])

  # format output and combine unresolved masses
  targets <-
    dplyr::slice(iso, target_isos) %>%
    dplyr::arrange(.data$mass)

  list(targets = targets, annot = l)

}
