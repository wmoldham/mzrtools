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
    dplyr::group_by(.data$shift) %>%
    dplyr::mutate(group = assign_groups(.data$mass, ...))
    # dplyr::mutate(group = assign_groups(.data$mass)) %>%
    dplyr::group_by(.data$shift, .data$group) %>%
    dplyr::summarise(mass = mean(.data$mass)) %>%
    dplyr::ungroup()

  list(targets = targets, annot = l)

}


# based on Orbitrap detector
delta_mz <-
  function(mz,
           nominal_mz = 200,
           nominal_resolution = 70000) {
    1.67 * mz * sqrt(mz / nominal_mz) / nominal_resolution
  }


assign_groups <- function(masses, ...) {
  if (length(masses) == 1) {
    "A"
  } else if (max(masses) - min(masses) < delta_mz(mean(masses), ...)) {
    "A"
  } else {
    groups <-
      dist(masses) %>%
      hclust() %>%
      cutree(h = delta_mz(mean(masses, na.rm = TRUE), ...))
    LETTERS[groups]
  }
}
