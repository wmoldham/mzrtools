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

  iso_list <- l$iso_list

  labels <- unlist(lapply(iso_info[tracer], "[[", "label"))

  ## TODO: Account for increased mass shift with oxygen (and other) labels ##

  shifts <- matrix(nrow = length(labels), ncol = 1, dimnames = list(labels, NA))
  for (nm in names(labels)) {
    idx <- which(iso_info[[nm]]$isotope == labels[[nm]])
    shifts[labels[[nm]], ] <- iso_info[[nm]]$shift[[idx]]
  }

  iso_shifts <- as.matrix(iso_list[, labels, drop = FALSE]) %*% shifts

  target_isos <- which(iso_shifts == iso_list[, "shift"])

  # format output and combine unresolved masses
  targets <-
    dplyr::slice(iso_list, target_isos) %>%
    dplyr::arrange(.data$mass) %>%
    dplyr::select(.data$mass, .data$shift, dplyr::everything())

  list(elements = l$elements,
       isotopes = l$isotopes,
       iso_list = iso_list,
       targets = targets)

}
