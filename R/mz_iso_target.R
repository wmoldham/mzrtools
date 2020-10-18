#' Identify LC-MS target isotopes
#'
#' \code{mz_iso_target} takes a molecular formula and tracer element(s) to
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
#'
#' @export
#'
#' @examples
#' mz_iso_target("C3H4O3", tracer = c("C", "H"))
#'
mz_iso_target <- function(molecule, tracer = "C", polarity = "negative", ...) {

  if (any(tracer %nin% names(iso_info))) {
    errs <- setdiff(tracer, names(iso_info))
    stop(stringr::str_c("Unsupported tracer(s) (",
                        stringr::str_c(errs, collapse = ", "),
                        ") submitted"))
  }

  l <- mz_iso_annotate(molecule = molecule, polarity = polarity)

  if (any(tracer %nin% names(l$elements))) {
    errs <- setdiff(tracer, names(l$elements))
    stop(stringr::str_c("Tracer(s) (",
                        stringr::str_c(errs, collapse = ", "),
                        ") are not elements in the molecule"))
  }

  iso_list <- l$iso_list

  labels <- unlist(lapply(iso_info[tracer], "[[", "label"))

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

  l[["tracer"]] <- tracer
  l[["targets"]] <- targets

  l

}
