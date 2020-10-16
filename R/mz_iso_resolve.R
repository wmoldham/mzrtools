#' Filter unlabeled resolved isotopes
#'
#' \code{mz_iso_resolve} will return the molecular isotopes of interest
#'     (\emph(e.g.), M0 or label-containing) and those isotopes that cannot be
#'     resolved from these targets by the mass spectrometer.
#'
#' @inheritParams mz_iso_target
#'
#' @return A tibble containing the unit mass shift (`m`), group, mass or
#'     m/\emph{z}, and isotope composition.
#'
#' @export
#'
#' @examples
#' mz_iso_resolve("C5H8O5")
#' mz_iso_resolve("C5H8O5", tracer = c("C", "H"))
#'
mz_iso_resolve <- function(mol, tracer, pol = "negative", ...) {

  # l <- mz_iso_target(mol = mol, tracer = tracer, pol = pol)
  l <- mz_iso_target(mol = mol, tracer = tracer, pol = pol, ...)
  iso_list <- l$iso_list
  targets <- l$targets

  resolved <-
    targets %>%
    dplyr::select(.data$shift, .data$mass) %>%
    dplyr::left_join(iso_list, by = "shift", suffix = c(".tar", ".iso")) %>%
    dplyr::mutate(diff = abs(.data$mass.tar - .data$mass.iso),
                  delta = delta_mz(.data$mass.tar, ...)) %>%
    dplyr::filter(.data$diff <= .data$delta) %>%
    dplyr::select(-c(.data$mass.iso, .data$diff, .data$delta)) %>%
    dplyr::rename(mass = .data$mass.tar)

  l[["resolved"]] <- resolved

  l

}


# based on Orbitrap detector
delta_mz <-
  function(mz,
           nominal_mz = 200,
           nominal_resolution = 70000) {
    1.67 * mz * sqrt(mz / nominal_mz) / nominal_resolution
  }
