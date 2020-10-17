#' Calculate the exact mass or m/\emph{z} of a molecule
#'
#' \code{mz_calculate} returns the exact mass of a molecule based on the
#'     molecular formula or it will return the expected m/\emph{z} based on the
#'     gain or loss of a proton if provided an argument for scan polarity.
#'     Other adducts should be specifically annotated in the molecular formula.
#'
#' @inheritParams mz_atomize
#' @param polarity Type of m/\emph{z} to return. Accepts one of "neutral", which
#'     returns the monoisotopic mass; "positive", which returns the \[M+H\]+ mass;
#'     or "negative", which returns the \[M-H\]- mass.
#'
#' @return A numeric vector length 1.
#'
#' @export
#'
#' @examples
#' mz_calculate("C5H8O5-")
#' mz_calculate("C2H7NO3S")
#' mz_calculate("C10H16N5O13P3-3")
#'
mz_calculate <- function(molecule, polarity = "neutral") {

  if (polarity %nin% c("neutral", "positive", "negative")) {
    warning("Unrecognized polarity, neutral mass returned.")
  }

  elements <- unlist(mz_atomize(molecule))
  masses <- atomic_mass[names(elements)]
  total <- sum(masses * elements)

  switch(polarity,
         "neutral" = total,
         "positive" = total + atomic_mass[["proton"]],
         "negative" = total - atomic_mass[["proton"]],
         total)

}
