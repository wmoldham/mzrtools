mz_iso_probs <- function(mol,
                         tracer,
                         purity = list(C = 1, N = 1, O = 1, S = 1, Se = 1)) {

  elements <- mz_atomize(mol)

  isotopes <- elements[names(elements) %in% names(iso_info)]

  mat <- mz_iso_list(mol)

  molecule_m <- vector(mode = "list", length = length(isotopes))
  names(molecule_m) <- names(isotopes)
  for (nm in names(isotopes)) {

    # get tracer info
    if (nm %in% tracer) {
      label_no <- elements[[nm]]
    } else {
      label_no <- 0
    }
    label_purity <- purity[[nm]]

    # get isotope abundances
    probs <- iso_info[[nm]]$abundance

    # extract columns from matrix
    f1 <- stringr::str_detect(colnames(mat), stringr::str_c(nm, "\\d"))


    # calculate probabilities for various label amounts
    element_m <- vector(mode = "list", length = label_no + 1)
    names(element_m) <- stringr::str_c(nm, 0:label_no)

    for (t in 0:label_no) {

      m <- mat[, f1]

      total <- rowSums(m) - t
      numerator <- factorial(total)

      f2 <- stringr::str_which(colnames(m), iso_info[[nm]]$label)
      m[, f2] <- m[, f2] - t
      denominator <- suppressWarnings(apply(factorial(m), 1, prod))
      ratio <- numerator/denominator
      abundance <- matrix(rep(probs, length(total)), nrow = length(total), byrow = TRUE)
      p <- ratio * apply(abundance ^ m, 1, prod)

      # change probabilities of impossible combinations to 0
      impossible <- which(is.nan(ratio))
      if (length(impossible) != 0) {
        p[impossible] <- 0
      }

      # tracer purity
      label_prob <- label_purity ^ t

      element_m[[t + 1]] <- p * label_prob

    }

    molecule_m[[nm]] <- data.frame(element_m)

  }

  combos <- expand.grid(lapply(molecule_m, names), stringsAsFactors = FALSE)
  choices <- do.call(cbind, molecule_m)
  names(choices) <- stringr::str_extract(names(choices), "[A-Z][a-z]?\\d+$")

  if (ncol(combos == 1)) {
    out <- choices
  } else {
    out <- apply(combos, 1, function(x) {
      f <- which(colnames(choices) %in% x)
      apply(choices[, f], 1, prod)
    })
    colnames(out) <- apply(combos, 1, paste, collapse = ".")
  }
  solve(out)
}
