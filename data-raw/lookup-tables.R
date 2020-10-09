### LOOKUP TABLES ###

# libraries ---------------------------------------------------------------

library(dplyr)
library(stringr)
library(tidyr)

conflicted::conflict_prefer("filter", "dplyr")


# generate atomic mass lookup table ---------------------------------------

raw_atomic_mass_table <-
  readr::read_delim(
    "data-raw/element-list.txt",
    delim = "=",
    col_names = FALSE,
    trim_ws = TRUE,
    skip = 7
  )

atomic_mass_table <-
  raw_atomic_mass_table %>%
  mutate(index = rep(seq(1, nrow(raw_atomic_mass_table) / 7), each = 7)) %>%
  pivot_wider(names_from = X1, values_from = X2) %>%
  rename(
    "element_number" = "Atomic Number",
    "symbol" = "Atomic Symbol",
    "enrichment" = "Isotopic Composition",
    "mass_number" = "Mass Number",
    "mass" = "Relative Atomic Mass",
    "rel_mass" = "Standard Atomic Weight"
  ) %>%
  select(-Notes) %>%
  separate(rel_mass, c("low_rel_mass", "high_rel_mass"), ",") %>%
  mutate(across(c(enrichment, mass, low_rel_mass, high_rel_mass), ~str_extract(., "\\d+\\.\\d+"))) %>%
  mutate(across(c(-symbol), as.numeric))

elements <-
  group_by(atomic_mass_table, element_number) %>%
  mutate(mass_delta = abs(mass - low_rel_mass)) %>%
  filter(mass_delta == min(mass_delta)) %>%
  ungroup() %>%
  select(symbol, mass)

atomic_mass <-
  elements$mass %>%
  c(., 5.48579909070e-4, 1.007276466879, 0, 0) %>%
  rlang::set_names(c(elements$symbol, "electron", "proton", "+", "-"))


# define isotope properties -----------------------------------------------

iso_info <-
  list(
    C = list(number = 2,
             isotope = c("C12", "C13"),
             mass = c(12, 13.003355),
             abundance = c(0.9893, 0.0107),
             shift = c(0, 1)),
    H = list(number = 2,
             isotope = c("H1", "H2"),
             masse = c(1.007825, 2.014102),
             abundance = c(0.999885, 0.00115),
             shift = c(0, 1)),
    O = list(number = 3,
             isotope = c("O16", "017", "O18"),
             mass = c(15.994915, 16.999132, 17.999160),
             abundance = c(0.99757, 0.00038, 0.00205),
             shift = c(0, 1, 2)),
    N = list(number = 2,
             isotope = c("N14", "N15"),
             mass = c(14.003074, 15.000109),
             abundance = c(0.99632, 0.00368),
             shift = c(0, 1))
  )

# isotope_properties <-
#   tibble::tribble(
#     ~element, ~mass, ~abundance, ~M,
#     "C12", 12,        98.93,   0,
#     "C13", 13.003355,  1.07,   1,
#     "H1",   1.007825, 99.9885, 0,
#     "H2",   2.014102,  0.0115, 1,
#     "O16", 15.994915, 99.757,  0,
#     "O17", 16.999132,  0.038,  1,
#     "O18", 17.999160,  0.205,  2,
#     "N14", 14.003074, 99.632,  0,
#     "N15", 15.000109,  0.368,  1,
#     "S32", 31.972071, 94.93,   0,
#     "S33", 32.971458,  0.76,   1,
#     "S34", 33.967867,  4.29,   2,
#     "S36", 35.967081,  0.02,   4
#   ) %>%
#   mutate(abundance = abundance/100)


# create tibble of relevant isotopes --------------------------------------

# isotope_table <-
#   tibble::tribble(
#     ~element, ~isotopes, ~names, ~heavy,
#     "C", 2, c("C12", "C13"), "C13",
#     "H", 2, c("H1", "H2"), "H2",
#     "O", 3, c("O16", "O17", "O18"), c("O17", "O18"),
#     "N", 2, c("N14", "N15"), "N15",
#     "S", 4, c("S32", "S33", "S34", "S36"), c("S33", "S34", "S36")
#   )


# update R/sysdata.rda ----------------------------------------------------

usethis::use_data(
  atomic_mass,
  iso_info,
  # isotope_properties,
  # isotope_table,
  internal = TRUE,
  overwrite = TRUE
)
