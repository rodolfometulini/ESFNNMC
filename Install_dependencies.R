required_packages <- c(
  "dplyr",
  "tidyr",
  "tibble",
  "sp",
  "sf",
  "mapview",
  "ggplot2",
  "spdep"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
} else {
  message("All required packages are already installed.")
}