#' Basilisk environment for CNAlign
#' @import basilisk
cnalign_env <- basilisk::BasiliskEnvironment(
  envname = "cnalign_env",
  pkgname = "CNAlign",
  pip = "true",
  packages = c(
    "gurobipy==12.0.2",
    "numpy==2.1.2"
  )
)
