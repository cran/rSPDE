## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

set.seed(1)

## ----inla_link, include = FALSE-----------------------------------------------
inla_link <- function() {
  sprintf("[%s](%s)", "`R-INLA`", "https://www.r-inla.org")
}

