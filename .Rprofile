myinstall <- function() {
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "devtools::install(quick = TRUE, keep_source = TRUE)"
    )
  } else {
    devtools::install(quick = TRUE, keep_source = TRUE)
  }
}
mydocument <- function() {
  rstudioapi::restartSession(
   "roxygen2::roxygenise(load_code = roxygen2::load_installed)" 
  )
}
