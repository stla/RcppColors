#' @useDynLib RcppColors, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @noRd
NULL


#' @title HSLuv color specification
#' @description Converts a color given in HSLuv coordinates to a hexadecimal 
#'   string or a RGB color specification  
#'
#' @param h the hue, a number between \code{0} and \code{360}
#' @param s the saturation, a number between \code{0} and \code{100}
#' @param l the lightness, a number between \code{0} and \code{100}
#'
#' @return The \code{hsluv} function returns a hexadecimal string representing 
#'   a color, and the \code{hsluv2rgb} returns the RGB coordinates of this 
#'   color, a named vector of three integers between \code{0} and \code{255}.
#' @export
#'
#' @examples
#' xx
hsluv <- function(h = 360, s = 100, l = 100){
  hsluv_cpp(h, s, l)
}

#' @rdname hsluv
#' @export
hsluv2rgb <- function(h = 360, s = 100, l = 100){
  `names<-`(hsluv2rgb_cpp(h, s, l), c("r", "g", "b"))
}