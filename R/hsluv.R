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
#' @param alpha opacity, a number between \code{0} and \code{1}, 
#'   or \code{NULL} 
#'
#' @return The \code{hsluv} function returns a hexadecimal string representing 
#'   a color, and the \code{hsluv2rgb} returns the RGB coordinates of this 
#'   color, a named vector of three integers between \code{0} and \code{255}.
#' @export
#'
#' @examples
#' saturation <- 100
#' f <- Vectorize(
#'   function(x, y){
#'     z <- complex(real = x, imaginary = y)
#'     modulus <- Mod(z)
#'     if(modulus > 1){
#'       return("#ffffff")
#'     }
#'     radians <- Arg(z)
#'     if(radians < 0){
#'       radians <- radians + 2*pi
#'     }
#'     degrees <- 360 * radians / 2 / pi
#'     hsluv(h = degrees, s = saturation, l = 100*modulus)
#'   }
#' )
#' 
#' x <- y <- seq(-1, 1, length.out = 200L)
#' image <- outer(x, y, f)
#' 
#' opar <- par(mar = c(0, 0, 0, 0))
#' plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
#' rasterImage(image, -1, -1, 1, 1)
#' par(opar)
hsluv <- function(h = 360, s = 100, l = 100, alpha = NULL){
  if(is.null(alpha)){
    hsluv_cpp(h, s, l)
  }else{
    hsluv_alpha_cpp(h, s, l, alpha)
  }
}

#' @rdname hsluv
#' @export
hsluv2rgb <- function(h = 360, s = 100, l = 100){
  `names<-`(hsluv2rgb_cpp(h, s, l), c("r", "g", "b"))
}

#' @title HSL color specification
#' @description Converts a color given in HSL coordinates to a hexadecimal 
#'   string  
#'
#' @param h the hue, a number between \code{0} and \code{360}
#' @param s the saturation, a number between \code{0} and \code{100}
#' @param l the lightness, a number between \code{0} and \code{100}
#' @param alpha opacity, a number between \code{0} and \code{1}, 
#'   or \code{NULL} 
#'
#' @return The \code{hsl} function returns a hexadecimal string representing 
#'   the corresponding color.
#' @export
#'
#' @examples
#' saturation <- 100
#' f <- Vectorize(
#'   function(x, y){
#'     z <- complex(real = x, imaginary = y)
#'     modulus <- Mod(z)
#'     if(modulus > 1){
#'       return("#ffffff")
#'     }
#'     radians <- Arg(z)
#'     if(radians < 0){
#'       radians <- radians + 2*pi
#'     }
#'     degrees <- 360 * radians / 2 / pi
#'     hsl(h = degrees, s = saturation, l = 100*modulus)
#'   }
#' )
#' 
#' x <- y <- seq(-1, 1, length.out = 200L)
#' image <- outer(x, y, f)
#' 
#' opar <- par(mar = c(0, 0, 0, 0))
#' plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
#' rasterImage(image, -1, -1, 1, 1)
#' par(opar)
hsl <- function(h = 360, s = 100, l = 100, alpha = NULL){
  if(is.null(alpha)){
    hsl_cpp(h, s, l)
  }else{
    hsl_alpha_cpp(h, s, l, alpha)
  }
}
