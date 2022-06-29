isComplex <- function(x){
  is.complex(x) || is.numeric(x)
}

isString <- function(x){
  is.character(x) && length(x) == 1L && !is.na(x)
}

isBooleanTriplet <- function(x){
  is.logical(x) && length(x) == 3L && !anyNA(x)
}

#' @title Color mappings functions
#' @description Functions mapping each complex number to a color.
#'
#' @param Z complex number, vector or matrix
#' @param bkgcolor background color; it is applied for the \code{NA} values 
#'   of \code{Z}  
#' @param nancolor color for infinite and \code{NaN} values 
#' @param reverse logical vector of length three; for each color component 
#'   (e.g. R, G, B), whether to reverse it (e.g. \code{R -> 255-R})
#'
#' @return A string or a character vector or a character matrix, 
#'   having the same size as \code{Z}. Each entry is a color given 
#'   by a hexadecimal string.
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{ModularForm.png}{options: style="max-width:60\%;"}\out{</div>}
#' }
#' \if{latex}{
#'   \out{\begin{center}}\figure{ModularForm.png}\out{\end{center}}
#' }
#' 
#' @export
#' 
#' @rdname colorMaps
#'
#' @examples
#' library(RcppColors)
#' 
#' iota <- function(z){
#'   (z + 1i) / (1i*z + 1)
#' }
#' f <- function(z){
#'   q <- exp(2i * pi * z)
#'   r <- q - 4*q^2 + 2*q^3 + 8*q^4 - 5*q^5 - 8*q^6 + 6*q^7 - 23*q^9
#'   r / Mod(r)
#' }
#' g <- function(z){
#'   ifelse(
#'     Mod(z) >= 1, 
#'     NA_complex_,
#'     f(iota(Conj(z)))
#'   )
#' }
#' 
#' x <- y <- seq(-1, 1, len = 1500)
#' W <- outer(y, x, function(x, y) complex(real = x, imaginary = y))
#' Z <- g(W)
#' image <- colorMap1(Z)
#' 
#' opar <- par(mar = c(0,0,0,0), bg = "#15191E")
#' plot(
#'   c(-100, 100), c(-100, 100), type = "n", 
#'   xlab = "", ylab = "", axes = FALSE, asp = 1
#' )
#' rasterImage(image, -100, -100, 100, 100)
#' par(opar)
colorMap1 <- function(
  Z, bkgcolor = "#15191e", nancolor = "#000000", 
  reverse = c(FALSE, FALSE, FALSE)
){
  stopifnot(isComplex(Z))
  stopifnot(isString(bkgcolor))
  stopifnot(isString(nancolor))
  stopifnot(isBooleanTriplet(reverse))
  ismatrix <- is.matrix(Z)
  storage.mode(Z) <- "complex"
  if(!ismatrix){
    Z <- cbind(Z)
  }
  P <- ColorMap1(Z, bkgcolor, nancolor, reverse[1], reverse[2], reverse[3])
  if(!ismatrix){
    P <- c(P)
  }
  P
}

#' @rdname colorMaps
#' @export
colorMap2 <- function(
    Z, bkgcolor = "#15191e", nancolor = "#000000", 
    reverse = c(FALSE, FALSE, FALSE)
){
  stopifnot(isComplex(Z))
  stopifnot(isString(bkgcolor))
  stopifnot(isString(nancolor))
  stopifnot(isBooleanTriplet(reverse))
  storage.mode(Z) <- "complex"
  ismatrix <- is.matrix(Z)
  if(!ismatrix){
    Z <- cbind(Z)
  }
  P <- ColorMap2(Z, bkgcolor, nancolor, reverse[1], reverse[2], reverse[3])
  if(!ismatrix){
    P <- c(P)
  }
  P
}