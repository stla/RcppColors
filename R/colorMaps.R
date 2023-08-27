isComplex <- function(x) {
  is.complex(x) || is.numeric(x)
}

isString <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

isBooleanTriplet <- function(x) {
  is.logical(x) && length(x) == 3L && !anyNA(x)
}

isNumber <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x)
}


#' @title Color mappings functions
#' @description Functions mapping each complex number to a color.
#'
#' @param Z complex number, vector or matrix
#' @param bkgcolor background color; it is applied for the \code{NA} values 
#'   of \code{Z}  
#' @param nancolor color for infinite and \code{NaN} values 
#' @param reverse logical vector of length three; for each component of the
#'   color space (R, G, B or H, S, L), whether to reverse it (e.g. 
#'   \code{R -> 255-R})
#' @param s saturation, a number between 0 and 100
#' @param n number of rays drawn in a cycle; it should be a positive integer 
#'   but any non-zero numeric value is accepted
#' @param nthreads number of threads used for parallel computation
#'
#' @return A string or a character vector or a character matrix, 
#'   having the same size as \code{Z}. Each entry is a color given 
#'   by a hexadecimal string.
#' 
#' \if{html}{
#'   \figure{ModularForm.png}{options: style="max-width: 60\%; text-align: center; display: block;"}
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
  reverse = c(FALSE, FALSE, FALSE), 
  nthreads = 1L
){
  stopifnot(isComplex(Z))
  stopifnot(isString(bkgcolor))
  stopifnot(isString(nancolor))
  stopifnot(isBooleanTriplet(reverse))
  nthreads <- as.integer(nthreads)
  stopifnot(nthreads >= 1L)
  ismatrix <- is.matrix(Z)
  storage.mode(Z) <- "complex"
  if(!ismatrix){
    Z <- cbind(Z)
  }
  P <- ColorMap1(
    Z, bkgcolor, nancolor, reverse[1], reverse[2], reverse[3], nthreads
  )
  if(!ismatrix){
    P <- c(P)
  }
  P
}

#' @rdname colorMaps
#' @export
colorMap2 <- function(
    Z, bkgcolor = "#15191e", nancolor = "#000000", 
    reverse = c(FALSE, FALSE, FALSE),
    nthreads = 1L
){
  stopifnot(isComplex(Z))
  stopifnot(isString(bkgcolor))
  stopifnot(isString(nancolor))
  stopifnot(isBooleanTriplet(reverse))
  nthreads <- as.integer(nthreads)
  stopifnot(nthreads >= 1L)
  storage.mode(Z) <- "complex"
  ismatrix <- is.matrix(Z)
  if(!ismatrix){
    Z <- cbind(Z)
  }
  P <- ColorMap2(
    Z, bkgcolor, nancolor, reverse[1], reverse[2], reverse[3], nthreads
  )
  if(!ismatrix){
    P <- c(P)
  }
  P
}

#' @rdname colorMaps
#' @export
colorMap3 <- function(
    Z, bkgcolor = "#15191e", nancolor = "#000000", 
    s = 80, n = 5, 
    nthreads = 1L
){
  stopifnot(isComplex(Z))
  stopifnot(isString(bkgcolor))
  stopifnot(isString(nancolor))
  stopifnot(isNumber(s))
  s <- as.double(s)
  stopifnot(s >= 0, s <= 100)
  stopifnot(isNumber(n))
  n <- as.double(n)
  stopifnot(n != 0)
  nthreads <- as.integer(nthreads)
  stopifnot(nthreads >= 1L)
  storage.mode(Z) <- "complex"
  ismatrix <- is.matrix(Z)
  if(!ismatrix){
    Z <- cbind(Z)
  }
  P <- ColorMap3(
    Z, bkgcolor, nancolor, s, n, nthreads
  )
  if(!ismatrix){
    P <- c(P)
  }
  P
}

#' @rdname colorMaps
#' @export
colorMap4 <- function(
    Z, bkgcolor = "#15191e", nancolor = "#000000", 
    reverse = c(FALSE, FALSE, FALSE),
    nthreads = 1L
){
  stopifnot(isComplex(Z))
  stopifnot(isString(bkgcolor))
  stopifnot(isString(nancolor))
  stopifnot(isBooleanTriplet(reverse))
  nthreads <- as.integer(nthreads)
  stopifnot(nthreads >= 1L)
  storage.mode(Z) <- "complex"
  ismatrix <- is.matrix(Z)
  if(!ismatrix){
    Z <- cbind(Z)
  }
  P <- ColorMap4(
    Z, bkgcolor, nancolor, reverse[1], reverse[2], reverse[3], nthreads
  )
  if(!ismatrix){
    P <- c(P)
  }
  P
}

#' @rdname colorMaps
#' @export
colorMap5 <- function(
    Z, bkgcolor = "#15191e", nancolor = "#000000", 
    reverse = c(FALSE, FALSE, FALSE),
    nthreads = 1L
){
  stopifnot(isComplex(Z))
  stopifnot(isString(bkgcolor))
  stopifnot(isString(nancolor))
  stopifnot(isBooleanTriplet(reverse))
  nthreads <- as.integer(nthreads)
  stopifnot(nthreads >= 1L)
  storage.mode(Z) <- "complex"
  ismatrix <- is.matrix(Z)
  if(!ismatrix){
    Z <- cbind(Z)
  }
  P <- ColorMap5(
    Z, bkgcolor, nancolor, reverse[1], reverse[2], reverse[3], nthreads
  )
  if(!ismatrix){
    P <- c(P)
  }
  P
}

