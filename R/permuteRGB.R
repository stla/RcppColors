#' @title RGB permutation
#' @description Permutes the R-G-B components of a color.
#'
#' @param hexcolor vector or matrix or array of hexadecimal colors
#' @param permutation a character string with three letters \code{"r"}, 
#'   \code{"g"} and \code{"b"}
#'
#' @return The colors after applying the permutation.
#' @export
#'
#' @examples
#' library(RcppColors)
#' x <- y <- seq(-1.7, 1.7, length.out = 512L)
#' zarray <- outer(y, x, function(x, y) {
#'   z <- x + 1i*y
#'   (1 + 1i) * log(sin((z^3 - 1)))
#' })
#' # image
#' img1 <- colorMap1(zarray)
#' # r -> b, g -> r, b -> g
#' img2 <- permuteRGB(img1, "brg")
#' # plot
#' opar <- par(mar = c(0,0,0,0), mfrow = c(1, 2), bg = "#002240")
#' plot(
#'   c(0, 1), c(0, 1), type = "n", asp = 1,
#'   xlab = NA, ylab = NA, axes = FALSE
#' )
#' rasterImage(img1, 0, 0, 1, 1, interpolate = TRUE)
#' plot(
#'   c(0, 1), c(0, 1), type = "n", asp = 1,
#'   xlab = NA, ylab = NA, axes = FALSE
#' )
#' rasterImage(img2, 0, 0, 1, 1, interpolate = TRUE)
#' par(opar)
permuteRGB <- function(hexcolor, permutation = "gbr") {
  stopifnot(nchar(permutation) == 3L)
  perm <- substring(permutation, 1L:3L, 1L:3L)
  tmp <- paste0(sort(perm), collapse = "")
  if(tmp != "bgr") {
    stop("Invalid `permutation` argument.")
  }
  perm <- c("r" = 1L, "g" = 2L, "b" = 3L)[perm]
  perm[perm] <- 1L:3L
  dims     <- dim(hexcolor)
  hexcolor <- c(hexcolor)
  R <- substring(hexcolor, 2L, 3L)
  G <- substring(hexcolor, 4L, 5L)
  B <- substring(hexcolor, 6L, 7L)
  hexcolor <- cbind(R, G, B)[, perm]
  hexcolor <- paste0("#", hexcolor[, 1L], hexcolor[, 2L], hexcolor[, 3L])
  dim(hexcolor) <- dims
  hexcolor
}
