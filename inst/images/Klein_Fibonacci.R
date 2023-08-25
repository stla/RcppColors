# Klein-Fibonacci map ####
library(jacobi)
library(RcppColors)

# the modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i*z) / (1i - z)
}
# background color
bkgcol <- "#ffffff"
# make the color mapping
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  w <- PhiInv(z)
  ifelse(
    Mod(z) > 0.96, 
    NA_complex_,
    ifelse(
      y < 0, 
      -1/w, w
    )
  )
}
x <- seq(-1, 1, length.out = 2048)
y <- seq(-1, 1, length.out = 2048)
Z <- outer(x, y, f)
K <- kleinj(Z) / 1728
G <- K / (1 - K - K*K)
image <- colorMap4(G, bkgcolor = bkgcol)

#
svg("x.svg", width = 15, height = 15)
opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(
  c(-1, 1), c(-1, 1), type = "n", 
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(image, -1, -1, 1, 1)


# now we add the Dedekind tessellation (the white lines)
library(PlaneGeometry)
isInteger <- function(x) abs(x - floor(x)) < x * 1e-6
abline(h = 0, col = "white", lwd = 2)
N <- 150L
for(n in 1L:N) {
  if(isInteger(n/2) && ((n/2L) %% 2L == 1L)) {
    next
  }
  for(p in 1:n) {
    q <- sqrt(n*n - p*p + 4L)
    cases <- (isInteger(q) && isInteger(q/2) && (n %% 2L == 1L)) ||
      (isInteger(q) && isInteger(q/4) && (n %% 4L == 0L))
    if(cases) {
      circ <- Circle$new(center = c(q, p)/n, radius = 2/n)
      draw(circ, border = "white", lwd = 2)
      circ <- Circle$new(center = c(-q, p)/n, radius = 2/n)
      draw(circ, border = "white", lwd = 2)
      circ <- Circle$new(center = c(q, -p)/n, radius = 2/n)
      draw(circ, border = "white", lwd = 2)
      circ <- Circle$new(center = c(-q, -p)/n, radius = 2/n)
      draw(circ, border = "white", lwd = 2)
    }
  }
}
par(opar)



dev.off()
pngf <- "KleinFibonacciDedekind.png"
rsvg::rsvg_png(
  "x.svg", pngf, width = 512, height = 512
)
