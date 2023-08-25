library(jacobi)
library(RcppColors)


# modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i*z) / (1i - z)
}

bkgcol <- "#15191e"

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
x <- seq(-1, 1, len = 2048)
y <- seq(-1, 1, len = 2048)
Z <- outer(x, y, f)
K <- kleinj(Z) / 1728
G <- K / (1 - K - K*K)
image <- colorMap4(G)
#
#svg("x.svg", width = 15, height = 15)
opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)


dev.off()
pngf <- sprintf("zzpic%03d.png", i)
rsvg::rsvg_png(
  "x.svg", pngf, width = 512, height = 512
)
