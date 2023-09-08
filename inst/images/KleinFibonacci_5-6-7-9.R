library(RcppColors)
library(jacobi)

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
image <- colorMap5(G, bkgcolor = bkgcol)
# plot
opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(
  c(-1, 1), c(-1, 1), type = "n", xaxs = "i", yaxs = "i", 
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(image, -1, -1, 1, 1)
par(opar)

# 
image5 <- colorMap5(G, bkgcolor = bkgcol)
image6 <- colorMap6(G, bkgcolor = bkgcol)
image7 <- colorMap7(G, bkgcolor = bkgcol)
image9 <- colorMap9(G, bkgcolor = bkgcol)

svglite::svglite("x.svg", width = 10, height = 10)
layout(rbind(c(1, 2), c(3, 4)))
#
opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(
  c(-1, 1), c(-1, 1), type = "n", xaxs = "i", yaxs = "i", 
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(image5, -1, -1, 1, 1)
#
plot(
  c(-1, 1), c(-1, 1), type = "n", xaxs = "i", yaxs = "i", 
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(image6, -1, -1, 1, 1)
#
plot(
  c(-1, 1), c(-1, 1), type = "n", xaxs = "i", yaxs = "i", 
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(image7, -1, -1, 1, 1)
#
plot(
  c(-1, 1), c(-1, 1), type = "n", xaxs = "i", yaxs = "i", 
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(image9, -1, -1, 1, 1)
par(opar)
dev.off()

rsvg::rsvg_png("x.svg", "Klein_5-6-7-9.png", width = 512, height = 512)
