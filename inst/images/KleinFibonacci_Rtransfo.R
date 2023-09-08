library(RcppColors)
library(jacobi)

# a modular transformation
R <- function(z, t) {
  a <- pi*t/3
  ((sqrt(3)*cos(a) - sin(a)) * z - 2*sin(a))/
    (2*sin(a) * z + sqrt(3)*cos(a) + sin(a))
}

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

t_ <- seq(0, 3, length.out= 91L)[-1L]
for(i in 1:90) {
  KRt <- R(K, t_[i])
  G <- KRt / (1 - KRt - KRt*KRt)
  image <- colorMap5(G, bkgcolor = bkgcol)
  svglite::svglite("x.svg", width = 10, height = 10)
  opar <- par(mar = c(0,0,0,0), bg = bkgcol)
  plot(
    c(-1, 1), c(-1, 1), type = "n", xaxs = "i", yaxs = "i", 
    xlab = NA, ylab = NA, axes = FALSE, asp = 1
  )
  rasterImage(image, -1, -1, 1, 1)
  par(opar)
  dev.off()
  rsvg::rsvg_png(
    "x.svg", sprintf("zzpic%03d.png", i), width = 512, height = 512
  )
}

library(gifski)
pngFiles <- Sys.glob("zzpic*.png")
gifski(
  pngFiles,
  "KleinFiboRt5.gif",
  width = 512, height = 512,
  delay = 1/10
)


file.remove(pngFiles)

