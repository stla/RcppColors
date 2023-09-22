library(RcppColors)

Klein <- function(z) {
  1728 * (z * (z**10 + 11 * z**5 - 1))**5 / 
    (-(z**20 + 1) + 228 * (z**15 - z**5) - 494 * z**10)**3
}

Mobius <- function(z) {
  w <- (z - 20)/(3*z + 1i)
  w / Mod(w) / 2
}

Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i*z) / (1i - z)
}

library(jacobi)
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  kleinj(PhiInv(Mobius(kleinj(z))))
}

x <- y <- seq(1, 5, len = 512)

Z <- outer(y, x, f)

image <- colorMap5(Z)

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", xaxs = "i", yaxs = "i",
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)





Mobius <- function(gamma, t, z){
  mgamma <- Mod(gamma)
  h <- sqrt(1.0 - mgamma * mgamma)
  z1 <- complex(real = cos(t * pi / 2.0), imag = sin(t * pi / 2.0))
  d2 <- h^t * z1
  d1 <- Conj(d2)
  a <- complex(real = Re(d1), imag = -Im(d1) / h)
  b <- Im(d2) * gamma / h
  c <- Conj(b)
  d <- Conj(a)
  (a * z + b)/(c * z + d)
}

f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  k <- kleinj(z / 1728)
  k <- k / Mod(k) / (Arg(k)+2*pi)/(2*pi)
  kleinj(PhiInv(Mobius(0.5+0.5i, 1, k))/ 1728)
}

x <- y <- seq(1, 5, len = 512)

Z <- outer(y, x, f)

image <- colorMap5(Z, reverse = c(T,T,F))

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", xaxs = "i", yaxs = "i",
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)



f <- function(z, t) {
  Klein(Mobius(0.5 + 0.5i, t, Klein(z)))
}

x <- y <- seq(-6, 6, len = 1024)

t_ <- seq(0, 2, length.out = 81)[-1L]

for(i in seq_along(t_)) {
  Z <- outer(y, x, function(x, y){
    f(complex(real = x, imaginary = y), t_[i])
  })
  image <- colorMap3(Z, bkgcolor = "#002240", s = 100, n = 3)
  svg("zpic.svg")
  opar <- par(mar = c(0,0,0,0), bg = "#002240")
  plot(c(-100, 100), c(-100, 100), type = "n", xaxs = "i", yaxs = "i", 
       xlab = "", ylab = "", axes = FALSE, asp = 1)
  rasterImage(image, -100, -100, 100, 100)
  par(opar)
  dev.off()
  fl <- sprintf("zpic%03d.png", i)
  rsvg::rsvg_png(
    "zpic.svg", fl, width = 512, height = 512
  )
}

pics <- Sys.glob("zpic*.png")
library(gifski)
gifski(
  png_files = pics, gif_file = "kaleidoscope.gif", 
  width = 512, height = 512, delay = 1/2)


opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", xaxs = "i", yaxs = "i",
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)


