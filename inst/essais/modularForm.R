library(RcppColors)

iota <- function(z){
  (z + 1i) / (1i*z + 1)
}
f <- function(z){
  q <- exp(2i * pi * z)
  r <- q - 4*q^2 + 2*q^3 + 8*q^4 - 5*q^5 - 8*q^6 + 6*q^7 - 23*q^9
  r / Mod(r)
}
g <- function(z){
  ifelse(
    Mod(z) >= 1, 
    NA_complex_,
    f(iota(Conj(z)))
  )
}

x <- y <- seq(-1, 1, len = 5000)
W <- outer(y, x, function(x, y) complex(real = x, imaginary = y))
Z <- g(W)
image <- colorMap1(Z, bkgcolor = "#002240", nthreads = 1L)

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)

svg("ModularForm.svg")
opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png(
  "ModularForm.svg", "ModularForm.png", width = 512, height = 512
)

