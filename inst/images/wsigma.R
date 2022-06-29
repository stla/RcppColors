library(RcppColors)
library(jacobi)

f <- Vectorize(function(z){
  wsigma(z, tau = 1i + Mod(z))
})

x <- y <- seq(-2, 2, len = 240)
Z <- outer(y, x, function(x, y){
  f(complex(real = x, imaginary = y)) 
})

image <- colorMap1(Z, bkgcolor = "#002240")

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)


