library(RcppColors)
library(pracma)

FresnelS <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  sqrt(pi/2) * (1+1i)/4 * (erfz((1+1i)/sqrt(2) * z) - 1i * erfz((1-1i)/sqrt(2) * z))
}

x <- seq(-7, 7, length.out = 256)
y <- seq(-7, 7, length.out = 256)
Z <- outer(x, y, FresnelS)


image <- colorMap14(Z/40, bkgcolor = "#002240")

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)

