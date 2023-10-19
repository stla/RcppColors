library(RcppColors)
library(pracma)

f <- function(x, y) {
  erfz(complex(real = x, imaginary = y))
}

x <- seq(-6, 6, length.out = 256)
y <- seq(-6, 6, length.out = 256)
Z <- outer(x, y, f)


image <- colorMap1(Z, bkgcolor = "#002240")

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)

