library(RcppColors)
library(jacobi)

f <- Vectorize(function(z){
  wsigma(z, omega = c(1, 0.25 + 1i))
})

x <- y <- seq(-5, 5, len = 512)
Z <- outer(y, x, function(x, y){
  f(complex(real = x, imaginary = y)) 
})

image <- colorMap5(Z, bkgcolor = "#002240")

image[which(nchar(image)==8,arr.ind=TRUE)] <- "#ff0000"


svg("x.svg")
opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = NA, ylab = NA, axes = FALSE, asp = 1, xaxs="i", yaxs="i"
)
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png("x.svg", "wsigma_cm5.png", width = 512, height = 512)
