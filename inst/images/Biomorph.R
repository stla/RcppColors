library(RcppColors)
library(jacobi)

f <- function(z){
  2 * (z^3-2) / 3
}

biomorph <- Vectorize(function(z){
  zz <- NA_complex_
  for(k in 1L:15L){
    w <- 0.8*(2*z^3/3 - 0.8) + (1 - 0.8)*z
    z <- 0.7*(2*w^3/3 - 0.8) + (1 - 0.7)*w
    if(Mod(z) > 30){
      zz <- z
      break;
    }
  }
  colorMap2(as.matrix(zz))
})

x <- y <- seq(-1.5, 1.5, len = 3000)
image <- outer(y, x, function(x, y){
  biomorph(complex(real = x, imaginary = y)) 
})

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)


svg("biomorph.svg")
opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png(
  "biomorph.svg", "biomorph.png", width = 512, height = 512
)



