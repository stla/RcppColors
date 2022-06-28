library(RcppColors)
library(jacobi)

f <- Vectorize(function(q){
  if(Mod(q) > 1  || (Im(q) == 0 && Re(q) <= 0)){
    z <- NA_complex_
  }else{
    z <- EisensteinE(6, q)
  }
  colorMap2(
    as.matrix(z), bkgcolor = "#002240"
  )
})

x <- y <- seq(-1, 1, len = 2000)
image <- outer(y, x, function(x, y){
  f(complex(real = x, imaginary = y)) 
})

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)



svg("E6.svg")
opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png(
  "E6.svg", "E6.png", width = 512, height = 512
)
