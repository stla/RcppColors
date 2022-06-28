library(RcppColors)
library(jacobi)

f <- Vectorize(function(q){
  if(Mod(q) >= 1){
    z <- NA_complex_
  }else{
    tau <- -1i * log(q) / pi
    if(Im(tau) <= 0){
      z <- NA_complex_
    }else{
      z <- kleinj(tau) / 1728
    } 
  }
  z
})

x <- y <- seq(-1, 1, len = 3000)
Z <- outer(y, x, function(x, y){
  f(complex(real = x, imaginary = y)) 
})
image <- colorMap2(1/Z, bkgcolor = "#002240", reverse = c(T,T,T))

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)


svg("kleinj.svg")
opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png(
  "kleinj.svg", "kleinj.png", width = 512, height = 512
)

