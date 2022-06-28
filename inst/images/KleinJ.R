library(RcppColors)
library(jacobi)

f <- Vectorize(function(q){
  if(Mod(q) > 1){
    z <- NA_complex_
  }else{
    tau <- -1i * log(q) / pi
    if(Im(tau) <= 0){
      z <- NA_complex_
    }else{
      z <- kleinj(tau)
    } 
  }
  colorMap2(
    as.matrix(z), bkgcolor = "#002240"
  )
})


x <- y <- seq(-1, 1, len = 1000)
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


