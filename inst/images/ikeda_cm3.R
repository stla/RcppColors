library(RcppColors)

ikeda <- Vectorize(function(x, y, tau0 = 0, gamma = 2.5){
  for(k in 1L:5L){
    tau <- tau0 - 6.0/(1.0 + x*x + y*y)
    newx <- 0.97 + gamma * (x*cos(tau) - y*sin(tau))
    y <- gamma * (x*sin(tau)+y*cos(tau))
    x <- newx
  }
  complex(real = x, imaginary = y)
})

x <- y <- seq(-3, 3, len = 2048)
Z <- outer(y, x, function(x, y) ikeda(x, y))
image <- colorMap3(Z, s = 80, n = 10)


opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n",
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)
