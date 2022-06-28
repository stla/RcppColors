library(RcppColors)
ikeda <- Vectorize(function(x, y, tau0 = 0, gamma = 2.5){
  for(k in 1L:5L){
    tau <- tau0 - 6.0/(1.0 + x*x + y*y)
    newx <- 0.97 + gamma * (x*cos(tau) - y*sin(tau))
    y <- gamma * (x*sin(tau)+y*cos(tau))
    x <- newx
  }
  z <- complex(real = x, imaginary = y)
  colorMap1(z, reverse = c(TRUE, FALSE, FALSE))
})


x <- y <- seq(-3, 3, len = 3000)
image <- outer(y, x, function(x, y) ikeda(x, y))

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)


svg("Ikeda.svg")
opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png(
  "Ikeda.svg", "Ikeda.png", width = 512, height = 512
)

