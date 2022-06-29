library(RcppColors)
library(jacobi)
library(raster)


x <- y <- seq(-1, 1, len = 2000)
W <- outer(y, x, function(x, y){
  complex(real = x, imaginary = y)
})
W[Mod(W) >= 1] <- NA_complex_
Tau <- -1i * log(W) / pi
Tau[Im(Tau) <= 0] <- NA_complex_
nas <- is.na(Tau)
Tau[nas] <- 1i
Z <- kleinj(Tau) / 1728
Z[nas] <- NA_complex_
image <- colorMap2(1/Z, bkgcolor = "maroon", reverse = c(T,T,T))

opar <- par(mar = c(0,0,0,0), bg = "maroon")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)

################################################################################
npixels <- 2500L
RASTER <- raster(nrow = npixels, ncol = npixels)
values(RASTER) <- 1L:ncell(RASTER)
extent(RASTER) <- c(-180, 180, -180, 180)
crs(RASTER) <- NULL

x <- y <- seq(-1, 1, len = npixels)
W0 <- outer(y, x, function(x, y){
  complex(real = x, imaginary = y)
})

eps <- 1e-6

Mobius <- function(gamma, t, W){
  mgamma <- Mod(gamma)
  h <- sqrt(1.0 - mgamma * mgamma)
  z1 <- complex(real = cos(t * pi / 2.0), imag = sin(t * pi / 2.0))
  d2 <- h^t * z1
  d1 <- Conj(d2)
  a <- complex(real = Re(d1), imag = -Im(d1) / h)
  b <- Im(d2) * gamma / h
  c <- Conj(b)
  d <- Conj(a)
  W <- (a * W + b)/(c * W + d)
  W[Mod(W) >= 1-eps] <- NA_complex_
  W[abs(Re(W)) < eps & Im(W) <= eps] <- NA_complex_
  W
}

t_ <- head(seq(0, 2, length.out = 361), -1L)
for(i in 295:length(t_)){
  W <- Mobius(gamma = 0.6 + 0.7i, t = t_[i], W0)
  Tau <- -1i * log(W) / pi
  Tau[Im(Tau) <= eps] <- NA_complex_
  Z <- kleinj(Tau, transfo = TRUE) / 1728
  image <- colorMap2(1/Z, bkgcolor = "pink", reverse = c(T,T,T))
  colortable(RASTER) <- c(image)
  svg(sprintf("zzpic%03d.svg", i))
  opar <- par(mar = c(0,0,0,0)+1, bg = "pink")
  plot(RASTER)
  par(opar)
  dev.off()
}


for(i in 1:npixels){
  for(j in 1:npixels){
    z <- kleinj(Tau[i, j])
  }
}


opar <- par(mar = c(0,0,0,0), bg = "pink")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)


colortable(RASTER) <- c(image)








svg("xkleinj.svg")
opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png(
  "xkleinj.svg", "xkleinj.png", width = 512, height = 512
)



svg("xkleinj.svg")
opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(RASTER)
par(opar)
dev.off()

