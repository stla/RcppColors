library(RcppColors)

#f[z_] := z^3; DensityPlot[Sin[20Pi Abs[f[x + I y]]], {x, -2.5, 2.5}, {y, -2.5, 2.5}

f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  sin(22 * pi / z^1.5)
}


x <- y <- seq(-15, 15, length.out = 1024)
Z <- outer(x, y, f)


image <- colorMap5(Z, bkgcolor = "#002240")

opar <- par(mar = c(0,0,0,0), bg = "#002240")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)

# anim 1: exponent 1->1.5
f <- function(alpha) {
  function(x, y) {
    z <- complex(real = x, imaginary = y)
    sin(22 * pi / z^alpha)
  }
}
x <- y <- seq(-15, 15, length.out = 1024L)
alpha_ <- seq(1, 1.5, length.out = 50L)

for(i in seq_along(alpha_)) {
  print(i)
  Z <- outer(x, y, f(alpha_[i]))
  image <- colorMap5(Z, bkgcolor = "#002240")
  svg("x.svg")
  opar <- par(mar = c(0,0,0,0), bg = "#002240")
  plot(c(-100, 100), c(-100, 100), type = "n", xaxs = "i", yaxs = "i", 
       xlab = NA, ylab = NA, axes = FALSE, asp = 1)
  rasterImage(image, -100, -100, 100, 100)
  par(opar)
  dev.off()
  rsvg::rsvg_png(
    "x.svg", sprintf("zzpic%03d.png", i), width = 512, height = 512
  )
}

library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  c(pngs, rev(pngs)),
  "sin_z_pow_alpha.gif",
  width = 512, height = 512,
  delay = 1/12
)

file.remove(pngs)

# anim 2: zoom exponent 1.5