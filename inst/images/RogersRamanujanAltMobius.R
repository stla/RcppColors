# Alternating Rogers-Ramanujan map ####
library(jacobi)
library(RcppColors)

# the modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i*z) / (1i - z)
}

# background color
bkgcol <- "#ffffff"

# Möbius transformation of order 3
Mob <- function(z, t) {
  a <- pi*t/3
  ((sqrt(3)*cos(a) - sin(a)) * z - 2*sin(a))/
    (2*sin(a) * z + sqrt(3)*cos(a) + sin(a))
}

# smooth stair function
xmsinx <- function(x) x - sin(x)
s <- function(x) {
  (xmsinx(xmsinx(2*x)) + 4*pi) / (2*pi)  
}

# make the color mapping
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  tau <- PhiInv(z)
  ifelse(
    Mod(z) > 0.98, 
    NA_complex_,
    RRa(exp(1i * pi * tau))
  )
}
x <- y <- seq(-1, 1, length.out = 1024L)
R <- outer(x, y, Vectorize(f))^5

# animation ####
t_ <- seq(-2*pi, pi, length.out = 121)[-1L]

for(i in seq_along(t_)) {
  print(i)
  img <- colorMap5(Mob(R, s(t_[i])), bkgcolor = bkgcol)
  img <- permuteRGB(img, "grb")
  svg("x.svg", width = 16, height = 16)
  opar <- par(mar = c(0, 0, 0, 0))
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
       axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
  rasterImage(img, 0, 0, 1, 1)
  par(opar)
  dev.off()
  rsvg::rsvg_png(
    "x.svg", file = sprintf("zzpic%03d.png", i), 
    width = 512, height = 512
  )
}

# mount animation
library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  pngs,
  "RogersRamanujanAltMobius.gif",
  width = 512, height = 512,
  delay = 1/11
)

file.remove(pngs)
