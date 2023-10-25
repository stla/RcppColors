library(jacobi)
library(RcppColors)

# the modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i*z) / (1i - z)
}

# background color
bkgcol <- "#ffffff"


# make the color mapping
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  tau <- PhiInv(z)
  ifelse(
    Mod(z) >= 1, 
    NA_complex_,
    wp(tau/2, omega=c(0.5, 0.5*tau))
  )
}
x <- y <- seq(-1, 1, length.out = 1024L)
W <- outer(x, y, Vectorize(f))

img <- colorMap5(W) # 7, 14 is nice as well (for Im(W)^2 / Mod(W)^2)

opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)



svg("x.svg", width = 16, height = 16)
opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)
dev.off()

rsvg::rsvg_png(
  "x.svg", "wp_Cayley_cm14.png", width = 512, height = 512
)


# animation Möbius ####

# Möbius transformation of order 2
Mob <- function(z, gamma = -0.5 + 0.5i, t){
  h <- sqrt(1-Mod(gamma)^2)
  d2 <- h^t * (cos(t*pi/2) + 1i*sin(t*pi/2))
  d1 <- Conj(d2)
  a <- Re(d1) - 1i*Im(d1)/h
  b <- gamma * Im(d2)/h
  c <- Conj(b)
  d <- Conj(a)
  (a*z + b) / (c*z + d)
}

f <- function(x, y, t) {
  z <- complex(real = x, imaginary = y)
  z <- Mob(z, t = t)
  tau <- PhiInv(z)
  ifelse(
    Mod(z) >= 1, 
    NA_complex_,
    wp(tau/2, omega=c(0.5, 0.5*tau))
  )
}


s <- function(x) x

t_ <- seq(0, 2, length.out = 61)[-1L]

for(i in seq_along(t_)) {
  print(i)
  #MobW <- Mob(W, t = s(t_[i])) # c pas ça qu'il faut faire... il faut refaire un outer
  W <- outer(x, y, Vectorize(f), t = s(t_[i]))
  img <- colorMap14(Im(W)^2 / Mod(W)^2)#, bkgcolor = bkgcol)
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
  "wp_Cayley.gif",
  width = 512, height = 512,
  delay = 1/8
)

file.remove(pngs)

