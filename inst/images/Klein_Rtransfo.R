library(RcppColors)
library(jacobi)

# a modular transformation
R <- function(z, t) {
  a <- pi*t/3
  ((sqrt(3)*cos(a) - sin(a)) * z - 2*sin(a))/
    (2*sin(a) * z + sqrt(3)*cos(a) + sin(a))
}

# the modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i*z) / (1i - z)
}
# background color
bkgcol <- "gold3"
# make the color mapping
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  w <- PhiInv(z)
  ifelse(
    Mod(z) > 0.96, 
    NA_complex_,
    ifelse(
      y < 0, 
      -1/w, w
    )
  )
}
x <- seq(-1, 1, length.out = 2048)
y <- seq(-1, 1, length.out = 2048)
Z <- outer(x, y, f)
K <- kleinj(Z) / 1728


xmsin <- function(x) x - sin(x)


s <- function(x) {
  (2*x - sin(2*x) + 4*pi) / (6*pi)  
}

xmsin <- function(x) x - sin(x)
s <- function(x) {
  (xmsin(xmsin(2*x)) + 4*pi) / (2*pi)  
}

t_ <- seq(-2*pi, pi, length.out= 91L)[-1L]
for(i in 60:60) {
  KRt <- R(K, s(t_[i]))
  image <- colorMap5(KRt, bkgcolor = bkgcol)
  svglite::svglite("x.svg", width = 10, height = 10)
  opar <- par(mar = c(0,0,0,0), bg = bkgcol)
  plot(
    c(-1, 1), c(-1, 1), type = "n", xaxs = "i", yaxs = "i", 
    xlab = NA, ylab = NA, axes = FALSE, asp = 1
  )
  rasterImage(image, -1, -1, 1, 1)
  par(opar)
  dev.off()
  rsvg::rsvg_png(
    "x.svg", sprintf("vvpic%03d.png", i), width = 1024, height = 1024
  )
}

library(gifski)
pngFiles <- Sys.glob("vvpic*.png")
gifski(
  pngFiles,
  "KleinRt5_smooth.gif",
  width = 512, height = 512,
  delay = 1/10
)


file.remove(pngFiles)

# with Dedekind ####
library(PlaneGeometry)
i <- 60L
KRt <- R(K, s(t_[i]))
image <- colorMap5(KRt, bkgcolor = "gold3")
svg("x.svg")
opar <- par(mar = c(0,0,0,0), bg = "gold3")
plot(
  c(-1, 1), c(-1, 1), type = "n", xaxs = "i", yaxs = "i", 
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(t(image), -1, -1, 1, 1)
# now we add the Dedekind tessellation (the white lines)
isInteger <- function(x) abs(x - floor(x)) < x * 1e-6
abline(h = 0, col = "gold3", lwd = 2)
N <- 50L
for(n in 1L:N) {
  if(isInteger(n/2) && ((n/2L) %% 2L == 1L)) {
    next
  }
  for(p in 1:n) {
    q <- sqrt(n*n - p*p + 4L)
    cases <- (isInteger(q) && isInteger(q/2) && (n %% 2L == 1L)) ||
      (isInteger(q) && isInteger(q/4) && (n %% 4L == 0L))
    if(cases) {
      circ <- Circle$new(center = c(q, p)/n, radius = 2/n)
      draw(circ, border = "gold3", lwd = 2)
      circ <- Circle$new(center = c(-q, p)/n, radius = 2/n)
      draw(circ, border = "gold3", lwd = 2)
      circ <- Circle$new(center = c(q, -p)/n, radius = 2/n)
      draw(circ, border = "gold3", lwd = 2)
      circ <- Circle$new(center = c(-q, -p)/n, radius = 2/n)
      draw(circ, border = "gold3", lwd = 2)
    }
  }
}
circ <- Circle$new(center = c(0, 0), radius = 0.96)
draw(circ, border ="black", lwd=2)
par(opar)
dev.off()

rsvg::rsvg_png(
  "x.svg", "KleinDedekind.png", width = 512, height = 512
)
