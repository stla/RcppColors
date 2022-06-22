saturation <- 100
f <- Vectorize(
  function(x, y){
    z <- complex(real = x, imaginary = y)
    modulus <- Mod(z)
    if(modulus > 1){
      return("#ffffff")
    }
    radians <- Arg(z)
    if(radians < 0){
      radians <- radians + 2*pi
    }
    degrees <- 360 * radians / 2 / pi
    hsluv(h = degrees, s = saturation, l = 100*modulus)
  }
)

x <- y <- seq(-1, 1, length.out = 200L)
image <- outer(x, y, f)

opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
rasterImage(image, -1, -1, 1, 1)
par(opar)

