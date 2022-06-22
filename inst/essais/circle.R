x <- y <- seq(-1, 1, length.out = 100L)
f <- Vectorize(
  function(x, y){
    z <- complex(real = x, imaginary = y)
    radians <- Arg(z)
    if(radians < 0){
      radians <- radians + 2*pi
    }
    degrees <- 360 * radians / 2 / pi
    modulus <- Mod(z)
    hsluv(h = degrees, s = modulus, l = 50)
  }
)
image <- outer(x, y, f)

opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1))
rasterImage(image, -1, -1, 1, 1)
par(opar)

