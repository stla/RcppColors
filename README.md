# RcppColors

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/RcppColors/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/RcppColors/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Six C++ functions are exposed by this package:

```cpp
std::string rgb2hex(double r, double g, double b);
std::string rgba2hex(double r, double g, double b, double a);
std::string hsluv2hex(double h, double s, double l);
std::string hsluv2hex(double h, double s, double l, double alpha);
std::string hsv2hex(double h, double s, double v);
std::string hsv2hex(double h, double s, double v, double alpha);
```

```
r, g, b ∈ [0, 255] (red, green, blue)
a, alpha ∈ [0, 1] (opacity)
h ∈ [0, 360] (hue)
s,l,v ∈ [0, 100] (saturation, lightness, value)
```

## Usage in a package with **Rcpp**

The **LinkingTo** field in the **DESCRIPTION** file should look like

```yaml
LinkingTo: 
    Rcpp,
    RcppColors
```

Then, in your **C++** file, you can call the above functions like this:

```cpp
#include <RcppColors.h>

std::string mycolor = RcppColors::rgb2hex(0.0, 128.0, 255.0);
```

## Color maps

```r
library(RcppColors)
library(Bessel)
x <- y <- seq(-4, 4, len = 1500)
# complex grid
W <- outer(y, x, function(x, y) complex(real = x, imaginary = y))
# computes Bessel values
Z <- matrix(BesselY(W, nu = 3), nrow = nrow(W), ncol = ncol(W))
# maps them to colors
image <- colorMap1(Z)
# plot
opar <- par(mar = c(0,0,0,0), bg = "#15191E")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)
```

![](https://raw.githubusercontent.com/stla/RcppColors/main/inst/images/BesselY.png)


```r
library(RcppColors)
library(Carlson)
library(rgl)
library(Rvcg)

mesh <- vcgSphere(subdivision = 8)

color <- apply(mesh$vb[-4L, ], 2L, function(xyz){
  if(sum(xyz == 0) >= 2){
    z <- as.matrix(NA_complex_)
  }else{
    a <- xyz[1]
    b <- xyz[2]
    c <- xyz[3]
    z <- as.matrix(Carlson_RJ(a, b, c, 1i, 1e-5))
  }
  colorMap1(z)
})

mesh$material <- list(color = color)

open3d(windowRect = c(50, 50, 562, 562), zoom = 0.75)
bg3d("whitesmoke")
shade3d(mesh)
```

![](https://raw.githubusercontent.com/stla/RcppColors/main/inst/images/CarlsonBall.gif)

