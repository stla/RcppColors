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
