# RcppColors

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/RcppColors/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/RcppColors/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Four C++ functions are exposed by this package:

```cpp
std::string rgb2hex(double r, double g, double b);
std::string rgba2hex(double r, double g, double b, double a);
std::string hsluv2hex(double h, double s, double l);
std::string hsluv2hex(double h, double s, double l, double alpha);
```

```
r, g, b ∈ [0, 255] (red, green, blue)
a, alpha ∈ [0, 1] (opacity)
h ∈ [0, 360] (hue)
s,l ∈ [0, 100] (saturation, lightness)
```