# RcppColors

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/RcppColors/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/RcppColors/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Two C++ functions are exposed by this package:

```cpp
std::string rgb2hex(double r, double g, double b);
std::string hsluv2hex(double h, double s, double l);
```