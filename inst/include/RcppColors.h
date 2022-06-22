#ifndef RCPPCOLORS_H
#define RCPPCOLORS_H

#include <Rcpp.h>

namespace RcppColors {
  std::string rgb2hex(double, double, double);
  Rcpp::IntegerVector hsluv2rgb(double, double, double);
  std::string hsluv2hex(double, double, double);
}  // namespace RcppColors

#endif