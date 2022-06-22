#include <RcppColors.h>

std::string _hsluv2hex_(double, double, double);
Rcpp::IntegerVector _hsluv2rgb_(double, double, double);

// [[Rcpp::export]]
Rcpp::IntegerVector hsluv2rgb_cpp(double h, double s, double l) {
  return _hsluv2rgb_(h, s, l);
}

// [[Rcpp::export]]
std::string hsluv_cpp(double h, double s, double l) {
  return _hsluv2hex_(h, s, l);
}
