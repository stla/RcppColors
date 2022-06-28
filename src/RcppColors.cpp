#include <RcppColors.h>

std::string _hsluv2hex_(double, double, double);
Rcpp::IntegerVector _hsluv2rgb_(double, double, double);
std::string opacity(double);
std::string rgb_to_hex(int, int, int);

// [[Rcpp::export]]
Rcpp::IntegerVector hsluv2rgb_cpp(double h, double s, double l) {
  return _hsluv2rgb_(h, s, l);
}

// [[Rcpp::export]]
std::string hsluv_cpp(double h, double s, double l) {
  return _hsluv2hex_(h, s, l);
}

// [[Rcpp::export]]
std::string hsluv_alpha_cpp(double h, double s, double l, double alpha) {
  return _hsluv2hex_(h, s, l) + opacity(alpha);
}

std::string rgb2hex(double r, double g, double b) {
  int ri = (int)round(r);
  int gi = (int)round(g);
  int bi = (int)round(b);
  return rgb_to_hex(ri, gi, bi);
}

cplx fromCplx(Rcomplex zr) {
  cplx z(zr.r, zr.i);
  return z;
}

double modulo2(double a, double p) {
  double i = a > 0 ? std::floor(a / p) : std::ceil(a / p);
  return a - i * p;
}

std::string colormap1(cplx z, std::string nancolor) {
  double x = z.real();
  double y = z.imag();
  if(std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) {
    return nancolor;
  }
  double a = std::arg(z);
  double r = modulo2(std::abs(z), 1.0);
  double g = fabs(modulo2(a, 0.5)) * 2.0;
  double b = fabs(modulo2(x * y, 1));
  if(std::isnan(b)) {
    return nancolor;
  }
  return rgb2hex((int)lround((1.0 - cos(r - 0.5)) * 8.0 * 255.0),
                 (int)lround((1.0 - cos(g - 0.5)) * 8.0 * 255.0),
                 (int)lround((1.0 - cos(b - 0.5)) * 8.0 * 255.0));
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix ColorMap1(Rcpp::ComplexMatrix Z,
                                std::string bkgcolor,
                                std::string nancolor) {
  const int m = Z.nrow();
  const int n = Z.ncol();
  Rcpp::CharacterMatrix P(m, n);
  for(int j = 0; j < n; j++) {
    Rcpp::CharacterVector Pj(m);
    const Rcpp::ComplexVector Zj = Z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      if(Rcpp::ComplexVector::is_na(Zj(i))) {
        Pj(i) = bkgcolor;
      } else {
        Pj(i) = colormap1(fromCplx(Zj(i)), nancolor);
      }
    }
    P(Rcpp::_, j) = Pj;
  }
  return P;
}