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

std::string colormap1(cplx z,
                      std::string nancolor,
                      bool revr,
                      bool revg,
                      bool revb) {
  double x = z.real();
  double y = z.imag();
  if(std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) {
    return nancolor;
  }
  double alpha = std::arg(z);
  double u = modulo2(std::abs(z), 1.0);
  double v = fabs(modulo2(alpha, 0.5)) * 2.0;
  double w = fabs(modulo2(x * y, 1));
  if(std::isnan(w)) {
    return nancolor;
  }
  double r = (1.0 - cos(u - 0.5)) * 8.0;
  double g = (1.0 - cos(v - 0.5)) * 8.0;
  double b = (1.0 - cos(w - 0.5)) * 8.0;
  if(revr) {
    r = 1.0 - r;
  }
  if(revg) {
    g = 1.0 - g;
  }
  if(revb) {
    b = 1.0 - b;
  }
  return rgb2hex((int)lround(r * 255.0), 
                 (int)lround(g * 255.0),
                 (int)lround(b * 255.0));
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix ColorMap1(Rcpp::ComplexMatrix Z,
                                std::string bkgcolor,
                                std::string nancolor,
                                bool revr,
                                bool revg,
                                bool revb) {
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
        Pj(i) = colormap1(fromCplx(Zj(i)), nancolor, revr, revg, revb);
      }
    }
    P(Rcpp::_, j) = Pj;
  }
  return P;
}

std::string colormap2(cplx z,
                      std::string nancolor,
                      bool revh,
                      bool revs,
                      bool revl) {
  double x = z.real();
  double y = z.imag();
  if(std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) {
    return nancolor;
  }
  double arg = std::arg(z);
  if(arg < 0) {
    arg += 2.0 * M_PI;
  }
  double h = arg * 57.29577951308232087680; /* (180 / pi) */
  double w = 2 * M_PI * log1p(std::abs(z));
  double s = 100.0 * sqrt((1.0 + sin(w)) / 2.0);
  double l = 100.0 * (1.0 + cos(w)) / 2.0;
  if(revh) {
    h = 360.0 - h;
  }
  if(revs) {
    s = 100.0 - s;
  }
  if(revl) {
    l = 100.0 - l;
  }
  return _hsluv2hex_(h, s, l);
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix ColorMap2(Rcpp::ComplexMatrix Z,
                                std::string bkgcolor,
                                std::string nancolor,
                                bool revh,
                                bool revs,
                                bool revl) {
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
        Pj(i) = colormap2(fromCplx(Zj(i)), nancolor, revh, revs, revl);
      }
    }
    P(Rcpp::_, j) = Pj;
  }
  return P;
}
