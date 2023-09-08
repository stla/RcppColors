#include <RcppColors.h>
#ifdef _OPENMP
#include <omp.h>
#endif

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


// [[Rcpp::export]]
std::string hsl_cpp(double h, double s, double l) {
  return RcppColors::hsl2hex(h, s, l);
}

// [[Rcpp::export]]
std::string hsl_alpha_cpp(double h, double s, double l, double alpha) {
  return RcppColors::hsl2hex(h, s, l) + opacity(alpha);
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


/* ------------------------------------------------------------------------- */
/* ~~~~~~ color map 1 ~~~~~~ */
/* ------------------------------------------------------------------------- */

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
                                bool revb,
                                const unsigned int nthreads) {
  const int m = Z.nrow();
  const int n = Z.ncol();
  Rcpp::CharacterMatrix P(m, n);
  
  if(nthreads == 1) {
    Rcpp::CharacterVector Pj(m);
    for(int j = 0; j < n; j++) {
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
  } else {
    Rcomplex zij;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) collapse(2) private(zij)
#endif
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < m; i++) {
        zij = Z(i,j); 
        if(Rcpp::ComplexVector::is_na(zij)) {
          P(i,j) = bkgcolor;
        } else {
          P(i,j) = colormap1(fromCplx(zij), nancolor, revr, revg, revb);
        }
      }
    }
  }

  return P;
}


/* ------------------------------------------------------------------------- */
/* ~~~~~~ color map 2 ~~~~~~ */
/* ------------------------------------------------------------------------- */

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
  double h = arg * 57.29577951308232087680; // (180 / pi) //
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
                                bool revl,
                                const unsigned int nthreads) {
  const int m = Z.nrow();
  const int n = Z.ncol();
  Rcpp::CharacterMatrix P(m, n);
  
  if(nthreads == 1) {
    Rcpp::CharacterVector Pj(m);
    for(int j = 0; j < n; j++) {
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
  } else {
    Rcomplex zij;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) collapse(2) private(zij)
#endif
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < m; i++) {
        zij = Z(i,j); 
        if(Rcpp::ComplexVector::is_na(zij)) {
          P(i,j) = bkgcolor;
        } else {
          P(i,j) = colormap2(fromCplx(zij), nancolor, revh, revs, revl);
        }
      }
    }
  }
  
  return P;
}


/* ------------------------------------------------------------------------- */
/* ~~~~~~ color map 3 ~~~~~~ */
/* ------------------------------------------------------------------------- */

double perFract(double x, double t, double m, double M) {
  x = x / t;
  return m + (M - m) * (x - std::floor(x));  
}


std::string colormap3(cplx z,
                      std::string nancolor,
                      double s,
                      double r) {
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
  double ph = perFract(h, 360.0 / r, 216.0, 360.0) / 360.0;
  double plogm = perFract(log1p(std::abs(z)), 2 * M_PI / r, 0.6, 1.0);
  double l = 100.0 * ph * plogm;
  return _hsluv2hex_(h, s, l);
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix ColorMap3(Rcpp::ComplexMatrix Z,
                                std::string bkgcolor,
                                std::string nancolor,
                                double s,
                                double r,
                                const unsigned int nthreads) {
  const int m = Z.nrow();
  const int n = Z.ncol();
  Rcpp::CharacterMatrix P(m, n);
  
  if(nthreads == 1) {
    Rcpp::CharacterVector Pj(m);
    for(int j = 0; j < n; j++) {
      const Rcpp::ComplexVector Zj = Z(Rcpp::_, j);
      for(int i = 0; i < m; i++) {
        if(Rcpp::ComplexVector::is_na(Zj(i))) {
          Pj(i) = bkgcolor;
        } else {
          Pj(i) = colormap3(fromCplx(Zj(i)), nancolor, s, r);
        }
      }
      P(Rcpp::_, j) = Pj;
    }
  } else {
    Rcomplex zij;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) collapse(2) private(zij)
#endif
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < m; i++) {
        zij = Z(i,j); 
        if(Rcpp::ComplexVector::is_na(zij)) {
          P(i,j) = bkgcolor;
        } else {
          P(i,j) = colormap3(fromCplx(zij), nancolor, s, r);
        }
      }
    }
  }
  
  return P;
}


/* ------------------------------------------------------------------------- */
/* ~~~~~~ color map 4 ~~~~~~ */
/* ------------------------------------------------------------------------- */

double b(double r) {
  double out;
  if(r == 0) {
    out = 0;
  } else if(std::isinf(r)) {
    out = 1;
  } else {
    out = atan(log(r)) * M_1_PI + 0.5;
  }
  return out;
}

std::string colormap4(cplx z,
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
  double bmz = b(std::abs(z));
  double s = 100 * (1 - bmz * bmz);
  double l = 100 * (1 - (1 - bmz)*(1 - bmz));
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
Rcpp::CharacterMatrix ColorMap4(Rcpp::ComplexMatrix Z,
                                std::string bkgcolor,
                                std::string nancolor,
                                bool revh,
                                bool revs,
                                bool revl,
                                const unsigned int nthreads) {
  const int m = Z.nrow();
  const int n = Z.ncol();
  Rcpp::CharacterMatrix P(m, n);
  
  if(nthreads == 1) {
    Rcpp::CharacterVector Pj(m);
    for(int j = 0; j < n; j++) {
      const Rcpp::ComplexVector Zj = Z(Rcpp::_, j);
      for(int i = 0; i < m; i++) {
        if(Rcpp::ComplexVector::is_na(Zj(i))) {
          Pj(i) = bkgcolor;
        } else {
          Pj(i) = colormap4(fromCplx(Zj(i)), nancolor, revh, revs, revl);
        }
      }
      P(Rcpp::_, j) = Pj;
    }
  } else {
    Rcomplex zij;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) collapse(2) private(zij)
#endif
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < m; i++) {
        zij = Z(i,j); 
        if(Rcpp::ComplexVector::is_na(zij)) {
          P(i,j) = bkgcolor;
        } else {
          P(i,j) = colormap4(fromCplx(zij), nancolor, revh, revs, revl);
        }
      }
    }
  }
  
  return P;
}


// -------------------------------------------------------------------------- //
// ~~~~~~ color map 5 ~~~~~~ //
// -------------------------------------------------------------------------- //

std::string colormap5(cplx z,
                      std::string nancolor,
                      bool revh,
                      bool revs,
                      bool revi) {
  double x = z.real();
  double y = z.imag();
  if(std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) {
    return nancolor;
  }
  double arg = std::arg(z);
  if(arg < 0) {
    arg += M_PI;
  }
  double h = std::max(0.0, std::min(360.0 * arg/2/M_PI, 359.999));
  double w = 2 * M_PI * std::log1p(fabs(arg));
  double s = 100 * sqrt(( 1.0 + sin(w) ) / 2.0);
  double i = 100 * ( 1.0 + cos(w) ) / 2.0;
  if(revh) {
    h = 360.0 - h;
  }
  if(revs) {
    s = 100.0 - s;
  }
  if(revi) {
    i = 100.0 - i;
  }
  return RcppColors::hsi2hex(h, s, i);
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix ColorMap5(Rcpp::ComplexMatrix Z,
                                std::string bkgcolor,
                                std::string nancolor,
                                bool revh,
                                bool revs,
                                bool revi,
                                const unsigned int nthreads) {
  const int m = Z.nrow();
  const int n = Z.ncol();
  Rcpp::CharacterMatrix P(m, n);
  
  if(nthreads == 1) {
    Rcpp::CharacterVector Pj(m);
    for(int j = 0; j < n; j++) {
      const Rcpp::ComplexVector Zj = Z(Rcpp::_, j);
      for(int i = 0; i < m; i++) {
        if(Rcpp::ComplexVector::is_na(Zj(i))) {
          Pj(i) = bkgcolor;
        } else {
          Pj(i) = colormap5(fromCplx(Zj(i)), nancolor, revh, revs, revi);
        }
      }
      P(Rcpp::_, j) = Pj;
    }
  } else {
    Rcomplex zij;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) collapse(2) private(zij)
#endif
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < m; i++) {
        zij = Z(i,j); 
        if(Rcpp::ComplexVector::is_na(zij)) {
          P(i,j) = bkgcolor;
        } else {
          P(i,j) = colormap5(fromCplx(zij), nancolor, revh, revs, revi);
        }
      }
    }
  }
  
  return P;
}


/* ------------------------------------------------------------------------- */
/* ~~~~~~ color map 6 ~~~~~~ */
/* ------------------------------------------------------------------------- */

std::string colormap6(cplx z,
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
  double h = arg * 57.29577951308232087680; // (180 / pi) //
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
  return RcppColors::hsl2hex(h, s, l);
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix ColorMap6(Rcpp::ComplexMatrix Z,
                                std::string bkgcolor,
                                std::string nancolor,
                                bool revh,
                                bool revs,
                                bool revl,
                                const unsigned int nthreads) {
  const int m = Z.nrow();
  const int n = Z.ncol();
  Rcpp::CharacterMatrix P(m, n);
  
  if(nthreads == 1) {
    Rcpp::CharacterVector Pj(m);
    for(int j = 0; j < n; j++) {
      const Rcpp::ComplexVector Zj = Z(Rcpp::_, j);
      for(int i = 0; i < m; i++) {
        if(Rcpp::ComplexVector::is_na(Zj(i))) {
          Pj(i) = bkgcolor;
        } else {
          Pj(i) = colormap6(fromCplx(Zj(i)), nancolor, revh, revs, revl);
        }
      }
      P(Rcpp::_, j) = Pj;
    }
  } else {
    Rcomplex zij;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) collapse(2) private(zij)
#endif
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < m; i++) {
        zij = Z(i,j); 
        if(Rcpp::ComplexVector::is_na(zij)) {
          P(i,j) = bkgcolor;
        } else {
          P(i,j) = colormap6(fromCplx(zij), nancolor, revh, revs, revl);
        }
      }
    }
  }
  
  return P;
}


/* ------------------------------------------------------------------------- */
/* ~~~~~~ color map 7 ~~~~~~ */
/* ------------------------------------------------------------------------- */

std::string colormap7(cplx z,
                      std::string nancolor,
                      bool revh,
                      bool revs,
                      bool revv) {
  double x = z.real();
  double y = z.imag();
  if(std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) {
    return nancolor;
  }
  double arg = std::arg(z);
  if(arg < 0) {
    arg += 2.0 * M_PI;
  }
  double h = arg * 57.29577951308232087680; // (180 / pi) //
  double w = 2 * M_PI * log1p(std::abs(z));
  double s = 100.0 * sqrt((1.0 + sin(w)) / 2.0);
  double v = 100.0 * (1.0 + cos(w)) / 2.0;
  if(revh) {
    h = 360.0 - h;
  }
  if(revs) {
    s = 100.0 - s;
  }
  if(revv) {
    v = 100.0 - v;
  }
  return RcppColors::hsv2hex(h, s, v);
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix ColorMap7(Rcpp::ComplexMatrix Z,
                                std::string bkgcolor,
                                std::string nancolor,
                                bool revh,
                                bool revs,
                                bool revv,
                                const unsigned int nthreads) {
  const int m = Z.nrow();
  const int n = Z.ncol();
  Rcpp::CharacterMatrix P(m, n);
  
  if(nthreads == 1) {
    Rcpp::CharacterVector Pj(m);
    for(int j = 0; j < n; j++) {
      const Rcpp::ComplexVector Zj = Z(Rcpp::_, j);
      for(int i = 0; i < m; i++) {
        if(Rcpp::ComplexVector::is_na(Zj(i))) {
          Pj(i) = bkgcolor;
        } else {
          Pj(i) = colormap7(fromCplx(Zj(i)), nancolor, revh, revs, revv);
        }
      }
      P(Rcpp::_, j) = Pj;
    }
  } else {
    Rcomplex zij;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) collapse(2) private(zij)
#endif
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < m; i++) {
        zij = Z(i,j); 
        if(Rcpp::ComplexVector::is_na(zij)) {
          P(i,j) = bkgcolor;
        } else {
          P(i,j) = colormap7(fromCplx(zij), nancolor, revh, revs, revv);
        }
      }
    }
  }
  
  return P;
}
