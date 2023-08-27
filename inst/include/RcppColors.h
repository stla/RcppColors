#ifndef RCPPCOLORS_H
#define RCPPCOLORS_H

#include <Rcpp.h>

typedef std::complex<double> cplx;

static std::string rgb_to_hex(int r, int g, int b) {
  std::stringstream ss;
  ss << "#";
  ss << std::setfill('0') << std::setw(6) << std::hex << (r << 16 | g << 8 | b);
  return ss.str();
}

typedef struct Triplet_tag Triplet;
struct Triplet_tag {
  double a;
  double b;
  double c;
};

/* for RGB */
static const Triplet m[3] = {
    {3.24096994190452134377, -1.53738317757009345794, -0.49861076029300328366},
    {-0.96924363628087982613, 1.87596750150772066772, 0.04155505740717561247},
    {0.05563007969699360846, -0.20397695888897656435, 1.05697151424287856072}};

/* for XYZ */
static const Triplet m_inv[3] = {
    {0.41239079926595948129, 0.35758433938387796373, 0.18048078840183428751},
    {0.21263900587151035754, 0.71516867876775592746, 0.07219231536073371500},
    {0.01933081871559185069, 0.11919477979462598791, 0.95053215224966058086}};

static const double ref_u = 0.19783000664283680764;
static const double ref_v = 0.46831999493879100370;

static const double kappa = 903.29629629629629629630;
static const double epsilon = 0.00885645167903563082;

static const double ratio_pi_180 = 0.01745329251994329;
static const double ratio_180_pi = 57.29577951308232286;

typedef struct Bounds_tag Bounds;
struct Bounds_tag {
  double a;
  double b;
};

static void get_bounds(double l, Bounds bounds[6]) {
  double tl = l + 16.0;
  double sub1 = (tl * tl * tl) / 1560896.0;
  double sub2 = (sub1 > epsilon ? sub1 : (l / kappa));
  int channel;
  int t;
  for(channel = 0; channel < 3; channel++) {
    double m1 = m[channel].a;
    double m2 = m[channel].b;
    double m3 = m[channel].c;
    for(t = 0; t < 2; t++) {
      double top1 = (284517.0 * m1 - 94839.0 * m3) * sub2;
      double top2 = (838422.0 * m3 + 769860.0 * m2 + 731718.0 * m1) * l * sub2 -
                    769860.0 * t * l;
      double bottom = (632260.0 * m3 - 126452.0 * m2) * sub2 + 126452.0 * t;
      bounds[channel * 2 + t].a = top1 / bottom;
      bounds[channel * 2 + t].b = top2 / bottom;
    }
  }
}

static double ray_length_until_intersect(double theta, const Bounds* line) {
  return line->b / (sin(theta) - line->a * cos(theta));
}

static double max_chroma_for_lh(double l, double h) {
  double min_len = DBL_MAX;
  double hrad = h * ratio_pi_180;
  Bounds bounds[6];
  int i;
  get_bounds(l, bounds);
  for(i = 0; i < 6; i++) {
    double len = ray_length_until_intersect(hrad, &bounds[i]);
    if(len >= 0 && len < min_len) {
      min_len = len;
    }
  }
  return min_len;
}

static double dot_product(const Triplet* t1, const Triplet* t2) {
  return (t1->a * t2->a + t1->b * t2->b + t1->c * t2->c);
}

/* Used for rgb conversions */
static double from_linear(double c) {
  if(c <= 0.0031308) {
    return 12.92 * c;
  } else {
    return 1.055 * pow(c, 1.0 / 2.4) - 0.055;
  }
}

static void xyz_to_rgb(Triplet* in_out) {
  double r = from_linear(dot_product(&m[0], in_out));
  double g = from_linear(dot_product(&m[1], in_out));
  double b = from_linear(dot_product(&m[2], in_out));
  in_out->a = r;
  in_out->b = g;
  in_out->c = b;
}

static double l2y(double l) {
  if(l <= 8.0) {
    return l / kappa;
  } else {
    double x = (l + 16.0) / 116.0;
    return (x * x * x);
  }
}

static void luv_to_xyz(Triplet* in_out) {
  if(in_out->a <= 0.00000001) {
    /* Black will create a divide-by-zero error. */
    in_out->a = 0.0;
    in_out->b = 0.0;
    in_out->c = 0.0;
    return;
  }
  double var_u = in_out->b / (13.0 * in_out->a) + ref_u;
  double var_v = in_out->c / (13.0 * in_out->a) + ref_v;
  double y = l2y(in_out->a);
  double x = -(9.0 * y * var_u) / ((var_u - 4.0) * var_v - var_u * var_v);
  double z = (9.0 * y - (15.0 * var_v * y) - (var_v * x)) / (3.0 * var_v);
  in_out->a = x;
  in_out->b = y;
  in_out->c = z;
}

static void lch_to_luv(Triplet* in_out) {
  double hrad = in_out->c * ratio_pi_180;
  double u = cos(hrad) * in_out->b;
  double v = sin(hrad) * in_out->b;
  in_out->b = u;
  in_out->c = v;
}

static void hsluv_to_lch(Triplet* in_out) {
  double h = in_out->a;
  double s = in_out->b;
  double l = in_out->c;
  double c;
  /* White and black: disambiguate chroma */
  if(l > 99.9999999 || l < 0.00000001) {
    c = 0.0;
  } else {
    c = max_chroma_for_lh(l, h) / 100.0 * s;
  }
  /* Grays: disambiguate hue */
  if(s < 0.00000001) {
    h = 0.0;
  }
  in_out->a = l;
  in_out->b = c;
  in_out->c = h;
}

static void hsluv_to_rgb(double h,
                         double s,
                         double l,
                         double* pr,
                         double* pg,
                         double* pb) {
  Triplet tmp = {h, s, l};
  hsluv_to_lch(&tmp);
  lch_to_luv(&tmp);
  luv_to_xyz(&tmp);
  xyz_to_rgb(&tmp);
  *pr = tmp.a;
  *pg = tmp.b;
  *pb = tmp.c;
}

static Rcpp::IntegerVector _hsluv2rgb_(double h, double s, double l) {
  if(h < 0.0 || h > 360.0) {
    Rcpp::stop("Invalid value of `h`.");
  }
  if(s < 0.0 || s > 100.0) {
    Rcpp::stop("Invalid value of `s`.");
  }
  if(l < 0.0 || l > 100.0) {
    Rcpp::stop("Invalid value of `l`.");
  }
  Rcpp::IntegerVector out(3);
  double r, g, b;
  hsluv_to_rgb(h, s, l, &r, &g, &b);
  out[0] = (int)round(255.0 * r);
  out[1] = (int)round(255.0 * g);
  out[2] = (int)round(255.0 * b);
  return out;
}

static std::string _hsluv2hex_(double h, double s, double l) {
  Rcpp::IntegerVector out = _hsluv2rgb_(h, s, l);
  std::string outhex = (out[0] == 0 && out[1] == 0 && out[2] == 0)
                           ? "#000000"
                           : rgb_to_hex(out[0], out[1], out[2]);
  return outhex;
}

static std::string opacity(double alpha) {
  if(alpha < 0.0 || alpha > 1.0) {
    Rcpp::stop("Invalid value of `alpha`.");
  }
  int alphai = (int)(round(255.0 * alpha));
  std::stringstream ss;
  ss << std::setfill('0') << std::setw(2) << std::hex << alphai;
  return ss.str();
}

static double modulo(double a, double p) {
  double i = std::floor(a / p);
  return a - i * p;
}

// used for hsv -> rgb
static double f(double x, double h, double s, double v) {
  double k = modulo(x + h / 60.0, 6.0);
  double mm = std::max(0.0, std::min(k, std::min(4.0 - k, 1.0)));
  v /= v / 100.0;
  s /= s / 100.0;
  return 255.0 * (v - v * s * mm);
}

static std::array<int, 3> _hsi2rgb_(double h, double s, double i) {
  if(h < 0.0 || h > 360.0) {
    Rcpp::stop("Invalid value of `h`.");
  }
  if(s < 0.0 || s > 100.0) {
    Rcpp::stop("Invalid value of `s`.");
  }
  if(i < 0.0 || i > 100.0) {
    Rcpp::stop("Invalid value of `i`.");
  }
  std::array<int, 3> out;
  double r, g, b;
  i = i / 100;
  double is = i*s / 100;
  double second = i - is;
  h = ratio_pi_180 * h;
  const double PIover3 = M_PI / 3;
  if(h < 2*PIover3) {
    r = i + is * cos(h) / cos(PIover3 - h);
    b = second;
    g = i + 2*is + b - r; 
  } else if(h < 4*PIover3) {
    g = i + is * cos(h - 2*PIover3) / cos(M_PI + h);
    r = second;
    b = i + 2*is + r - g;
  } else {
    b = i + is * cos(h - 4*PIover3) / cos(5*PIover3 - h);
    g = second;
    r = i + 2*is + g - b;
  }
  out[0] = (int)round(85.0 * r);
  out[1] = (int)round(85.0 * g);
  out[2] = (int)round(85.0 * b);
  return out;
}


namespace RcppColors {

  static inline std::string rgb2hex(double r, double g, double b) {
    int ri = (int)round(r);
    int gi = (int)round(g);
    int bi = (int)round(b);
    if(ri < 0 || ri > 255) {
      Rcpp::stop("Invalid value of `r`.");
    }
    if(gi < 0 || gi > 255) {
      Rcpp::stop("Invalid value of `g`.");
    }
    if(bi < 0 || bi > 255) {
      Rcpp::stop("Invalid value of `b`.");
    }
    return rgb_to_hex(ri, gi, bi);
  }

  static inline std::string rgba2hex(double r, double g, double b, double a) {
    const std::string ahex = opacity(a);
    return rgb2hex(r, g, b) + ahex;
  }

  static inline std::string hsluv2hex(double h, double s, double l) {
    return _hsluv2hex_(h, s, l);
  }

  static inline std::string hsluv2hex(double h,
                                      double s,
                                      double l,
                                      double alpha) {
    const std::string alphahex = opacity(alpha);
    return _hsluv2hex_(h, s, l) + alphahex;
  }

  static inline std::string hsv2hex(double h, double s, double v) {
    if(h < 0.0 || h > 360.0) {
      Rcpp::stop("Invalid value of `h`.");
    }
    if(s < 0.0 || s > 100.0) {
      Rcpp::stop("Invalid value of `s`.");
    }
    if(v < 0.0 || v > 100.0) {
      Rcpp::stop("Invalid value of `v`.");
    }
    int r = (int)(f(5.0, h, s, v));
    int g = (int)(f(3.0, h, s, v));
    int b = (int)(f(1.0, h, s, v));
    return rgb_to_hex(r, g, b);
  }

  static inline std::string hsv2hex(double h,
                                    double s,
                                    double v,
                                    double alpha) {
    const std::string alphahex = opacity(alpha);
    return hsv2hex(h, s, v) + alphahex;
  }

  static inline std::string hsi2hex(double h, double s, double i) {
    std::array<int, 3> rgb = _hsi2rgb_(h, s, i);
    return rgb_to_hex(rgb[0], rgb[1], rgb[2]); 
  }

}  // namespace RcppColors

#endif