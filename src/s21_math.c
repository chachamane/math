// created by pizpotli, tokenpah, lauriege, muquele
#include "s21_math.h"

long double s21_sin(double x) {
  while (x > s21_PI || x < -s21_PI) {
    if (x > 0) {
      x = x - 2 * s21_PI;
    } else {
      x = x + 2 * s21_PI;
    }
  }
  long double y = x, rezult = x, fac = 1.;
  while (s21_fabs(y) > s21_EPS) {
    y = -1 * y * x * x / (2 * fac * (2 * fac + 1));
    fac += 1.;
    rezult += y;
  }
  return rezult;
}

long double s21_cos(double x) {
  long double rezult;
  while (x > s21_PI || x < -s21_PI) {
    if (x > 0) {
      x = x - 2 * s21_PI;
    } else {
      x = x + 2 * s21_PI;
    }
  }
  rezult = s21_sin(s21_PI / 2 - x);
  return rezult;
}

long double s21_tan(double x) {
  long double rezult;
  if (s21_cos(x) != 0) {
  rezult = s21_sin(x) / s21_cos(x);
  } else {
    rezult = s21_NAN;
  }
  return rezult;
}

long double s21_asin(double x) {
  long double rezult;
  if (x > 1 || x < -1) {
    rezult = s21_NAN;
  } else {
  rezult = s21_atan(x / (s21_sqrt(1 - x * x)));
  }
  return rezult;
}

long double s21_acos(double x) {
  long double rezult;
  if (x > 1 || x < -1) {
    rezult = s21_NAN;
  } else {
    if (x > 0 && x < 1) {
      rezult = s21_atan(s21_sqrt(1 - x * x) / x);
    }
    if (x < 0 && x > -1) {
      rezult = s21_PI + s21_atan(s21_sqrt(1 - x * x) / x);
    }
  }
  if (x == -1) {
    rezult = 3;
  }
  if (x == 0) {
    rezult = 1;
  }
  if (x == 1) {
    rezult = 0;
  }
  return rezult;
}

long double s21_atan(double x) {
  double a, n, marctg, marctg_sravnenie;
  int b, znak, result_znak;
  result_znak = b = znak = n = 1;
  marctg = 0;
  a = x;
  if (x == 1) {
    marctg = 0.785398;
  } else if (x == -1) {
    marctg = -0.785398;
  } else if (x > -1 && x < 1) {
    while (n > s21_EPS) {
      marctg_sravnenie = marctg;
      marctg = marctg + znak * (a / b);
      a = a * x * x;
      b = b + 2;
      znak = -znak;
      n = abs_sravnenie(marctg, marctg_sravnenie);
    }
  } else
  {
    if (x < 0) {
      a = -x;
      result_znak = -result_znak;
    }
    marctg = s21_PI / 2;
    while (n > s21_EPS) {
      marctg_sravnenie = marctg;
      marctg = marctg - znak * (1 / (a * b));
      znak = -znak;
      a = a * x * x;
      b = b + 2;
      n = abs_sravnenie(marctg, marctg_sravnenie);
    }
  }
  return marctg * result_znak;
}

double abs_sravnenie(double a1, double a2) {
  double n_sravnenie;
  if (a1 < 0) {
    a1 = -a1;
  }

  if (a2 < 0) {
    a2 = -a2;
  }

  n_sravnenie = a1 - a2;

  if (n_sravnenie < 0) {
    n_sravnenie = -n_sravnenie;
  }
  return n_sravnenie;
}

long double s21_pow(double base, double exp) {
  long double rezult;
  if(base != s21_INF && exp == s21_NAN) {
  rezult = s21_NAN;
} else {
    if (base == 0 && exp < 0){
        rezult = s21_INF;
    } else {
        if (exp == s21_INF && (base < -1 || base > 1)) {
            rezult = s21_INF;
        } else {
            if (exp == s21_INF && (base > -1 && base < 1)) {
                rezult = 0;
            } else {
                if (exp == s21_INF && (base == -1 || base == 1)) {
                    rezult = 1;
                } else {
                    if (base == s21_INF) {
                        if (exp < 0){
                            rezult = 0;
                        }
                        if (exp == 0){
                            rezult = 1;
                        }
                        if (exp > 0){
                            rezult = s21_INF;
                        }
                    } else {
  int i = exp;
  if (base < 0 && (exp - i) != 0) {
    rezult = s21_NAN;
  } else {
    if (exp != 0) {
      rezult = s21_exp(exp * s21_log(base));
    } else {
      rezult = 1;}
    }
  }}}}}
    }
  return rezult;
}

long double s21_sqrt(double x) {
  long double result = 4, y = 0;
  while (s21_fabs(result - y) > s21_EPS) {
    if (x < 0) {
      result = s21_NAN;
      break;
    }
    y = result;
    result = (y + x / y) / 2;
  }
  return result;
}

long double s21_fabs(long double x) {
  if (x < 0) {
    x *= -1.0;
  }
  return x;
}

int s21_abs(int x) {
  long int rezult = x;
  if (rezult < 0) {
    rezult *= -1;
  }
  return rezult;
}

long double s21_floor(double x) {
  int i;
  double y;
  if (x < 0) {
    y = x * -1;
  } else {
    y = x;
  }
  for (i = 0; i < y; i++) {
  }
  if (x < 0) {
    i *= -1;
  }
  if (x >= 0) {
    i = x;
  }
  return i;
}

long double s21_fmod(double x, double y) {
  int znak = 0;
  if (x < 0) {
    znak++;
    x *= -1;
  }
  if (y < 0) {
    y *= -1;
  }
  if (y != 0) {
    while (x >= y) {
      x = x - y;
    }
    while (znak != 0) {
      x *= -1;
      znak--;
    }
  } else {
    x = s21_NAN;
  }
  return x;
}

long double s21_exp(double x) {
  long double exp_f = 1, n = 1, p = 1;
  if (x > -24) {
    while (s21_fabs(p) > s21_EPS) {
      p *= x / n;
      exp_f += p;
      n++;
    }
  } else {
    exp_f = 0;
  }
  return exp_f;
}

long double s21_log(double x) {
  long double ans = 0.0;

  if (x < -s21_EPS) {
    ans = 0 / (x - x);
  } else if (x < s21_EPS) {
    ans = -1 / (x - x);
  } else if (x < 0.1 || (x > 1.9 && x < 2 - s21_EPS)) {
    ans = 2 * s21_log(s21_sqrt(x));
  } else if (x > s21_EPS && x < 2 - s21_EPS) {
    long double mul = x - 1.0;
    long double cntr = 1.0;

    long double cur = mul;

    while (cur > s21_EPS || cur < -s21_EPS) {
      ans += cur / cntr;

      cntr += 1.0;
      cur = (-cur) * mul;
    }

  } else {
    ans = s21_log((2 * x) / (x + 1)) - s21_log(2 / (x + 1));
  }

  return ans;
}

long double s21_ceil(double x) {
  int i = 0;
  long double rezult;
  i = x;
  if (i < x) {
    rezult = i + 1;
  }
  if (i >= x) {
    rezult = i;
  }
  return rezult;
}
