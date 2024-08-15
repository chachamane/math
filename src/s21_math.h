// created by pizpotli, tokenpah, lauriege, muquele
#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#include <stdio.h>

#define s21_PI 3.141592653589793238462643
#define s21_EPS 1e-20
#define s21_EXP 2.7182818284590452353602874713526624
#define s21_MAX_DOUBLE 1.7976931348623158e+308
#define s21_MAX_INT 2147483647
#define s21_INF 1.0 / 0.0
#define s21_NAN 0.0 / 0.0

long double s21_tan(double x);
long double s21_cos(double x);
long double s21_sin(double x);
long double s21_ceil(double x);
long double s21_floor(double x);
long double s21_atan(double x);
long double s21_acos(double x);
long double s21_asin(double x);
long double s21_fmod(double x, double y);
long double s21_exp(double x);
long double s21_log(double x);
long double s21_sqrt(double x);
long double s21_pow(double base, double exp);
long double s21_fabs(long double x);
int s21_abs(int x);
double abs_sravnenie(double a1, double a2);

#endif  // SRC_S21_MATH_H_
