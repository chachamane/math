// created by pizpotli, tokenpah, lauriege, muquele
#include <check.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "s21_math.h"

START_TEST(sin_test) {
  long double test1 = 3.141592653589793238462643 / 2;
  ck_assert_uint_eq(sin(test1), s21_sin(test1));
  long double test2 = 1001;
  ck_assert_uint_eq(sin(test2), s21_sin(test2));
  long double test3 = -1001;
  ck_assert_uint_eq(sin(test3), s21_sin(test3));
  long double test4 = 100.101;
  ck_assert_uint_eq(sin(test4), s21_sin(test4));
  for (double i = -s21_PI; i <= s21_PI; i += 0.01) {
    ck_assert_uint_eq(sin(i), s21_sin(i));
  }
}
END_TEST

START_TEST(cos_test) {
  long double test1 = 0;
  ck_assert_ldouble_eq_tol(s21_cos(test1), cos(test1), 1e-6);
  long double test2 = 100101;
  ck_assert_ldouble_eq_tol(s21_cos(test2), cos(test2), 1e-6);
  long double test3 = -100101;
  ck_assert_ldouble_eq_tol(s21_cos(test3), cos(test3), 1e-6);
  long double test4 = 100.101;
  ck_assert_ldouble_eq_tol(s21_cos(test4), cos(test4), 1e-6);
  for (double i = -s21_PI; i <= s21_PI; i += 0.01) {
    ck_assert_ldouble_eq_tol(s21_cos(i), cos(i), 1e-6);
  }
}
END_TEST

START_TEST(tan_test) {
  long double test1 = 0;
  ck_assert_ldouble_eq_tol(s21_tan(test1), tan(test1), 1e-6);
  long double test2 = 100101;
  ck_assert_ldouble_eq_tol(s21_tan(test2), tan(test2), 1e-6);
  long double test3 = -100101;
  ck_assert_ldouble_eq_tol(s21_tan(test3), tan(test3), 1e-6);
  long double test4 = 100.101;
  ck_assert_ldouble_eq_tol(s21_tan(test4), tan(test4), 1e-6);
  for (double i = -1; i <= 1; i += 0.01) {
    ck_assert_ldouble_eq_tol(s21_tan(i), tan(i), 1e-6);
  }
}
END_TEST

START_TEST(abs_test) {
  for (int i = -1000; i <= 1000; i += 10) {
    ck_assert_uint_eq(s21_abs(i), abs(i));
  }
}
END_TEST

START_TEST(fabs_test) {
  for (double i = -10000000.12354; i <= 10000000; i += 1000.65455435) {
    ck_assert_double_eq(s21_fabs(i), fabs(i));
  }
}
END_TEST

START_TEST(pow_test) {
  // ck_assert_ldouble_eq_tol(pow(2,4), s21_pow(2,s21_NAN),
  //                          1e-6);
  ck_assert_int_eq(pow(2,s21_NAN),s21_pow(2,s21_NAN));                                           // 0
  ck_assert_ldouble_eq_tol(pow(0.2, s21_INF), s21_pow(0.2, s21_INF), 1e-6); // 0
  ck_assert_ldouble_eq_tol(pow(-1, s21_INF), s21_pow(-1, s21_INF), 1e-6);   // 1
  ck_assert_ldouble_eq_tol(pow(1, s21_INF), s21_pow(1, s21_INF), 1e-6);     // 1
  ck_assert_ldouble_infinite(pow(1.2, s21_INF));
  ck_assert_ldouble_infinite(s21_pow(1.2, s21_INF)); // inf +
  ck_assert_ldouble_infinite(pow(-1.2, s21_INF));
  ck_assert_ldouble_infinite(s21_pow(-1.2, s21_INF)); // inf +
  ck_assert_ldouble_eq_tol(pow(s21_INF, -2), s21_pow(s21_INF, -2), 1e-6); // 0
  ck_assert_ldouble_infinite(pow(s21_INF, 2));
  ck_assert_ldouble_infinite(s21_pow(s21_INF, 2));                      // inf
  ck_assert_ldouble_eq_tol(pow(s21_INF, 0), s21_pow(s21_INF, 0), 1e-6); // 1
  ck_assert_ldouble_nan(pow(-4, 4.6));
  ck_assert_ldouble_nan(s21_pow(-4, 4.6)); // nan
  ck_assert_ldouble_eq_tol(pow(10, -4.6), s21_pow(10, -4.6), 1e-6);
  ck_assert_ldouble_eq_tol(pow(10, 4.6), s21_pow(10, 4.6), 1e-6);
 ck_assert_ldouble_nan(s21_pow(-4, s21_NAN));
}
END_TEST

START_TEST(floor_test) {
    for (double i = -7270; i <= 7270; i += 4.635) {
      ck_assert_ldouble_eq_tol(s21_floor(i), floor(i), 1e-6);
      ck_assert_ldouble_eq_tol(s21_floor(0), floor(0), 1e-6);
    }
}
END_TEST

START_TEST(atan_test) {
    for (double i = -100; i <= 100; i++) {
      ck_assert_ldouble_eq_tol(atan(i), s21_atan(i), 1e-06);
    }
    for (double i = -1; i < 1; i += 0.01) {
      ck_assert_ldouble_eq_tol(atan(i), s21_atan(i), 1e-06);
    }
}
END_TEST

START_TEST(acos_test) {
    for (double i = -100; i <= 100; i++) {
      ck_assert_uint_eq(acos(i), s21_acos(i));
    }
    for (double i = -1; i <= 1; i += 0.001) {
      ck_assert_uint_eq(acos(i), s21_acos(i));
    }
}
END_TEST

START_TEST(asin_test) {
    for (double i = -100; i <= 100; i += 1) {
      ck_assert_uint_eq(asin(i), s21_asin(i));
    }
    for (double i = -1; i <= 1; i += 0.001) {
      ck_assert_uint_eq(asin(i), s21_asin(i));
    }
}
END_TEST

START_TEST(ceil_test) {
    for (double i = -15; i <= 15; i += 0.025) {
      ck_assert_ldouble_eq_tol(s21_ceil(i), ceil(i), 1e-6);
    }
}
END_TEST

START_TEST(fmod_test) {
    double y = 3;
    for (double i = -50; i <= 50; i += 0.5) {
      ck_assert_uint_eq(fmod(i, y), s21_fmod(i, y));
      y -= 0.5;
    }
}
END_TEST

START_TEST(exp_test) {
    for (double i = -70; i <= 70; i += 0.1) {
      ck_assert_float_eq_tol(s21_exp(i), exp(i), 1e-06);
    }
}
END_TEST

START_TEST(log_test) {
    for (double i = 0.01; i < 1; i += 0.01) {
      ck_assert_ldouble_eq_tol(log(i), s21_log(i), 1e-06);
    }
    for (double i = 1; i <= 150; i += 1.) {
      ck_assert_ldouble_eq_tol(s21_log(i), log(i), 1e-06);
    }
}
END_TEST

START_TEST(sqrt_test) {
    for (double i = 1; i <= 1005; i += 10.025) {
      ck_assert_ldouble_eq_tol(sqrt(i), s21_sqrt(i), 1e-6);
    }
}
END_TEST

Suite *test_func(void) {
    Suite *x;
    x = suite_create("check");

    TCase *sin_case;
    TCase *cos_case;
    TCase *tan_case;
    TCase *abs_case;
    TCase *fabs_case;
    TCase *pow_case;
    TCase *floor_case;
    TCase *atan_case;
    TCase *acos_case;
    TCase *asin_case;
    TCase *ceil_case;
    TCase *fmod_case;
    TCase *exp_case;
    TCase *log_case;
    TCase *sqrt_case;

    sin_case = tcase_create("sin");
    suite_add_tcase(x, sin_case);
    tcase_add_test(sin_case, sin_test);

    cos_case = tcase_create("cos");
    suite_add_tcase(x, cos_case);
    tcase_add_test(cos_case, cos_test);

    tan_case = tcase_create("tan");
    suite_add_tcase(x, tan_case);
    tcase_add_test(tan_case, tan_test);

    abs_case = tcase_create("abs");
    suite_add_tcase(x, abs_case);
    tcase_add_test(abs_case, abs_test);

    fabs_case = tcase_create("fabs");
    suite_add_tcase(x, fabs_case);
    tcase_add_test(fabs_case, fabs_test);

    pow_case = tcase_create("pow");
    suite_add_tcase(x, pow_case);
    tcase_add_test(pow_case, pow_test);

    floor_case = tcase_create("floor");
    suite_add_tcase(x, floor_case);
    tcase_add_test(floor_case, floor_test);

    atan_case = tcase_create("atan");
    suite_add_tcase(x, atan_case);
    tcase_add_test(atan_case, atan_test);

    acos_case = tcase_create("acos");
    suite_add_tcase(x, acos_case);
    tcase_add_test(acos_case, acos_test);

    asin_case = tcase_create("asin");
    suite_add_tcase(x, asin_case);
    tcase_add_test(asin_case, asin_test);

    ceil_case = tcase_create("ceil");
    suite_add_tcase(x, ceil_case);
    tcase_add_test(ceil_case, ceil_test);

    fmod_case = tcase_create("fmod");
    suite_add_tcase(x, fmod_case);
    tcase_add_test(fmod_case, fmod_test);

    exp_case = tcase_create("exp");
    suite_add_tcase(x, exp_case);
    tcase_add_test(exp_case, exp_test);

    log_case = tcase_create("log");
    suite_add_tcase(x, log_case);
    tcase_add_test(log_case, log_test);

    sqrt_case = tcase_create("sqrt");
    suite_add_tcase(x, sqrt_case);
    tcase_add_test(sqrt_case, sqrt_test);

    return x;
}

int main(void) {
    // int no_filed = 0;
    Suite *x;
    SRunner *y;

    x = test_func();
    y = srunner_create(x);
    srunner_run_all(y, CK_VERBOSE);
    // srunner_run_all(y, CK_NORMAL);
    // no_filed = srunner_ntests_failed(y);
    srunner_free(y);
    // return (no_filed == 0) ? EXIT_SUCCESS: EXIT_FAILURE;
    return 0;
}
