/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * All material in this file is covered by the                     *
 * GNU General Public License.                                     *
 *                                                                 *
 * Copyright (C) 2024-2025 Ramanakumar Sankar, except as noted.    *
 * Copyright (C) 1998-2023 Timothy Dowling, except as noted.       *
 *                                                                 *
 * The fcmp() function is derived from the fcmp() function that is *
 * Copyright (c) 1998-2000 Theodore C. Belding                     *
 * University of Michigan Center for the Study of Complex Systems. *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License     *
 * as published by the Free Software Foundation; either version 2  *
 * of the License, or (at your option) any later version.          *
 * A copy of this License is in the file:                          *
 *   $EPIC_PATH/License.txt                                        *
 *                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
 *                                                                 *
 * You should have received a copy of the GNU General Public       *
 * License along with this program; if not, write to the Free      *
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     *
 * Boston, MA 02110-1301, USA.                                     *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef EPIC_FUNCS_UTIL_H
#define EPIC_FUNCS_UTIL_H

/* * * * * * * * * epic_funcs_util.h * * * * * * * * * * * * 
 *                                                         *
 * Header file for utility functions that are independent  *
 * of the EPIC Model; the source code is in                *
 *   $EPIC_PATH/src/shared/epic_funcs_util.c.              *
 *                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h> 
#include <stdlib.h>   
#include <errno.h>    
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#include "epic_datatypes.h"

/*
 * Parameters.
 */
#define N_STR 256

/*
 * Logical:
 */
#define TRUE    1
#define FALSE   0
#define SUCCESS TRUE
#define FAILURE FALSE

#define SETUP_UTIL   0
#define RUN_UTIL     1
#define CLEANUP_UTIL 2

/*
 * Mathematical:
 */
#if !defined(M_E)
#  define M_E         2.7182818284590452354
#  define M_LOG2E     1.4426950408889634074
#  define M_LOG10E    0.43429448190325182765
#  define M_LN2       0.69314718055994530942
#  define M_LN10      2.30258509299404568402
#  define M_PI        3.14159265358979323846
#  define M_PI_2      1.57079632679489661923
#  define M_PI_4      0.78539816339744830962
#  define M_1_PI      0.31830988618379067154
#  define M_2_PI      0.63661977236758134308
#  define M_2_SQRTPI  1.12837916709551257390
#  define M_SQRT2     1.41421356237309504880
#  define M_SQRT1_2   0.70710678118654752440
#endif

#define DEG (M_PI/180.)

#undef MIN
#define MIN(x,y) ({ \
         const double _x = (double)(x); \
         const double _y = (double)(y); \
         _x < _y ? _x : _y; })

#undef MAX
#define MAX(x,y) ({ \
         const double _x = (double)(x); \
         const double _y = (double)(y); \
         _x > _y ? _x : _y; })

#undef LIMIT_RANGE
#define LIMIT_RANGE(min,x,max) ({ \
         const double _min = (double)(min); \
         const double _x   = (double)(x);   \
         const double _max = (double)(max); \
         _x < _min ? _min : ( _x > _max ? _max : _x ); })

#undef IMIN
#define IMIN(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i < _j ? _i : _j; })

#undef IMAX
#define IMAX(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i > _j ? _i : _j; })

#undef NINT
#define NINT(x) ({ \
         const double _x = (double)(x); \
         _x > 0. ? (int)(_x+.5) : (int)(_x-.5); })

#undef SIGN
#define SIGN(x) ({ \
         const double _x = (double)(x); \
         _x == 0. ? 0. : (_x > 0. ? 1. : -1.); })

#undef NR_SIGN
#define NR_SIGN(a,b) ((b) > 0. ? fabs(a) : -fabs(a))

#undef SQR
#define SQR(x) ({ \
          const double _x = (double)(x); \
          _x*_x; })

#undef GET_DIGIT
#define GET_DIGIT(ival,div) (((ival)-(ival)%(div))/(div))

/* 
 * DS stands for Dennis and Schnabel (1996).
 */
#define DS_EPS      machine_epsilon()
#define DS_SQRT_EPS sqrt(DS_EPS)
#define DS_TOL_X    pow(DS_EPS,2./3.)
#define DS_TOL_F    pow(DS_EPS,1./3.)
#define DS_TOL_MIN  pow(DS_EPS,2./3.)
#define DS_ALPHA    1.e-4
#define DS_MAX_STEP 100.

#define DS_LINE_STEP   0
#define DS_HOOK_STEP   1    /* Not yet implemented. */
#define DS_DOGLEG_STEP 2
/*
 * Values for status (retcode in DS96):
 */
#define DS_X_ACCEPTED        0
#define DS_X_NO_PROGRESS     1
#define DS_REDUCE_DELTA      2
#define DS_INCREASE_DELTA    3
#define DS_MAX_IT_EXCEEDED   4
#define DS_MAX_TAKEN_5       5
#define DS_INITIAL           6
#define DS_SINGULAR_JACOBIAN 7

/*
 * Shift macros:
 */
#define QT(i,j)       qt[j+(i)*n]
#define R(i,j)         r[j+(i)*n]
#define A_EIG(r,c) a_eig[(c)+klen*(r)]
#define AA_EIG(r,c) aa_eig[(c)+klen*(r)]

/* * * * * * * * * * * * *
 *                       *
 * Function prototypes.  *
 *                       *
 * * * * * * * * * * * * */

/* 
 * Memory allocation. 
 */
int *ivector(int   nl,
             int   nh,
             char *calling_func);

void free_ivector(int  *m,
                  int   nl,
                  int   nh,
                  char *calling_func);

double *fvector(int   nl,
                int   nh,
                char *calling_func);

void free_fvector(double *m,
                  int     nl,
                  int     nh,
                  char   *calling_func);

double *dvector(int   nl,
                int   nh,
                char *calling_func);

void free_dvector(double *m,
                  int     nl,
                  int     nh,
                  char   *calling_func);

double_triplet *dtriplet(int   nl,
                         int   nh,
                         char *calling_func);

void free_dtriplet(double_triplet *m,
                   int             nl,
                   int             nh,
                   char           *calling_func);

/*
 *  The following functions are adapted from 
 *  Numerical Recipes in C:
 */
void spline(int             n,
            double_triplet *table,
            double          yp0, 
            double          ypn);

double splint(register double          xx, 
              register double_triplet *table,
              register double          dx);

double linint(register double          xx,
              register double_triplet *table,
              register double          dx);

void spline_pchip(int             n,
                  double_triplet *table);

double splint_pchip(double          xx,
                    double_triplet *table,
                    double          h);

void periodic_spline_pchip(int             n,
                           double_triplet *table);

double pchst(double arg1,
             double arg2);

double lagrange_interp(register double *f,
                       register double *x,
                       register int     order);

double inc_gamma_nr(double aa,
                    double xx);

double gamma_nr(double xx);

double sech2(double);

void exp_integral_setup(double_triplet *exp3table,
                        double_triplet *exp4table,
                        int             numdatapoints);

double normed_legendre(int    l,
                       int    m,
                       double x);

double machine_epsilon(void);

int find_root(double  x1,
              double  x2,
              double  xacc,
              double *x_root,
              double  (*func)(double));

int broyden_root(int     n,
                 double *x,
                 void  (*vecfunc)(int,double *,double *),
                 double tol_f,
                 int    max_it);

int global_step(int     n,
                double *x_old,
                double  f_old,
                double *g,
                double *r,
                double *sn,
                double  max_step,
                double *delta,
                int     step_type,
                int    *status,
                double *x,
                double *f,
                double *fvec,
                void   (*vecfunc)(int,double *,double *));

int line_search(int     n,
                double *x_old,
                double  f_old,
                double *g,
                double *sn,
                double  max_step,
                int    *status,
                double *x,
                double *f,
                double *fvec,
                void  (*vecfunc)(int,double *,double *));

int dogleg_driver(int     n,
                  double *x_old,
                  double  f_old,
                  double *g,
                  double *r,
                  double *sn,
                  double  max_step,
                  double *delta,
                  int    *status,
                  double *x,
                  double *f,
                  double *fvec,
                  void  (*vecfunc)(int,double *,double *));

int dogleg_step(int     n,
                double *g,
                double *r,
                double *sn,
                double  newt_length,
                double  max_step,
                double *delta,
                int    *first_dog,
                double *s_hat,
                double *nu_hat,
                double *s);

int trust_region(int     n,
                 double *x_old,
                 double  f_old,
                 double *g,
                 double *s,
                 int     newt_taken,
                 double  max_step,
                 int     step_type,
                 double *r,
                 double *delta,
                 int    *status,
                 double *x_prev,
                 double *f_prev,
                 double *x,
                 double *f,
                 double *fvec,
                 void  (*vecfunc)(int,double *,double *));

int qr_decompose(int     n,
                 double *r,
                 double *c,
                 double *d);

void qr_update(int     n,
               double *r,
               double *qt,
               double *u,
               double *v);

void qr_rotate(int     n,
               double *r,
               double *qt,
               int     i,
               double  a,
               double  b);

void lu_decompose(int    n,
                 double *a,
                 int    *index,
                 double *d);

void lu_backsub(int     n,
                double *a,
                int    *index,
                double *b);

void lu_improve(int     n,
                double *a,
                double *alu,
                int    *index,
                double *b,
                double *x);

int find_place_in_table(int             n,
                        double_triplet *table,
                        double          x,
                        double         *dx);

int hunt_place_in_table(int             n,
                        double_triplet *table,
                        double          x,
                        double         *dx,
                        int             il);

/*
 * Choices for pivot_type in tridiag():
 */
#define WITHOUT_PIVOTING 0
#define WITH_PIVOTING    1

void tridiag(int     n,
             double *a,
             double *b,
             double *c,
             double *r,
             double *u,
             int    pivot_type);

void band_decomp(int     n,
                 int     m1,
                 int     m2,
                 double *a,
                 double *al,
                 int    *index,
                 double *d);

void band_back_sub(int     n,
                   int     m1,
                   int     m2,
                   double *a,
                   double *al,
                   int    *index,
                   double *b);

void band_multiply(int     n,
                   int     m1,
                   int     m2,
                   double *a,
                   double *x,
                   double *b);

void band_improve(int     n,
                  int     m1,
                  int     m2,
                  double *aorig,
                  double *a,
                  double *al,
                  int    *index,
                  double *b,
                  double *x);

double poly_interp(int    n,
                  double *xa,
                  double *ya,
                  double  x,
                  double *dy);

void poly_coeff(int     n,
                double *x,
                double *y,
                double *coeff);

double nth_trapezoidal(int   n,
                      double (*func)(double),
                      double a,
                      double b);

double romberg_integral(double (*func)(double),
                        double a,
                        double b,
                        double tol);

void crank_nicolson(int     n,
                    double  dt,
                    double *z,
                    double *A,
                    double *mu,
                    double *rho,
                    double *ANS);

void compact_differentiation(int     action,
                             int     n,
                             double *z,
                             double *f,
                             double *dfdz);

void compact_integration(int     n,
                         double *z,
                         double *dfdz,
                         double *f);

void hqr(double *a,
         int     n,
         double *wr,
         double *wi);

void quicksort(double *mag, 
               int     left, 
               int     right);

void swap(double *mag,
          int     i,
          int     j);

void four1(double *,
           unsigned long, 
           int              );

void realft(double *,
            unsigned long,
            int             );

/*
 * NOTE: For LINUX, with -D_BSD_SOURCE cabs() is defined, such that
 * a type-mismatch error occurs if we name the function c_abs() as cabs().
 */
complex c_num(double x, double y);
complex c_mult(complex z1,complex z2);
complex c_add(complex z1,complex z2);
complex c_sub(complex z1,complex z2);
complex c_exp(complex z);
double  c_abs(complex z);
double  c_real(complex z);
double  c_imag(complex z);

int fcmp(double x1,
         double x2);

void least_squares(double *x,
                   double *y,
                   int     n,
                   double *a);

void savitzky_golay(double *c,
                    int     np,
                    int     nl,
                    int     nr,
                    int     ld,
                    int     m);

double random_number(long *idum);

double lat_centric_to_graphic(double lat,
                              double rerp);

double lat_graphic_to_centric(double lat,
                              double rerp);

double surface_area_oblate(double a,
                           double c);

void util_error(char *calling_function,
                char *Message);

/* * * * * * * * * *  end of epic.h  * * * * * * * * * * * * * * * * * * * * */ 
#endif

