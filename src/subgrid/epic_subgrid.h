/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2024-2025 Ramanakumar Sankar                      *
 * Copyright (C) 2013-2023 Timothy Dowling                         *
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

#include "epic_datatypes.h"

#ifndef EPIC_SUBGRID_H
#define EPIC_SUBGRID_H

/* * * * * * * * * * * * * epic_subgrid.h  * * * * * * * * * * * * *
 *                                                                 *
 *       TE Dowling, VK Parimi, RP LeBeau,                         *
 *                                                                 *
 *       Header file for epic_subgrid.c                            *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Defines.
 */
#define NU_TURB_EPSILON (1.e-3)
#define NU_TURB_MIN (1.e-3 * planet->kinvisc)

/*
 * Shift macros.
 */
#define D_WALL(k, j, i) d_wall[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define DIFF_COEF(k, j, i) diff_coef[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define TAU11(j, i) tau11[i + (j) * Iadim - Shift2d]
#define TAU22(j, i) tau22[i + (j) * Iadim - Shift2d]
#define TAU33(j, i) tau33[i + (j) * Iadim - Shift2d]
#define TAU12(j, i) tau12[i + (j) * Iadim - Shift2d]
#define DIE(j, i) die[i + (j) * Iadim - Shift2d]
#define KIE(j, i) kie[i + (j) * Iadim - Shift2d]
#define TAU_WALL(j, i) tau_wall[i + (j) * Iadim - Shift2d]

#define LPHH(k, j, i) lphh[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define LPUU(k, j, i) lpuu[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define LPVV(k, j, i) lpvv[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

#define BUFF1(j, i) buff1[i + (j) * Iadim - Shift2d]
#define BUFF2(j, i) buff2[i + (j) * Iadim - Shift2d]

#define COEF(itmp, j, i) coef[itmp][i + (j) * Iadim - Shift2d]

/*
 * Function prototypes.
 */
double max_nu_nondim(int order);

void set_max_nu(double *max_nu_horizontal);

void set_hyperviscosity(boolean modify);

void scalar_horizontal_subgrid(double **Buff2D);

void scalar_horizontal_diffusion(double **Buff2D);

void scalar_hyperviscosity(int nu_order, double nu_hyper, double **Buff2D, int kstart, int kend,
                           double *h);

void adiabatic_adjustment(void);

void laplacian_h(int kstart, int kend, double *hh, double *diff_coeff, double *lph);

void scalar_vertical_subgrid(double **Buff2D);

void scalar_vertical_diffusion(double **Buff2D);

void uv_horizontal_subgrid(double **Buff2D);

void uv_horizontal_diffusion(double **Buff2D);

void divergence_damping(double nudiv_nondim, double **Buff2D);

void uv_hyperviscosity(int nu_order, double nu_hyper, double **Buff2D);

void laplacian_uv(int kstart, int kend, double *uu, double *vv, double viscosity, double *lpuu,
                  double *lpvv);

void uv_vertical_subgrid(double **Buff2D);

void uv_vertical_diffusion(double **Buff2D);

void make_arrays_subgrid(void);

void free_arrays_subgrid(void);

void init_subgrid(void);

void set_diffusion_coef(void);

void source_sink_turb(double **Buff2D);

void source_sink_SA(double **Buff2D);

void dwall_SA(double *d_wall);

void fp_init_prof(void);

double delta_SA(int K, int J, int I);

void tau_surface(int index, double *tau_wall, double *buffji);

double law_of_the_wall(int K, int J, int I, int index, double u_tan);

double func_utau(double u_tau, double u_tan, double dwall);

double invert_fv1(double t_vis);

/* * * * * * * * * * * * * * * end of epic_subgrid.h * * * * * * * */
#endif
