/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2024-2025 Ramanakumar Sankar                      *
 * Copyright (C) 1998-2023 Timothy E. Dowling                      *
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

#ifndef EPIC_H
#define EPIC_H
/* * * * * * * * * * * * * epic.h  * * * * * * * * * * * * * * * * * * * * * *
 *                                                                           *
 *       Timothy E. Dowling                                                  *
 *                                                                           *
 *       Header file for the EPIC atmospheric model                          *
 *                                                                           *
 *       This file provides a glossary of terms for the compiler.            *
 *                                                                           *
 *       Global declarations are made in epic_globals.c for dynamical-core   *
 *       values, and in the physics packages as needed, and are referenced   *
 *       here with the "extern" modifier.                                    *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#define _XOPEN_SOURCE_EXTENDED 1
#ifdef sun4
/* To get definitions of u_short, etc. */
#define __EXTENSIONS__
#endif
#endif

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <time.h>
#include <unistd.h>

#include <netcdf.h>

#include "epic_datatypes.h"
#include "epic_funcs_util.h"
#include "epic_io_macros.h"
#include "epic_microphysics.h"
#include "epic_subgrid.h"

/*
 * The following macro is useful to insert into code while trying to corral a problem.
 */
#define DEBUG_MILESTONE(comment)                                                                   \
    fprintf(stderr, "node %2d, %s(), %2d: " #comment "\n", IAMNODE, dbmsname, ++idbms);            \
    fflush(stderr);

/*
 * There are two implemented ways to average potential vorticity onto the u and v grids.
 * Define PV_SCHEME to be SADOURNY_1975 or ARAKAWA_LAMB_1981:
 *
 * NOTE: See comments in get_kin() related to this choice and the Hollingsworth-Kallberg
 * instability.
 */
#define SADOURNY_1975 1
#define ARAKAWA_LAMB_1981 2

#define PV_SCHEME ARAKAWA_LAMB_1981

/*
 * Logical:
 */
#define TRUE 1
#define FALSE 0
#define YES TRUE
#define NO FALSE
#define SUCCESS TRUE
#define FAILURE FALSE
#define NOT_SET -100
#define OFF 0
#define ON 1
#define ACTIVE 1
#define PASSIVE 2

#define NOT_LISTED -1
#define LISTED_AND_OFF 0
#define LISTED_AND_ON 1

#define HAS_STANDARD_NAME TRUE
#define NEEDS_STANDARD_NAME FALSE

#define APPLY_TO_UV 1
#define APPLY_TO_UV_TENDENCIES 2

#define EPIC_ALLOC 1
#define EPIC_APPLY 2
#define EPIC_FREE 3

/*
 * Mathematical:
 */
#if !defined(M_E)
#define M_E 2.7182818284590452354
#define M_LOG2E 1.4426950408889634074
#define M_LOG10E 0.43429448190325182765
#define M_LN2 0.69314718055994530942
#define M_LN10 2.30258509299404568402
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#define M_PI_4 0.78539816339744830962
#define M_1_PI 0.31830988618379067154
#define M_2_PI 0.63661977236758134308
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2 1.41421356237309504880
#define M_SQRT1_2 0.70710678118654752440
#endif

#define SQRT3 1.7320508075688772
#define SQRT3_2 0.8660254037844386
#define SQRT5 2.2360679774997898
#define ONE_SQRT5 0.4472135954999579
#define SQRT15 3.8729833462074170
#define ONE_3 0.3333333333333333

#undef  DEG
#define DEG (M_PI / 180.)

#define DT ((double)grid.dt)

/*
 * NOTE: Updates to these macros should also be done to any
 *       corresponding ones in epic_funcs_util.h
 */
#undef MIN
#define MIN(x, y)                                                                                  \
    ({                                                                                             \
        const double _x = (double)(x);                                                             \
        const double _y = (double)(y);                                                             \
        _x < _y ? _x : _y;                                                                         \
    })

#undef MAX
#define MAX(x, y)                                                                                  \
    ({                                                                                             \
        const double _x = (double)(x);                                                             \
        const double _y = (double)(y);                                                             \
        _x > _y ? _x : _y;                                                                         \
    })

#undef LIMIT_RANGE
#define LIMIT_RANGE(min, x, max)                                                                   \
    ({                                                                                             \
        const double _min = (double)(min);                                                         \
        const double _x = (double)(x);                                                             \
        const double _max = (double)(max);                                                         \
        _x < _min ? _min : (_x > _max ? _max : _x);                                                \
    })

#undef ILIMIT_RANGE
#define ILIMIT_RANGE(min, i, max)                                                                  \
    ({                                                                                             \
        const int _min = (int)(min);                                                               \
        const int _i = (int)(i);                                                                   \
        const int _max = (int)(max);                                                               \
        _i < _min ? _min : (_i > _max ? _max : _i);                                                \
    })

#undef IMIN
#define IMIN(i, j)                                                                                 \
    ({                                                                                             \
        const int _i = (int)(i);                                                                   \
        const int _j = (int)(j);                                                                   \
        _i < _j ? _i : _j;                                                                         \
    })

#undef IMAX
#define IMAX(i, j)                                                                                 \
    ({                                                                                             \
        const int _i = (int)(i);                                                                   \
        const int _j = (int)(j);                                                                   \
        _i > _j ? _i : _j;                                                                         \
    })

#undef NINT
#define NINT(x)                                                                                    \
    ({                                                                                             \
        const double _x = (double)(x);                                                             \
        _x > 0. ? (int)(_x + .5) : (int)(_x - .5);                                                 \
    })

#undef SIGN
#define SIGN(x)                                                                                    \
    ({                                                                                             \
        const double _x = (double)(x);                                                             \
        _x == 0. ? 0. : (_x > 0. ? 1. : -1.);                                                      \
    })

#undef NR_SIGN
#define NR_SIGN(a, b) ((b) > 0. ? fabs(a) : -fabs(a))

#undef GET_DIGIT
#define GET_DIGIT(ival, div) (((ival) - (ival) % (div)) / (div))

#if defined(EPIC_MPI)
/* Defined in mpispec_init(): */
extern MPI_Datatype EPIC_MPI_COMPLEX;
#endif

/*
 * Physical constants.
 * Main reference: http://physics.nist.gov/cuu/Constants/index.html.
 */
#define N_A 6.02214199e+23               /* Avogadro constant                            */
#define K_B 1.3806503e-23                /* Boltzmann constant, J/K                      */
#define H_PLANCK 6.62606876e-34          /* Planck constant, J s                         */
#define M_PROTON 1.67262158e-27          /* proton mass, kg                              */
#define G_NEWTON 6.67428e-11             /* Newton's constant of gravitation, m^3/kg/s^2 */
#define R_GAS 8.314472e+3                /* molar gas constant, J/K/kmol                 */
#define YEAR 31536000                    /* secs/year, must be an int                    */
#define DAY 86400                        /* secs/day,  must be an int                    */
#define PA_FROM_TORR (1.01325e+5 / 760.) /* pressure unit conversion                     */
#define PA_FROM_CGS 0.1                  /* pressure unit conversion                     */

/*
 * Vertical coordinate types.
 *
 * NOTE: These need to be in the same order as in the string array vert_coord_type[] in
 * epic_initial.c.
 */
#define COORD_FIRST 1
#define COORD_ISENTROPIC 1
#define COORD_ISOBARIC 2
#define COORD_HYBRID 3
#define COORD_LAST 3

/*
 * Minimum Q.
 */
#define Q_MIN EPSILON

/*
 * Bookkeeping:
 */
#define NODE0 0

#define INPUT 1
#define OUTPUT 2

#define SKIP_BC_JI 1
#define APPLY_BC_JI 2

#define NUM_WORKING_BUFFERS 8 /* JI-plane buffers, Buff2D */

/* Filter is applied poleward of LAT0 [deg]: */
#define LAT0 45.

#define MASS 0
#define MOLAR 1
#define TEND 2

/*
 * Vorticity types:
 */
#define POTENTIAL 1
#define ABSOLUTE 2
#define RELATIVE 3

/*
 * High-latitude filter options:
 */
#define APPLY_HIGHLAT_FILTER TRUE
#define DONT_APPLY_HIGHLAT_FILTER FALSE

/*
 * H types:
 */
#define DRY 1
#define TOTAL 2
#define MOIST TOTAL

/*
 * Evaluation on sigmatheta surfaces or theta surfaces.
 */
#define ON_SIGMATHETA 1
#define ON_THETA 2

/*
 * get_jlohi, get_ilohi flags:
 */
#define SETUP_GET_JLOHI (-INT_MAX)
#define SETUP_GET_ILOHI (-INT_MAX)

/*
 * portion types (must be positive):
 */
#define SIZE_DATA 1
#define POST_SIZE_DATA 2
#define HEADER_DATA 3
#define EXTRACT_HEADER_DATA 4
#define POST_HEADER_DATA 5
#define VAR_DATA 6
#define EXTRACT_DATA 7
#define ALL_DATA 8

/*
 * Layer spacing choices.
 */
#define SPACING_LOGP 1
#define SPACING_P 2
#define SPACING_THETA 3
#define SPACING_FROM_FILE 4

/*
 * Misc argument types.
 */
#define PASSING_THETA 1
#define PASSING_T 2

#define DONT_UPDATE_THETA 1
#define UPDATE_THETA 2

#define DONT_APPLY_HYPERVISCOSITY 1
#define APPLY_HYPERVISCOSITY 2

#define DONT_APPLY_HIGH_LAT_FILTER 1
#define APPLY_HIGH_LAT_FILTER 2

#define USE_DEFAULTS 1
#define USE_PROMPTS 2

#define INITMODE 0
#define CHANGEMODE 1

#define HORIZONTAL_AND_VERTICAL 1
#define JUST_HORIZONTAL 2
#define JUST_VERTICAL 3

#define STEP_PROGS_AND_SYNC_DIAGS 1
#define SYNC_DIAGS_ONLY 2

#define CALC_PHI3NK 0
#define PASSING_PHI3NK 1

#define DONT_ADJUST_AMPLITUDE FALSE
#define ADJUST_AMPLITUDE TRUE

#define STOP_SHORT_BC 1
#define MIRROR_SYMMETRIC_BC 2
#define PERIODIC_BC 3
#define ZERO_CURVATURE_BC 4
#define ZERO_BC 5
#define LINEAR_BC 6
#define LEAST_SQRS_BC 7

/*
 * Time-index pointers:
 */
extern int IT_ZERO, IT_MINUS1, IT_MINUS2;

/*
 *  Potential vorticity averaging schemes:
 */
#define AL81 1
#define HA90 2

extern const chem_element ElementAG89[LAST_ATOMIC_NUMBER + 1], ElementGS01[LAST_ATOMIC_NUMBER + 1],
    *Element;

extern char Message[N_STR];
/*
 *  Structures:
 */
extern planetspec *planet, venus, earth, mars, jupiter, saturn, titan, uranus, neptune, triton,
    pluto, hot_jupiter, held_suarez, goldmine, venus_llr05;
extern gridspec grid;
extern thermospec thermo;
extern variablespec var;

#define TIME (double)difftime(var.model_time, var.start_time)

#if defined(EPIC_MPI)
#define IAMNODE (para.iamnode)
#define ILO (para.mylo[0])
#define JLO (para.mylo[1])
#define KLO (para.mylo[2])
#define IHI (para.myhi[0])
#define JHI (para.myhi[1])
#define KHI (para.myhi[2])
#define IPAD (para.npad[0])
#define JPAD (para.npad[1])
#define KPAD (para.npad[2])
#define JFIRST (para.jfirst)
#define JLAST (para.jlast)
#define IS_SPOLE (para.is_spole)
#define IS_NPOLE (para.is_npole)
#else
#define IAMNODE 0
#define ILO (grid.ilo)
#define JLO (grid.jlo)
#define KLO (1)
#define IHI (grid.ni)
#define JHI (grid.nj)
#define KHI (grid.nk)
#define IPAD (grid.pad[0])
#define JPAD (grid.pad[1])
#define KPAD (grid.pad[2])
#define JFIRST (grid.jfirst)
#define JLAST (grid.jlast)
#define IS_NPOLE (grid.is_npole)
#define IS_SPOLE (grid.is_spole)
#endif

#define KLOPAD (KLO - KPAD)
#define KHIPAD (KHI + KPAD)
#define JLOPAD (JFIRST - JPAD)
#define JHIPAD (JLAST + JPAD)
#define JHIPADPV (JHI + JPAD)
#define ILOPAD (ILO - IPAD)
#define IHIPAD (IHI + IPAD)

#define IADIM (IHI - ILO + 1 + 2 * IPAD)
#define JADIM (JHI - JLO + 1 + 2 * JPAD)
#define KADIM (KHI - KLO + 1 + 2 * KPAD)
#define NELEM2D (IADIM * JADIM)
#define NELEM3D (NELEM2D * KADIM)

/*
 * Declare globally the following shift integers to speed up multidimensional
 * array referencing. Set Ishift = ILO-IPAD, etc. in make_arrays().
 */
extern int Ishift, Jshift, Kshift, Iadim, Jadim, Kadim, Nelem2d, Nelem3d, Shift2d, Shift3d, Shiftkj;

/*
 * Solar longitude, L_s.
 */
#define L_s var.l_s.value

#define U(it, k, j, i) var.u.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d + (it) * Nelem3d]
#define V(it, k, j, i) var.v.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d + (it) * Nelem3d]

#define DUDT(it, k, j, i) var.u.tendency[i + (j) * Iadim + (k) * Nelem2d - Shift3d + (it) * Nelem3d]
#define DVDT(it, k, j, i) var.v.tendency[i + (j) * Iadim + (k) * Nelem2d - Shift3d + (it) * Nelem3d]

#define H(k, j, i) var.h.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define THETA(k, j, i) var.theta.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define FPARA(k, j, i) var.fpara.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

/*
 * Q is mass mixing ratio [density_i/density_dry_air]
 *
 * NOTE: Some textbooks use "q" for specific humidity [density_i/density_total]
 *       and "w" for mass mixing ratio, but w is also vertical velocity.
 */
#define Q(is, ip, k, j, i) var.species[is].phase[ip].q[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

/*
 * The following arrays are related to Q, but are defined locally for a given species
 * and phase.
 *
 * X = n_i/n_tot is mole fraction or number fraction, which for an ideal gases,
 * is also the volume mixing ratio.
 *
 * PARTIAL_H = Q*HDRY is partial hybrid density.
 */
#define X(k, j, i) x[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define PARTIAL_H(k, j, i) partial_h[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

#define NU_TURB(k, j, i) var.nu_turb.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

#define H3(k, j, i) var.h3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define HDRY2(k, j, i) var.hdry2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define HDRY3(k, j, i) var.hdry3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define P2(k, j, i) var.p2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define P3(k, j, i) var.p3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define PDRY3(k, j, i) var.pdry3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define THETA2(k, j, i) var.theta2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define T2(k, j, i) var.t2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define T3(k, j, i) var.t3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define RHO2(k, j, i) var.rho2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define RHO3(k, j, i) var.rho3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define EXNER2(k, j, i) var.exner2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define EXNER3(k, j, i) var.exner3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define FGIBB2(k, j, i) var.fgibb2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define PHI2(k, j, i) var.phi2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define PHI3(k, j, i) var.phi3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define MONT2(k, j, i) var.mont2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

#define HEAT3(k, j, i) var.heat3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define PV2(k, j, i) var.pv2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define EDDY_PV2(k, j, i) var.eddy_pv2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define MOLAR_MASS3(k, j, i) var.molar_mass3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define RI2(k, j, i) var.ri2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define EDDY_REL_VORT2(k, j, i) var.eddy_rel_vort2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define DIV_UV2(k, j, i) var.div_uv2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define W3(k, j, i) var.w3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define Z2(k, j, i) var.z2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define Z3(k, j, i) var.z3.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define DZDT2(k, j, i) var.dzdt2.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define DIFFUSION_COEF_UV(k, j, i)                                                                 \
    var.diffusion_coef_uv.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define DIFFUSION_COEF_THETA(k, j, i)                                                              \
    var.diffusion_coef_theta.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define DIFFUSION_COEF_MASS(k, j, i)                                                               \
    var.diffusion_coef_mass.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define HEAT_MC(k, j, i) var.heat_mc.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define CLOUD_BASE(k, j, i) var.cloud_base.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

#define DQDT_MC(is, ip, k, j, i)                                                                   \
    var.species[is].phase[ip].dqdt[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define LAMBDA_MC(is, k, j, i)                                                                     \
    var.species[is].lambda_mc.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define MB_MC(is, k, j, i) var.species[is].mb_mc.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define CWF(is, k, j, i) var.species[is].cwf.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define DADT(is, k, j, i) var.species[is].dAdt.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define PBASE_MC(is, j, i) var.species[is].pbase_mc.value[i + (j) * Iadim - Shift2d]

#define PHI_IN(k, j, i) phi_in[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define VARPHI(k, j, i) varphi[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

#define PHI_SURFACE(j, i) var.phi_surface.value[i + (j) * Iadim - Shift2d]
#define GRAVITY2(k, j) var.gravity2.value[j + (k) * Jadim - Shiftkj]
#define U_SPINUP(k, j, i) var.u_spinup.value[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

#define PBOT(j, i) var.pbot.value[i + (j) * Iadim - Shift2d]

#define UH(j, i) uh[i + (j) * Iadim - Shift2d]
#define VH(j, i) vh[i + (j) * Iadim - Shift2d]
#define KIN(j, i) kin[i + (j) * Iadim - Shift2d]
#define MONT2D(j, i) mont2d[i + (j) * Iadim - Shift2d]
#define PHI2D(j, i) phi2d[i + (j) * Iadim - Shift2d]
#define MONT2D(j, i) mont2d[i + (j) * Iadim - Shift2d]
#define BERN2D(j, i) bern2d[i + (j) * Iadim - Shift2d]
#define P2D(j, i) p2d[i + (j) * Iadim - Shift2d]
#define RHO2D(j, i) rho2d[i + (j) * Iadim - Shift2d]
#define PV2D(j, i) pv2d[i + (j) * Iadim - Shift2d]
#define H_PT(j, i) h_pt[i + (j) * Iadim - Shift2d]
#define DH_PTDT(j, i) dh_ptdt[i + (j) * Iadim - Shift2d]
#define DBX(j, i) dbx[i + (j) * Iadim - Shift2d]
#define DBY(j, i) dby[i + (j) * Iadim - Shift2d]
#define RHS(j, i) rhs[i + (j) * Iadim - Shift2d]

#define U1D(j) u1d[j - Jshift]
#define P1D(j) p1d[j - Jshift]
#define RHO1D(j) rho1d[j - Jshift]
#define PHI1D(j) phi1d[j - Jshift]
#define MONT1D(j) mont1d[j - Jshift]

#define UM(j, i) um[i + (j) * Iadim - Shift2d]
#define UP(j, i) up[i + (j) * Iadim - Shift2d]
#define U2D(j, i) u2d[i + (j) * Iadim - Shift2d]
#define VM(j, i) vm[i + (j) * Iadim - Shift2d]
#define VP(j, i) vp[i + (j) * Iadim - Shift2d]
#define V2D(j, i) v2d[i + (j) * Iadim - Shift2d]
#define WM(j, i) wm[i + (j) * Iadim - Shift2d]
#define WP(j, i) wp[i + (j) * Iadim - Shift2d]

#undef A
#define A(j, i) a[i + (j) * Iadim - Shift2d]
#define FF(j, i) ff[i + (j) * Iadim - Shift2d]
#define A_PRED(j, i) a_pred[i + (j) * Iadim - Shift2d]
#define A_DIFF1(j, i) a_diff1[i + (j) * Iadim - Shift2d]
#define A_DIFF2(j, i) a_diff2[i + (j) * Iadim - Shift2d]
#define AA1(j, i) aa1[i + (j) * Iadim - Shift2d]
#define AAA(j, i) aaa[i + (j) * Iadim - Shift2d]

#define DELTAA(j, i) deltaa[i + (j) * Iadim - Shift2d]
#define DELA(j, i) dela[i + (j) * Iadim - Shift2d]

#define WHTOP(j, i) whtop[i + (j) * Iadim - Shift2d]
#define WHBOT(j, i) whbot[i + (j) * Iadim - Shift2d]

#undef UU
#define UU(j, i) uu[i + (j) * Iadim - Shift2d]
#define VV(j, i) vv[i + (j) * Iadim - Shift2d]
#define HH(j, i) hh[i + (j) * Iadim - Shift2d]
#define ZE(j, i) ze[i + (j) * Iadim - Shift2d]
#define DI(j, i) di[i + (j) * Iadim - Shift2d]
#define DIV(j, i) div[i + (j) * Iadim - Shift2d]
#define GH1(j, i) gh1[i + (j) * Iadim - Shift2d]
#define GH2(j, i) gh2[i + (j) * Iadim - Shift2d]
#define LPH(j, i) lph[i + (j) * Iadim - Shift2d]
#define LAPH(j, i) laph[i + (j) * Iadim - Shift2d]
#define DDXVAR(j, i) ddxvar[i + (j) * Iadim - Shift2d]
#define D2DX2VAR(j, i) d2dx2var[i + (j) * Iadim - Shift2d]

#define BUFF3D(k, j, i) buff3d[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define BUFFJI(j, i) buffji[i + (j) * Iadim - Shift2d]
#define BUFFJ(j) buffj[j - Jshift]
#define BUFFI(i) buffi[i - Ishift]

#define RH(j, i) rh[i + (j) * Iadim - Shift2d]

#define UNITY(k, j, i) unity[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define OLD(k, j, i) old[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

#define WIND(u, it, k, j, i) (u)->value[i + (j) * Iadim + (k) * Nelem2d - Shift3d + (it) * Nelem3d]
#define DWINDDT(u, it, k, j, i)                                                                    \
    (u)->tendency[i + (j) * Iadim + (k) * Nelem2d - Shift3d + (it) * Nelem3d]

/*
 * Shift macro used in read_array() and write_array()
 */
#define BUFF_SUBARRAY(n, k, j, i)                                                                  \
    buff_subarray[(i) + (j) * i_len + (k) * ji_len + (n) * kji_len - shift_buff]

/*
 * Set BSWAP, for byte-order swapping.
 */
#define DO_SWAP TRUE
#define DONT_SWAP FALSE
#if defined(decalpha) || defined(decmips) || defined(ncube2) || defined(LINUX)
#define BSWAP DO_SWAP
#else
#define BSWAP DONT_SWAP
#endif

/* macros for cloud top pressure */
#define CONSTANT_CLOUD_BASE 0
#define USE_BELT_ZONE 1
#define USE_IF_889NM 2

/* macros for vapor initialization */
#define VMR_COLD_TRAP 0
#define VMR_SLOPE 1

/*
 * Function prototypes:
 */

void make_arrays(void);
void free_arrays(void);

void set_var_props(void);
void free_var_props(void);

int setup_read_array(void);
int setup_write_array(void);

void read_array(int node, int dim, int *start, int *end, char *name, int index, void *array,
                int array_type, int nc_id);

void write_array(int node, int dim, int *start, int *end, int stretch_ni, char *name, int index,
                 void *array, int array_type, int nc_id);

void var_read(char *infile, int portion, unsigned int time_index);

void var_write(char *outfile, int portion, unsigned int time_index, int stretch_ni);

int lookup_netcdf(char *infile, int *nc_id, int *ngatts, char ***gattname, int *num_vars,
                  char ***varname);

void define_netcdf(int portion, int stretch_ni, int nc_id);

void handle_file_compression(char *file_name);

void get_jlohi(int node, int num_nodes, int *jlo, int *jhi);

void get_ilohi(int node, int num_nodes, int *ilo, int *ihi);

void init_phi_surface(void);

int is_phi_surface_file(const struct dirent *entry);

void setup_mu_p(void);

double mu_p(double pressure);

void init_with_u(int floor_tp, int ceiling_tp, init_defaultspec *def);

void init_with_ref(void);

void init_fpara_as_fpe(void);

void init_species(  init_defaultspec   *def, int prompt_mode);
void change_species(change_defaultspec *def, int prompt_mode, int reinit_mode);

void init_vmr_via_deep_value(int is, double *x, double *mole_fraction,
                             double *mole_fraction_over_solar, double *rh_max, 
                             double *vmr_slope, double *vmr_pcrit, int prompt_mode);

void init_vmr_via_obs(int is, double *x, double *mole_fraction, int prompt_mode);

void init_vmr_with_slope(int is, double mole_frac0, double slope);

void trim_species(int is, double *rh_max, int prompt_mode);

double CIRS_data(int index, double p, double lat);

void set_p2_etc(int theta_switch);

void timestep(int action, double **Buff2D);

void adams_bashforth_step(int index);

void leapfrog_step(int index);

void uv_core(double **Buff2D);

void uv_pgrad(double **Buff2D);

void uv_sponge(void);

void uv_vertical_advection(void);

void uv_vert_upwind_3rd_order(void);

double get_sigma(double pbot, double p);

double get_p_sigma(double pbot, double sigma);

void calc_h(int jj, double *p, double *h);

double get_p(int index, int kk, int J, int I);

double get_kin(double *u2d, double *v2d, int kk, int J, int I);

int get_index(char *name);

double molar_mixing_ratio(int is, int ip, int kk, int J, int I);

double get_var(int species_index, int phase_index, int IT, int kk, int J, int I);

double get_var_mean2d(double *a, int index, int kk);

double onto_kk(int index, double topval, double botval, int kk);

double get_brunt2(int kk, int J, int I);

double get_richardson(int kk, int J, int I);

void stability_factor(int J, int I, double *stab_factor);

void get_sounding(double pressure, char *output_name, double *pt_output_value);

double return_cp(double fp, double p, double temp);

double isobaric_specific_heat(double temp, int index);

void store_pgrad_vars(int action, int phi3nk_status);

void store_diag(void);

int cfl_dt(void);

void prompt_extract_on(char *def_extract_str, int *def_extract_species_fraction_type,
                       int *def_mc_diag_extract_sum, boolean modify);

int prompt_species_on(char *def_species_str, boolean modify);

double solar_longitude(time_t date);

double solar_declination(double l_s);

double equation_of_time(time_t date, double l_s);

double east_longitude_of_solar_noon(time_t date, double l_s);

void season_string(double l_s, char *outstring);

double radius_vector(time_t date);

double eccentric_anomaly(double M, double e);

double true_anomaly(double E, double e);

void bcast_char(int node, char *str, int num);

void bcast_int(int node, int *val, int num);

void bcast_double(int node, double *val, int num);

void print_model_description(void);

void print_zonal_info(void);

void print_vertical_column(int J, int I, char *filename);

void node0_barrier_print(char *calling_function, char *Message, int display_time);

void read_spacing_file(init_defaultspec *def, int mode);

int read_t_vs_p(int portion);

void read_meridional_plane(char *infile, int portion, int *np, int *nlat, double *logp, double *lat, double *value);

double t_yp(double p, double latr, int mode, init_defaultspec *def);

double fp_yp(double p, double latr, int mode);

double p_phi(int J, double phi);

void advection(double **Buff2D);

void hsu_advection(int is, int ip, double *buff3d, int direction, double **Buff2D);

void upwind_3rd_order(int is, int ip, double *buff3d, int direction, double **Buff2D);

void relax_vapor();

void phi_from_u(int kk, double *u1d, double *phi1d, int j0, double phi0);

void mont_from_u(int kk, double *u1d, double *mont1d, int j0, double mont0);

void parse_species_name(char *species_name, int *num_elements, char ***symbols, int **counts);

double solar_fraction(char *species, int type, char *min_element);

void relative_humidity(int species_index, double *buffji, int K);

void vertical_modes(int J, int I);

void divergence(int kk, double *uu, double *vv, double *di);

void vorticity(int surface_type, int type, int kk, double *uu, double *vv, double *hh,
               double *pv2d);

double t_rad(int K, int J, int I);

double temp_eq(double lat, double press);

double b_vir(char *chem_a, char *chem_b, double temperature);

double b1_vir(char *chem_a, char *chem_b, double temperature);

double b2_vir(char *chem_a, char *chem_b, double temperature);

double sum_xx(double (*bfunc)(char *, char *, double), double temperature);

double avg_molar_mass(int kk, int J, int I);

void sync_x_to_q(int is, int ip, double *x);

double molar_mass(int index);

double mass_diffusivity(int vapor_index, double temperature, double pressure);

double thermal_conductivity(double temperature);

double u_venus(double p, double lat), u_earth(double p, double lat), u_mars(double p, double lat),
    u_jupiter(double p, double lat), u_saturn(double p, double lat), u_titan(double p, double lat),
    u_uranus(double p, double lat), u_neptune(double p, double lat), u_triton(double p, double lat),
    u_pluto(double p, double lat), u_hot_jupiter(double p, double lat),
    u_null(double p, double lat), u_amp(double p, double lat), cloud_top_jupiter(double lat);

double *extract_scalar(int index, int k, double *buff);

void insert_scalar(int index, int k, double *buff);

void scdswap(char *arr, int nbytes, int cnt);

void thermo_setup(double *cpr);

double return_temp(double fp, double p, double theta);

double alt_return_temp(double fp, double p, double mu, double density);

double return_density(double fp, double p, double theta, double mu, int temp_type);

double p_from_t_rho_mu(double temperature, double rho, double mu);

double return_theta(double fp, double p, double temperature, double *theta_ortho, double *theta_para);

double return_press(double fp, double temperature, double theta);

double rho_minus_rho(double temperature);

void turn_on_phases(void);

/*
 * Global variables used to communicate with rho_minus_rho().
 */
extern double RHOMRHO_fp, RHOMRHO_p, RHOMRHO_mu, RHOMRHO_density;

double th_minus_th_p(double p);
double th_minus_th_t(double temperature);
/*
 * Global variables used to communicate with th_minus_th_p() and
 * th_minus_th_t().
 */
extern double THMTH_fp, THMTH_temperature, THMTH_theta, THMTH_p;

double fpe_minus_fpe(double fpe);
/*
 * Global variables used to communicate with fpe_minus_fpe().
 */
extern double FPEMFPE_p, FPEMFPE_theta;

double sgth_minus_sgth(double sigma);
/*
 * Global variables used to communicate with sgth_minus_sgth().
 */
extern double SGTHMSGTH_theta, SGTHMSGTH_sigmatheta;

double polar_radius(double a, double J2, double GM, double omega);
double polar_radius_cubic(double ep);
/*
 * Global variables used to communicate with polar_radius_cubic().
 */
extern double PRC_J2, PRC_a3w2_GM;

double equatorial_radius(double c, double J2, double GM, double omega);
double equatorial_radius_cubic(double ep);
/*
 * Global variables used to communicate with equatorial_radius_cubic().
 */
extern double ERC_J2, ERC_c3w2_GM;

double p_sigmatheta(double theta, double sigmatheta, double pbot, double hi_p, double lo_p);

double return_enthalpy(double fp, double pressure, double temperature, double *fgibb, double *fpe, double *uoup);

double return_fpe(double temperature);

double blackbody_fraction(double wavelength, double temperature);

/*
 * Use double precision to increase accuracy of diagnostic theta calculations.
 */
double return_sigmatheta(register double theta, register double p, register double pbot);
double f_sigma(register double sigma);
double g_sigma(register double sigma);

/* Values for theta_flag: */
#define THETA_VIA_VAR 0
#define THETA_VIA_REF 1

void set_lonlat(void);
void set_re_rp(void);
void set_fmn(void);
void set_gravity(void);
void set_dsgth(void);

double flux_top(int is, int ip, int J, int I, double *wtop);
double flux_bot(int is, int ip, int J, int I, double *wbot);

void calc_heating(void);

void calc_w(int action);

double pioneer_venus_u(double pressure);

double galileo_u(double pressure, double lat);

double cassini_cirs_u(double pressure);

void source_sink(void);

int restore_mass(int species_index, int phase_index);

void zonal_filter(int index, double *buffji);

void timeplane_bookkeeping(void);

double input_double(char *prompt, double def, boolean modify);

int input_int(char *prompt, int def, boolean modify);

void input_string(char *prompt, char *def, char *ans, boolean modify);

int time_mod(int *time, int step);

void check_nan(char *message);

void bc_lateral(double *pt, int dim);

void epic_error(char *calling_function, char *Message);

void epic_warning(char *challing_function, char *Message);

#if defined(EPIC_MPI)
void mpispec_init(void);
#endif

void declare_copyright(void);

void inquire_radiation_scheme(boolean modify);

void radiative_heating(int action);
void rt_heating(int action);
void newtonian_cooling(int action);
void heating_from_file(int action);
void global_heating_cooling(int action);

void perturbation_heating(void);

/*
 * OpenMARS (formerly MACDA)
 * Running "change -openmars" converts an OpenMARS data file
 * into various EPIC extract.nc files.
 */
void openmars_make_arrays(openmars_gridspec *openmars_grid);
void openmars_free_arrays(openmars_gridspec *openmars_grid);
void openmars_var_read(openmars_gridspec *openmars_grid, char *openmars_infile, int portion, int itime);

void openmars_epic_nc(void);

void openmars_conversion(openmars_gridspec *openmars_grid, char *openmars_infile, char *openmars_outfile_qb,
                         char *openmars_outfile_uvpt);

#define OPENMARS_PS(j, i) openmars_grid->ps[i + openmars_grid->ni * (j)]
#define OPENMARS_TSURF(j, i) openmars_grid->tsurf[i + openmars_grid->ni * (j)]
#define OPENMARS_THETA(k, j, i)                                                                    \
    openmars_grid->theta[i + openmars_grid->ni * (j + openmars_grid->nj * (k))]
#define OPENMARS_TEMP(k, j, i)                                                                     \
    openmars_grid->temp[i + openmars_grid->ni * (j + openmars_grid->nj * (k))]
#define OPENMARS_U(k, j, i) openmars_grid->u[i + openmars_grid->ni * (j + openmars_grid->nj * (k))]
#define OPENMARS_V(k, j, i) openmars_grid->v[i + openmars_grid->ni * (j + openmars_grid->nj * (k))]

#define SHORTWAVE_HEATING(k, j, i) shortwave_heating[i + (j) * Iadim + (k) * Nelem2d - Shift3d]
#define LONGWAVE_HEATING(k, j, i) longwave_heating[i + (j) * Iadim + (k) * Nelem2d - Shift3d]

/*
 * EMARS
 * Running "change -emars" converts an EMARS data file (i.e. emars_v1.0_back_*.nc)
 * into various EPIC extract.nc files.
 */
void emars_make_arrays(emars_gridspec *emars_grid);
void emars_free_arrays(emars_gridspec *emars_grid);
void emars_var_read(emars_gridspec *emars_grid, char *emars_infile, int portion, int itime);

void emars_epic_nc(void);

void emars_conversion(emars_gridspec *emars_grid, char *emars_infile,
                      char *emars_outfile_qb, char *emars_outfile_uvpt);

#define EMARS_PS(j, i) emars_grid->ps[i + emars_grid->ni * (j)]
#define EMARS_TSURF(j, i) emars_grid->tsurf[i + emars_grid->ni * (j)]
#define EMARS_PHI_SURFACE(j, i) emars_grid->phi_surface[i + emars_grid->ni * (j)]
#define EMARS_P(k, j, i) emars_grid->p[i + emars_grid->ni * (j + emars_grid->nj * (k))]
#define EMARS_THETA(k, j, i) emars_grid->theta[i + emars_grid->ni * (j + emars_grid->nj * (k))]
#define EMARS_TEMP(k, j, i) emars_grid->temp[i + emars_grid->ni * (j + emars_grid->nj * (k))]
#define EMARS_U(k, j, i) emars_grid->u[i + emars_grid->ni * (j + emars_grid->nj * (k))]
#define EMARS_V(k, j, i) emars_grid->v[i + emars_grid->ni * (j + emars_grid->nj * (k))]

/*
 * Cloud microphysics macros that reference EPIC model data types.
 */
#define N_0R(is) planet->cloud[is - FIRST_SPECIES].n_0r
#define E_R(is) planet->cloud[is - FIRST_SPECIES].e_r
#define D_Icrit(is) planet->cloud[is - FIRST_SPECIES].d_icrit

#define f1r(is) planet->cloud[is - FIRST_SPECIES].f1r
#define f2r(is) planet->cloud[is - FIRST_SPECIES].f2r
#define f1s(is) planet->cloud[is - FIRST_SPECIES].f1s
#define f2s(is) planet->cloud[is - FIRST_SPECIES].f2s

#define P_EXP_LIQ(is) planet->cloud[is - FIRST_SPECIES].p_exp_liq
#define P_EXP_ICE(is) planet->cloud[is - FIRST_SPECIES].p_exp_ice
#define P_EXP_SNOW(is) planet->cloud[is - FIRST_SPECIES].p_exp_snow

#define ALPHA_RAUT(is) planet->cloud[is - FIRST_SPECIES].alpha_raut
#define Q_LIQ_0(is) planet->cloud[is - FIRST_SPECIES].q_liq_0

/* thermodynamics variables macros */
#define T_triple_pt(is) var.species[is].triple_pt_t     /* Triple point temperature   */
#define p_triple_pt(is) var.species[is].triple_pt_p     /* Triple point pressure      */
#define T_critical_pt(is) var.species[is].critical_pt_t /* Critical point temperature */
#define p_critical_pt(is) var.species[is].critical_pt_p /* Critical point pressure    */

/* microphysics variables macros  */
#define LIQUID_00(is) planet->cloud[is - FIRST_SPECIES].rain_threshold
#define SOLID_00(is) planet->cloud[is - FIRST_SPECIES].snow_threshold
#define RHO_RAIN(is) planet->cloud[is - FIRST_SPECIES].rain_density
#define RHO_SNOW(is) planet->cloud[is - FIRST_SPECIES].snow_density
#define T_00(is) planet->cloud[is - FIRST_SPECIES].t_00
#define GAMMA_SNOW(is) planet->cloud[is - FIRST_SPECIES].gamma_s
#define GAMMA_RAIN(is) planet->cloud[is - FIRST_SPECIES].gamma_r
#define GAMMA_ICE(is) planet->cloud[is - FIRST_SPECIES].gamma_i

/* reference pressure [Pa] */
#define P_REF (1.e+5)

/* Latent heat macros  */
#define Lf(is) var.species[is].Lf                                        /* fusion       */
#define Lv(is) var.species[is].Lv                                        /* vaporization */
#define Ls(is) var.species[is].Ls                                        /* sublimation  */
#define Lc(is, ip0, ip1, t) var.species[is].enthalpy_change(ip0, ip1, t) /* condensation */

/* Coefficients in the power-laws */
#define COEFF_XI(is) planet->cloud[is - FIRST_SPECIES].powerlaw_x_i
#define COEFF_YI(is) planet->cloud[is - FIRST_SPECIES].powerlaw_y_i
#define COEFF_XS(is) planet->cloud[is - FIRST_SPECIES].powerlaw_x_s
#define COEFF_YS(is) planet->cloud[is - FIRST_SPECIES].powerlaw_y_s
#define COEFF_XR(is) planet->cloud[is - FIRST_SPECIES].powerlaw_x_r
#define COEFF_YR(is) planet->cloud[is - FIRST_SPECIES].powerlaw_y_r
#define COEFF_A(is) planet->cloud[is - FIRST_SPECIES].powerlaw_a
#define COEFF_B(is) planet->cloud[is - FIRST_SPECIES].powerlaw_b
#define COEFF_C(is) planet->cloud[is - FIRST_SPECIES].powerlaw_c
#define COEFF_D(is) planet->cloud[is - FIRST_SPECIES].powerlaw_d
#define COEFF_M(is) planet->cloud[is - FIRST_SPECIES].powerlaw_m /* alpha */
#define COEFF_N(is) planet->cloud[is - FIRST_SPECIES].powerlaw_n /* beta  */

#define ANT_AS(is) planet->cloud[is - FIRST_SPECIES].a_s
#define ANT_BS(is) planet->cloud[is - FIRST_SPECIES].b_s
#define ANT_AL(is) planet->cloud[is - FIRST_SPECIES].a_l
#define ANT_BL(is) planet->cloud[is - FIRST_SPECIES].b_l

/* Schmidt parameter D89 p3103 eq.B13 */
#define SC(k, j, i) (dynvis / (rho * mass_diffusivity(is, T3(k, j, i), P3(k, j, i))))

/* Dry air density */
#define RHO_DRY(k, j, i) (PDRY(k, j, i) / (T(k, j, i) * planet->rgas))

/* Gas constant macro for vapor of all the species */
#define GAS_R(is) (R_GAS / var.species[is].molar_mass)

/*
 * Cloud microphysics functions that reference EPIC model data types like
 * planetspec and double_triplet (which are not defined in epic_microphysics.h).
 */
void cloud_microphysics(void);

/*
 * RS 02/24/2020 Note that moist_convection() also references datatypes in epic.h.
 */
void moist_convection(void);

void calculate_base_pressure(int **kbase);

void RAS_flux(int is, int kbase_ji, int J, int I, double *cwf_is);

double terminal_velocity(int is, int ip, double temperature, double pressure, double partial_rho);

int read_enthalpy_change_data(int species_index, double_triplet **pt_hi, double_triplet **pt_hf,
                              double_triplet **pt_hg);

double enthalpy_change(int species_index, int init_phase, int final_phase, double temperature,
                       int ndat, double_triplet *hi, double_triplet *hf, double_triplet *hg);
double enthalpy_change_H_2O(int init_phase, int final_phase, double temperature);
double sat_vapor_p_H_2O(double temperature);

double enthalpy_change_NH_3(int init_phase, int final_phase, double temperature);
double sat_vapor_p_NH_3(double temperature);

double enthalpy_change_H_2S(int init_phase, int final_phase, double temperature);
double sat_vapor_p_H_2S(double temperature);

double enthalpy_change_CH_4(int init_phase, int final_phase, double temperature);
double sat_vapor_p_CH_4(double temperature);

double enthalpy_change_C_2H_2(int init_phase, int final_phase, double temperature);
double sat_vapor_p_C_2H_2(double temperature);

double enthalpy_change_C_2H_4(int init_phase, int final_phase, double temperature);
double sat_vapor_p_C_2H_4(double temperature);

double enthalpy_change_C_2H_6(int init_phase, int final_phase, double temperature);
double sat_vapor_p_C_2H_6(double temperature);

double enthalpy_change_CO_2(int init_phase, int final_phase, double temperature);
double sat_vapor_p_CO_2(double temperature);

double enthalpy_change_NH_4SH(int init_phase, int final_phase, double temperature);
double sat_vapor_p_NH_4SH(double temperature);

double enthalpy_change_O_3(int init_phase, int final_phase, double temperature);
double sat_vapor_p_O_3(double temperature);

double enthalpy_change_N_2(int init_phase, int final_phase, double temperature);
double sat_vapor_p_N_2(double temperature);

double enthalpy_change_PH_3(int init_phase, int final_phase, double temperature);
double sat_vapor_p_PH_3(double temperature);

/* * * * * * * * * *  end of epic.h  * * * * * * * * * * * * * * * * * * * * */
#endif
