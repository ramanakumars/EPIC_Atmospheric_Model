/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2002-2023 Timothy E. Dowling                      *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License     *
 * as published by the Free Software Foundation; either version 2  *
 * of the License, or (at your option) any later version.          *
 * A copy of this License is in the file:                          *
 *   $EPIC_PATH/GNU_General_Public_License.                        *
 *                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
 *                                                                 *
 * You should have received a copy of the GNU General Public       *
 * License along with this program; if not, write to the Free      *
 * Software Foundation, Inc., 59 Temple Place - Suite 330,         *
 * Boston, MA  02111-1307, USA.                                    *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * gen_t_vs_p.c  * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Smooth an input t_vs_p.Uranus with a Savitzky-Golay filter,     *
 * and record these details in the file header.                    *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* 
 * Compile with:                                             
 * cc -I$EPIC_PATH/include -DEPIC_PATH=\"$EPIC_PATH\" -lm -o gen_t_vs_p gen_t_vs_p.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*
 * Global
 */
char
  Message[128];

/*
 * Set minimum N^2 value [s^-2].
 * Set to -1.e+20 for no effect.
 *
 * NOTE: The resulting values of N^2 may go below N2_floor
 *       because here we are assuming representative constant parameters.
 *       The idea is to experiment with N2_floor and the smoothing to get the desired result.
 */
#define N2_floor 0.00009

/*
 * Function prototypes:
 */
void savitzky_golay(double *c,
                    int     np,
                    int     nl,
                    int     nr,
                    int     ld,
                    int     m);
void lu_decompose(int     n,
                  double *a,
                 int     *index,
                 double  *d);
void lu_backsub(int     n,
                double *a,
                int    *index,
                double *b);
int *ivector(int  nl, 
             int  nh);
void free_ivector(int  *m, 
                  int  nl, 
                  int  nh);
double *dvector(int nl, 
                int nh);
void  free_dvector(double *m, 
                   int     nl, 
                   int     nh);
void util_error(char *calling_function,
                char *Message);

void input_string(char   *prompt, 
                  char   *def, 
                  char   *ans);

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

#undef C
#define C(i,j) c[j+i*(np)]

#define N_STR 128

int main(int   argc,
         char *argv[]) 
{
  char   
    input_file[N_STR],
    header[N_STR];
  int
    kk,k,ntp;
  int
    np,nhw,nl,nr;
  double
    *tdat,*pdat,*dtdat,
    *smooth_t,
    *count,
    *c,
     pup,pdn,
     a,b,
     tmp;
  FILE 
    *infile,
    *t_vs_p;

  /*
   * Prompt for input file, and open it.
   */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="main";

  sprintf(header,"");
  input_string("Input t_vs_p file:",header,input_file);
  infile = fopen(input_file,"r");
  if (!infile) {
    sprintf(Message,"Unable to open input file %s",input_file);
    util_error(dbmsname,Message);
  }
  
  /* Read 6-line header */
  for (kk = 0; kk < 6; kk++) {
    fgets(header,N_STR,infile);
  }
  /* Read number of data points */
  fscanf(infile,"%d",&ntp);

  /* Allocate memory */
  pdat     = dvector(0,ntp-1);
  tdat     = dvector(0,ntp-1);
  dtdat    = dvector(0,ntp-1);
  smooth_t = dvector(0,ntp-1);
  count    = dvector(0,ntp-1);

  /* 
   * Store in order of increasing altitude.
   */
  for (kk = ntp-1; kk >= 0; kk--) { 
    fscanf(infile,"%lf %lf %lf",pdat+kk,tdat+kk,dtdat+kk);
  }

  fclose(infile);

  if (N2_floor > -1.e+20) {
    /*
     * Modify T(p) such that N^2 >= N2_floor.
     */
    double
      g     = 8.69,  /* Uranus gravity [m s-2] */
      cp    = 9280.,  /* c_p [J kg-1 K-1]        */
      kappa = .375;   /* R/c_p                   */
    double
      T1,Tinfty,p_p0_kappa;

    Tinfty = g*g/(N2_floor*cp);

    for (k = 1; k < ntp; k++) {
      p_p0_kappa = pow(pdat[k]/pdat[k-1],kappa);
      T1         = Tinfty*p_p0_kappa/(Tinfty/tdat[k-1]+p_p0_kappa-1.);
      tdat[k]    = MIN(tdat[k],T1);
    }
  }
  
  /* 
   * Smooth temperature data using a Savitzky-Golay filter.
   *
   * Set half-width, nhw.
   */
  nhw = 15;

  np = 2*nhw+1;
  c  = dvector(0,np*np-1);
  /*
   * Calculate Savitzky-Golay coefficients.
   */
  for (nl = 0; nl < np; nl++) {
    nr = np-nl-1;
    savitzky_golay(c+nl*np,np,nl,nr,0,2);
  }

  for (k = 0; k < ntp; k++) {
    if (k < nhw) {
      /* Skewed by left boundary */
      nl = k;
      nr = np-nl-1;
    }
    else if (k > ntp-1-nhw) {
      /* Skewed by right boundary */
      nr = ntp-1-k;
      nl = np-nr-1;
    }
    else {
      /* Centered */
      nl = nhw;
      nr = nhw;
    }
    smooth_t[k] = C(nl,0)*tdat[k];
    for (kk = 1; kk <= nl; kk++) {
      smooth_t[k] += C(nl,kk)*tdat[k-kk];
    }
    for (kk = 1; kk <= nr; kk++) {
      smooth_t[k] += C(nl,np-kk)*tdat[k+kk];
    } 
  }

  /* Open t_vs_p.Uranus.smoothed */
  t_vs_p = fopen(EPIC_PATH"/data/Uranus/t_vs_p.Uranus.smoothed","w");
  if (!t_vs_p) {
    fprintf(stderr,"Error: cannot open %s\n",EPIC_PATH"/data/Uranus/t_vs_p.Uranus.smoothed");
    exit(1);
  }

  /* Write preamble: */
  fprintf(t_vs_p," Temperature versus pressure for Uranus, generated by epic/data/Uranus/gen_t_vs_p.c.\n");
  fprintf(t_vs_p," Input: %s\n",input_file);
  fprintf(t_vs_p,"\n");
  fprintf(t_vs_p," Data have been smoothed with a Savitzky-Golay filter, nhw = %d. N2_floor = %9.7f.\n",nhw,N2_floor);
  fprintf(t_vs_p," NOTE: The third column is needed for the file's format, but used rarely in some Newtonian cooling schemes.\n");
  fprintf(t_vs_p,"#   p[hPa]      T[K]     dT[K]\n");

  /* 
   * Write number of data points.
   */
  fprintf(t_vs_p,"%d\n",ntp);

  /* Write data: */
  for (kk = ntp-1; kk >= 0; kk--) {
    fprintf(t_vs_p,"   %10.3e  %7.2f     0.\n",pdat[kk],smooth_t[kk]);
  }
  fclose(t_vs_p);

  /* Free allocated memory: */
  free_dvector(c,       0,np*np-1);
  free_dvector(smooth_t,0,ntp-1  );
  free_dvector(pdat,    0,ntp-1  );
  free_dvector(tdat,    0,ntp-1  );
  free_dvector(dtdat,   0,ntp-1  );

  return 0;
}

/*======================= savitzky_golay() ===================================*/

/*
 * Calculate Savitzky-Golay weighting coefficients to compute a smooth value
 * or smooth derivative of noisy data.
 *
 * Based on Numerical Recipes in C, p. 652.
 * Note the wrap-around ordering of the coefficients returned in c[].
 *
 * Assumes zero-based arrays.
 */

#undef  A
#define A(i,j) a[j+(m+1)*i]
 
void savitzky_golay(double *c,
                    int     np,
                    int     nl,
                    int     nr,
                    int     ld,
                    int     m)
{
  int
    imj,ipj,j,k,kk,mm;
  int
   *index;
  double
    d,fac,sum;
  double
    *a,
    *b;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="savitzky_golay";

  /*
   * Screen for inconsistent arguments.
   */
  if (np < nl+nr+1) {
    sprintf(Message,"np=%d < nl+nr+1=%d",np,nl+nr+1);
    util_error(dbmsname,Message);
  }
  else if (nl < 0) {
    sprintf(Message,"nl=%d < 0",nl);
    util_error(dbmsname,Message);
  }
  else if (nr < 0) {
    sprintf(Message,"nr=%d < 0",nr);
    util_error(dbmsname,Message);
  }
  else if (ld > m) {
    sprintf(Message,"ld=%d > m=%d",ld,m);
    util_error(dbmsname,Message);
  }
  else if (nl+nr < m) {
    sprintf(Message,"nl+nr=%d < m=%d",nl+nr,m);
    util_error(dbmsname,Message);
  }

  /* Allocate memory */
  index = ivector(0,m);
  a     = dvector(0,(m+1)*(m+1)-1);
  b     = dvector(0,m);

  for (ipj = 0; ipj <= (m << 1); ipj++) {
    sum = (ipj ? 0. : 1.);
    for (k = 1; k <= nr; k++) {
      sum += pow((double)k,(double)ipj);
    }
    for (k = 1; k <= nl; k++) {
      sum += pow((double)-k,(double)ipj);
    }
    mm = IMIN(ipj,2*m-ipj);
    for (imj = -mm; imj <= mm; imj+=2) {
      A((ipj+imj)/2,(ipj-imj)/2) = sum;
    }
  }

  lu_decompose(m+1,a,index,&d);

  for (j = 0; j < m+1; j++) {
    b[j] = 0.;
  }
  b[ld] = 1.0;

  lu_backsub(m+1,a,index,b);

  for (kk = 0; kk < np; kk++) {
    c[kk] = 0.;
  }

  for (k = -nl; k <= nr; k++) {
    sum = b[0];
    fac = 1.;
    for (mm = 0; mm < m; mm++) {
      sum += b[mm+1]*(fac *= k);
    }
    kk = ((np-k)%np);
    c[kk] = sum;
  }

  /* Free allocated memory. */
  free_ivector(index,0,m);
  free_dvector(a,0,(m+1)*(m+1)-1);
  free_dvector(b,0,m);

  return;
}

/*======================= end of savitzky_golay() ============================*/

/*======================= lu_decompose() =====================================*/

/*
 * Adapted from Numerical Recipes in C, pp. 46-47.
 * Assumes zero-based indexing, with (i,j) referring to (row,column).
 */

#undef  A
#define A(i,j) a[(j)+n*(i)]

void lu_decompose(int     n,
                  double *a,
                  int    *index,
                  double *d)
{
  int
    i,imax,j,k;
  static int
    n_max = 0;
  double
    tiny = 1.e-20,
    big,dum,sum,temp;
  static double
   *vv;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lu_decompose";

  /*
   * Allocate memory.
   */
  if (n_max == 0) {
    n_max = n;
    vv    = dvector(0,n_max-1);
  }
  else if (n > n_max) {
    free_dvector(vv,0,n_max-1);
    n_max = n;
    vv    = dvector(0,n_max-1);
  }

  *d = 1.;
  for (i = 0; i < n; i++) {
    big = 0.;
    for (j = 0; j < n; j++) {
      temp = fabs(A(i,j));
      if (temp > big) {
        big = temp;
      }
    }
    if (big == 0.) {
      sprintf(Message,"matrix singular");
      util_error(dbmsname,Message);
    }
    vv[i] = 1./big;
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = A(i,j);
      for (k = 0; k < i; k++) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
    }
    big = 0.;
    for (i = j; i < n; i++) {
      sum = A(i,j);
      for (k = 0; k < j; k++) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
 
      dum = vv[i]*fabs(sum);
      if (dum >= big) {
        big  = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum       = A(imax,k);
        A(imax,k) = A(j,   k);
        A(j,   k) = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    index[j] = imax;

    if (A(j,j) == 0.) {
      A(j,j) = tiny;
      fprintf(stderr,"Warning, %s(), matrix is singular\n",dbmsname);
    }

    if (j != n-1) {
      dum = 1./A(j,j);
      for (i = j+1; i < n; i++) {
        A(i,j) *= dum;
      }
    }
  }
 
  return;
}

/*======================= end of lu_decompose() ==============================*/

/*======================= lu_backsub() =======================================*/

/*
 * Adapted from Numerical Recipes in C, p. 47.
 * Assumes zero-based indexing, with (i,j) referring to (row,column).
 */

#undef  A
#define A(i,j) a[(j)+n*(i)]

void lu_backsub(int     n,
                double *a,
                int    *index,
                double *b)
{
  int
    ii = -1,
    i,ip,j; 
  double
    sum;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lu_backsub";

  for (i = 0; i < n; i++) {
    ip    = index[i];
    sum   = b[ip];
    b[ip] = b[i];
    if (ii != -1) {
      for (j = ii; j <= i-1; j++) {
        sum -= A(i,j)*b[j];
      }
    }
    else if (sum != 0.) {
      ii = i;
    }
    b[i] = sum;
  } 

  for (i = n-1; i >= 0; i--) {
    sum = b[i];
    for (j = i+1; j < n; j++) {
      sum -= A(i,j)*b[j];
    }
    b[i] = sum/A(i,i);
  }

  return;
}

/*======================= end of lu_backsub() ================================*/

/*======================= ivector() ============================================*/
      
/*
 *  Allocates memory for a 1D int array 
 *  with range [nl..nh].
 */

#undef  DEBUG

int *ivector(int  nl, 
             int  nh)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  int         
    *m;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="ivector";

  if (nh < nl) {
    sprintf(Message,"range (%d,%d)",nl,nh);
    util_error(dbmsname,Message);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (int *)calloc(len_safe,sizeof(int));
  if (!m) {
    sprintf(Message,"nl=%d,nh=%d,len_safe=%d",nl,nh,len_safe);
    util_error(dbmsname,Message);
  }
  m -= nl_safe;

  return m;
}

/*======================= end of ivector() ====================================*/

/*======================= free_ivector() ======================================*/

#undef  DEBUG

/*
 *  Frees memory allocated by ivector().
 */

void free_ivector(int  *m, 
                  int  nl, 
                  int  nh)
{
  int  
    nl_safe;

  if (m) {
    nl_safe = (nl < 0) ? nl : 0;
    m += nl_safe;
    free(m);
  }

#if defined(DEBUG)
  if (m) {
    fprintf(stderr,"error in free_ivector()\n");
  }
  else {
    fprintf(stderr,"error in free_ivector()\n");
  }
#endif

  return;
}

/*======================= end of free_ivector() ================================*/

/*======================= dvector() ==========================================*/

double *dvector(int nl, 
                int nh)
      /*
       *  Allocates memory for a 1D double array 
       *  with range [nl..nh].
       */
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  double         
    *m;

  if (nh < nl) {
    fprintf(stderr,"called dvector(%d,%d) \n",nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe - nl_safe + 1);

  m = (double *)calloc(len_safe, sizeof(double));
  if (!m) {
    fprintf(stderr, "calloc error in dvector \n");
    exit(1);
  }
  m -= nl_safe;
  return m;
}

/*======================= end of dvector() ====================================*/

/*======================= free_dvector() ======================================*/

void  free_dvector(double *m, 
                   int     nl, 
                   int     nh)
      /*
       *  Frees memory allocated by dvector().
       */
{
  int  
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m += nl_safe;
  free(m);
}

/*======================= end of free_dvector() ===============================*/

/*======================= util_error() =======================================*/
/*
 * Prints calling function name, and Message to stderr, then aborts.
 */
void util_error(char *calling_function,
                char *Message)
{
  fprintf(stderr,"\n** Error: %s(): %s\n",calling_function,Message);
  fflush(stderr);

  exit(1);
}

/*======================= end of util_error() ================================*/

/*====================== input_string() ======================================*/

/* 
 * Read in a string, or set to default if input is a return ('\n').
 *
 * C.Santori, T.Dowling 
 */

void input_string(char *prompt, 
                  char *def, 
                  char *ans) 
{
  char 
    c,
    buffer[N_STR];
  int  
    len;

  fprintf(stdout,"%s[%s]: ",prompt,def);
  for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
    buffer[len]=c;
  }
  buffer[len] = '\0';
  if (len == 0) {
    strcpy(ans,def);
  }
  else {
    strcpy(ans,buffer);
    strcpy(def,buffer);
  }
}

/*====================== end input_string() ==================================*/

/* * * * * * * * * * end of gen_t_vs_p.c * * * * * * * * * * * * * */
