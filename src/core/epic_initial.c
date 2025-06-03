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

/* * * * * * * * * * epic_initial.c  * * * * * * * * * * * * * * * * * * * * 
 *                                                                         *
 *  Timothy E. Dowling                                                     *
 *                                                                         *
 *  Creates zonally-symmetric initial epic.nc input file for the           *
 *  EPIC atmospheric model.                                                *
 *                                                                         *
 *  The horizontal differencing scheme uses the staggered C grid.          *
 *                                                                         *
 *  For the full-globe geometry, the numbering is:                         *
 *                                                                         *
 *                            -180 deg             +180 deg                *
 *                               |                   |                     *
 *       j = nj+1   -- pv---0----pv---0----pv---0---pv -- north pole       *
 *                     |    :    |    :    |    :    |                     *
 *                     |    :    |    :    |    :    |                     *
 *       j = nj+1/2 -- u....h....u....h....u....h....u                     *
 *                     |    :    |    :    |    :    |                     *
 *       j = nj     -- pv---v----pv---v----pv---v---pv                     *
 *                     |    :    |    :    |    :    |                     *
 *                     u....h....u11..h11..u....h....u -- equator          *
 *                     |    :    |    :    |    :    |                     *
 *       j = 1      -- pv---v---pv11--v11--pv---v---pv                     *
 *                     |    :    |    :    |    :    |                     *
 *       j = 1/2    -- u....h....u....h....u....h....u                     *
 *                     |    :    |    :    |    :    |                     *
 *                     |    :    |    :    |    :    |                     *
 *       j = 0      -- pv---0----pv---0----pv---0---pv -- south pole       *
 *                     |         |         |         |                     *
 *                   i = 0     i = 1     i = ni    i = ni+1                *
 *                                                                         *
 *   In the EPIC model, the vertical coordinate is called "sigmatheta"     *
 *   and sometimes abbreviated "sgth", to indicate its flexible nature.    *
 *                                                                         *
 *   The vertical staggering of prognostic variables is shown below. The   *
 *   variable q stands for the mass mixing ratios of optional species.     *
 *                                                                         *
 *     k = 1/2    .......theta,q,fpara..... sigmatheta[1]      = sgth_top  *
 *                                                                         *
 *     k = 1             u,v,h,nu_turb      sigmatheta[2]                  *
 *                                                                         *
 *     k = 3/2    .......theta,q,fpara..... sigmatheta[3]                  *
 *                                                                         *
 *     k = nk            u,v,h,nu_turb      sigmatheta[2*nk]               * 
 *                                                                         *
 *     k = nk+1/2 .......theta,q,fpara..... sigmatheta[2*nk+1] = sgth_bot  *
 *                                                                         *
 *   The vertical staggering of key diagnostic variables is shown below.   *
 *                                                                         *
 *     k = 1/2    .......w3,heat3,p3....... sigmatheta[1]      = sgth_top  *
 *                                                                         *
 *     k = 1             phi2,div_uv2       sigmatheta[2]                  *
 *                                                                         *
 *     k = 3/2    .......w3,heat3,p3....... sigmatheta[3]                  *
 *                                                                         *
 *     k = nk            phi2,div_uv2       sigmatheta[2*nk]               * 
 *                                                                         *
 *     k = nk+1/2 .......w3,heat3,p3....... sigmatheta[2*nk+1] = sgth_bot  *
 *                                                                         *
 *     (gas giants have an abyssal layer, k = nk+1)                        *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <epic_datatypes.h>
#include <epic.h>
int MODIFY = 1;

/*
 *  Data from the books "Venus," "Mars," "Saturn," "Uranus," 
 *  University of Arizona Press; Lindal et al (1983, Icarus) for Titan. 
 *
 *  NOTE: planet->cp and planet->kappa are reassigned according to the value of cpr 
 *        returned from the initialization routine thermo_setup(), which is the
 *        low-temperature-limit value of cp/rgas, in order to maintain self-consistency.
 *        Use the function return_cp() rather than planet->cp unless you actually want
 *        the low-temperature reference value.  
 */

init_defaultspec 
  defaults;

/*
 * Function prototypes:
 */
void read_defaults(init_defaultspec *def);

void write_defaults(init_defaultspec *def);

/*
 * The growth factor 1.+ALPHA_BOUNDARY is used
 * for geometrically decreasing boundary-layer spacing 
 * in the top sponge and bottom planetary boundary layer (PBL).
 */
#undef  ALPHA_BOUNDARY
#define ALPHA_BOUNDARY 0.50

/*
 * Make sure these indices line up with the string array advection_scheme[] below.
 */
#define HSU_ADVECTION          0
#define THIRD_ORDER_UPWIND     1
   
/*======================= main() ============================================*/

int main(int   argc,
         char *argv[])
{
  char   
    header[N_STR],       /*  character string buffer                     */
    system_id[32],       /*  Name of planet or system to study           */
    outfile[FILE_STR],   /*  used to output data                         */
    sflag[80],           /*  string to hold command-line flags           */
    out_str[128];        /*  string to hold user inquiries               */ 
  register int    
    spacing_type,        /*  choice of layer spacing                     */
    ki,ii,index,ip,      /*  utility indices                             */
    K,J,I,kk,            /*  counters                                    */
    floor_tp,            /*  bottom index of shortened t_vs_p data       */
    ceiling_tp,          /*  top index of shortened t_vs_p data          */
    is,                  /*  species index                               */
    count,
    nk_reg;
  unsigned int
    time_index = 0;
  int
    itmp,
    num_advection_schemes      = 2,
    num_uv_timestep_schemes    = 1,
    num_turbulence_schemes     = 3;
  const char
    *radiation_scheme[]
      = {"off",
         "Correlated k",
         "Newtonian",
         "Heating from file",
         "Global heating-cooling"};
  /*
   * NOTE: The name for advection schemes in non-flux form (aka non-divergence form, not h-weighted)
   *       need to lead with the 13 characters "Non-flux form".
   */
  const char
    *advection_scheme[]
      = {"Predictor-corrector (Hsu, Konor, and Arakawa)",
         "Non-flux form, 3rd-order upwind"},
    *uv_timestep_scheme[]
      = {"3rd-order Adams-Bashforth"},
    *turbulence_scheme[]
      = {"on",
         "on_vertical_only",
         "off"},
    *vert_coord_type[]
      = {"isentropic","isobaric","hybrid"};
  double
    fgibb,fpe,uoup,
    theta,theta_ortho,theta_para,
    temperature,
    ptop,pbot,avg;
  static double
    *Buff2D[NUM_WORKING_BUFFERS];
  register double
    tmp,tmp2,            /*  temporary storage                          */
    fpara,pressure,mu,
    sigma,sg1,sg2,slope,theta_knee,sigma_knee,
    sgth,log_sgth,neglogp,negp,
    dlnp;
  double  
    *neglogpdat,         /*  -log(pdat)                                 */
    *thetadat,           /*  theta corresponding to t_vs_p data         */
     p_d,sgth_d;
  double_triplet
    *buff_triplet;
  struct tm
    date_start;
  time_t
    current_time;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="epic_initial";

#if defined(EPIC_MPI)
  sprintf(Message,"not designed to be run on multiple processors");
  epic_error(dbmsname,Message);
#endif
  
  /* EPIC Model version number: */
  grid.epic_version = 5.30;

  declare_copyright();

  /* 
   * Interpret command-line arguments: 
   */
  if (argc > 1) {
    for (count = 1; count < argc; count++) {
      sscanf(argv[count],"%s",sflag);
      if (strcmp(sflag,"-help") == 0 ||
          strcmp(sflag,"-h")    == 0) {
        /* Print help, exit: */
        system("more "EPIC_PATH"/help/epic_initial.help");
        exit(1);
      }
      else {
        fprintf(stderr,"Unrecognized epic_initial command-line flag: %s \n",sflag);
        exit(1);
      }
    }
  }

  /*
   *  Read default parameter settings:
   */
  read_defaults(&defaults);

  if (IAMNODE == NODE0) {
    int
      current_year;

    current_time = time(NULL);
    strftime(header,N_STR,"%Y",localtime(&current_time));
    current_year = atoi(header);

    /*  Print welcome statement: */
    fprintf(stdout,"\n");    
    fprintf(stdout,"        ______/      ___       __   ___/        __     \n");
    fprintf(stdout,"       /            /     /       /          /   __/ \n");  
    fprintf(stdout,"      ______/      ______/       /          /       \n");
    fprintf(stdout,"     /            /             /          (       \n");    
    fprintf(stdout," _________/   ___/       _________/    \\_______/\n");
    fprintf(stdout,"\n\n            WELCOME TO THE EPIC MODEL\n");
    fprintf(stdout,"       Celebrating %1d years as open source\n\n",current_year-1998);
    fprintf(stdout,"             Version: %4.2f\n",grid.epic_version);
    fprintf(stdout,"      Floating-point: double precision\n");
    fprintf(stdout,"          Geometries: globe, f-plane \n");
    fprintf(stdout,"         Atmospheres: Venus, Earth, Jupiter, Saturn, Titan, Uranus, Neptune \n");
    fprintf(stdout,"          Benchmarks: Held_Suarez, Venus_LLR05\n");
    fprintf(stdout,"  Under construction: Mars, Triton, Pluto, Brown Dwarfs, Extrasolar Planets \n");
    fprintf(stdout,"      Google Scholar: https://scholar.google.com/citations?user=6LJV7-UAAAAJ&hl=en\n\n");
  }

  /*
   * Choose atmosphere or benchmark to initialize:
   */
  ii = FALSE;
  while (ii == FALSE) {
    char
      *ptr;

    input_string("Choose atmosphere or benchmark to initialize \n",defaults.system_id,system_id,MODIFY);

    /* 
     * Set up system.
     *
     * NOTE: The easiest way to add a new planet is to start with a similar existing system
     * and edit it appropriately. 
     */
    ii = TRUE;
    if      (strcmp(system_id,"Venus")   == 0) {
      planet                  = &venus;
    }
    else if (strcmp(system_id,"Earth") == 0) {
      planet                  = &earth;
    }
    else if (strcmp(system_id,"Mars") == 0) {
      planet                  = &mars;
    }
    else if (strcmp(system_id,"Jupiter") == 0) {
      planet                  = &jupiter;
    }
    else if (strcmp(system_id,"Saturn")  == 0) {
      planet                  = &saturn;
    }
    else if (strcmp(system_id,"Titan") == 0) {
      planet                  = &titan;
    }
    else if (strcmp(system_id,"Uranus")  == 0) {
      planet                  = &uranus;
    }
    else if (strcmp(system_id,"Neptune") == 0) {
      planet                  = &neptune;
    }
    else if (strcmp(system_id,"Pluto") == 0) {
      planet                  = &pluto;
    }
    else if (strcmp(system_id,"Hot_Jupiter") == 0) {
      planet                  = &hot_jupiter;
    }
    else if (strcmp(system_id,"Held_Suarez") == 0) {
      planet                  = &held_suarez;
    }
    else if (strcmp(system_id,"Goldmine") == 0) {
      planet                  = &goldmine;
    }
    else if (strcmp(system_id,"Venus_LLR05") == 0) {
      planet                  = &venus_llr05;
    }
    else {
      /* unimplemented system */
      sprintf(Message,"\n\"%s\" not defined (Note: system names are capitalized).\n\n",system_id);
      epic_error(dbmsname,Message);
    }
  }

  /*
   * Set vertical coordinate type, both as a readable string, grid.vertical_coordinate, and
   * as an integer, grid.coord_type, the latter for use in switch() branches (which are faster
   * than strcmp() calls).
   */
  ii = COORD_FIRST;
  fprintf(stdout,"Input vertical-coordinate type: %d => %s\n",ii,vert_coord_type[ii-1]);
  ii++;
  sprintf(header,"                                %d => %s\n",ii,vert_coord_type[ii-1]);
  for (ii = COORD_FIRST+2; ii <= COORD_LAST; ii++) {
    sprintf(Message,"                                %d => %s\n",ii,vert_coord_type[ii-1]);
    strcat(header,Message);
  }

  grid.coord_type = defaults.coord_type = input_int(header,defaults.coord_type,MODIFY);
  if (grid.coord_type < COORD_FIRST ||
      grid.coord_type > COORD_LAST    ) {
    fprintf(stdout,"%d is out of valid range [%d,%d]:  %d => %s\n",
                    grid.coord_type,COORD_FIRST,COORD_LAST,COORD_FIRST,vert_coord_type[COORD_FIRST-1]);
  }
  else {
    sprintf(grid.vertical_coordinate,"%s",vert_coord_type[grid.coord_type-1]);
  }

  /*
   *  Initialize time.
   */
  strftime(header,N_STR,"%Y_%m_%d_%H:%M:%S (UTC)",gmtime(&defaults.start_time));
  if (defaults.start_date_input_type == NOT_SET) {
    defaults.start_date_input_type = 0;
    sprintf(Message,"Input starting date and time for model's time dimension:\n"
                    " 0 => default, %s",header);
  }
  else {
    sprintf(Message,"Input starting date and time for model's time dimension:\n"
                    " 0 => previous, %s",header);
  }
  fprintf(stdout,"%s\n",Message);
  current_time = time(NULL);
  strftime(header,N_STR,"%Y %m %d %H:%M:%S (UTC)",gmtime(&current_time));
  sprintf(Message," 1 => current date and time, %s",header);
  fprintf(stdout,"%s\n",Message);
  defaults.start_date_input_type = 
    input_int(" 2 => prompt for values (UTC)\n",defaults.start_date_input_type,MODIFY);
  if (defaults.start_date_input_type == 0) {
    var.start_time = defaults.start_time;
  }
  else if (defaults.start_date_input_type == 1) {
    var.start_time = time(NULL);
  }
  else if (defaults.start_date_input_type == 2) {
    /* Trigger timezone evaluation */
    date_start.tm_year  =  0;
    date_start.tm_mon   =  0;
    date_start.tm_mday  =  0;
    date_start.tm_hour  =  0;
    date_start.tm_min   =  0;
    date_start.tm_sec   =  0;
    date_start.tm_isdst = -1;
    mktime(&date_start);

    defaults.UTC_start = *gmtime(&defaults.start_time);
    date_start.tm_year = input_int("  Input year (YYYY):",defaults.UTC_start.tm_year+1900,MODIFY)-1900;
    date_start.tm_mon  = input_int("   Input month (MM):",defaults.UTC_start.tm_mon+1,MODIFY)-1;
    date_start.tm_mday = input_int("     Input day (DD):",defaults.UTC_start.tm_mday,MODIFY);
    date_start.tm_hour = input_int("  Input hour (0-23):",defaults.UTC_start.tm_hour,MODIFY);
    date_start.tm_min  = input_int("Input minute (0-59):",defaults.UTC_start.tm_min,MODIFY);
    date_start.tm_sec  = input_int("Input second (0-59):",defaults.UTC_start.tm_sec,MODIFY);
    /*
     * Shift UTC input into local time (because "struct tm" does not keep track of timezone and assumes time is local).
     */
    date_start.tm_hour -= timezone/3600;  /* timezone is defined in time.h */

    /*
     * Make and store calendar time for starting point.
     */
    var.start_time = mktime(&date_start);
  }
  else {
    sprintf(Message,"unrecognized defaults.start_date_input_type=%d",defaults.start_date_input_type);
    epic_error(dbmsname,Message);
  }
  var.model_time = defaults.start_time = var.start_time;

  /*
   * Set solar longitude, L_s [deg], which is a function of time.
   */
  L_s = solar_longitude(var.model_time);

 /*
  * Prompt for initial winds.
  *
  * u_scale - initial wind scaling factor
  * du_vert - a flag: 1 = use vert. wind profile; 0 = du is constant w/ height
  */
  switch(planet->index) {
    case EARTH_INDEX:
      defaults.u_scale = 0.;
      grid.du_vert     = defaults.du_vert = 0.;
    break;
    case HELD_SUAREZ_INDEX:
      defaults.u_scale = 0.;
      grid.du_vert     = defaults.du_vert = 0.;
    break;
    case GOLDMINE_INDEX:
      defaults.u_scale = 0.;
      grid.du_vert     = defaults.du_vert = 0.;
    break;
    case VENUS_LLR05_INDEX:
      defaults.u_scale = 0.;
      grid.du_vert     = defaults.du_vert = 0.;
    break;
    case TITAN_INDEX:
      defaults.u_scale = 0.;
      grid.du_vert     = defaults.du_vert = 0.;
    break;
    default:
      defaults.u_scale = input_double("Input initial-wind scaling factor [zero=0., full strength=1.]\n",defaults.u_scale,MODIFY);
      if (defaults.u_scale != 0.) {
        /*
         * Prompt for vertical variation of zonal wind.
         */

          if (strcmp(planet->name,"jupiter") == 0) {
              fprintf(stdout,"Vertical wind shear mode: %d => use probe profile \n", WIND_SHEAR_PROBE);
              fprintf(stdout,"                          %d => use functional form (Garcia-Melendo 2005) \n", WIND_SHEAR_FUNCTION);
              grid.wind_shear_mode = defaults.grid_wind_shear_mode = input_int("",defaults.grid_wind_shear_mode,MODIFY);
              
              fprintf(stdout,"Vertical wind decay slope: %d => constant \n", WIND_DECAY_SLOPE_CONSTANT);
              fprintf(stdout,"                           %d => varying with latitude \n", WIND_DECAY_SLOPE_VARYING);
              grid.wind_decay_m_const = defaults.grid_wind_decay_m_const = input_double("",defaults.grid_wind_decay_m_const,MODIFY);

              if(grid.wind_decay_m_const == WIND_DECAY_SLOPE_CONSTANT) {
                fprintf(stdout,"Vertical wind decay slope: \n");
                grid.du_vert_m = defaults.du_vert_m = input_double("",defaults.du_vert_m,MODIFY);
              } else if (grid.wind_decay_m_const == WIND_DECAY_SLOPE_VARYING) {
                /* maybe will need to read in a file name at some point */
              } else {
                  sprintf(Message,"%d not a valid option", grid.wind_decay_m_const);
                  epic_error(dbmsname,Message);
              }
                

              if(grid.wind_shear_mode == WIND_SHEAR_FUNCTION) {
                fprintf(stdout,"Transition pressure for zero shear [mbar]: \n");
                grid.du_vert_pend = defaults.du_vert_pend = input_double("",defaults.du_vert_pend / 100.,1) * 100.;
              } else {
                fprintf(stdout,"Vertical wind profile: 1.0 => use full probe profile \n");
                fprintf(stdout,"                       0.0 => constant with height \n");
                grid.du_vert = defaults.du_vert = input_double("",defaults.du_vert,MODIFY);
                grid.du_vert_m = defaults.du_vert_m = 0.;
                grid.du_vert_pend = defaults.du_vert_pend = 0.;
              }

              fprintf(stdout,"Cloud top pressure for wind field: %d => constant (680mb)\n", CONSTANT_CLOUD_BASE);
              fprintf(stdout,"                                   %d => Use belt/zone difference\n", USE_BELT_ZONE);
              fprintf(stdout,"                                   %d => Use I/F value\n", USE_IF_889NM);

              grid.cloud_top_mode = defaults.grid_cloud_top_mode = input_int("",defaults.grid_cloud_top_mode,MODIFY);
              if(grid.cloud_top_mode == USE_IF_889NM) {
                  fprintf(stdout, "I/F at 400 mb: \n");
                  grid.if_400mb = defaults.grid_if_400mb = input_double("", defaults.grid_if_400mb,MODIFY);
                  fprintf(stdout, "I/F at 1000 mb: \n");
                  grid.if_1000mb = defaults.grid_if_1000mb = input_double("", defaults.grid_if_1000mb,MODIFY);
              } else {
                  grid.if_400mb = grid.if_1000mb = 0.;
              }
          } else if(strcmp(planet->name,"venus")   == 0) {
              fprintf(stdout,"Vertical wind profile: 1.0 => use full probe profile \n");
              fprintf(stdout,"                       0.0 => constant with height \n");
              grid.du_vert = defaults.du_vert = input_double("",defaults.du_vert,MODIFY);
              grid.du_vert_m = defaults.du_vert_m = 0.;
              grid.du_vert_pend = defaults.du_vert_pend = 0.;
          } else if (strcmp(planet->name,"saturn") == 0) {
                /*
                 * Cassini CIRS profile: see function cassini_cirs_u().
                 */
                fprintf(stdout,"Vertical wind profile: 1.0 => use Cassini CIRS profile \n");
                fprintf(stdout,"                       0.0 => constant with height \n");
                grid.du_vert = defaults.du_vert = input_double("",defaults.du_vert,MODIFY);
          } else {
            grid.du_vert = defaults.du_vert = 0.;
            grid.du_vert_m = defaults.du_vert_m = 0.;
            grid.du_vert_pend = defaults.du_vert_pend = 0.;
          }
      }
      else {
        grid.du_vert = defaults.du_vert = 0.;
        grid.du_vert_m = defaults.du_vert_m = 0.;
        grid.du_vert_pend = defaults.du_vert_pend = 0.;
      }
    break;
  }

  /*
   * Inquire about radiation scheme.
   */
  grid.radiation_index  = defaults.radiation_index;
  strcpy(grid.radiation_scheme,radiation_scheme[defaults.radiation_index]);
  grid.newt_cool_adjust = defaults.newt_cool_adjust;
  grid.zonal_average_rt = defaults.zonal_average_rt;

  inquire_radiation_scheme(MODIFY);

  if (strcmp(grid.radiation_scheme, "Global heating-cooling") == 0) {
    defaults.heat_top_pressure = grid.heat_top_pressure = input_double("Top of the heating region [hPa]\n", defaults.heat_top_pressure,MODIFY);
    defaults.cool_bot_pressure = grid.cool_bot_pressure = input_double("Bottom of the cooling region [hPa]\n", defaults.cool_bot_pressure,MODIFY);
    defaults.heat_rate = grid.heat_rate = input_double("Heating rate at the bottom [K/day]\n", defaults.heat_rate,MODIFY);
    defaults.cool_rate = grid.cool_rate = input_double("Cooling rate at the top [K/day]\n", defaults.cool_rate,MODIFY);
    
    // convert to appropriate units
    grid.heat_top_pressure = grid.heat_top_pressure * 100.;
    grid.cool_bot_pressure = grid.cool_bot_pressure * 100.;
    grid.heat_rate = grid.heat_rate / 86400.;
    grid.cool_rate = grid.cool_rate / 86400.;
  }

  defaults.radiation_index  = grid.radiation_index;
  defaults.newt_cool_adjust = grid.newt_cool_adjust;
  defaults.zonal_average_rt = grid.zonal_average_rt;
   
  /* 
   * Set velocity timestep scheme to 3rd-order Adams-Bashforth.
   */
  defaults.uv_timestep_scheme = 0;
  strcpy(grid.uv_timestep_scheme,uv_timestep_scheme[defaults.uv_timestep_scheme]);

  /* 
   * Prompt for turbulence scheme.
   */
  fprintf(stdout,"\n");
  fprintf(stdout,"Turbulence schemes:\n");
  for (ii = 0; ii < num_turbulence_schemes; ii++) {
    fprintf(stdout,"  %s\n",turbulence_scheme[ii]);
  }
  input_string("Input turbulence scheme\n",defaults.turbulence_scheme,grid.turbulence_scheme,MODIFY);

  /*
   * Grab the vertical u profile, if it exists.
   * Malloc the string space, and name the indices, and units for each 
   * prognostic and diagnostic into var.(var_name).info[0].index,name,units.
   * Also set enthalpy change and saturation vapor pressure for species.
   * Also, set molar mass and triple-point values for species.
   *
   * NOTE: Need to set grid.uv_timestep_scheme before calling set_var_props().
   */
  set_var_props();

  /* 
   * Set advection scheme for THETA, FPARA.
   *
   * NOTE: Avoid using HSU_ADVECTION for THETA, it tends to be numerically unstable, for reasons unknown.
   *       THIRD_ORDER_UPWIND has proven reliable for THETA.
   *
   * NOTE: FPARA is like a Q, so HSU_ADVECTION is recommended.
   */
  strcpy(var.theta.advection_scheme,advection_scheme[THIRD_ORDER_UPWIND]);
  strcpy(var.fpara.advection_scheme,advection_scheme[HSU_ADVECTION]);

  /*
   * Set advection schemes for mass variables (H, Qs).
   *
   * NOTE: These need to be flux-form (h-weighted), e.g. HSU_ADVECTION but not THIRD_ORDER_UPWIND.
   * NOTE: This must come before make_arrays().
   */
  strcpy(var.h.advection_scheme,advection_scheme[HSU_ADVECTION]);
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    strcpy(var.species[is].advection_scheme,advection_scheme[HSU_ADVECTION]);
  }

  /*
   * Set advection scheme for NU_TURB.
   * THIRD_ORDER_UPWIND is recommended. 
   */
  strcpy(var.nu_turb.advection_scheme,advection_scheme[THIRD_ORDER_UPWIND]);

  /*
   * Set equation of state to ideal.
   * NOTE: We have code available in the model for handling the virial equation of state,
   *       but we have not run across any case in the solar system that benefits signficantly from it,
   *       so we hardwire the model to use the ideal eos.
   */
  sprintf(grid.eos,"ideal");
  sprintf(defaults.eos,"ideal");

  /*
   * grid.k_sponge: number of top sponge layers
   */
  if (strcmp(planet->name,"Held_Suarez") == 0) {
    defaults.k_sponge = grid.k_sponge = -1;
  }
  else if (defaults.k_sponge == NOT_SET) {
    defaults.k_sponge     = -1;
    defaults.k_sponge     = grid.k_sponge = input_int("Input k_sponge (-1 = no effect)\n",defaults.k_sponge,MODIFY);
    defaults.spacing_type = SPACING_LOGP;
  }
  else {
    defaults.k_sponge = grid.k_sponge = input_int("Input k_sponge (-1 = no effect)\n",defaults.k_sponge,MODIFY);
  }

  /*
   * grid.n_bot_drag - number of layers at bottom with transitional Rayleigh drag towards abyssal wind profile.
   */
  if (strcmp(planet->type,"terrestrial") == 0) {
    defaults.n_bot_drag = grid.n_bot_drag = -1;
  }
  else if (defaults.n_bot_drag == NOT_SET) {
    defaults.n_bot_drag = -1;
    defaults.n_bot_drag = grid.n_bot_drag = input_int("Input number of bottom layers with transitional drag towards abyssal wind profile (-1 = no effect,MODIFY):\n",
                                                      defaults.n_bot_drag,MODIFY);
  }
  else {
    defaults.n_bot_drag = grid.n_bot_drag = input_int("Input number of bottom layers with transitional drag towards abyssal wind profile (-1 = no effect,MODIFY):\n",
                                                      defaults.n_bot_drag,MODIFY);
  }

  /*
   * grid.j_sponge: number of lateral sponge layers counting inwards from the southern and northern edges
   */
  if (defaults.j_sponge == NOT_SET) {
    defaults.j_sponge     = -1;
    defaults.j_sponge     = grid.j_sponge = input_int("Input j_sponge (-1 = no effect; 3 is typical)\n",defaults.j_sponge,MODIFY);
  }
  else {
    defaults.j_sponge = grid.j_sponge = input_int("Input j_sponge (-1 = no effect; 3 is typical)\n",defaults.j_sponge,MODIFY);
  }

  /*
   * Inquire about vertical spacing between layers.
   *
   * NOTE: The option SPACING_LOGP_W_BOUNDARIES has been removed, because the variable spacing
   *       was found to generate noise in the vertical dimension.
   */
  fprintf(stdout,"\n");
  sprintf(header,"Layer spacing initially even in: %1d => log pressure\n"
                 "                                 %1d => pressure \n"
                 "                                 %1d => theta\n"
                 "                                 %1d => from file \n",
                 SPACING_LOGP,SPACING_P,SPACING_THETA,SPACING_FROM_FILE);
  defaults.spacing_type = spacing_type = input_int(header,defaults.spacing_type,MODIFY);

  /*
   * Inquire about the vertical range and number of vertical layers.
   *
   * Grab the file, if the user wants to input one.
   * read_spacing_file() is a simple template for reading ascii data files.
   * grid.ptop = p at k=1/2 [Pa, internally]
   * grid.pbot = p at k=nk+1/2 [Pa, internally]
   * grid.nk = user's number of vert. layers
   */
  switch(defaults.spacing_type) {
    case SPACING_FROM_FILE:
      input_string("File containing layer spacing data (type 'initial -h' to see an example,MODIFY)\n",
                   defaults.layer_spacing_dat,defaults.layer_spacing_dat,MODIFY);
      read_spacing_file(&defaults,SIZE_DATA);
    break;
    case SPACING_P:
    case SPACING_LOGP:
      /* 
       * Prompt for ptop (external units are hPa, internal are Pa):
       */
      defaults.ptop = grid.ptop = 100.*input_double("Pressure at K = 1/2 (model's top,MODIFY) [hPa]\n",
                                                    defaults.ptop/100.,MODIFY);

      if (strcmp(planet->type,"gas-giant") == 0) {
        /*
         * Prompt for pbot for gas giant.
         */
        defaults.pbot = grid.pbot = 100.*
                                    input_double("Pressure at K = nk+1/2 (model's bottom,MODIFY) [hPa]\n",
                                                 defaults.pbot/100.,MODIFY);
      }
      else {
        defaults.pbot = grid.pbot = planet->p0;
      }

      if (defaults.nk == NOT_SET) {
        /* 
         * Prompt nk ~ 2*(model range in scale heights).
         */
        if (strcmp(planet->type,"gas-giant") == 0) {
          defaults.nk = 2*NINT(log(defaults.pbot/defaults.ptop));
        }
        else {
          defaults.nk = 20;
        }
      }
      defaults.nk = grid.nk = input_int("\nInput the number of vertical layers, nk\n",defaults.nk,MODIFY);
    break;
    case SPACING_THETA:
      defaults.thetatop = grid.thetatop = input_double("Potential temperature at K = 1/2 (model's top,MODIFY) [K]\n",
                                                       defaults.thetatop,MODIFY);
      defaults.thetabot = grid.thetabot = input_double("Potential temperature at K = nk+1/2 (model's bottom,MODIFY) [K]\n",
                                                       defaults.thetabot,MODIFY);
      if (defaults.nk == NOT_SET) {
        defaults.nk = 20;
      }
      defaults.nk = grid.nk = input_int("\nInput the number of vertical layers, nk\n",defaults.nk,MODIFY);
    break;
    default:
      sprintf(Message,"unrecognized defaults.spacing_type=%d",defaults.spacing_type);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * Choose geometry.
   *
   *           grid.geometry: name of geometry
   *             grid.wrap[]: periodicity flags for each dim. (globe: 0,0,1)
   *              grid.pad[]: boundary pad widths for each dim. (globe: 1,1,1)
   *                grid.jlo: low-end index for j (globe: jlo = 0)  
   *                grid.ilo: low-end index for i (typically 1)
   *       grid.globe_latbot: lowest latitude value
   *       grid.globe_lattop: highest latitude value
   *       grid.globe_lonbot: lowest longitude value
   *       grid.globe_lontop: highest longitude value
   *       grid.f_plane_lat0: central latitude of f_plane
   * grid.f_plane_half_width: half of the f_plane width
   */
  if (strcmp(system_id,"held_suarez") == 0) {
    sprintf(grid.geometry,    "globe");
    sprintf(defaults.geometry,"globe");
  }
  else if (strcmp(system_id,"venus_llr05") == 0) {
    sprintf(grid.geometry,    "globe");
    sprintf(defaults.geometry,"globe");
  }
  else {
    input_string("Choose geometry \n",defaults.geometry,grid.geometry,MODIFY);
  }

  if (strcmp(grid.geometry,"f-plane") == 0) {
    grid.wrap[2] = 0;
    grid.wrap[0] = 1;
    grid.pad[2]  = 1;
    grid.pad[0]  = 1;
    grid.f_plane_lat0 = defaults.f_plane_lat0
                      = input_double("Latitude of f-plane [deg] \n",
                                     defaults.f_plane_lat0,MODIFY);
    input_string("Choose mapping: cartesian or polar \n",
                 defaults.f_plane_map,grid.f_plane_map,MODIFY);
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      grid.wrap[1] = 1;
      grid.pad[1]  = 1;
      grid.jlo     = 1;
      grid.ilo     = 1;
      grid.f_plane_half_width = defaults.f_plane_half_width
                              = input_double("Half-width [km] \n",
                                             defaults.f_plane_half_width,MODIFY);
      /* convert km to m */
      grid.f_plane_half_width *= 1000.;
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      grid.wrap[1] = 0;
      grid.pad[1]  = 1;
      grid.jlo     = 0;
      grid.ilo     = 1;
      grid.f_plane_half_width = defaults.f_plane_half_width
                              = input_double("Radius [km] \n",
                                            defaults.f_plane_half_width,MODIFY);
      /* convert km to m */
      grid.f_plane_half_width *= 1000.;
    }
    else {
      fprintf(stderr,"Unrecognized f-plane mapping. \n");
      exit(1);
    }
    grid.include_nontrad_accel = FALSE;
  }
  else if (strcmp(grid.geometry,"globe") == 0) {
    grid.wrap[2] = 0;
    grid.wrap[1] = 0;
    grid.wrap[0] = 1;
    grid.pad[2]  = 1;
    grid.pad[1]  = 1;
    grid.pad[0]  = 1;
    grid.jlo     = 0;
    grid.ilo     = 1;

    grid.globe_latbot=defaults.globe_latbot
      = input_double("Lowest latitude [deg] \n",defaults.globe_latbot,MODIFY);
    if (grid.globe_latbot < -90.) {
      fprintf(stderr,"latbot must be >= -90.  Setting latbot to -90. \n");
      grid.globe_latbot = defaults.globe_latbot = -90.;
    }
    grid.globe_lattop=defaults.globe_lattop
      = input_double("Highest latitude [deg] \n",defaults.globe_lattop,MODIFY);
    /* Sanity check: */
    if (grid.globe_lattop < grid.globe_latbot) {
      sprintf(Message,"lattop=%f < latbot=%f",grid.globe_lattop,grid.globe_latbot);
      epic_error(dbmsname,Message);
    }
    if (grid.globe_lattop > 90.) {
      fprintf(stderr,"lattop must be <= 90. Setting lattop to 90. \n");
      grid.globe_lattop = defaults.globe_lattop = 90.;
    }

    grid.globe_lonbot=defaults.globe_lonbot
      = input_double("Lowest longitude [deg] \n",defaults.globe_lonbot,MODIFY);
    grid.globe_lontop=defaults.globe_lontop
      = input_double("Highest longitude [deg] \n",defaults.globe_lontop,MODIFY);
    /* Sanity check: */
    if (grid.globe_lontop < grid.globe_lonbot) {
      sprintf(Message,"lontop=%f < lonbot=%f",grid.globe_lontop,grid.globe_lonbot);
      epic_error(dbmsname,Message);
    }
    grid.include_nontrad_accel = TRUE;
  }
  else {
    fprintf(stderr,"initial: unrecognized geometry: %s \n",grid.geometry);
    exit(1);
  }

  /*
   * Set up thermodynamics subroutines.
   *
   * This fills in the 'thermo' struct with planet-relevant values.
   * It returns planet->cpr, which is the planet's reference cp.
   * The (->) operator allows us to reference the member of a pointed-to struct.
   * Here, we are passing the address of planet->cpr, so that its value can be set
   * by the subroutine. 
   * 
   *   planet->cpr: nondimensional ref. cp
   *    planet->cp: spec. heat at const. p
   *  planet->rgas: gas constant
   * planet->kappa: rgas/cp
   */
  thermo_setup(&planet->cpr);
  /* Assign thermodynamics function's reference cpr to planet->cp */
  planet->cp    = planet->cpr*planet->rgas;
  planet->kappa = 1./planet->cpr;

  /*
   * Turn on parameters that are present in all systems.
   */
  var.gravity2.on = TRUE;
  var.on_list[GRAVITY2_INDEX] = NOT_LISTED;

  /* 
   * Turn on core prognostic variables, u,v,h,theta.
   *
   * var.(prog_name).on: a boolean flag to tell us it is invoked or not
   * var.on_list[INDEX]: an array of tristate flags (NOT_LISTED,LISTED_AND_OFF,LISTED_AND_ON)
   */
  var.u.on     = TRUE;
  var.v.on     = TRUE; 
  var.h.on     = TRUE;
  var.theta.on = TRUE;
  var.on_list[    U_INDEX] = LISTED_AND_ON;
  var.on_list[    V_INDEX] = LISTED_AND_ON;
  var.on_list[    H_INDEX] = LISTED_AND_ON;
  var.on_list[THETA_INDEX] = LISTED_AND_ON;

  /*
   * Turn on turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
      strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
    var.nu_turb.on = TRUE;
    var.on_list[NU_TURB_INDEX] = LISTED_AND_ON;
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  /*
   * Turn on diagnostic variables that are used
   * in the model.  Variables not turned on here can
   * still be written to extract.nc.
   */
  var.h3.on         = TRUE;
  var.hdry2.on      = TRUE;
  var.hdry3.on      = TRUE;
  var.p2.on         = TRUE;
  var.p3.on         = TRUE;
  var.theta2.on     = TRUE;
  var.t2.on         = TRUE;
  var.t3.on         = TRUE;
  var.rho2.on       = TRUE;
  var.rho3.on       = TRUE;
  var.exner2.on     = TRUE;
  var.exner3.on     = TRUE;
  var.phi2.on       = TRUE;
  var.phi3.on       = TRUE;
  var.mont2.on      = TRUE;
  var.heat3.on      = TRUE;
  var.pv2.on        = TRUE;
  var.w3.on         = TRUE;
  var.z2.on         = TRUE;
  var.dzdt2.on      = TRUE;

  var.on_list[H3_INDEX      ] = LISTED_AND_ON;
  var.on_list[HDRY2_INDEX   ] = LISTED_AND_ON;
  var.on_list[HDRY3_INDEX   ] = LISTED_AND_ON;
  var.on_list[P2_INDEX      ] = LISTED_AND_ON;
  var.on_list[P3_INDEX      ] = LISTED_AND_ON;
  var.on_list[THETA2_INDEX  ] = LISTED_AND_ON;
  var.on_list[T2_INDEX      ] = LISTED_AND_ON;
  var.on_list[T3_INDEX      ] = LISTED_AND_ON;
  var.on_list[RHO2_INDEX    ] = LISTED_AND_ON;
  var.on_list[RHO3_INDEX    ] = LISTED_AND_ON;
  var.on_list[EXNER2_INDEX  ] = LISTED_AND_ON;
  var.on_list[EXNER3_INDEX  ] = LISTED_AND_ON;
  var.on_list[PHI2_INDEX    ] = LISTED_AND_ON;
  var.on_list[PHI3_INDEX    ] = LISTED_AND_ON;
  var.on_list[MONT2_INDEX   ] = LISTED_AND_ON;
  var.on_list[HEAT3_INDEX   ] = LISTED_AND_ON;
  var.on_list[PV2_INDEX     ] = LISTED_AND_ON;
  var.on_list[W3_INDEX      ] = LISTED_AND_ON;
  var.on_list[Z2_INDEX      ] = LISTED_AND_ON;
  var.on_list[DZDT2_INDEX   ] = LISTED_AND_ON;

  if (strcmp(planet->type,"gas-giant") == 0) {
    var.pbot.on             = TRUE;
    var.on_list[PBOT_INDEX] = NOT_LISTED;

    var.phi_surface.on             = FALSE;
    var.on_list[PHI_SURFACE_INDEX] = NOT_LISTED;
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    var.pbot.on             = FALSE;
    var.on_list[PBOT_INDEX] = NOT_LISTED;

    var.phi_surface.on             = TRUE;
    var.on_list[PHI_SURFACE_INDEX] = LISTED_AND_ON;
  }
  else {
    sprintf(Message,"planet->type=%s not yet implemented\n",planet->type);
    epic_error(dbmsname,Message);
  }

  fprintf(stdout,"\n Core prognostic variables: \n");

  fprintf(stdout,"  %2d  ><  %-s \n",U_INDEX,    var.u.info[    0].name);
  fprintf(stdout,"  %2d  ><  %-s \n",V_INDEX,    var.v.info[    0].name);
  fprintf(stdout,"  %2d  ><  %-s \n",H_INDEX,    var.h.info[    0].name);
  fprintf(stdout,"  %2d  ><  %-s \n",THETA_INDEX,var.theta.info[0].name); 

  /*
   * Inquire about optional prognostic variables.
   *
   * prompt_species_on() : a function to inquire about the species, and turn 
   *                        them on, as above.
   */
  if (strcmp(system_id,"held_suarez") == 0) {
    ;
  }
  else if (strcmp(system_id,"goldmine") == 0) {
    ;
  }
  else {
    prompt_species_on(defaults.species_str,MODIFY);
  }

  if (var.on_list[H_2O_INDEX] == LISTED_AND_ON ||
      var.on_list[NH_3_INDEX] == LISTED_AND_ON ||
      var.on_list[CH_4_INDEX] == LISTED_AND_ON   ) {
    if (defaults.cloud_microphysics == NOT_SET ||
        defaults.cloud_microphysics == OFF) {
      defaults.cloud_microphysics = ACTIVE;
    }
    sprintf(Message,"Cloud microphysics: %2d => off \n"
                    "                    %2d => active  (latent heat, phase changes, precipitation) \n"
                    "                    %2d => passive (advection only)\n",
                    OFF,ACTIVE,PASSIVE);
    defaults.cloud_microphysics = grid.cloud_microphysics = input_int(Message,defaults.cloud_microphysics,MODIFY);
  }
  else {
    defaults.cloud_microphysics = grid.cloud_microphysics = OFF;
  }

  /*
   * Turn on phases appropriate to cloud microphysics package.
   */
  turn_on_phases();

  /*
   * Prompt for moist convection.
   */
  if (grid.cloud_microphysics != OFF) {
    if (defaults.moist_convection == NOT_SET ||
        defaults.moist_convection == OFF) {
      defaults.moist_convection = OFF;
    }
    sprintf(Message,"Moist convection scheme: %2d => off \n"
                    "                         %2d => on (Sud & Walker, 1999) \n"
                    "                         %2d => passive (calculate diagnostic variables only) \n",
                    OFF,ACTIVE,PASSIVE);
    defaults.moist_convection = grid.moist_convection = input_int(Message,defaults.moist_convection,MODIFY);

    if (grid.moist_convection == ACTIVE) {
      sprintf(Message,"Maximum number of iterations: \n");
      defaults.max_mc_it = grid.max_mc_it = input_int(Message,defaults.max_mc_it,MODIFY);

      sprintf(Message,"Relaxation timescale (sec): \n");
      defaults.tau_relax = grid.tau_relax = input_double(Message,defaults.tau_relax,MODIFY);

      /*
       * Set the status to TRUE that relaxed Arakawa-Schubert (RAS)
       * will be calculated for the first time.
       */
      grid.first_RAS_upd = TRUE;
    }
  }
  else {
    defaults.moist_convection = grid.moist_convection = OFF;
    defaults.max_mc_it        = grid.max_mc_it        = 0;
    defaults.tau_relax        = grid.tau_relax        = 0.;
  }
  
  if(grid.cloud_microphysics != OFF) {
    if (defaults.grid_relax_vapor == NOT_SET ||
        defaults.grid_relax_vapor == OFF) {
      defaults.grid_relax_vapor = OFF;
    }
    sprintf(Message,"Relax vapor profile to initial state : %2d => on, or \n"
                    "                                       %2d => off\n", ON, OFF);
    defaults.grid_relax_vapor = grid.relax_vapor = input_int(Message,defaults.grid_relax_vapor,1);

    if(grid.relax_vapor == ACTIVE) {
        sprintf(Message,"Timescale [days] : \n");
        defaults.grid_relax_vapor_timescale = grid.relax_vapor_timescale = input_int(Message,defaults.grid_relax_vapor_timescale,1);
        grid.relax_vapor_timescale *= 24. * 3600.;
    }
  } else {
      defaults.grid_relax_vapor = OFF;
      grid.relax_vapor_timescale *= 0.;
  }

  /* 
   * Source-sink parameters.
   */
  if (var.fpara.on) {
    defaults.fpara_rate_scaling = var.fpara_rate_scaling = 
      input_double("Ortho-para conversion-rate scaling [nominal = 1.0]?\n",defaults.fpara_rate_scaling,MODIFY);
  }

  /*
   * Inquire about number of latitude grid points.
   *
   * NOTE: nj equals the number of interior v or q points, counting from 1 to nj,
   *       whereas nj+1 equals the number of u or h points, counting from 0 to nj
   *       for channel models and globes. 
   *       One consequence is that for a 1D model, nj = 0.
   *
   * grid.jlo   : low-end index for j (globe: jlo = 0)
   * grid.ilo   : low-end index for i (typically 1)
   * grid.nj    : number of lat. grid points
   * grid.dlt   : lat. grid spacing (deg)
   */
  if (strcmp(grid.geometry,"globe") == 0) {
    sprintf(Message,"\nLatitude indexing: nj   = No. V or PV points interior to the boundaries, from 1 to nj, whereas"
                    "\n                   nj+1 = No. U or H points (they have no boundary pts), from 0 to nj."
                    "\nInput the number of latitude gridpoints (nj, 0 => one u,h point)\n");
    defaults.nj = grid.nj = input_int(Message,defaults.nj,MODIFY);
    if (grid.globe_latbot == -90. &&
        grid.globe_lattop ==  90.) {
      /* 
       * Poles are offset by extra sqrt(.5)*dlt.
       * NOTE: This is more accurate than the .5*dlt prescribed by Arakawa and Lamb.
       */
      grid.dlt = 180./((grid.nj+1)+sqrt(.5)+sqrt(.5));
    }
    else if (grid.globe_latbot == -90.) {
      grid.dlt = (grid.globe_lattop+90.)/((grid.nj+1)+sqrt(.5));
    }  
    else if (grid.globe_lattop ==  90.) {
      grid.dlt = (90.-grid.globe_latbot)/((grid.nj+1)+sqrt(.5));
    }
    else {
      grid.dlt = (grid.globe_lattop-grid.globe_latbot)/(grid.nj+1);
    }
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      defaults.nj = grid.nj =
        input_int("\nInput the number of gridpoints on a side\n",defaults.nj,MODIFY);
      grid.dlt = 180./(grid.nj);
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      /* 
       * "latitude" runs from 0 at the edge (j = 1) 
       *  to 90 in the center (j = nj+1) */
      defaults.nj = grid.nj =
        input_int("\nInput the number of radial gridpoints\n",defaults.nj,MODIFY);
      grid.dlt   = 90./((grid.nj+1-1)+sqrt(.5));
    }
  }
  else {
    sprintf(Message,"unrecognized grid.geometry=%s",grid.geometry);
    epic_error(dbmsname,Message);
  }
  /*
   * To handle staggered C-grid.
   */
  if (strcmp(grid.geometry,   "f-plane")   == 0 &&
      strcmp(grid.f_plane_map,"cartesian") == 0   ) {
    /* Periodic in J direction */
    grid.jlast  = grid.nj;
    grid.jfirst = grid.jlo;
  }
  else {
    /* Not periodic in J direction */
    grid.jlast  = grid.nj -1;
    grid.jfirst = grid.jlo+1;
  }

  /*
   * Inquire about number of longitude gridpoints.
   *
   * frexp() : a standard c math.h function.
   * grid.ni : number of long. grid points (must be power of 2)
   * grid.dln: long. grid spacing (deg)
   */
  if (strcmp(grid.geometry,"globe") == 0) {
    if (IAMNODE == NODE0) {
      fprintf(stdout,"\nInput the number of longitude gridpoints (ni), \n");
    }
    defaults.ni = grid.ni = 
      input_int("which must be an integer power of two\n",defaults.ni,MODIFY);
    /*
     * Verify that ni is a power of 2.
     */
    if (frexp((double)grid.ni,&itmp) != 0.5) {
      sprintf(Message,"ni=%d is not a power of 2",grid.ni);
      epic_warning(dbmsname,Message);
    }
    grid.dln = (grid.globe_lontop-grid.globe_lonbot)/(grid.ni);
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      defaults.ni = grid.ni = grid.nj;
      grid.dln    = grid.dlt;
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      if (IAMNODE == NODE0) {
        fprintf(stdout,"\nInput the number of longitude gridpoints, \n");
      }
      defaults.ni = grid.ni = 
        input_int("which must be an integer power of two\n",defaults.ni,MODIFY);
      /*
       * Verify that ni is a power of 2.
       */
      if (frexp((double)grid.ni,&itmp) != 0.5) {
        sprintf(Message,"ni=%d is not a power of 2",grid.ni);
        epic_warning(dbmsname,Message);
      }
      grid.dln = 360./(grid.ni);
    }
  }

  /*
   * Input dt.
   */
  defaults.dt = grid.dt = 
    input_int("Input timestep, dt [s, int to prevent roundoff]\n",defaults.dt,MODIFY);

  var.ntp = read_t_vs_p(SIZE_DATA);

  /* * * * * * * * * * * * * * * * * * *
   *                                   *
   *  Allocate memory for variables.   *
   *                                   *
   * * * * * * * * * * * * * * * * * * */

  make_arrays();

  for (I = 0; I < NUM_WORKING_BUFFERS; I++) {
    Buff2D[I] = dvector(0,Nelem2d-1,dbmsname);
  }

  /*
   * NOTE: these must come after make_arrays().
   */
  read_t_vs_p(POST_SIZE_DATA);

  if (spacing_type == SPACING_FROM_FILE) {
    read_spacing_file(&defaults,ALL_DATA);
  }

  if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
    /*
     * Flag DUDT and DVDT for timeplanes IT_MINUS2 and IT_MINUS1 with DBL_MAX
     * to indicate this is the initial timestep.
     */
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          DUDT(IT_MINUS2,K,J,I) = DBL_MAX;
          DUDT(IT_MINUS1,K,J,I) = DBL_MAX;
          DVDT(IT_MINUS2,K,J,I) = DBL_MAX;
          DVDT(IT_MINUS1,K,J,I) = DBL_MAX;
        }
      }
    }
  }
  else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
    /*
     * Flag U and V for timeplane IT_MINUS1 with DBL_MAX
     * to indicate this is the initial timestep.
     */
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          U(IT_MINUS1,K,J,I) = DBL_MAX;
          V(IT_MINUS1,K,J,I) = DBL_MAX;
        }
      }
    }
  }
  else {
    sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
    epic_error(dbmsname,Message);
  }

  /*
   * NOTE: timeplane_bookkeeping() should come after setting the initial-step
   *       DBL_MAX flags in U or DUDT.
   */
  timeplane_bookkeeping();

  /*
   * Set lon, lat.
   * set_lonlat():
   * this propogates grid.lon, grid.lat (the dim arrays) with proper values 
   * according to ni, nj, and latbot, lattop, dlt, dln, etc.
   */
  set_lonlat();

  /* Call set_dsgth() below after sigmatheta is defined. */

  /* 
   * Allocate memory: 
   */
  thetadat     = dvector( 0,var.ntp-1,dbmsname);
  neglogpdat   = dvector( 0,var.ntp-1,dbmsname);
  buff_triplet = dtriplet(0,var.ntp-1,dbmsname);

  /*
   * Interpolate on -log p:
   */
  for (ki = 0; ki < var.ntp; ki++) {
    neglogpdat[ki] = -log(var.pdat[ki]);
  }

  for (ki = 0; ki < var.ntp; ki++) {
    fpara        = return_fpe(var.tdat[ki]);
    thetadat[ki] = return_theta(fpara,var.pdat[ki],var.tdat[ki],&theta_ortho,&theta_para);
  }

  /*
   * Top and bottom reference values:
   */
  switch (spacing_type) {
    case SPACING_P:
    case SPACING_LOGP:
    case SPACING_FROM_FILE:
      grid.ptop = grid.p_ref[      1] = defaults.ptop;
      grid.pbot = grid.p_ref[2*KHI+1] = defaults.pbot;

      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = -log(var.pdat[ki]);
        buff_triplet[ki].y = thetadat[ki];
      }
      spline_pchip(var.ntp,buff_triplet);

      neglogp           = -log(defaults.ptop);
      ki                = find_place_in_table(var.ntp,buff_triplet,neglogp,&p_d);
      defaults.thetatop = grid.thetatop = grid.theta_ref[1] = splint_pchip(neglogp,buff_triplet+ki,p_d);

      neglogp                 = -log(defaults.pbot);
      ki                      = find_place_in_table(var.ntp,buff_triplet,neglogp,&p_d);
      defaults.thetabot = grid.thetabot = grid.theta_ref[2*KHI+1] = splint_pchip(neglogp,buff_triplet+ki,p_d);
    break;
    case SPACING_THETA:
      grid.thetatop = grid.theta_ref[      1] = defaults.thetatop;
      grid.thetabot = grid.theta_ref[2*KHI+1] = defaults.thetabot;

      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = thetadat[ki];
        buff_triplet[ki].y = -log(var.pdat[ki]);
      }
      spline_pchip(var.ntp,buff_triplet);

      tmp           = defaults.thetatop; 
      ki            = find_place_in_table(var.ntp,buff_triplet,tmp,&sgth_d);
      defaults.ptop = grid.ptop = grid.p_ref[1] = exp(-splint_pchip(tmp,buff_triplet+ki,sgth_d));

      tmp                 = defaults.thetabot;
      ki                  = find_place_in_table(var.ntp,buff_triplet,tmp,&sgth_d);
      defaults.pbot = grid.pbot = grid.p_ref[2*KHI+1] = exp(-splint_pchip(tmp,buff_triplet+ki,sgth_d));
    break;
    default:
      sprintf(Message,"unrecognized spacing_type=%d",spacing_type);
      epic_error(dbmsname,Message);
    break;
  }

  /************************************************
   *                                              *
   * Set the pressure at the bottom of the model. *
   *                                              *
   ************************************************/

  /*
   * Set the J index (latitude) to apply T(p) sounding data.
   */

  /* Limit to valid range */
  defaults.lat_tp = LIMIT_RANGE(grid.globe_latbot,defaults.lat_tp,grid.globe_lattop);

  if (strcmp(planet->name,"Venus") == 0) {
    defaults.lat_tp = 4.;
  }
  else if (strcmp(planet->name,"Jupiter") == 0) {
    defaults.lat_tp = input_double("Latitude to apply T(p) sounding data [deg]\n",defaults.lat_tp,MODIFY);
  }
  else {
    defaults.lat_tp = 0.;
  }
  /* Limit to valid range (again) */
  defaults.lat_tp = LIMIT_RANGE(grid.globe_latbot,defaults.lat_tp,grid.globe_lattop);

  grid.jtp = 0;
  for (J = JLO; J <= JHI; J++) {
    if (fabs(grid.lat[2*J+1]-defaults.lat_tp) <= .5*grid.dlt) {
      grid.jtp = J;
      break;
    }
  }

  if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * For gas giants let the bottom of the model start out as a constant pressure surface.
     *
     * NOTE: The array PBOT(J,I) could be used to specify a horizontally varying bottom pressure
     *       if so desired. 
     */
    K = grid.nk;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        P3(K,J,I) = PBOT(J,I) = grid.pbot;
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    /*
     * Set the surface geopotential.
     */
    init_phi_surface();

    if (grid.coord_type == COORD_ISENTROPIC) {
      /*
       * NOTE: a horizontally varying bottom pressure is not yet implemented here. 
       */
      K = grid.nk;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          P3(K,J,I) = grid.pbot;
        }
      }
      /* No need to apply bc_lateral() here. */
    }
    else {
      /*
       * Set P3 for K=nk.
       * Reset defaults.pbot, grid.pbot to maximum bottom pressure.
       */
      K = grid.nk;
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          P3(K,J,I) = p_phi(J,PHI_SURFACE(J,I));
        }
      }
      /* No need to apply bc_lateral() here. */
    }

    defaults.pbot = -DBL_MAX;
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        defaults.pbot = MAX(P3(K,J,I),defaults.pbot);
      }
    }
    grid.pbot = defaults.pbot;
    /* No need to apply bc_lateral() here. */
  }
  else {
    sprintf(Message,"unrecognized planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }

  switch(grid.coord_type) {
    case COORD_ISOBARIC:
      defaults.p_sigma = 0.;
    break;
    case COORD_ISENTROPIC:
      defaults.p_sigma = FLT_MAX;
    break;
    case COORD_HYBRID:
      /*
       * Prompt for the target pressure for the transition region
       * between sigma and theta coordinates.
       */
      if (defaults.p_sigma == (double)NOT_SET) {
        switch(planet->index) {
          case VENUS_INDEX:
          case VENUS_LLR05_INDEX:
            defaults.p_sigma = 500.*100.;
          break;
          case EARTH_INDEX:
          case HELD_SUAREZ_INDEX:
            defaults.p_sigma = 900.*100.;
          break;
          case GOLDMINE_INDEX:
            defaults.p_sigma = 300.*100.;
          break;
          case MARS_INDEX:
            defaults.p_sigma = 0.2*100.;
          break;
          case TITAN_INDEX:
            defaults.p_sigma = 100.*100.;
          break;
          case URANUS_INDEX:
            defaults.p_sigma = 3000.*100.;
          break;
          default:
            defaults.p_sigma = sqrt(defaults.pbot*defaults.ptop);
          break;
        }
      }
      defaults.p_sigma = 100.*input_double("Pressure of transition between sigma and theta vertical coordinate [hPa]\n",
                                           defaults.p_sigma/100.,MODIFY);
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * Determine t_vs_p data indices corresponding to pbot, ptop.
   */
  floor_tp = 0;
  tmp      = DBL_MAX;
  for (ki = 0; ki < var.ntp; ki++) {
    tmp2 = fabs(var.pdat[ki]-defaults.pbot);
    /* 
     * Need var.pdat[floor_tp] <= pbot.
     */
    if (tmp2 <= tmp && var.pdat[ki] <= defaults.pbot) {
      floor_tp = ki;
      tmp      = tmp2;
    }
  }
  ceiling_tp = var.ntp-1;
  tmp        = DBL_MAX;
  for (ki = var.ntp-1; ki >= 0; ki--) {
    tmp2 = fabs(var.pdat[ki]-defaults.ptop);
    if (tmp2 < tmp && var.pdat[ki] <= defaults.ptop) {
      ceiling_tp = ki;
      tmp        = tmp2;
    }
  }
  
  /*
   * Set grid.zeta0, grid.zeta1, grid.hybrid_alpha
   *
   * NOTE: In the hybrid-coordinate case, these values may need to be adjusted to yield a
   *       hybrid coordinate, zeta, that is monotonic and reasonable.
   *
   * These parameters control the behavior of the function f_sigma().
   *
   * The value of the coordinate at the bottom of the model atmosphere, sigma = 0., is grid.zeta0.
   * This is a free parameter that should be chosen such that zeta is monotonic and smooth.
   */
  switch(grid.coord_type) {
    case COORD_ISOBARIC:
    case COORD_ISENTROPIC:
      grid.zeta0 = grid.theta_ref[2*KHI+1];
      grid.zeta1 = grid.theta_ref[2*0  +1];
    break;
    case COORD_HYBRID:
      switch(planet->index) {
        case VENUS_INDEX:
          /*
           * Grace Lee 3/22/05
           *
           * Recommended p_sigma = 500 hPa (for pbot > 500 hPa).
           */
          grid.zeta0        = 420.-53.333*log10(defaults.pbot/100.);
          grid.zeta1        = 420.-53.333*log10(defaults.ptop/100.);
          grid.hybrid_alpha = 20.;
        break;
        case VENUS_LLR05_INDEX:
          /*
           * Aaron Herrnstein 6/26/06
           */
          grid.zeta0        = 165.;
          grid.zeta1        = 480.;
          grid.hybrid_alpha = 20.;
        break;
        case EARTH_INDEX:
        case HELD_SUAREZ_INDEX:
          /*
           * Recommended p_sigma = 900 hPa.
           */
          grid.zeta0        = 260.;
          grid.zeta1        = 700.;
          grid.hybrid_alpha = 20.;
        break;
        case GOLDMINE_INDEX:
          /*
           * Recommended p_sigma = 300 hPa.
           */
          grid.zeta0        = 601.-116.*log10(defaults.pbot/100.);
          grid.zeta1        = 601.-116.*log10(defaults.ptop/100.);
          grid.hybrid_alpha = 20.;
        break;
        case MARS_INDEX:
          /*
           * TD 5/4/09
           * Recommended p_sigma = 0.2 hPa.
           */
          grid.zeta0        = 700.-300.*log10(defaults.pbot/100.);
          grid.zeta1        = 700.-300.*log10(defaults.ptop/100.);
          grid.hybrid_alpha = 20.;
        break;
        case JUPITER_INDEX:
          /*
           * CJP 2/9/05
           * The tangent line is written in the form of
           *     zeta = a0-a1*log10(p)
           * grid.zeta0 and grid.zeta1 are being determined from this equation
           * substituting p_bot and p_top respectively
           *
           * Recommended p_sigma = 100 hPa (for pbot > 100 hPa).
           */
          grid.zeta0        = 1036.-362.*log10(defaults.pbot/100.);
          grid.zeta1        = 1036.-362.*log10(defaults.ptop/100.);
          grid.hybrid_alpha = 20.;
          /*
           * NOTE: A check that zeta does not cross the theta curve might be useful.
           */
        break;
        case TITAN_INDEX:
          /*
           * Recommended p_sigma = 100 hPa.
           */
          grid.zeta0        = 30.;
          grid.zeta1        = grid.zeta0+(133.-grid.zeta0)*log(defaults.ptop/defaults.pbot)/log(10000./defaults.pbot);
          grid.hybrid_alpha = 20.;
        break;
        case URANUS_INDEX:
          /*
           *  Recommended p_sigma = 100 hPa.
           */
          grid.zeta0        = 341.0-106.0*log10(defaults.pbot/100.);
          grid.zeta1        = 341.0-106.0*log10(defaults.ptop/100.);
          grid.hybrid_alpha = 20.;
        break;
        default:
          tmp = DBL_MAX;
          for (ki = floor_tp; ki <= ceiling_tp; ki++) {
            tmp = MIN(tmp,thetadat[ki]);
          }
          grid.zeta0 = tmp*0.7;

          /*
           * Find the highest tangent point from grid.zeta0 to the theta profile,
           * call it theta_knee, sigma_knee.
           */
          tmp = DBL_MAX;
          for (ki = floor_tp; ki <= ceiling_tp; ki++) {
            if (var.pdat[ki] < defaults.ptop) {
              break;
            }
            tmp2 = fabs(var.pdat[ki]-defaults.ptop);
            if (tmp2 < tmp) {
              tmp        = tmp2;
              theta_knee = thetadat[ki];
              sigma_knee = get_sigma(defaults.pbot,var.pdat[ki]);
            }
          }
          for (ki = ceiling_tp; ki > floor_tp; ki--) {
            if (var.pdat[ki] < defaults.ptop) {
              continue;
            }
            sigma = get_sigma(defaults.pbot,var.pdat[ki]);
            slope = (thetadat[ki]-grid.zeta0)/(sigma-0.);
            tmp   = get_sigma(defaults.pbot,var.pdat[ki+1]);
            tmp2  = grid.zeta0+slope*tmp;
            if (tmp2 <= thetadat[ki+1]) {
              if (thetadat[ki] < theta_knee) {
                theta_knee = thetadat[ki];
                sigma_knee = sigma;
              }
            }
            else {
              break;
            }
          }
          if (sigma_knee > 0.) {
            grid.zeta1  = (theta_knee-grid.zeta0*(1.-sigma_knee))/sigma_knee;
          }
          else {
            sprintf(Message,"sigma_knee=%g",sigma_knee);
            epic_error(dbmsname,Message);
          }
          grid.hybrid_alpha = 20.;
          fprintf(stdout,"hybrid-coordinate parameters: grid.zeta0=%g, grid.zeta1=%g, grid.hybrid_alpha=%g\n",
                          grid.zeta0,grid.zeta1,grid.hybrid_alpha);
        break;
      }
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   *                                                       *
   * Set grid.p_ref[kk], grid.theta_ref[kk], and           *
   * the vertical coordinate grid.sigmatheta[kk].          *
   *                                                       *
   * grid.p_ref[kk]     :  typical pressure values         *
   * grid.theta_ref[kk] :  typical theta values            *
   * grid.sigmatheta[kk]:  vertical-coordinate array       *
   * grid.sgth_bot      :  sgth at bottom of model         *
   * grid.sgth_top      :  sgth at top of model            *
   *                                                       *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  grid.sgth_top = return_sigmatheta(grid.theta_ref[2*0  +1],defaults.ptop,defaults.pbot);
  grid.sgth_bot = return_sigmatheta(grid.theta_ref[2*KHI+1],defaults.pbot,defaults.pbot);

  grid.sigmatheta[1] = grid.sgth_top;

  /* 
   * Set interior points for grid.p_ref[kk] and grid.theta_ref[kk].
   */
  switch(spacing_type) {
    case SPACING_P:
      for (kk = 2; kk < 2*KHI+1; kk++) {
        grid.p_ref[kk] = defaults.ptop+(double)(kk-1)/(double)(2*KHI)*(defaults.pbot-defaults.ptop);
      }

      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = -var.pdat[ki];
        buff_triplet[ki].y = thetadat[ki];
      }
      spline_pchip(var.ntp,buff_triplet);
      ki = -2;
      for (K = KLO; K <= KHI; K++) {
        kk                   = 2*K+1;
        negp                 = -grid.p_ref[kk];
        ki                   = hunt_place_in_table(var.ntp,buff_triplet,negp,&p_d,ki);
        grid.theta_ref[kk]   = splint_pchip(negp,buff_triplet+ki,p_d);
        /*
         * For consistency, use the same averaging scheme as in the calculation of THETA2 from THETA.
         */
        grid.theta_ref[kk-1] = .5*(grid.theta_ref[kk]+grid.theta_ref[kk-2]);
      }
    break;
    case SPACING_LOGP:
      for (K = KLO; K < KHI; K++) {
        neglogp        = -log(defaults.ptop)
                         +(double)(K)/(double)(KHI)*(-log(defaults.pbot)+log(defaults.ptop));
        grid.p_ref[2*K+1] = exp(-neglogp);
      }
      /*
       * Fill in layer values.
       */
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K;
        grid.p_ref[kk] = onto_kk(P2_INDEX,grid.p_ref[kk-1],grid.p_ref[kk+1],kk);
      }

      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = neglogpdat[ki];
        buff_triplet[ki].y = thetadat[ki];
      }
      spline_pchip(var.ntp,buff_triplet);
      ki = -2;
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K+1;
        negp               = -log(grid.p_ref[kk]);
        ki                 = hunt_place_in_table(var.ntp,buff_triplet,negp,&p_d,ki);
        grid.theta_ref[kk] = splint_pchip(negp,buff_triplet+ki,p_d);
        /*
         * For consistency, use the same averaging scheme as in the calculation of THETA2 from THETA.
         */
        grid.theta_ref[kk-1] = .5*(grid.theta_ref[kk]+grid.theta_ref[kk-2]);
      }
    break;
    case SPACING_THETA:
      for (kk = 2; kk < 2*KHI+1; kk++) {
        grid.theta_ref[kk] = defaults.thetatop+(double)(kk-1)/(double)(2*KHI)*(defaults.thetabot-defaults.thetatop);
      }

      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = thetadat[ki];
        buff_triplet[ki].y = neglogpdat[ki];
      }
      spline_pchip(var.ntp,buff_triplet);
      ki = -2;
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K+1;
        tmp            = grid.theta_ref[kk];
        ki             = hunt_place_in_table(var.ntp,buff_triplet,tmp,&sgth_d,ki);
        grid.p_ref[kk] = exp(-splint_pchip(tmp,buff_triplet+ki,sgth_d));
        /*
         * For consistency, use the same averaging scheme as in the calculation of P2 from P.
         */
        grid.p_ref[kk-1] = sqrt(grid.p_ref[kk-2]*grid.p_ref[kk]);
      }
    break;
    case SPACING_FROM_FILE:
      /*
       * Have grid.p_ref[K] from a file, just need grid.theta_ref[K].
       */
      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = neglogpdat[ki];
        buff_triplet[ki].y = thetadat[ki];
      }
      spline_pchip(var.ntp,buff_triplet);
      ki = -2;
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K+1;
        negp               = -log(grid.p_ref[kk]);
        ki                 = hunt_place_in_table(var.ntp,buff_triplet,negp,&p_d,ki);
        grid.theta_ref[kk] = splint_pchip(negp,buff_triplet+ki,p_d);
        /*
         * For consistency, use the same averaging scheme as in the calculation of THETA2 from THETA.
         */
        grid.theta_ref[kk-1] = .5*(grid.theta_ref[kk]+grid.theta_ref[kk-2]);
      }
    break;
    default:
      sprintf(Message,"unrecognized spacing_type=%d",spacing_type);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * Set above-model grid.p_ref[kk]
   */
  grid.p_ref[0] = SQR(grid.p_ref[1])/grid.p_ref[2];

  /* 
   * Assign values for abyssal layer, used for gas giants.
   */
  grid.p_ref[2*KHI+2] = SQR(grid.p_ref[2*KHI+1])/grid.p_ref[2*KHI  ];
  grid.p_ref[2*KHI+3] = SQR(grid.p_ref[2*KHI+2])/grid.p_ref[2*KHI+1];

  /*
   * Set above-model grid.theta_ref[kk]
   */
  grid.theta_ref[0] = 2.*grid.theta_ref[1]-grid.theta_ref[2];

  /* 
   * Assign values for abyssal layer, used for gas giants.
   */
  grid.theta_ref[2*KHI+2] = grid.theta_ref[2*KHI+1];
  grid.theta_ref[2*KHI+3] = grid.theta_ref[2*KHI+2];

  /*
   * Set grid.t_ref[kk], grid.rho_ref[kk].
   *
   * NOTE: Need to set grid.rho_ref[kk] before calling set_re_rp().
   */
  switch(spacing_type) {
    case SPACING_P:
      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = -var.pdat[ki];
        buff_triplet[ki].y =  var.tdat[ki];
      }
      spline_pchip(var.ntp,buff_triplet);
      ki = -2;
      for (kk = 1; kk <= 2*KHI+1; kk++) {
        negp           = -grid.p_ref[kk];
        ki             = hunt_place_in_table(var.ntp,buff_triplet,negp,&p_d,ki);
        grid.t_ref[kk] = splint_pchip(negp,buff_triplet+ki,p_d);
      }

      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = -var.pdat[ki];
        fpara              = return_fpe(var.tdat[ki]);

        /* NOTE: This is only a crude estimate of the molar mass, mu. */
        mu                 = R_GAS/planet->rgas;
        buff_triplet[ki].y = return_density(fpara,var.pdat[ki],var.tdat[ki],mu,PASSING_T);
      }
      spline_pchip(var.ntp,buff_triplet);
      ki = -2;
      for (kk = 1; kk <= 2*KHI+1; kk++) {
        negp             = -grid.p_ref[kk];
        ki               = hunt_place_in_table(var.ntp,buff_triplet,negp,&p_d,ki);
        grid.rho_ref[kk] = splint_pchip(negp,buff_triplet+ki,p_d);
      }
    break;
    case SPACING_LOGP:
    case SPACING_THETA:
    case SPACING_FROM_FILE:
      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = neglogpdat[ki];
        buff_triplet[ki].y = var.tdat[ki];
      }
      spline_pchip(var.ntp,buff_triplet);
      ki = -2;
      for (kk = 1; kk <= 2*KHI+1; kk++) {
        negp           = -log(grid.p_ref[kk]);
        ki             = hunt_place_in_table(var.ntp,buff_triplet,negp,&p_d,ki);
        grid.t_ref[kk] = splint_pchip(negp,buff_triplet+ki,p_d);
      }

      for (ki = 0; ki < var.ntp; ki++) {
        buff_triplet[ki].x = neglogpdat[ki];
        fpara              = return_fpe(var.tdat[ki]);

        /* NOTE: This is only a crude estimate of the molar mass, mu. */
        mu                 = R_GAS/planet->rgas;
        buff_triplet[ki].y = return_density(fpara,var.pdat[ki],var.tdat[ki],mu,PASSING_T);
      }
      spline_pchip(var.ntp,buff_triplet);
      ki = -2;
      for (kk = 1; kk <= 2*KHI+1; kk++) {
        negp             = -log(grid.p_ref[kk]);
        ki               = hunt_place_in_table(var.ntp,buff_triplet,negp,&p_d,ki);
        grid.rho_ref[kk] = splint_pchip(negp,buff_triplet+ki,p_d);
      }
    break;
    default:
      sprintf(Message,"unrecognized spacing_type=%d",spacing_type);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * set_re_rp(): fills the arrays grid.re[K] and grid.rp[K], assuming a uniformly rotating fluid
   *
   *  set_fmn() : Coriolis parameter and geometric map factors according to geometry.
   *  grid.rln  : Longitude map factor, r
   *  grid.rlt  : Latitude map factor, R
   *  grid.f    : Coriolis parameter
   *  grid.m    : Longitude map factor 1/(r*dln*DEG) = 1./dx
   *  grid.n    : Latitude map factor 1/(R*dlt*DEG)  = 1./dy
   *  grid.mn   : 1/(grid-box area)
   *
   * set_gravity(): must be called after specifying grid.p_ref[kk], grid.re[K], and grid.rp[K].
   */
  set_re_rp();
  set_fmn();
  set_gravity();

  switch(grid.coord_type) {
    case COORD_HYBRID:
      /*
       * Set grid.sigma_sigma.
       *
       * NOTE: grid.sigmatheta[kk] is not yet set at this point.
       */
      grid.sigma_sigma = get_sigma(grid.p_ref[2*KHI+1],defaults.p_sigma);

      /*
       * Set grid.k_sigma.
       * The model's vertical coordinate is approximately a sigma coordinate (in temperature units)
       * for K = KHI up to the top of layer K = grid.k_sigma.
       */
      grid.k_sigma = -1;
      for (K = 0; K <= KHI; K++) {
        sigma = get_sigma(grid.p_ref[2*KHI+1],grid.p_ref[2*K+1]);
        if (sigma <= grid.sigma_sigma) {
          grid.k_sigma = K;
          break;
        }
      }
      if (grid.k_sigma == -1) {
        sprintf(Message,"Error calculating grid.k_sigma");
        epic_error(dbmsname,Message);
      }
      else if (grid.k_sigma < 2) {
        sprintf(Message,"Need grid.k_sigma=%d >= 2",grid.k_sigma);
        epic_error(dbmsname,Message);
      }
      grid.k_sigma = IMIN(KHI-1,grid.k_sigma);
    break;
    case COORD_ISOBARIC:
      grid.sigma_sigma =  FLT_MAX;
      grid.k_sigma     = -INT_MAX;
    break;
    case COORD_ISENTROPIC:
      grid.sigma_sigma = -FLT_MAX;
      grid.k_sigma     =  INT_MAX;
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * Set grid.sigmatheta[kk] for interior points.
   */
  if (grid.coord_type == COORD_ISENTROPIC) {
    for (kk = 2; kk < 2*KHI+1; kk++) {
      grid.sigmatheta[kk] = grid.theta_ref[kk];
    }
  }
  else {
    pbot = grid.p_ref[2*KHI+1];
    for (kk = 2; kk < 2*KHI+1; kk++) {
      pressure            = grid.p_ref[kk];
      theta               = grid.theta_ref[kk];
      grid.sigmatheta[kk] = return_sigmatheta(theta,pressure,pbot);
    }
  }

  /*
   * Set above-model grid.sigmatheta[kk]
   */
  grid.sigmatheta[0] = 2.*grid.sigmatheta[1]-grid.sigmatheta[2];

  /* 
   * Assign values for abyssal layer, used for gas giants.
   */
  grid.sigmatheta[2*KHI+1] = grid.sgth_bot;
  grid.sigmatheta[2*KHI+2] = 2.*grid.sigmatheta[2*KHI+1]-grid.sigmatheta[2*KHI  ];
  grid.sigmatheta[2*KHI+3] = 2.*grid.sigmatheta[2*KHI+2]-grid.sigmatheta[2*KHI+1];

  /*
   * Set grid.dsgth and grid.dsgth_inv.
   */
  set_dsgth();

  /*
   * Set grid.h_min[K], the minimum thickness (a.k.a. hybrid density) for
   * each layer. Konor and Arakawa use delta p_min = 1 hPa to find h_min, 
   * but this will not work for planets that have small atmospheric pressures
   * (for example, the surface pressure of Mars is ~6 hPa). 
   * Instead, we set h_min equal to a fraction of the reference value of h in each layer.
   */  
  tmp = 0.03;
  for (K = KLO; K <= KHI; K++) {
    kk            = 2*K;
    grid.h_min[K] = (tmp/grid.g[kk][2*grid.jtp+1])*(grid.p_ref[kk+1]-grid.p_ref[kk-1])*grid.dsgth_inv[kk];
  }
  grid.h_min[0] = grid.h_min[KLO]*grid.h_min[KLO]/grid.h_min[KLO+1];

  if (var.fpara.on) {
   /*
    * Preliminary initialization of para-hydrogen fraction, fpara.
    * Set to fpe(T), for representative T for layer.
    *
    * NOTE: The function init_fpara_as_fpe() cannot be used here because
    *       P2 and THETA2 are not known yet.
    */
    for (K = KLO; K <= KHI+1; K++) {
      pressure = grid.p_ref[2*K+1];
      get_sounding(pressure,"temperature",&temperature);
      fpe = return_fpe(temperature);
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FPARA(K,J,I) = fpe;
        }
      }
    }
    /* Assume fpara doesn't change with height above top of model. */
    K = 0;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        FPARA(K,J,I) = FPARA(K+1,J,I);
      }
    }
    bc_lateral(var.fpara.value,THREEDIM);
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   *                                                                     *
   * Initialize u, v, h, theta, and mixing ratios of optional species.   *
   *                                                                     *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  init_with_ref();
  init_species(&defaults,USE_PROMPTS);

  setup_mu_p();

  if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * NOTE: If gradient-balance in the meridional plane, which is the purpose of init_with_u(),
     *       is not desired, then add an appropriate initialization scheme.
     */
    init_with_u(floor_tp,ceiling_tp,&defaults);
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    init_with_ref();
  }
  else {
    sprintf(Message,"planet->type=%s, need an initialization scheme for this case",planet->type);
    epic_error(dbmsname,Message);
  }

  /*
   * Set U_SPINUP(K,J,I), the target zonal wind for Rayleigh drag.
   * Use the zonal average of U.
   */
  for (K = 0; K <= grid.nk+1; K++) {
    for (J = JLO; J <= JHI; J++) {
      avg = 0.;
      for (I = ILO; I <= IHI; I++) {
        avg += U(grid.it_uv,K,J,I);
      }
      avg /= grid.ni;
      for (I = ILO; I <= IHI; I++) {
        U_SPINUP(K,J,I) = avg;
      }
    }
  }

  /*
   * Initialize DZDT2.
   * Currently we are just setting it to zero.
   */
  memset(var.dzdt2.value,0,Nelem3d*sizeof(double));

  /*
   * Update fpara as fpe:
   */
  if (var.fpara.on) {
    init_fpara_as_fpe();
  }

  init_species(&defaults,USE_DEFAULTS);

  /* 
   * Set hyperviscosity coefficients.
   */
  grid.nudiv_nondim = defaults.nudiv_nondim;
  grid.nu_nondim    = defaults.nu_nondim;
  grid.nu_order     = defaults.nu_order;
  set_hyperviscosity(MODIFY);
  defaults.nudiv_nondim = grid.nudiv_nondim;
  defaults.nu_nondim    = grid.nu_nondim;
  defaults.nu_order     = grid.nu_order;

  /*
   * Initialize turbulence-model variables.
   */
  if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
      strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
    /*
     * NOTE: init_subgrid() must be called here before synching diagnostic variables.
     */
    init_subgrid();
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  grid.phi0 = 0.;

  /*
   * Calculate most commonly used diagnostic variables, in case they are needed below.
   *
   * NOTE: For planet->type "terrestrial" and grid.coord_type == COORD_ISENTROPIC,
   *       need to calculate PHI3(KHI,J,I) on the grid.thetabot isentropic surface
   *       and call store_pgrad_vars with PASSING_PHI3NK; this is not yet implemented.
   */
  set_p2_etc(UPDATE_THETA);
  store_pgrad_vars(SYNC_DIAGS_ONLY,CALC_PHI3NK);
  store_diag();

  /*
   * Prompt for which variables to write to extract.nc.
   */
  prompt_extract_on(defaults.extract_str,&defaults.extract_species_fraction_type,&defaults.mc_diag_extract_sum,MODIFY);

  grid.extract_species_fraction_type = defaults.extract_species_fraction_type;
  grid.mc_diag_extract_sum           = defaults.mc_diag_extract_sum;

  /*
   * Write defaults file:
   */
  write_defaults(&defaults);

  /* 
   * Print out zonal-wind information: 
   */
  print_zonal_info();

  /* 
   * Print out vertical information: 
   */
  print_vertical_column(grid.jtp,ILO,"vertical.dat");

  /*
   * Call vertical_modes() to write vertical eigenvalues and eigenvectors to a file.
   */
  if (KHI >= 4) {
    vertical_modes(grid.jtp,ILO);
  }

  /* 
   * Output epic.nc file (netCDF = Network Common Data Form):
   */
  sprintf(outfile,"epic.nc");
  time_index = 0;
  var_write(outfile,ALL_DATA,time_index,0);

  /* 
   * Free allocated memory: 
   */
  free_arrays(); 
  free_dtriplet(buff_triplet,0,var.ntp-1,         dbmsname);
  free_dvector(neglogpdat,   0,var.ntp-1,         dbmsname);
  free_dvector(thetadat,     0,var.ntp-1,         dbmsname);
  free_var_props();
  for (I = 0; I < NUM_WORKING_BUFFERS; I++) {
    free_dvector(Buff2D[I],0,Nelem2d-1,dbmsname);
  }

  /*
   * NOTE: In initial (but not change or epic), the structure memory that the pointer "planet" points to 
   *       is not dynamically allocated, but is the memory of one of the statically allocated planets
   *       in epic_globals.c, therefore we do not call "free(planet)" here.
   */

  return 0;
}

/*======================= end of main() =====================================*/

/*======================= read_defaults() ===================================*/

/*
 * NOTE: To use the READ* macros in epic_io_macros.h, dummy
 *       io and fd variables are declared.
 */

void read_defaults(init_defaultspec *def) 
{
  int
    index,
    nc_id,nc_err;
  double
    solar;
  char
    min_element[4];
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts   =0,
    num_vars =0;
  nc_type
    the_nc_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_defaults";

  nc_err = lookup_netcdf("init_defaults.nc",
                         &nc_id,&ngatts,&gattname,&num_vars,&varname);

  if (nc_err == NC_NOERR) {
    READI(&def->start_date_input_type,def_start_date_input_type,1);
    READTIME(&def->start_time,def_start_time);  

    READC(def->geometry,def_geometry,GEOM_STR);
    READC(def->system_id,def_system_id,32);
    READC(def->f_plane_map,def_f_plane_map,GEOM_STR);
    READC(def->eos,def_eos,8);
    READC(def->extract_str,def_extract_str,N_STR);
    READC(def->species_str,def_species_str,N_STR);
    READC(def->layer_spacing_dat,def_layer_spacing_dat,N_STR);
    READC(def->turbulence_scheme,def_turbulence_scheme,N_STR);

    READI(&def->nk,def_nk,1);
    READI(&def->nj,def_nj,1);
    READI(&def->ni,def_ni,1);
    READI(&def->dt,def_dt,1);
    READI(&def->newt_cool_adjust,def_newt_cool_adjust,1);
    READI(&def->cloud_microphysics,def_cloud_microphysics,1);
    READI(&def->moist_convection,def_moist_convection,1);
    READI(&def->max_mc_it,def_max_mc_it,1);
    READD(&def->tau_relax,def_tau_relax,1);
    READI(&def->grid_wind_shear_mode,grid_wind_shear_mode,1);
    READI(&def->grid_wind_decay_m_const,grid_wind_decay_m_const,1);

    READI(def->on,def_on,LAST_SPECIES+1);

    READD(&def->du_vert,def_du_vert,1);
    READD(&def->du_vert_m,du_vert_m,1);
    READD(&def->du_vert_pend,du_vert_pend,1);

    READI(&def->grid_cloud_top_mode,grid_cloud_top_mode,1);
    READD(&def->grid_if_400mb,grid_if_400mb,1);
    READD(&def->grid_if_1000mb,grid_if_1000mb,1);

    READI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
    READI(&def->radiation_index,def_radiation_index,1);
    READI(&def->zonal_average_rt,def_zonal_average_rt,1);
    READI(&def->extract_species_fraction_type,def_extract_species_fraction_type,1);
    READI(&def->mc_diag_extract_sum,def_mc_diag_extract_sum,1);
    READI(&def->spacing_type,def_spacing_type,1);
    READI(&def->coord_type,def_coord_type,1);
    READI(&def->k_sponge,def_k_sponge,1);
    READI(&def->n_bot_drag,def_n_bot_drag,1);
    READI(&def->j_sponge,def_j_sponge,1);
    

    READD(&def->heat_rate,grid_heat_rate,1);
    READD(&def->cool_rate,grid_cool_rate,1);
    READD(&def->heat_top_pressure,grid_heat_top_pressure,1);
    READD(&def->cool_bot_pressure,grid_cool_bot_pressure,1);
    
    READI(&def->grid_relax_vapor,grid_relax_vapor,1);
    READD(&def->grid_relax_vapor_timescale,grid_relax_vapor_timescale,1);

    READD(&def->globe_lonbot,def_globe_lonbot,1);
    READD(&def->globe_lontop,def_globe_lontop,1);
    READD(&def->globe_latbot,def_globe_latbot,1);
    READD(&def->globe_lattop,def_globe_lattop,1);
    READD(&def->lat_tp,def_lat_tp,1);
    READD(&def->f_plane_half_width,def_f_plane_half_width,1);
    READD(&def->f_plane_lat0,def_f_plane_lat0,1);
    READD(&def->ptop,def_ptop,1);
    READD(&def->pbot,def_pbot,1);
    READD(&def->thetatop,def_thetatop,1);
    READD(&def->thetabot,def_thetabot,1);
    READD(&def->p_sigma,def_p_sigma,1);
    READD(&def->nudiv_nondim,def_nudiv_nondim,1);
    READD(&def->u_scale,def_u_scale,1);
    READI(&def->nu_order,def_nu_order,1);
    READD(&def->nu_nondim,def_nu_nondim,1);
    READD(&def->fpara_rate_scaling,def_fpara_rate_scaling,1);

    READD(def->mole_fraction,def_mole_fraction,LAST_SPECIES+1);
    READD(def->mole_fraction_over_solar,def_mole_fraction_over_solar,LAST_SPECIES+1);
    READD(def->rh_max,def_rh_max,LAST_SPECIES+1);
  }
  else {
    /*
     * If the file is not readable, use standard defaults.
     */
    strcpy(def->geometry,"globe");
    strcpy(def->system_id,  "Jupiter");
    strcpy(def->f_plane_map,"polar");
    strcpy(def->eos,        "ideal");
    strcpy(def->extract_str,"none");
    strcpy(def->species_str,"none");
    strcpy(def->layer_spacing_dat,"layer_spacing.dat");
    strcpy(def->turbulence_scheme,"off");

    def->start_date_input_type   = NOT_SET;
    def->start_time              = (time_t)(-2145916800);  /* 1902_01_01_00:00:00 (UTC) */

    def->nk                      =  NOT_SET;
    def->nj                      =  64;
    def->ni                      =  128;
    def->dt                      =  120;
    def->newt_cool_adjust        =  FALSE;
    def->cloud_microphysics      =  NOT_SET;

    def->on[    U_INDEX]    =  TRUE;
    def->on[    V_INDEX]    =  TRUE;
    def->on[    H_INDEX]    =  TRUE;
    def->on[THETA_INDEX]    =  TRUE;
    for (index = THETA_INDEX+1; index <= LAST_SPECIES; index++) {
      def->on[index]   = FALSE;
    }
    
    def->grid_wind_shear_mode           = WIND_SHEAR_FUNCTION;
    def->grid_wind_decay_m_const        = WIND_DECAY_SLOPE_CONSTANT;
    def->du_vert                        = .0;
    def->du_vert_m                      = 0.;
    def->du_vert_pend                   = 30000.e2;

    def->grid_cloud_top_mode = CONSTANT_CLOUD_BASE;
    def->grid_if_400mb = 0.09;
    def->grid_if_1000mb = 0.05;

    def->uv_timestep_scheme             =  0;
    def->radiation_index                =  0;
    def->zonal_average_rt               =  TRUE;
    def->extract_species_fraction_type  =  MASS;
    def->mc_diag_extract_sum            =  0000;
    def->spacing_type                   =  SPACING_LOGP;
    def->coord_type                     =  NOT_SET;
    def->k_sponge                       =  NOT_SET;
    def->n_bot_drag                     =  NOT_SET;
    def->j_sponge                       =  NOT_SET;

    def->heat_top_pressure  = 20000.; // hPa
    def->cool_bot_pressure  = 50.;    // hPa
    def->heat_rate = 0.008; // K/day
    def->cool_rate = 0.01;  // K/day
  
    def->grid_relax_vapor = OFF;
    def->grid_relax_vapor_timescale = 5;

    def->globe_lonbot       = -180.;
    def->globe_lontop       =  180.;
    def->globe_latbot       = -90.;
    def->globe_lattop       =  90.;
    def->lat_tp             =  0.;
    def->f_plane_half_width =  90.;
    def->f_plane_lat0       =  0.;
    def->ptop               =  .1*100.;
    def->pbot               =  10000.*100.;
    def->thetatop           = 820.;
    def->thetabot           = 320.;
    def->p_sigma            = (double)NOT_SET;
    def->nudiv_nondim       =  0.;
    def->u_scale            =  1.;
    def->nu_order           =  6;
    def->nu_nondim          =  0.1;
    def->fpara_rate_scaling =  1.0;
    for (index = FIRST_SPECIES; index <= LAST_SPECIES; index++) {
      solar = solar_fraction(var.species[index].info[0].name,MOLAR,min_element);

      def->mole_fraction[index]            = solar;
      def->mole_fraction_over_solar[index] = 1.;
      def->rh_max[index]                   = 1.;
    }
  }

  return;
}

/*======================= end of read_defaults() ============================*/

/*======================= write_defaults() ==================================*/

/*
 * NOTE: To use the WRITE* macros in epic_io_macros.h, dummy
 *       io and fd variables are declared.
 */

void write_defaults(init_defaultspec *def)
{
  int
    nc_id,nc_err;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="write_defaults";

  nc_err = nc_create("init_defaults.nc",NC_CLOBBER,&nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  WRITEI(&def->start_date_input_type,def_start_date_input_type,1);
  WRITETIME(&def->start_time,def_start_time);  

  WRITEC(def->geometry,def_geometry,GEOM_STR);
  WRITEC(def->system_id,def_system_id,32);
  WRITEC(def->f_plane_map,def_f_plane_map,GEOM_STR);
  WRITEC(def->eos,def_eos,8);
  WRITEC(def->extract_str,def_extract_str,N_STR);
  WRITEC(def->species_str,def_species_str,N_STR);
  WRITEC(def->layer_spacing_dat,def_layer_spacing_dat,N_STR);
  WRITEC(def->turbulence_scheme,def_turbulence_scheme,N_STR);

  WRITEI(&def->nk,def_nk,1);
  WRITEI(&def->nj,def_nj,1);
  WRITEI(&def->ni,def_ni,1);
  WRITEI(&def->dt,def_dt,1);
  WRITEI(&def->newt_cool_adjust,def_newt_cool_adjust,1);
  WRITEI(&def->cloud_microphysics,def_cloud_microphysics,1);
  WRITEI(&def->moist_convection,def_moist_convection,1);
  WRITEI(&def->max_mc_it,def_max_mc_it,1);
  WRITED(&def->tau_relax,def_tau_relax,1);

  WRITEI(def->on,def_on,LAST_SPECIES+1);
  WRITEI(&def->grid_wind_shear_mode,grid_wind_shear_mode,1);
  WRITEI(&def->grid_wind_decay_m_const,grid_wind_decay_m_const,1);
  WRITED(&def->du_vert_m,du_vert_m,1);
  WRITED(&def->du_vert_pend,du_vert_pend,1);
  WRITED(&def->du_vert,def_du_vert,1);
  WRITEI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
  WRITEI(&def->radiation_index,def_radiation_index,1);
  WRITEI(&def->zonal_average_rt,def_zonal_average_rt,1);
  WRITEI(&def->extract_species_fraction_type,def_extract_species_fraction_type,1);
  WRITEI(&def->mc_diag_extract_sum,def_mc_diag_extract_sum,1);
  WRITEI(&def->spacing_type,def_spacing_type,1);
  WRITEI(&def->coord_type,def_coord_type,1);
  WRITEI(&def->k_sponge,def_k_sponge,1);
  WRITEI(&def->n_bot_drag,def_n_bot_drag,1);
  WRITEI(&def->j_sponge,def_j_sponge,1);
  
  WRITEI(&def->grid_cloud_top_mode,grid_cloud_top_mode,1);
  WRITED(&def->grid_if_400mb,grid_if_400mb,1);
  WRITED(&def->grid_if_1000mb,grid_if_1000mb,1);
    
  WRITED(&def->heat_rate,grid_heat_rate,1);
  WRITED(&def->cool_rate,grid_cool_rate,1);
  WRITED(&def->heat_top_pressure,grid_heat_top_pressure,1);
  WRITED(&def->cool_bot_pressure,grid_cool_bot_pressure,1);
  
  WRITEI(&def->grid_relax_vapor,grid_relax_vapor,1);
  WRITED(&def->grid_relax_vapor_timescale,grid_relax_vapor_timescale,1);

  WRITED(&def->globe_lonbot,def_globe_lonbot,1);
  WRITED(&def->globe_lontop,def_globe_lontop,1);
  WRITED(&def->globe_latbot,def_globe_latbot,1);
  WRITED(&def->globe_lattop,def_globe_lattop,1);
  WRITED(&def->lat_tp,def_lat_tp,1);
  WRITED(&def->f_plane_half_width,def_f_plane_half_width,1);
  WRITED(&def->f_plane_lat0,def_f_plane_lat0,1);
  WRITED(&def->ptop,def_ptop,1);
  WRITED(&def->pbot,def_pbot,1);
  WRITED(&def->thetatop,def_thetatop,1);
  WRITED(&def->thetabot,def_thetabot,1);
  WRITED(&def->p_sigma,def_p_sigma,1);
  WRITED(&def->nudiv_nondim,def_nudiv_nondim,1);
  WRITED(&def->u_scale,def_u_scale,1);
  WRITEI(&def->nu_order,def_nu_order,1);
  WRITED(&def->nu_nondim,def_nu_nondim,1);
  WRITED(&def->fpara_rate_scaling,def_fpara_rate_scaling,1);
  WRITED(def->mole_fraction,def_mole_fraction,LAST_SPECIES+1);
  WRITED(def->mole_fraction_over_solar,def_mole_fraction_over_solar,LAST_SPECIES+1);
  WRITED(def->vmr_slope,def_vmr_slope,LAST_SPECIES+1);
  WRITED(def->vmr_pcrit,def_vmr_pcrit,LAST_SPECIES+1);
  WRITED(def->rh_max,def_rh_max,LAST_SPECIES+1);

  nc_close(nc_id);

  return;
}

/*======================= end of write_defaults() ===========================*/

/* * * * * * * * * * * *  end of epic_initial.c  * * * * * * * * * * * * * * */







