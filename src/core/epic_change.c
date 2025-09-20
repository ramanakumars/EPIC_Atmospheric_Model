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

/* * * * * * * * * * epic_change.c * * * * * * * * * * * * * * * * * 
 *                                                                 *
 *  Makes changes to parameters in epic.nc                         *
 *                                                                 *
 *  Options for introducing perturbations:                         *
 *                                                                 *
 *    Use -spots spots.dat to add vortices                         *
 *        -waves waves.dat to add waves                            *
 *        -heat  spots.dat to add thermal perturbation             *
 *    These can all be done at the same time.                      *
 *                                                                 *
 *  Options for converting data to isentropic coordinates:         *
 *                                                                 *
 *    Use -openmars to process an OpenMARS file                    *
 *    Use -emars to process an EMARS file                          *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <epic_datatypes.h>
#include <epic.h>

// boolean value on whether the we are loading from a defaults file
// or if the user is allowed to modify the defaults option
int MODIFY = 1;

/*
 * Local function prototypes:
 */
void add_spots(char    *spots_file,
               double  *pert);

void add_noise(char *noise_file);

void add_waves(char     *waves_file,
               double   *pert);

void add_heat(char *heat_file);

void modify_progs_from_pert(double *pert);

int number_objects_in_file(char *objects_file);

void read_spots_file(char   *spots_file,
                     double *ampspot,
                     double *lonspot,
                     double *latspot,
                     double *pspot,
                     double *aspot,
                     double *bspot,
                     double *cspot_up,
                     double *cspot_down,
                     int     adjust_amplitude);

void read_waves_file(char   *waves_file,
                     double *latwave,
                     double *ampwave,
                     double *wnwave,
                     double *pwave,
                     double *cwave,
                     double *fwhmwave);

void read_defaults(change_defaultspec  *def, char *defaults_file);

void write_defaults(change_defaultspec *def);

/*
 * Conversion of non-EPIC input into isentropic-coordinate output.
 *
 *   Input options:
 *     OpenMars Mars data
 *     EMARS Mars data (emars_v1.0_back_*.nc format)
 *
 *   The customized data-processing functions are in epic_funcs_init.c,
 *   the associated function prototypes and shift macros are in epic.h,
 *   and the associated data types are in epic_datatypes.h.
 */

#undef  PERT
#define PERT(k,j,i) pert[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

/*======================= main() =====================================*/

/*
 *  NOTE: structures planet, grid, and var are declared globally in epic.h.
 */

int main(int   argc,
         char *argv[])
{
  char   
    spots_file[FILE_STR], /*  added-spot locations and sizes            */
    noise_file[FILE_STR], /*  added-windfield noise locations and sizes */
    waves_file[FILE_STR], /*  wave perturbation parameters              */
    heat_file[FILE_STR],  /*  file for thermal perturbation             */
    sflag[80],            /*  string to hold command-line flags         */
    infile[ FILE_STR],
    outfile[FILE_STR],
    openmars_infile[FILE_STR],
    openmars_outfile_qb[FILE_STR],
    openmars_outfile_uvpt[FILE_STR],
    emars_infile[FILE_STR],
    emars_outfile_qb[FILE_STR],
    emars_outfile_uvpt[FILE_STR],
    defaults_file[FILE_STR],
    buffer[16];
  int    
    time_index,
    openmars_itime,
    emars_itime,
    K,J,I,
    kk,jj,
    is,itmp,
    count,index,ii;
  int
    spots      = FALSE,
    noise      = FALSE,
    waves      = FALSE,
    heat       = FALSE,  /* RS 04/16/2020 adding thermal source */
    stretch_ni = FALSE,
    openmars   = FALSE,
    emars      = FALSE,
    recloud    = FALSE;
  openmars_gridspec
    *openmars_grid;
  emars_gridspec
    *emars_grid;
  double  
    dx0,dt;
  double
    *p,*h;
  static double
    *Buff2D[NUM_WORKING_BUFFERS];
  change_defaultspec
    defaults;
  init_defaultspec
    defaults_init;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="epic_change";

#if defined(EPIC_MPI)
  sprintf(Message,"not designed to be run on multiple processors");
  epic_error(dbmsname,Message);
#endif

  declare_copyright();

  /* 
   * Interpret command-line arguments: 
   */
  /* Start with defaults: */
  sprintf(spots_file,"none");
  sprintf(waves_file,"none");
  sprintf( heat_file,"none");
  sprintf(defaults_file,"change_defaults.nc");
  if (argc > 1) {
    /* Read flags: */
    for (count = 1; count < argc; count++) {
      sscanf(argv[count],"%s",sflag);
      if (strcmp(sflag,"-spots") == 0) {
        sscanf(argv[++count],"%s",spots_file);
        spots = TRUE;
      }
      else if (strcmp(sflag,"-noise") == 0) {
        sscanf(argv[++count],"%s",noise_file);
        noise = TRUE;
      }
      else if (strcmp(sflag,"-waves") == 0) {
        sscanf(argv[++count],"%s",waves_file);
        waves = TRUE;
      }
      else if (strcmp(sflag,"-heat") == 0) {
        sscanf(argv[++count],"%s",heat_file);
        heat = TRUE;
      }
      else if (strcmp(sflag,"-stretch_ni") == 0) {
        sscanf(argv[++count],"%d",&stretch_ni);
        /*
         * Verify that stretch_ni is a power of 2.
         */
        if (frexp((double)stretch_ni,&itmp) != 0.5) {
          sprintf(Message,"-stretch_ni %d is not a power of 2",stretch_ni);
          epic_error(dbmsname,Message);
        }
      }
      else if (strcmp(sflag,"-openmars") == 0) {
        /*
         * OpenMARS (formerly MACDA)
         * Input an OpenMARS .nc file and convert to isentropic-coordinate QB* and UVPT* files.
         */
        openmars_grid = (openmars_gridspec *)calloc(1,sizeof(openmars_gridspec));

        openmars = TRUE;
      }
      else if (strcmp(sflag,"-emars") == 0) {
        /*
         * EMARS (Ensemble Mars Atmosphere Reanalysis System)
         * Input an EMARS .nc file and convert to isentropic-coordinate QB* and UVPT* files.
         */
        emars_grid = (emars_gridspec *)calloc(1,sizeof(emars_gridspec));

        emars = TRUE;
      }
      else if (strcmp(sflag, "-from_defaults") == 0) {
        sscanf(argv[++count],"%s",defaults_file);

        // if we load in from a defaults file, then turn off the modify flag
        // since we are not asking for input from the user
        MODIFY = 0;
      }
      else if (strcmp(sflag,"-help") == 0 ||
               strcmp(sflag,"-h")    == 0) {
        /* Print help, exit: */
        system("more "EPIC_PATH"/help/epic_change.help");
        exit(1);
      }
      else {
        sprintf(Message,"Unrecognized change command-line flag: %s \n",sflag);
        epic_error(dbmsname,Message);
      }
    }
  }

  if (openmars || emars) {
    planet = &mars;
  }
  else {
    /* Allocate memory */
    if((planet=( planetspec *)malloc(sizeof(planetspec))) == 0) {
      sprintf(Message,"allocating space for planetspec \n");
      epic_error(dbmsname,Message);
    }
  }

#if defined(EPIC_MPI)
  MPI_Init(&argc,&argv);
  para.comm = MPI_COMM_WORLD;
  MPI_Comm_set_errhandler(para.comm,MPI_ERRORS_RETURN);
  MPI_Comm_rank(para.comm,&para.iamnode);
  MPI_Comm_size(para.comm,&para.nproc);
  para.ndim = NINT(log((double)para.nproc)/log(2.));
  /*
   * epic_change.c is not set up to be run across multiple processors.
   */
  if (para.nproc > 1) {
    sprintf(Message,"Not set up to be run across multiple processors");
    epic_error(dbmsname,Message);
  }
#endif

  /*
   *  Read in default parameter settings:
   */
  read_defaults(&defaults,defaults_file);

  /* Time-plane index. */
  time_index = 0;

  /*
   * Determine model size and allocate memory for arrays.
   */
  if (openmars) {
    int
      nc_id,nc_err;

    input_string("Input OpenMARS data file [netCDF format]\n",defaults.openmars_infile,openmars_infile,MODIFY);
    openmars_itime = 0;
    openmars_var_read(openmars_grid,openmars_infile,SIZE_DATA,openmars_itime);
    openmars_make_arrays(openmars_grid);
    /*
     * NOTE: For OpenMARS, the portion POST_SIZE_DATA refers to the 1D constant
     *       arrays for the dimensions, but not the 2D and 3D variable fields.  The latter
     *       portion is referred to by VAR_DATA. The reason for the two-part size-data
     *       partition for OpenMARS is that grid.nj depends on the latitude spacing in OpenMARS,
     *       which comes from the POST_SIZE_DATA segment.
     */
    openmars_var_read(openmars_grid,openmars_infile,POST_SIZE_DATA,openmars_itime);

     /*
      * Use an EPIC file, openmars_epic.nc, to establish an appropriate EPIC Mars environment.
      * Check whether openmars_epic.nc already exists, and if not, create it.
      */
    nc_err = nc_open("./openmars_epic.nc",NC_NOWRITE,&nc_id);
    if (nc_err) {
      openmars_epic_nc();
    }
    else {
      /*
       * The file openmars_epic.nc already exists.
       */
      nc_close(nc_id);
    }

    sprintf(defaults.infile,"./openmars_epic.nc");
    sprintf(infile,"%s",defaults.infile);
  }
  else if (emars) {
    int
      nc_id,nc_err;

    input_string("Input EMARS data file [netCDF format]\n",defaults.emars_infile,emars_infile,MODIFY);
    emars_itime = 0;
    emars_var_read(emars_grid,emars_infile,SIZE_DATA,emars_itime);
    emars_make_arrays(emars_grid);
    /*
     * NOTE: For EMARS, the portion POST_SIZE_DATA refers to the 1D constant
     *       arrays for the dimensions, but not the 2D and 3D variable fields.  The latter
     *       portion is referred to by VAR_DATA. The reason for the two-part size-data
     *       partition for EMARS is that grid.nj depends on the latitude spacing in EMARS,
     *       which comes from the POST_SIZE_DATA segment.
     */
    emars_var_read(emars_grid,emars_infile,POST_SIZE_DATA,emars_itime);

     /*
      * Use an EPIC file, emars_epic.nc, to establish an appropriate EPIC Mars environment.
      * Check whether emars_epic.nc already exists, and if not, create it.
      */
    nc_err = nc_open("./emars_epic.nc",NC_NOWRITE,&nc_id);
    if (nc_err) {
      emars_epic_nc();
    }
    else {
      /*
       * The file emars_epic.nc already exists.
       */
      nc_close(nc_id);
    }

    sprintf(defaults.infile,"./emars_epic.nc");
    sprintf(infile,"%s",defaults.infile);
  }
  else {
    input_string("Input EPIC file [netCDF format]\n",defaults.infile,infile,MODIFY);
  }

  /* NOTE: time_index is not used for SIZE_DATA */
  var_read(infile,SIZE_DATA,time_index);

  set_var_props();
  make_arrays();

  for (I = 0; I < NUM_WORKING_BUFFERS; I++) {
    Buff2D[I] = dvector(0,Nelem2d-1,dbmsname);
  }

  /*
   * Read in rest of input data.
   */
  var_read(infile,POST_SIZE_DATA,time_index);

  /* timeplane_bookkeeping() must come after reading in variables. */
  timeplane_bookkeeping();

  /* 
   * Set lon, lat, etc. 
   */
  set_lonlat();
  set_fmn();
  set_gravity();
  set_dsgth();

  /*
   * Set up thermodynamics subroutines.
   *
   * NOTE: The return value cpr is a low-temperature reference value, and
   *       should not be used otherwise.  Use return_cp() for a given
   *       thermodynamical state.
   */
  thermo_setup(&planet->cpr);

  /* 
   * Store diagnostic variables. 
   */
  if (!openmars && !emars) {
    fprintf(stdout,"\nCalculating and storing diagnostic variables...");fflush(stdout);
  }

  /*
   * Reconstitute the prognostic variable H from input P.
   */

  /* Allocate memory */
  p = dvector(0,2*grid.nk+1,dbmsname);
  h = dvector(0,2*grid.nk+1,dbmsname);

  for (J = JLO; J <= JHI; J++) {
    jj = 2*J+1;
    for (I = ILO; I <= IHI; I++) {
      for (kk = 1; kk <= 2*KHI+1; kk++) {
        p[kk] = get_p(P2_INDEX,kk,J,I);
      }
      calc_h(jj,p,h);

      for (K = KLO; K <= KHI; K++) {
        H(K,J,I) = h[2*K];
      }
      K = 0;
      H(K,J,I) = SQR(H(K+1,J,I))/H(K+2,J,I);
      K = KHI+1;
      H(K,J,I) = SQR(H(K-1,J,I))/H(K-2,J,I);
    }
  }
  bc_lateral(var.h.value,THREEDIM);

  if (openmars) {
    openmars_conversion(openmars_grid,openmars_infile,openmars_outfile_qb,openmars_outfile_uvpt);

    goto cleanup;
  }
  else if (emars) {
    emars_conversion(emars_grid,emars_infile,emars_outfile_qb,emars_outfile_uvpt);

    goto cleanup;
  }
  else {
    /*
     * No long jump (goto).
     */
    ;
  }

  /*
   * NOTE: For planet->type "terrestrial" and grid.coord_type == COORD_ISENTROPIC,
   *       need to calculate PHI3(KHI,J,I) on the grid.thetabot isentropic surface
   *       and call store_pgrad_vars with PASSING_PHI3NK; this is not yet implemented.
   */
  set_p2_etc(UPDATE_THETA);
  store_pgrad_vars(SYNC_DIAGS_ONLY,CALC_PHI3NK);
  store_diag();

  /*
   * Print out a listing of important model parameters.
   */
  print_model_description();

  /* 
   *  Print out vertical information:
   */
  print_vertical_column(JLO,ILO,"vertical.dat");

  /*
   *  Change parameters as instructed.
   */
  grid.dt = input_int("\nInput timestep\n", grid.dt,MODIFY);

  /*
   * Inquire about radiation scheme.
   */
  inquire_radiation_scheme(MODIFY);

  if (strcmp(grid.radiation_scheme, "Global heating-cooling") == 0) {
    defaults.heat_top_pressure = grid.heat_top_pressure = input_double("Top of the heating region [hPa]\n", defaults.heat_top_pressure, MODIFY);
    defaults.cool_bot_pressure = grid.cool_bot_pressure = input_double("Bottom of the cooling region [hPa]\n", defaults.cool_bot_pressure, MODIFY);
    defaults.heat_rate = grid.heat_rate = input_double("Heating rate at the bottom [K/day]\n", defaults.heat_rate, MODIFY);
    defaults.cool_rate = grid.cool_rate = input_double("Cooling rate at the top [K/day]\n", defaults.cool_rate, MODIFY);
    
    // convert to appropriate units
    grid.heat_top_pressure = grid.heat_top_pressure * 100.;
    grid.cool_bot_pressure = grid.cool_bot_pressure * 100.;
    grid.heat_rate = grid.heat_rate / 86400.;
    grid.cool_rate = grid.cool_rate / 86400.;
  }


  if (var.fpara.on) {
    /*
     * Inquire about fpara_rate_scaling.
     */
    sprintf(Message,"Input ortho-para H2 conversion rate scaling [nominal is 1.0]:\n");
    var.fpara_rate_scaling = input_double(Message,var.fpara_rate_scaling,MODIFY);
  }

  /*
   * Inquire whether to change the status of cloud microphysics.
   */
  if (grid.cloud_microphysics != OFF) {
    sprintf(Message,"Cloud microphysics: %2d => off \n"
                    "                    %2d => active  (latent heat, phase changes, precipitation) \n"
                    "                    %2d => passive (advection only)\n",
                    OFF,ACTIVE,PASSIVE);
    grid.cloud_microphysics = input_int(Message,grid.cloud_microphysics,MODIFY);

    /* 
     * Check if vapor needs to be trimmed or re-initialized.
     */
    defaults.reinit_cloud = input_int("Re-initialize vapor? (0 = no, 1 = trim excess, 2 = re-initialize):\n",recloud,MODIFY);
    switch(defaults.reinit_cloud) {
      case 0:
        fprintf(stdout,"Not modifying vapor.\n");
      break;

      case 1:
        fprintf(stdout,"Trimming supersaturated vapor.\n");
        change_species(&defaults,USE_PROMPTS,CHANGEMODE); 
      break;

      case 2:
        fprintf(stdout,"Re-initializing vapor. \n");
        change_species(&defaults,USE_PROMPTS,INITMODE); 
      break;

      default:
        sprintf(Message,"Unrecognized recloud input, defaulting to 0 = no.");
        epic_warning(dbmsname,Message);
        fprintf(stdout,"Not modifying vapor.\n");
      break;
    }
    
    if (grid.moist_convection == NOT_SET || grid.moist_convection == OFF) {
      grid.moist_convection = OFF;
    }
    sprintf(Message,"Moist convective scheme: %2d => off \n"
                    "                         %2d => on (Sud & Walker, 1999) \n"
		    "                         %2d => passive (calculate diagnostic variables only) \n",
                    OFF,ACTIVE,PASSIVE);
    grid.moist_convection = input_int(Message,grid.moist_convection,MODIFY);

    if( grid.moist_convection == ACTIVE) {
      sprintf(Message,"Max number of iterations: \n");
      grid.max_mc_it = input_int(Message,grid.max_mc_it,MODIFY);

      sprintf(Message,"Relaxation timescale: \n");
      grid.tau_relax = input_double(Message,grid.tau_relax,MODIFY);

      /* 
       * Reset the counter to state that RAS will be calculated for the first time.
       */
      grid.first_RAS_upd = 1;
    }
    else {
      grid.moist_convection = OFF;
      grid.max_mc_it        = 0;
      grid.tau_relax        = 0.;
    }
  }

  /* for relaxing the vapor back to the initial state */
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
   * Set sponges, drag layers:
   */
  grid.k_sponge   = input_int("Input k_sponge (-1 = no effect):\n",
                              grid.k_sponge,MODIFY);
  grid.n_bot_drag = input_int("Input number of bottom layers with transitional drag towards abyssal wind profile (-1 = no effect):\n",
                              grid.n_bot_drag,MODIFY);
  grid.j_sponge = input_int("Input j_sponge (-1 = no effect; 3 is typical):\n",
                            grid.j_sponge,MODIFY);

  /*
   * Revisit hyperviscosity coefficients (since dt may have changed).
   */
  set_hyperviscosity(MODIFY);
  
  /*
   * Add perturbations if requested, in the form of spots (vortices) and/or waves.
   */
  if (spots || waves) {
    double
      *pert;

    /* 
     * Allocate memory for perturbation streamfunction
     */
    pert = dvector(0,Nelem3d-1,dbmsname);

    if (spots) {
      add_spots(spots_file,pert);
    }

    if (waves) {
      add_waves(waves_file,pert);
    }

    modify_progs_from_pert(pert);

    /*
     * Free allocated memory.
     */
    free_dvector(pert,0,Nelem3d-1,dbmsname);

    /*
     * Update most commonly used diagnostic variables, in case they are needed.
     *
     * NOTE: For planet->type "terrestrial" and grid.coord_type == COORD_ISENTROPIC,
     *       need to calculate PHI3(KHI,J,I) on the grid.thetabot isentropic surface
     *       and call store_pgrad_vars with PASSING_PHI3NK; this is not yet implemented.
     */
    set_p2_etc(UPDATE_THETA);
    store_pgrad_vars(SYNC_DIAGS_ONLY,CALC_PHI3NK);
    store_diag();
  }

  if (noise) {
    add_noise(noise_file);
  }

  /*
   * Add thermal perturbation if requested:
   */
  if (heat) {
    add_heat(heat_file);
  }

  /*
   * Prompt for which variables to write to extract.nc.
   */
  defaults.extract_str[0] = '\0';
  for (index = FIRST_INDEX; index <= LAST_INDEX; index++) {
    if (var.extract_on_list[index] == LISTED_AND_ON) {
      sprintf(buffer," %d",index);
      strcat(defaults.extract_str,buffer);
    }
  }
  prompt_extract_on(defaults.extract_str,&grid.extract_species_fraction_type,&grid.mc_diag_extract_sum,MODIFY);

  /*
   * Write epic.nc file.
   */
  input_string("Name of output file\n",defaults.outfile,outfile,MODIFY);
  var_write(outfile,ALL_DATA,time_index,stretch_ni);

  /*---------------*
   * Cleanup block *
   *---------------*/
  cleanup:

  /* Write defaults file if we are not using a defaults input: */
  if(MODIFY) {
    write_defaults(&defaults);
  }

  /*
   * Free dynamically allocated memory.
   */
  free_dvector(p,0,2*grid.nk+1,dbmsname);
  free_dvector(h,0,2*grid.nk+1,dbmsname);

  free_arrays();
  free_var_props();

  if (openmars) {
    openmars_free_arrays(openmars_grid);
    free(openmars_grid);
  }
  else if (emars) {
    emars_free_arrays(emars_grid);
    free(emars_grid);
  }

  if (openmars || emars) {
    /*
     * No need to free planet structure, since it
     * was not dynamically allocated.
     */
    ;
  }
  else {
    free(planet);
  }

  for (I = 0; I < NUM_WORKING_BUFFERS; I++) {
    free_dvector(Buff2D[I],0,Nelem2d-1,dbmsname);
  }

  return 0;
}

/*======================= end of main() =====================================*/

/*======================= add_spots() =======================================*/

/*
 * Add vortices via a perturbation streamfunction.
 *
 * Example spots.dat file:
 * ---------------------------------------------------------------------------------------
 *  Number of vortices: 2
 *  lon[deg] lat[deg] press[hPa]  a[deg] b[deg] c_up[scale_hts] c_down[scale_hts] amp[m/s]
 *   30.      -33.     680.       3.0    2.5        2.5               3.0          100.
 *   60.      -33.5    680.       3.0    2.5        2.5               3.0          100.
 * ---------------------------------------------------------------------------------------
 *
 * NOTE: Not MPI ready. 
 */

/*
 * Implemented styles for the vortex perturbation streamfunction:
 */
#define GAUSSIAN_ELLIPSOID              0
#define POLYNOMIAL_GAUSSIAN_ELLIPSOID   1
#define POLYNOMIAL_GAUSSIAN_ORDER       2.0

/*
 * Choose streamfunction style from those listed above:
 */
#define SPOT_PERT  GAUSSIAN_ELLIPSOID

void add_spots(char    *spots_file,
               double  *pert)
{
  register int
    K,J,I,
    kk,jj,
    ispot;
  int
    nspots = 0;
  double
    rr,xspot,yspot,zspot,
    lon_width,lon_half_width,
    pressure;
  double
    *lonspot,*latspot,*pspot,
    *aspot,*bspot,*cspot_up,*cspot_down,*ampspot;

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="add_spots";

  if (strcmp(spots_file,"none") == 0) {
    /* Return if there is nothing to do: */
    return;
  }

  /* Read vortex description file: */
  lon_width      = grid.globe_lontop-grid.globe_lonbot;
  lon_half_width = .5*lon_width;

  nspots = number_objects_in_file(spots_file);

  if (nspots == 1) {
    fprintf(stdout,"Generating spot..."); fflush(stdout);
  }
  else {
    fprintf(stdout,"Generating spots..."); fflush(stdout);
  }

  /* 
   * Allocate memory:
   */
  lonspot      = dvector(0,nspots-1,dbmsname);
  latspot      = dvector(0,nspots-1,dbmsname);
  pspot        = dvector(0,nspots-1,dbmsname);
  aspot        = dvector(0,nspots-1,dbmsname);
  bspot        = dvector(0,nspots-1,dbmsname);
  cspot_up     = dvector(0,nspots-1,dbmsname);
  cspot_down   = dvector(0,nspots-1,dbmsname);
  ampspot      = dvector(0,nspots-1,dbmsname);

  /* 
   * Read in vortex information: 
   */
  read_spots_file(spots_file,ampspot,lonspot,latspot,pspot,
                  aspot,bspot,cspot_up,cspot_down,ADJUST_AMPLITUDE);

  /* 
   * Calculate streamfunction for vortices, and store in PERT.
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        pressure = P2(K,J,I);
        for (ispot = 0; ispot < nspots; ispot++) {
          /* Account for periodicity in x-direction: */
          xspot  = (grid.lon[2*I+1]-lonspot[ispot]); 
          if (xspot > lon_half_width) { 
            xspot -= lon_width; 
          } else if (xspot < -lon_half_width) { 
            xspot += lon_width; 
          } 
          xspot /= aspot[ispot]; 

          yspot  = (grid.lat[2*J+1]-latspot[ispot])/bspot[ispot];

	  rr     = xspot*xspot+yspot*yspot;
          if (pressure <= pspot[ispot]){
            zspot = -log(pressure/pspot[ispot])/cspot_up[ispot];
	  }
          else{
            zspot =  log(pressure/pspot[ispot])/cspot_down[ispot];
          }
          rr += zspot*zspot; 
      
#if (SPOT_PERT == GAUSSIAN_ELLIPSOID)
          PERT(K,J,I) += ampspot[ispot]*exp(-rr);
#elif (SPOT_PERT == POLYNOMIAL_GAUSSIAN_ELLIPSOID)
          rr = sqrt(rr);
          if ( rr <= pow(2.0,1./POLYNOMIAL_GAUSSIAN_ORDER) ) {
            rr = pow(rr,POLYNOMIAL_GAUSSIAN_ORDER);
            PERT(K,J,I) += ampspot[ispot]*(1.+(1.-rr)*exp(-rr+2.))/(1.+exp(2.));
          }
#endif
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(pert,THREEDIM);

  /* Free allocated memory: */
  free_dvector(ampspot,   0,nspots-1,dbmsname);
  free_dvector(cspot_down,0,nspots-1,dbmsname);
  free_dvector(cspot_up,  0,nspots-1,dbmsname);
  free_dvector(bspot,     0,nspots-1,dbmsname);
  free_dvector(aspot,     0,nspots-1,dbmsname);
  free_dvector(pspot,     0,nspots-1,dbmsname);
  free_dvector(latspot,   0,nspots-1,dbmsname);
  free_dvector(lonspot,   0,nspots-1,dbmsname);

  return;
}

#undef GAUSSIAN_ELLIPSOID
#undef POLYNOMIAL_GAUSSIAN_ELLIPSOID
#undef POLYNOMIAL_GAUSSIAN_ORDER

/*======================= end of add_spots() ================================*/

/*======================= add_noise() =======================================*/

/*
* Add noise to the wind field.
* Kunio Sayanagi, 11-07-07
* This function perturbs the wind velocity field.
* It takes the same input file as add_spots, and is based on add_spots
*
* NOTE: Not MPI ready.
*/
#undef  PERT_U
#define PERT_U(k,j,i) pert_u[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  PERT_V
#define PERT_V(k,j,i) pert_v[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

void add_noise(char       *noise_file)
{
 register int
   K,J,I,jj,ispot,i;
 int
   nspots=0;
 double 
   rr,xspot,yspot,zspot,fspot,
   pressure,
   *ampspot,
   *lonspot,*latspot,*pspot,
   *aspot,*bspot,*cspot_up,*cspot_down,
   *pert_u,
   *pert_v,
    factor,
   lon_width,
   lon_half_width;
 char
   *char_pt,
   buffer[FILE_STR];
 FILE
   *spots;
 /*
  * The following are part of DEBUG_MILESTONE(.) statements:
  */
 int
   idbms=0;
 static char
   dbmsname[]="add_spots";

 int
  kount1, kount2, kount3;

 if (strcmp(noise_file,"none") == 0) {
   /* Return if there is nothing to do: */
   return;
 }

 /* Read vortex description file: */
 lon_width      = grid.globe_lontop - grid.globe_lonbot;
 lon_half_width = 0.5 * lon_width;
 nspots = number_objects_in_file(noise_file);

 fprintf(stdout,"add_noise: Perturbing the wind field ...\n"); fflush(stdout);

 /*
  * Allocate memory:
  */
 lonspot    = dvector(0,nspots-1, dbmsname);
 latspot    = dvector(0,nspots-1, dbmsname);
 pspot      = dvector(0,nspots-1, dbmsname);
 aspot      = dvector(0,nspots-1, dbmsname);
 bspot      = dvector(0,nspots-1, dbmsname);
 cspot_up   = dvector(0,nspots-1, dbmsname);
 cspot_down = dvector(0,nspots-1, dbmsname);
 ampspot    = dvector(0,nspots-1, dbmsname);
 pert_u     = dvector(0,Nelem3d-1,dbmsname);
 pert_v     = dvector(0,Nelem3d-1,dbmsname);

 /*
  * Read in vortex information:
  */
 read_spots_file(noise_file,
                 ampspot,
                 lonspot,
                 latspot,
                 pspot,
                 aspot,
                 bspot,
                 cspot_up,
                 cspot_down,
                 DONT_ADJUST_AMPLITUDE);

 /* Clear PERT memory. */
 memset(pert_u,0,Nelem3d*sizeof(double));
 memset(pert_v,0,Nelem3d*sizeof(double));

 /*
  * Calculate the perturbation fields, and store
  * in PERT memory.
  */
 for (K = KLO; K <= KHI; K++) {
   for (J = JLO; J <= JHI; J++) {
     for (I = ILO; I <= IHI; I++) {
       pressure = P3(K,J,I);
       for (ispot = 0; ispot < nspots; ispot++) {
         /* Account for periodicity in x-direction: */
         xspot  = (grid.lon[2*I+1]-lonspot[ispot]);
         if (xspot > lon_half_width) {
           xspot -= lon_width;
         } else if (xspot < -lon_half_width) {
           xspot += lon_width;
         }
         xspot /= aspot[ispot];

         yspot  = (grid.lat[2*J+1]-latspot[ispot])/bspot[ispot];

         rr     = xspot*xspot+yspot*yspot;
         if (pressure <= pspot[ispot]){
           zspot = -log(pressure/pspot[ispot])/cspot_up[ispot];
         }
         else{
           zspot =  log(pressure/pspot[ispot])/cspot_down[ispot];
         }
         rr += zspot*zspot;

         /*
          * Kind of kludgie, but when ispot = even, perturb U and when odd, perturb V
          */
         if ((ispot%2) == 0) {
           PERT_U(K,J,I) += ampspot[ispot]*exp(-rr);
         }
         else {
           PERT_V(K,J,I) += ampspot[ispot]*exp(-rr);
         }
       }
     }
   }
 }

 /* Need to apply bc_lateral() here. */
 bc_lateral(pert_u,THREEDIM);
 bc_lateral(pert_v,THREEDIM);

 /*
  * Modify U and V.
  *
  * NOTE: Not MPI ready.
  */
 printf("add_noise: Modifing U and V ... \n");
 for (K = KLO; K <= KHI; K++) {
   for (J = JLO; J <= JHI; J++) {
     jj = 2*J+1;
     for (I = ILO; I <= IHI; I++) {
       U(grid.it_uv,K,J,I) += PERT_U(K,J,I);
     }
   }
   /* Need to apply bc_lateral() here. */
   bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);

   for (J = JFIRST; J <= JHI; J++) {
     jj = 2*J;
     for (I = ILO; I <= IHI; I++) {
       V(grid.it_uv,K,J,I) += PERT_V(K,J,I);
     }
   }

   /* Need to apply bc_lateral() here. */
   bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

 } /* (end loop over K) */


 fprintf(stdout,"done.\n"); fflush(stdout);

 /* Free allocated memory: */
 free_dvector(pert_u,    0,Nelem3d-1,dbmsname);
 free_dvector(pert_v,    0,Nelem3d-1,dbmsname);
 free_dvector(ampspot,   0,nspots-1, dbmsname);
 free_dvector(cspot_down,0,nspots-1, dbmsname);
 free_dvector(cspot_up,  0,nspots-1, dbmsname);
 free_dvector(bspot,     0,nspots-1, dbmsname);
 free_dvector(aspot,     0,nspots-1, dbmsname);
 free_dvector(pspot,     0,nspots-1, dbmsname);
 free_dvector(latspot,   0,nspots-1, dbmsname);
 free_dvector(lonspot,   0,nspots-1, dbmsname);

 return;
}

/*======================= end of add_noise() ================================*/

/*======================= add_waves() =======================================*/

/*
 * Add sinusoidal perturbation based on parameters in waves_file.
 * The vertical size c is in scale heights [sc_ht].
 * The amplitude [m2s-2] refers to the perturbation streamfunction, i.e., 
 * PERT = MONT for theta coordinates or PERT = PHI for pressure coordinates.
 *
 * Example waves_file:
 * ------------------------------------------------------------------
 * Number of waves: 2
 * lat[deg]   amp[m2s-2]  wn[]  p[hPa]  c[sc_ht]   fwhm[deg]
 *   30.0      0.02       1.0   700.0      1.5        1.50
 *   40.0      0.10       5.0   700.0      1.5        3.00
 * ------------------------------------------------------------------
 *
 * NOTE: wn[] is domain wavenumber, not planetary wavenumber.
 *
 * Tim Dowling and Raul Morales-Juberias, June 2022
 *
 * NOTE: Not MPI ready. 
 */

void add_waves(char   *waves_file,
               double *pert)
{
  int
    K,J,I,
    kkbot,jj,
    iwave;
  int
    nwaves = 0;
  long
    seed = -1;
  double
    xwave,ywave,zwave,kwave,
    rlt,rln,rr,
    phase_offset,
    coef;
  double
   *latwave,
   *ampwave,
   *wnwave,
   *pwave,
   *cwave,
   *fwhmwave;

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="add_waves";

  if (strcmp(waves_file,"none") == 0) {
    /* Return if there is nothing to do: */
    return;
  }
		
  /* 
   * Read number of waves:
   */
  nwaves = number_objects_in_file(waves_file);

  if (nwaves == 1) {
    fprintf(stdout,"Generating wave..."); fflush(stdout);
  }
  else {
    fprintf(stdout,"Generating waves..."); fflush(stdout);
  }

  /* 
   * Allocate memory:
   */
  latwave  = dvector(0,nwaves-1,dbmsname);
  ampwave  = dvector(0,nwaves-1,dbmsname);
  wnwave   = dvector(0,nwaves-1,dbmsname);
  pwave    = dvector(0,nwaves-1,dbmsname);
  cwave    = dvector(0,nwaves-1,dbmsname);
  fwhmwave = dvector(0,nwaves-1,dbmsname);

  /*
   * Read wave information:
   */
  read_waves_file(waves_file,latwave,ampwave,wnwave,pwave,cwave,fwhmwave);

  /* 
   * Calculate streamfunction for waves, and store in PERT.
   */

  if (strcmp(grid.geometry,"globe") != 0) {
    sprintf(Message,"grid.geometry = %s not yet implemented",grid.geometry);
    epic_error(dbmsname,Message);
  }

  /* Set kk to be kkbot = 2*grid.nk */
  kkbot = 2*grid.nk;

  /* Set meridional radius of curvature to be mid-channel value; note: nj = 2*nj/2 */
  rlt = grid.rlt[kkbot][grid.nj];

  for (iwave = 0; iwave < nwaves; iwave++) {
    phase_offset = 2.*M_PI*random_number(&seed);

    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;

      /* Normalized y coordinate, centered on wave */
      ywave  = rlt*(grid.lat[jj]-latwave[iwave])*DEG;
      ywave /= fwhmwave[iwave];

      /* Zonal radius of curvature; varies with latitude */
      rln = grid.rln[kkbot][jj];

      /* Wavenumber k = 2pi/wavelength [1/m] */
      kwave = 2.*M_PI/( rln*(grid.globe_lontop-grid.globe_lonbot)*DEG/wnwave[iwave] );

      for (I = ILO; I <= IHI; I++) {
        /* x coordinate [m], origin at west end of channel */
        xwave = rln*(grid.lon[2*I+1]-grid.globe_lonbot)*DEG;
        for (K = KLO; K < KHI; K++) {
          /* NOTE: bottom layer is not perturbed. */
          /* Normalized z coordinate, centered on wave */
          zwave        = -log(P2(K,J,I)/pwave[iwave])/cwave[iwave];
          rr           = ywave*ywave+zwave*zwave;
          PERT(K,J,I) += ampwave[iwave]*exp(-rr)*sin(kwave*xwave+phase_offset);
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(pert,THREEDIM);
	
  /* 
   * Free allocated memory:
   */
  free_dvector(latwave, 0,nwaves-1,dbmsname);
  free_dvector(ampwave, 0,nwaves-1,dbmsname);
  free_dvector(wnwave,  0,nwaves-1,dbmsname);
  free_dvector(pwave,   0,nwaves-1,dbmsname);
  free_dvector(cwave,   0,nwaves-1,dbmsname);
  free_dvector(fwhmwave,0,nwaves-1,dbmsname);
	
  return;
}

/*======================= end of add_waves() ================================*/

/*======================= modify_progs_from_pert() ==========================*/

/*
 * Use gradient-wind balance to modify the prognostic variables
 * given a streamfunction perturbation.
 */

#undef  UG
#define UG(k,j,i) ug[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  VG
#define VG(k,j,i) vg[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  ZETAG
#define ZETAG(j,i) zetag[i+(j)*Iadim-Shift2d]

#undef  KING
#define KING(j,i) king[i+(j)*Iadim-Shift2d]

#undef  DTEMP
#define DTEMP(k,j,i) dtemp[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  LAT_MIN
#define LAT_MIN 5.0

void modify_progs_from_pert(double *pert)
{
  double
   *ug,*vg,
   *zetag,*king,
   *dtemp,
   *p_hybrid,*theta_hybrid,*theta_sigma;
  double
    fpara,rgas,sigma,
    pbot,gsg,xi2,xi4,
    theta_ortho,theta_para;
  int
    K,J,I,
    kk,jj;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="modify_progs_from_pert";

  /* Allocate memory */
  ug    = dvector(0,Nelem3d-1,dbmsname);
  vg    = dvector(0,Nelem3d-1,dbmsname);
  dtemp = dvector(0,Nelem3d-1,dbmsname);

  zetag = dvector(0,Nelem2d-1,dbmsname);
  king  = dvector(0,Nelem2d-1,dbmsname);

  p_hybrid     = dvector(0,KHI,dbmsname);
  theta_hybrid = dvector(0,KHI,dbmsname);
  theta_sigma  = dvector(0,KHI,dbmsname);

  /*
   * Estimate the geostrophic wind (UG,VG) of the perturbed system.
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = JLOPAD; J <= JHIPADPV; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        UG(K,J,I) = U(grid.it_uv,K,J,I);
        VG(K,J,I) = V(grid.it_uv,K,J,I);
      }
    }
  }
  for (J = JLO+1; J <= JHI-1; J++) {
    jj = 2*J+1;
    if (fabs(grid.lat[jj]) > LAT_MIN) {
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K <= KHI; K++) {
          UG(K,J,I) -= (PERT(K,J+1,I-1)+PERT(K,J+1,I)-PERT(K,J-1,I-1)-PERT(K,J-1,I))*.25*grid.n[2*K][jj]/grid.f[jj];
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(ug,THREEDIM);

  for (J = JFIRST; J <= JHI; J++) {
    jj = 2*J;
    if (fabs(grid.lat[jj]) > LAT_MIN) {
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K <= KHI; K++) {
          VG(K,J,I) += (PERT(K,J,I+1)+PERT(K,J-1,I+1)-PERT(K,J,I-1)-PERT(K,J-1,I-1))*.25*grid.m[2*K][jj]/grid.f[jj];
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(vg,THREEDIM);

  /* 
   * Modify P, THETA, FPARA as appropriate.
   */

  switch(grid.coord_type) {
    case COORD_HYBRID:
      /*
       * For the hybrid coordinate, this is a challenging task, and the current approach is
       * approximate (see Dowling's notes, 12/2/10).
       *
       * Note that hydrostatic balance in isentropic and isobaric coordinates can be written:
       *
       *      ____isentropic____         ____isobaric____
       *
       *     dM/dlog(theta) = cp T       dPhi/dlog(p) = - R T
       *
       * These are quite similar in form.  Let <T> be the unperturbed temperature.  Then,
       *
       *     T = <T> + dPERT/dxi ,
       *
       * where PERT is the perturbation to either the Montgomery potential or the geopotential, as set above,
       * and xi is either cp log(theta) or -R log(p), respectively.  Our strategy is to define xi to be a blend,
       * using the transition function g(sigma):
       *
       *     xi = g cp log(theta) + (1.-g)(-R log(p)) ,
       *
       * where sigma, theta, and p are all the unperturbed values.
       */
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          pbot = P3(KHI,J,I);
          for (K = KLO; K < KHI; K++) {
            /*
             * Form the independent variable, xi, by blending from the pressure region to the theta region.
             */
            rgas         = R_GAS/avg_molar_mass(2*K+1,J,I);

            sigma        = get_sigma(pbot,P2(K,J,I));
            gsg          = (double)g_sigma(sigma);
            xi2          = gsg*planet->cp*log(THETA2(K,J,I))+(1.-gsg)*(-rgas*log(P2(K,J,I)));

            sigma        = get_sigma(pbot,P2(K+1,J,I));
            gsg          = (double)g_sigma(sigma);
            xi4          = gsg*planet->cp*log(THETA2(K+1,J,I))+(1.-gsg)*(-rgas*log(P2(K+1,J,I)));

            /* 
             * DTEMP(K,J,I) is on the p3-grid.
             */
            DTEMP(K,J,I) = (PERT(K,J,I)-PERT(K+1,J,I))/(xi2-xi4);
          }
        }
      }

      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          pbot = P3(KHI,J,I);
          for (K = KLO; K < grid.k_sigma; K++) {
            T3(K,J,I) += DTEMP(K,J,I);

            /*
             * Estimate para-hydrogen fraction with equilibrium value.
             */
            fpara = return_fpe(T3(K,J,I));
            if (var.fpara.on) {
              FPARA(K,J,I) = fpara;
            }

            /*
             * NOTE: Assuming the density doesn't change.
             */
            P3(K,J,I) = p_from_t_rho_mu(T3(K,J,I),RHO3(K,J,I),avg_molar_mass(2*K+1,J,I));

            /*
             * Use diagnostic value for THETA.
             */
            sigma        = get_sigma(pbot,P3(K,J,I));
            gsg          = (double)g_sigma(sigma);
            THETA(K,J,I) = (double)(((double)grid.sigmatheta[2*K+1]-(double)f_sigma(sigma))/gsg);
          }
          for (K = grid.k_sigma; K < KHI; K++) {
            T3(K,J,I) += DTEMP(K,J,I);

            /*
             * Estimate para-hydrogen fraction with equilibrium value.
             */
            fpara = return_fpe(T3(K,J,I));
            if (var.fpara.on) {
              FPARA(K,J,I) = fpara;
            }

            /*
             * NOTE: Assuming the density doesn't change.
             */
            P3(K,J,I) = p_from_t_rho_mu(T3(K,J,I),RHO3(K,J,I),avg_molar_mass(2*K+1,J,I));

            /*
             * Use thermodynamic value for THETA.
             */
            THETA(K,J,I) = return_theta(fpara,P3(K,J,I),T3(K,J,I),&theta_ortho,&theta_para);
          }
        }
      }
    break;
    case COORD_ISOBARIC:
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          for (K = KLO; K < KHI; K++) {
            rgas          = R_GAS/avg_molar_mass(2*K+1,J,I);

            DTEMP(K,J,I)  = (PERT(K,J,I)-PERT(K+1,J,I))/(rgas*log(P2(K+1,J,I)/P2(K,J,I)));
            T3(K,J,I)    += DTEMP(K,J,I);

            /*
             * Estimate para-hydrogen fraction with equilibrium value.
             */
            fpara = return_fpe(T3(K,J,I));
            if (var.fpara.on) {
              FPARA(K,J,I) = return_fpe(T3(K,J,I));
            }

            THETA(K,J,I)  = return_theta(fpara,P3(K,J,I),T3(K,J,I),&theta_ortho,&theta_para);
          }
        }
      }
    break;
    case COORD_ISENTROPIC:
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          for (K = KLO; K < KHI; K++) {
            DTEMP(K,J,I)  = (PERT(K,J,I)-PERT(K+1,J,I))/(planet->cp*log(THETA(K,J,I)/THETA(K+1,J,I)));
            T3(K,J,I)    += DTEMP(K,J,I);

            /*
             * Estimate para-hydrogen fraction with equilibrium value.
             */
            fpara = return_fpe(T3(K,J,I));
            if (var.fpara.on) {
              FPARA(K,J,I) = fpara;
            }

            /*
             * NOTE: Assuming the density doesn't change.
             */
            P3(K,J,I) = p_from_t_rho_mu(T3(K,J,I),RHO3(K,J,I),avg_molar_mass(2*K+1,J,I));
          }
        }
      }
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  for (K = KLO; K < KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      if (fabs(grid.lat[jj]) > LAT_MIN) {
        for (I = ILO; I <= IHI; I++) {
          U(grid.it_uv,K,J,I) = UG(K,J,I);
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;
      if (fabs(grid.lat[jj]) > LAT_MIN) {
        for (I = ILO; I <= IHI; I++) {
          V(grid.it_uv,K,J,I) = VG(K,J,I);
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);
  bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

  /*
   * The gradient-balance correction is based on the paper:
   *     McIntyre and Roulstone, 2002, Large-Scale Atmosphere-Ocean Dynamics. II Geometric Methods
   *         and Models. Cambridge Univ. Press, Ch. 8.
   */
  for (K = KLO; K < KHI; K++) {
    kk = 2*K;
    /*
     * Calculate geostrophic relative vorticity and kinetic energy per mass.
     */
    vorticity(ON_SIGMATHETA,RELATIVE,2*K,ug+(K-Kshift)*Nelem2d,vg+(K-Kshift)*Nelem2d,NULL,zetag);
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        KING(J,I) = get_kin(ug+(K-Kshift)*Nelem2d,vg+(K-Kshift)*Nelem2d,kk,J,I);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(king,TWODIM);
  
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      if (fabs(grid.lat[jj]) > LAT_MIN) {
        for (I = ILO; I <= IHI; I++) {
          UG(K,J,I) += (-.5*(ZETAG(J,I)+ZETAG(J+1,I))*UG(K,J,I)
                        -.25*grid.n[kk][jj]*(KING(J+1,I)+KING(J+1,I-1)-KING(J-1,I)-KING(J-1,I-1)))/grid.f[jj];
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;
      if (fabs(grid.lat[jj]) > LAT_MIN) {
        for (I = ILO; I <= IHI; I++) {
          VG(K,J,I) += (-.5*(ZETAG(J,I)+ZETAG(J,I+1))*VG(K,J,I)
                        +.25*grid.m[kk][jj]*(KING(J,I+1)+KING(J-1,I+1)-KING(J,I-1)-KING(J-1,I-1)))/grid.f[jj];
        }
      }
    }
  }

  /* Free allocated memory. */
  free_dvector(ug,   0,Nelem3d-1,dbmsname);
  free_dvector(vg,   0,Nelem3d-1,dbmsname);
  free_dvector(dtemp,0,Nelem3d-1,dbmsname);

  free_dvector(zetag,0,Nelem2d-1,dbmsname);
  free_dvector(king, 0,Nelem2d-1,dbmsname);

  free_dvector(p_hybrid,    0,KHI,dbmsname);
  free_dvector(theta_hybrid,0,KHI,dbmsname);
  free_dvector(theta_sigma, 0,KHI,dbmsname);

  fprintf(stdout,"done \n"); fflush(stdout);

  return;
}

/*======================= end of modify_progs_from_pert() ===================*/

/*======================= add_heat() ========================================*/

/*
 * Use the same file format as add_spots() to add a gaussian perturbation
 * to potential temperature, with ampspot interpreted as [K].
 *
 * Example heat.dat file:
 * ---------------------------------------------------------------------------------------
 *  Number of spots: 2
 *  lon[deg] lat[deg] press[hPa]  a[deg] b[deg] c_up[scale_hts] c_down[scale_hts] amp[K]
 *   30.      -33.     680.       3.0    2.5        2.5               3.0          -1.
 *   60.      -33.5    680.       3.0    2.5        2.5               3.0           2.
 * ---------------------------------------------------------------------------------------
 */

void add_heat(char *heat_file)
{
  register int
    K,J,I,
    kk,jj,
    ispot;
  int
    nspots = 0;
  double
    rr,xspot,yspot,zspot,
    lon_width,lon_half_width,
    pressure;
  double
    *lonspot,*latspot,*pspot,
    *aspot,*bspot,*cspot_up,*cspot_down,*ampspot;

  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="add_heat";

  if (strcmp(heat_file,"none") == 0) {
    /* Return if there is nothing to do: */
    return;
  }

  /*
   * Error if called for the pure isentropic-coordinate model.
   */
  if (grid.coord_type == COORD_ISENTROPIC) {
    sprintf(Message,"not implemented for grid.coord_type == COORD_ISENTROPIC");
    epic_error(dbmsname,Message);
  }

  /* Read perturbation description file: */
  lon_width      = grid.globe_lontop-grid.globe_lonbot;
  lon_half_width = .5*lon_width;

  nspots = number_objects_in_file(heat_file);

  if (nspots == 1) {
    fprintf(stdout,"Generating heat spot..."); fflush(stdout);
  }
  else {
    fprintf(stdout,"Generating heat spots..."); fflush(stdout);
  }

  /* 
   * Allocate memory:
   */
  lonspot      = dvector(0,nspots-1,dbmsname);
  latspot      = dvector(0,nspots-1,dbmsname);
  pspot        = dvector(0,nspots-1,dbmsname);
  aspot        = dvector(0,nspots-1,dbmsname);
  bspot        = dvector(0,nspots-1,dbmsname);
  cspot_up     = dvector(0,nspots-1,dbmsname);
  cspot_down   = dvector(0,nspots-1,dbmsname);
  ampspot      = dvector(0,nspots-1,dbmsname);

  /* 
   * Read in heat spot information: 
   */
  read_spots_file(heat_file,ampspot,lonspot,latspot,pspot,
                  aspot,bspot,cspot_up,cspot_down,DONT_ADJUST_AMPLITUDE);

  /*
   * Apply perturbation to THETA in the non-isentropic-coordinate part of the model.
   */
  for (K = grid.k_sigma; K <= KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        pressure = P3(K,J,I);
        for (ispot = 0; ispot < nspots; ispot++) {
           /* Account for periodicity in x-direction: */
           xspot  = (grid.lon[2*I+1]-lonspot[ispot]);
          if (xspot > lon_half_width) {
            xspot -= lon_width;
          } else if (xspot < -lon_half_width) {
            xspot += lon_width;
          }
          xspot /= aspot[ispot];

          yspot  = (grid.lat[2*J+1]-latspot[ispot])/bspot[ispot];

          rr = xspot*xspot+yspot*yspot;
          if (pressure <= pspot[ispot]) {
            zspot = -log(pressure/pspot[ispot])/cspot_up[ispot];
          }
          else{
            zspot =  log(pressure/pspot[ispot])/cspot_down[ispot];
          }
          rr += zspot*zspot;

          THETA(K,J,I) += ampspot[ispot]*exp(-rr);
        }
      }
    }
  }

  /*
   * NOTE: Only the prognostic variables are important in the output from change,
   *       but in case any diagnostic variables closely associated with THETA may
   *       affect the prognostic variables after add_heat():
   *
   * Update P2, THETA2, etc.
   */
  set_p2_etc(UPDATE_THETA);

  /* Free allocated memory: */
  free_dvector(ampspot,   0,nspots-1,dbmsname);
  free_dvector(cspot_down,0,nspots-1,dbmsname);
  free_dvector(cspot_up,  0,nspots-1,dbmsname);
  free_dvector(bspot,     0,nspots-1,dbmsname);
  free_dvector(aspot,     0,nspots-1,dbmsname);
  free_dvector(pspot,     0,nspots-1,dbmsname);
  free_dvector(latspot,   0,nspots-1,dbmsname);
  free_dvector(lonspot,   0,nspots-1,dbmsname);
 
  fprintf(stdout,"done \n"); fflush(stdout);

  return;
}
/*======================= end of add_heat() =================================*/

/*====================== number_objects_in_file() ===========================*/

/*
 * Retrieve the integer after the first colon, ':', which
 * is interpreted as the number of objects (spots, waves, etc).
 * in the given parameter file.
 */

int number_objects_in_file(char *objects_file)
{
  int 
    n_objects = 0;
  char
    *char_pt,
    buffer[FILE_STR];
  FILE
    *objects;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="number_objects_in_file";

  objects = fopen(objects_file,"r");
  if (!objects) {
    sprintf(Message,"unable to open %s \n",objects_file);
    epic_error(dbmsname,Message);
  }
  while (n_objects == 0) {
    fgets(buffer,FILE_STR,objects);
    char_pt = strchr(buffer,':');
    if (char_pt) {
      sscanf(char_pt+1,"%d",&n_objects);
    }
  }
  fclose(objects);

  return n_objects;
}

/*====================== end of number_objects_in_file() ====================*/

/*====================== read_spots_file() ==================================*/

void read_spots_file(char   *spots_file,
                     double *ampspot,
                     double *lonspot,
                     double *latspot,
                     double *pspot,
                     double *aspot,
                     double *bspot,
                     double *cspot_up,
                     double *cspot_down,
                     int     adjust_amplitude)
{
  register int 
    ispot;
  int
    nspots=0;
  double
    fspot,
    factor;
  char
    *char_pt,
    buffer[FILE_STR];
  FILE
    *spots;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_spots_file";

  spots = fopen(spots_file,"r");
  while (nspots == 0) {
    fgets(buffer,FILE_STR,spots);
    char_pt = strchr(buffer,':');
    if (char_pt) {
      sscanf(char_pt+1,"%d",&nspots);
    }
  }

  fgets(buffer,FILE_STR,spots);
  for (ispot = 0; ispot < nspots; ispot++) {
    fscanf(spots,"%lf %lf %lf %lf %lf %lf %lf %lf",
           lonspot+ispot,latspot+ispot,pspot+ispot,
           aspot+ispot,bspot+ispot,cspot_up+ispot,cspot_down+ispot,ampspot+ispot);

    /* Convert pspot from hPa to Pa: */
    pspot[ispot] *= 100.;

    /* 
     * Convert ampspot to amp for mont.  The parameter 'factor' 
     * should make the maximum spot velocity close to the input 
     * ampspot[m/s] for a gaussian spot perturbation streamfunction.
     */
    if (adjust_amplitude == ADJUST_AMPLITUDE) {
      fspot           = 2.*planet->omega_sidereal*sin(latspot[ispot]*DEG);
      factor          = 1.166;
      ampspot[ispot] *= factor*bspot[ispot]*DEG*planet->re*fabs(fspot);
    }

    fgets(buffer,FILE_STR,spots);
  }
  fclose(spots);

  return;

}

/*====================== end of read_spots_file() ===========================*/

/*====================== read_waves_file() ==================================*/

void read_waves_file(char   *waves_file,
                     double *latwave,
                     double *ampwave,
                     double *wnwave,
                     double *pwave,
                     double *cwave,
                     double *fwhmwave)
{
  register int 
    iwave;
  int
    kkbot,
    nwaves = 0;
  double
    coef,rlt;
  char
    *char_pt,
    buffer[FILE_STR];
  FILE
    *waves;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_waves_file";

  /* 
   * Coefficient for using FWHM (full-width at half-maximum) to specify the width of a 1D gaussian.
   */
  coef = 2.*sqrt(log(2.));

  /* Set kk to be kkbot = 2*grid.nk */
  kkbot = 2*grid.nk;

  /* Set meridional radius of curvature to be mid-channel value; note: nj = 2*nj/2 */
  rlt = grid.rlt[kkbot][grid.nj];

  waves = fopen(waves_file,"r");
  while (nwaves == 0) {
    fgets(buffer,FILE_STR,waves);
    char_pt = strchr(buffer,':');
    if (char_pt) {
      sscanf(char_pt+1,"%d",&nwaves);
    }
  }

  fgets(buffer,FILE_STR,waves);
  for (iwave = 0; iwave < nwaves; iwave++) {
    fscanf(waves,"%lf %lf %lf %lf %lf %lf",
                 latwave+iwave,ampwave+iwave,wnwave+iwave,pwave+iwave,cwave+iwave,fwhmwave+iwave);
		
    /* Convert pwave from hPa (mbar) to Pa: */
    pwave[iwave] *= 100.;

    /* Convert FWHM for use in 1D gaussian. */
    fwhmwave[iwave] *= rlt*DEG/coef;

    fgets(buffer,FILE_STR,waves);
  }

  fclose(waves);

  return;
}

/*====================== end of read_waves_file() ===========================*/

/*======================= read_defaults() ===================================*/

void read_defaults(change_defaultspec *def, char *defaults_file) 
{
  int
    nc_id,
    nc_err,
    index;
  char
    min_element[4];
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts    =0,
    num_progs =0;
  nc_type
    the_nc_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_defaults";

  nc_err = lookup_netcdf(defaults_file,
                         &nc_id,&ngatts,&gattname,&num_progs,&varname);
  if (nc_err == NC_NOERR) {
    READI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
    READC(def->infile,def_infile,N_STR);
    READC(def->outfile,def_outfile,N_STR);
    READC(def->openmars_infile,def_openmars_infile,N_STR);
    READC(def->emars_infile,def_emars_infile,N_STR);
    READC(def->extract_str,def_extract_str,N_STR);

    READI(&def->moist_convection,def_moist_convection,1);
    READI(&def->cloud_microphysics,def_cloud_microphysics,1);
    READI(&def->reinit_cloud,def_reinit_cloud,1);
    
    READI(&def->grid_relax_vapor,grid_relax_vapor,1);
    READD(&def->grid_relax_vapor_timescale,grid_relax_vapor_timescale,1);
    READD(def->rh_max,def_rh_max,LAST_SPECIES+1);

    READI(&def->max_mc_it,grid_max_mc_it,1);
    READD(&def->tau_relax,grid_tau_relax,1);
    
    READD(&def->heat_rate,grid_heat_rate,1);
    READD(&def->cool_rate,grid_cool_rate,1);
    READD(&def->heat_top_pressure,grid_heat_top_pressure,1);
    READD(&def->cool_bot_pressure,grid_cool_bot_pressure,1);
  }
  else {
    /*
     * If the file is not readable, use standard defaults.
     */
    def->uv_timestep_scheme = 0;

    strcpy(def->infile,"epic.nc");
    strcpy(def->outfile,"epic.nc");
    strcpy(def->openmars_infile,"");
    strcpy(def->emars_infile,"");
    strcpy(def->extract_str,"");

    def->radiation_scheme               = 0;

    def->moist_convection = 0;
    def->cloud_microphysics = 0;
    def->reinit_cloud = 0;
    def->max_mc_it = 0;
    def->tau_relax = 0.;
    
    def->grid_relax_vapor = OFF;
    def->grid_relax_vapor_timescale = 5;
    
    def->heat_top_pressure  = 20000.; // hPa
    def->cool_bot_pressure  = 50.;    // hPa
    def->heat_rate = 0.008; // K/day
    def->cool_rate = 0.01;  // K/day
  }

  return;
}

/*======================= end of read_defaults() ============================*/

/*======================= write_defaults() ==================================*/

void write_defaults(change_defaultspec *def)
{
  int
    nc_id,
    nc_err;
  nc_type
    the_nc_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="write_defaults";

  nc_err = nc_create("change_defaults.nc",NC_CLOBBER,&nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  WRITEI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
  WRITEC(def->infile,def_infile,N_STR);
  WRITEC(def->outfile,def_outfile,N_STR);
  WRITEC(def->openmars_infile,def_openmars_infile,N_STR);
  WRITEC(def->emars_infile,def_emars_infile,N_STR);
  WRITEC(def->extract_str,def_extract_str,N_STR);
  
  WRITEI(&def->radiation_scheme,def_radiation_scheme,1);

  WRITEI(&def->cloud_microphysics,def_cloud_microphysics,1);
  WRITEI(&def->moist_convection,def_moist_convection,1);
  WRITEI(&def->reinit_cloud,def_reinit_cloud,1);
  WRITED(def->rh_max,def_rh_max,LAST_SPECIES+1);

  WRITEC(def->infile,def_infile,N_STR);
  WRITEC(def->outfile,def_outfile,N_STR);
  WRITEC(def->extract_str,def_extract_str,N_STR);
  
  WRITED(&def->heat_rate,grid_heat_rate,1);
  WRITED(&def->cool_rate,grid_cool_rate,1);
  WRITED(&def->heat_top_pressure,grid_heat_top_pressure,1);
  WRITED(&def->cool_bot_pressure,grid_cool_bot_pressure,1);

  WRITEI(&def->max_mc_it,grid_max_mc_it,1);
  WRITED(&def->tau_relax,grid_tau_relax,1);
  
  WRITEI(&def->grid_relax_vapor,grid_relax_vapor,1);
  WRITED(&def->grid_relax_vapor_timescale,grid_relax_vapor_timescale,1);

  nc_close(nc_id);

  return;
}

/*======================= end of write_defaults() ===========================*/

/* * * * * * * * * * * * end of epic_change.c * * * * * * * * * * * * * * * * */

