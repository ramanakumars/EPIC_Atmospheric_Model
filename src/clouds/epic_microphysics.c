/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2024-2025 Csaba Palotai *A*, Ramanakumar Sankar   *
 * Copyright (C) 2002-2023 Csaba Palotai *A*, Timothy E. Dowling   *
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

/* * * * * * * * * *  epic_microphysics.c  * * * * * * * * * * * * *
 *                                                                 *
 *       Functions governing the hydrological cycle.               *
 *       This file includes the following functions:               *
 *                                                                 *
 *           cloud_microphysics()                                  *
 *           instantaneous_processes()                             *
 *           finite_rate_processes()                               *
 *           terminal_velocity()                                   *
 *           moist_convection()                                    *
 *           RAS_flux()                                            *
 *           calculate_base_pressure()                             *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>


/*========================== cloud_microphysics() ============================*/

  /*
   * Add latent heating to HEAT and transfer mass between the phases of each species,
   * as appropriate.
   */

void cloud_microphysics(void)
{
  register int
    K,J,I,
    is,ip,iq;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
 static int
    warned_once = FALSE,
    initialized = FALSE;
  static char
    dbmsname[]="cloud_microphysics";

#if EPIC_CHECK == TRUE 
  /*
   * Check that the phases are turned on for every species.
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        if (!var.species[is].phase[ip].on) {
          sprintf(Message,"%s is not ON",var.species[is].phase[ip].info[0].name);
          epic_error(dbmsname,Message);
        }
      }
    }
  }
#endif 

  /*
   * Initialize gamma values for terminal velocity calculations.
   */
  if (!initialized) {
    for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
      if (var.species[is].on) {
        GAMMA_ICE( is) = gamma_nr(3.+COEFF_YS(is));
        GAMMA_SNOW(is) = gamma_nr(4.+COEFF_YS(is));
        GAMMA_RAIN(is) = gamma_nr(4.+COEFF_YR(is));
        /* 
         * Make an initial call to Lc() for each species, to initialize
         * the associated enthalpy change data.
         */
        Lc(is,VAPOR,LIQUID,T_triple_pt(is));
      }
    }
    initialized = TRUE;
  }

  for (iq = 0; iq < grid.nq; iq++) {
    restore_mass(grid.is[iq],grid.ip[iq]);
  }

  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      for (K = KLO; K < KHI; K++) {
	for (J = JLO; J <= JHI; J++) {
          /*
           * The Qs at K = 0 and grid.nk are treated as boundary conditions.
           */
          for (I = ILO; I <= IHI; I++) {
	    /* * * * * * * * * * * * * * * * * * * * * * * * * *
	     *   Instantaneous processes: melting, freezing    *
	     * * * * * * * * * * * * * * * * * * * * * * * * * */
            instantaneous_processes(is,K,J,I);

            /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	     *   Finite-rate processes                                               *
	     *   Condensation/sublimation/evaporation, autoconversion, collection    *
	     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
            finite_rate_processes(is,K,J,I);
          }
        }
      }
      /*
       * Need to apply bc_lateral() here.
       */
      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        bc_lateral(var.species[is].phase[ip].q,THREEDIM);
      }
    }
  }

  for (iq = 0; iq < grid.nq; iq++) {
    restore_mass(grid.is[iq],grid.ip[iq]);
  }

  return;
}

/*======================== end of cloud_microphysics() ===========================*/

/*========================= instantaneous_processes() ============================*/

void instantaneous_processes(int is,
                             int K,
                             int J,
                             int I)
{
  double
    psmlti,psmlts;
  double
    q_ice,q_liquid,q_snow,q_rain;  
  boolean
    warm;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="instantaneous_processes";

  psmlti   = 0.;
  psmlts   = 0.;
  q_ice    = Q(is,ICE,   K,J,I);
  q_liquid = Q(is,LIQUID,K,J,I);
  q_snow   = Q(is,SNOW,  K,J,I);
  q_rain   = Q(is,RAIN,  K,J,I);

  /* 
   * Set useful booleans.
   */
  warm = (fcmp(T3(K,J,I),T_triple_pt(is)) >= 0);

  if (warm && q_ice > Q_MIN) { 
    /*  Phase change: melting: ice cloud -> liquid cloud  */
    psmlti = q_ice/DT;
    if (psmlti < Q_MIN) {
      psmlti = 0.;
    }
  }
  else if (!warm && q_liquid > Q_MIN) {
    /* Phase change: freezing: liquid cloud -> ice cloud */
    psmlti = q_liquid/DT;
    if (psmlti < Q_MIN) {
      psmlti = 0.;
    }
    else {
      psmlti *= -1.;
    }  
  }

  if (warm && q_snow > Q_MIN) {
    /*  Phase change: melting: snow -> rain  */
    psmlts = q_snow/DT;
    if (psmlts < Q_MIN) {
      psmlts = 0.;
    }
  }
  else if (!warm && q_rain > Q_MIN) {
    /* Phase change: freezing: rain -> snow */
    psmlts = q_rain/DT;
    if (psmlts < Q_MIN) {
      psmlts = 0.;
    }
    else {
      psmlts *= -1.;
    }	
  }
  Q(is,LIQUID,K,J,I) += psmlti*DT;
  Q(is,   ICE,K,J,I) -= psmlti*DT;
  Q(is,  RAIN,K,J,I) += psmlts*DT;
  Q(is,  SNOW,K,J,I) -= psmlts*DT;
         
  HEAT3(K,J,I) -=(psmlti+psmlts)*Lf(is);

  return;
}
/*==================== end of instantaneous_processes() =======================*/

/*========================= finite_rate_processes() ===========================*/

void finite_rate_processes(int is,
                           int K,
                           int J,
                           int I)
{
  register double
    tmp,rho,rho_inv,t3;
  static double
    c_inv[LAST_SPECIES+1],d_inv[LAST_SPECIES+1],
    qi_crit[LAST_SPECIES+1];
  double
    pcond,pint,pdepi,praut,psaut,prevap,psevap,
    psaci,psacw,pracw,psacr,
    required,d_mass,av_mass,
    subsat,supersat,av_ice,av_liquid,
    q_vapor,q_liquid,q_ice,q_rain,q_snow,
    q_sat,rh;       
  double
    N_I0,Q_I0,Q_C0,E_c,Q_Icrit,E_SI,N_0S,N_I,dynvis,
    lambda_r,lambda_s,A_S,B_S,A_R,B_R,M_Imax,
    Cp,fpara;
  boolean
    saturated,warm,subcritical;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  static int
    initialized = FALSE;
  int
    temp,sp,
    idbms=0;
  static char
    dbmsname[]="finite_rate_processes";

  /*
   * Note: 3/25/2011 CJP
   * Variables subsat,supersat,av_ice,av_liquid are introduced so we wouldn't have to update the mixing ratios 
   * instantaneously after every finite-rate process, but we would still avoid using more than the available amount
   * of moisture by these processes.   
   *
   * These variables might get updated in the last process and never be used again during a function call. 
   * The reason for this is that should additional processes be added to the model, these moisture variables would give
   * the actual available amount and the user would not have to modify the old code.
   */
  
  if (!initialized) {
    for (sp = FIRST_SPECIES; sp <= LAST_SPECIES; sp++) {
      if (var.species[sp].on) {
        qi_crit[sp] = pow(COEFF_C(sp)*COEFF_M(sp)*pow(D_Icrit(is),COEFF_N(sp)),1./(1.-COEFF_D(sp))); 

	c_inv[sp]   = 1./COEFF_C(sp);
	d_inv[sp]   = 1./COEFF_D(sp);
      }
    }
    initialized = TRUE;
  }

  rho       = PDRY3(K,J,I)/(T3(K,J,I)*planet->rgas);
  rho_inv   = 1./rho;
  t3        = T3(K,J,I);
  q_vapor   = Q(is,VAPOR, K,J,I);
  q_liquid  = Q(is,LIQUID,K,J,I);
  q_ice     = Q(is,ICE,   K,J,I);
  q_rain    = Q(is,RAIN,  K,J,I);
  q_snow    = Q(is,SNOW,  K,J,I);
  pcond     = 0.;
  pint      = 0.;
  pdepi     = 0.;
  praut     = 0.;
  psaut     = 0.;
  pracw     = 0.;
  psaci     = 0.;
  prevap    = 0.;
  psevap    = 0.;
  av_ice    = q_ice; 
  av_liquid = q_liquid;

  /*
  * Calculate saturation mixing ratio, supersaturation, subsaturation, and relative humidity (WMO definition).
   */
  q_sat    = (var.species[is].molar_mass*var.species[is].sat_vapor_p(T3(K,J,I)))/
             (molar_mass(HDRY3_INDEX)   *PDRY3(K,J,I)                          );
  supersat = MAX(0.,q_vapor-q_sat);
  subsat   = MAX(0.,q_sat-q_vapor);
  rh       = q_vapor/q_sat;
  
  /* 
   * Set useful booleans.
   */
  saturated   = (fcmp(rh,1.) > 0);
  warm        = (fcmp(t3,T_triple_pt(is))   >= 0);
  subcritical = (fcmp(t3,T_critical_pt(is)) <  0) ;

  if (var.fpara.on) {
    fpara = get_var(FPARA_INDEX,NO_PHASE,grid.it_h,2*K+1,J,I);
  }
  else {
    fpara = return_fpe(t3);
  }
  
  Cp = return_cp(fpara,P3(K,J,I),t3); 

  /* PCOND liquid condensation/evaporation */
  if (warm && subcritical) {    
    if (saturated) {             /* Phase change: condensation */
      pcond = supersat/(1.+ SQR(Lc(is,VAPOR,LIQUID,t3))*q_sat/(Cp*GAS_R(is)*SQR(t3)))/DT;
      if (pcond < Q_MIN) {
        pcond = 0.;
      }
      else { 
        supersat -= pcond*DT;  
      }
    }
    else if (!saturated && av_liquid > Q_MIN) {  
      /* Phase change: evaporation */
      pcond   = subsat/(1.+ SQR(Lc(is,VAPOR,LIQUID,t3))*q_sat/(Cp*GAS_R(is)*SQR(t3)))/DT;
      if (pcond < Q_MIN) {
        pcond = 0.;
      }
      else {
        pcond      = -MIN(pcond,av_liquid/DT);
        subsat    += pcond*DT;     /* pcond is now negative thus the + sign */
        av_liquid += pcond*DT; 
      }
    }
  }  /*====== End of PCOND ======*/

  /* PINIT + PDEPI: ice condensation/sublimation */
  A_S = (SQR(Ls(is))/(GAS_R(is)*t3)-1.)/(conductivity(planet->name,t3)*t3);
  B_S = 1./(mass_diffusivity(is,T3(K,J,I),P3(K,J,I))*q_sat*rho); 
  N_I = COEFF_C(is)*pow(rho*q_ice,COEFF_D(is)); 
  
  if (!warm) {     
    if (saturated) {   /* Phase change: vapor condensation to ice cloud */
      if (q_ice <= Q_MIN) {
        N_I0 = MIN(1.e+8,10000.*exp(0.1*(T_triple_pt(is)-t3)));
	Q_I0 = MIN(pow(c_inv[is]*N_I0,d_inv[is]),0.0015)*rho_inv;
	pint = MIN(Q_I0/DT,supersat/DT);

        if (pint < Q_MIN) {
          pint = 0.;
        }
        else {
	  supersat -= pint*DT;
	}  
      }
      else if (supersat > Q_MIN) {    /* PDEPI: vapor deposition of a small ice crystal */
        tmp   = 1./COEFF_N(is);  
        pdepi = 4.0*pow(COEFF_M(is),-tmp)*(rh-1.)*pow(rho*q_ice,tmp)*pow(N_I,1.-tmp)/(A_S+B_S);
        pdepi     = MIN(pdepi,supersat/DT);
 	if (pdepi < Q_MIN) {
          pdepi = 0.;
        }
        else {
          supersat -= pdepi*DT;
        }
      }
    }
    else if (!saturated && av_ice > Q_MIN) {   
      /* Phase change: ice cloud evaporation */
      tmp   = 1./COEFF_N(is);  
      pdepi = 4.0*pow(COEFF_M(is),-tmp)*(1.-rh)*pow(rho*q_ice,tmp)*pow(N_I,1.-tmp)/(A_S+B_S);
      pdepi = MIN(pdepi,av_ice/DT);
      if (pdepi < Q_MIN) {
        pdepi = 0.;
      }
      subsat -= pdepi*DT;
      av_ice -= pdepi*DT; 
      pdepi   *= -1.;
    }
  }
  /*====== End of PINT + PDEPI ======*/
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * *
   *          Precipitation related processes          *
   * * * * * * * * * * * * * * * * * * * * * * * * * * */
    
  dynvis   = dynvisc(planet->name,t3);
  N_0S     = MIN(2.e+8,2.0e+6* exp(0.12*(T_triple_pt(is)-t3)));
  lambda_s = pow(M_PI*RHO_SNOW(is)*N_0S/(rho*q_snow),0.25);  
       
  /* PRAUT: Phase change: autoconversion: liquid cloud -> rain */
  if (fcmp(av_liquid,Q_LIQ_0(is)) > 0) {
    praut = MAX(ALPHA_RAUT(is)*(av_liquid-Q_LIQ_0(is)),0.); 
    if (praut < Q_MIN) {
      praut = 0.;
    }
    else {
      av_liquid -= praut*DT;
    }
  }   /*============= End of PRAUT ============*/

  /* PSAUT: Phase change: autoconversion: ice crystals -> snow */
  Q_Icrit = qi_crit[is]*rho_inv; 
  if (fcmp(av_ice,Q_Icrit) > 0) {
    psaut   = MAX((av_ice-Q_Icrit)/DT,0.);   
    if (psaut < Q_MIN) {
      psaut = 0.;
    }
    else {
      av_ice -= psaut*DT;
    }
  }   /*============= End of PSAUT ============*/

  /* PRACW: Accretion of cloud liquid by rain  */
  lambda_r = pow(M_PI*RHO_RAIN(is)*N_0R(is)/(rho*q_rain),0.25);  

  if (q_rain > Q_MIN && av_liquid > Q_MIN) {
    pracw = M_PI*COEFF_XR(is)*q_liquid*E_R(is)*N_0R(is)*pow(P_REF/P3(K,J,I),P_EXP_LIQ(is))*gamma_nr(COEFF_YR(is)+3.)
           /(4.0*pow(lambda_r,COEFF_YR(is)+3.0));     
    pracw = MIN(av_liquid/DT,pracw);
     if (pracw < Q_MIN) {
        pracw = 0.;
    }   /*============= End of PRACW ============*/
    else {
      av_liquid -= pracw*DT;
    }	    
  }
  /* PSACI: Accretion of cloud ice by snow */ 
  if (q_snow > Q_MIN && av_ice > Q_MIN) {
    E_SI  = exp(0.05*(t3-T_triple_pt(is)));         
    psaci = M_PI*COEFF_XS(is)*av_ice*E_SI*N_0S*pow(P_REF/P3(K,J,I),P_EXP_ICE(is))*gamma_nr(COEFF_YS(is)+3.)
            /(4.0*pow(lambda_s,COEFF_YS(is)+3.0));
    psaci = MIN(av_ice/DT,psaci);
    if (psaci < Q_MIN) {
      psaci = 0.;
    }
    else {
      av_ice -= psaci*DT;
    }	    
  }   /*============== End of PSACI ============*/

  /* Evaporation of precipitation takes place if gridbox is subsaturated after cloud evaporation */
  if (q_rain > Q_MIN && warm && subcritical) {
    A_R    = (SQR(Lc(is,VAPOR,LIQUID,t3))/(GAS_R(is)*t3)-1.)
            /(conductivity(planet->name,t3)*t3);
    B_R    = 1./(mass_diffusivity(is,T3(K,J,I),P3(K,J,I))*q_sat*rho);
    tmp    = (f1r(is)/SQR(lambda_r))+f2r(is)*pow(SC(K,J,I),ONE_3)*sqrt(COEFF_XR(is)*rho/dynvis)*pow(P_REF/P3(K,J,I),P_EXP_LIQ(is)/2.)
             *gamma_nr((COEFF_YR(is)+5.)/2.)/pow(lambda_r,(COEFF_YR(is)+5.)/2.);
    if (subsat > Q_MIN) {
      /* 
       * Evaporation of rain
       */ 
      prevap = 2.*M_PI*N_0R(is)*subsat*tmp/(A_R+B_R);
      prevap = MIN(prevap,subsat/DT);
      prevap = MIN(prevap,q_rain/DT);
      if (prevap < Q_MIN) {
        prevap = 0.;
      }
      else {
        subsat -= prevap*DT;  /* reduce subsaturation by prevap */
        prevap  = -prevap;
      }
    }
    else if (supersat > Q_MIN) {
     /* 
      * Depositional growth of rain
      */ 
      prevap = 2.*M_PI*N_0R(is)*supersat*tmp/(A_R+B_R);
      prevap = MIN(prevap,supersat/DT);
      if (prevap < Q_MIN) {
        prevap = 0.;
      }
      else {
        supersat -= prevap*DT;
      }
    }
  }   /*====== End of PREVAP ======*/
  else if (q_snow > Q_MIN && !warm) { 
    A_S    = Ls(is)*rho*(Ls(is)-GAS_R(is)*t3)/(conductivity(planet->name,t3)*GAS_R(is)*SQR(t3));
    B_S    = 1./(mass_diffusivity(is,T3(K,J,I),P3(K,J,I))*q_sat);
    tmp    = (f1s(is)/SQR(lambda_s))+f2s(is)*pow(SC(K,J,I),ONE_3)*sqrt(COEFF_XS(is)*rho/dynvis)*pow(P_REF/P3(K,J,I),P_EXP_ICE(is)/2.)
             *gamma_nr((COEFF_YS(is)+5.)/2.)/pow(lambda_s,(COEFF_YS(is)+5.)/2.);

    if (subsat > Q_MIN) {
      /* 
       * Evaporation of snow
       */ 
      psevap  = 4.0*N_0S*subsat*tmp/(A_S+B_S);
      if (psevap < Q_MIN) {
        psevap = 0.;
      }
      else {
        psevap  = MIN(psevap,subsat/DT);
        psevap  = MIN(psevap,q_snow/DT);  
        subsat -= psevap*DT;
        psevap  = -psevap;
      } 
    }
    else if (supersat > Q_MIN) {
      /* 
       * Depositional growth of snow
       */ 
      psevap = 4.0*N_0S*supersat*tmp/(A_S+B_S);
      psevap = MIN(psevap,supersat/DT);
      if (psevap < Q_MIN) {
        psevap = 0.;
      }
      else {
        supersat -= psevap*DT;
      }
    }
  }  /*====== End of PSEVAP ======*/

  /* * * * * * * * * * * * * * * * * * * * * * * * * * *
   *      End of precipitation related processes       *
   * * * * * * * * * * * * * * * * * * * * * * * * * * */
	   
  Q(is, VAPOR,K,J,I) -= pcond*DT;
  Q(is,LIQUID,K,J,I) += pcond*DT;

  if (subcritical && fabs(pcond) > Q_MIN) {
    HEAT3(K,J,I) -= pcond*Lc(is,VAPOR,LIQUID,t3);
  }
  
  Q(is,VAPOR,K,J,I)  -= (pint+pdepi)*DT;
  Q(is,SOLID,K,J,I)  += (pint+pdepi)*DT; 

  Q(is,LIQUID,K,J,I) -= (praut+pracw)*DT;   /* autoconversion of rain */
  Q(is,  RAIN,K,J,I) += (praut+pracw)*DT;


  Q(is, ICE,K,J,I)   -= (psaut+psaci)*DT;   /* autoconversion of snow */
  Q(is,SNOW,K,J,I)   += (psaut+psaci)*DT;

  Q(is, RAIN,K,J,I)  += prevap*DT;
  Q(is,VAPOR,K,J,I)  -= prevap*DT;

  if (subcritical && fabs(prevap) > Q_MIN) {
    HEAT3(K,J,I) -= prevap*Lc(is,VAPOR,LIQUID,t3);
  }

  Q(is, SNOW,K,J,I)  += psevap*DT;   /* psevap is negative for evaporation */
  Q(is,VAPOR,K,J,I)  -= psevap*DT;
         
  HEAT3(K,J,I) += (pint+pdepi+psevap)*Ls(is);

  return;
}
/*======================= end of finite_rate_processes() =========================*/

/*========================== terminal_velocity() =================================*/

/*
 * Csaba J. Palotai 5/10/04  *A*
 * Returns terminal velocity [m/s] as a positive number.
 *
 * T. Dowling 7/12/24
 * Switched input arguments from (is,ip,kk,J,I) to (is,ip,temperature,pressure,partial_rho),
 * where partial_rho = rhodry*q, for mass-mixing ratio, q. This moves the responsibility
 * for removing any h-weighting on q up to the calling function.
 *
 * Coefficients are set in set_microphysics_params() in epic_microphysics_funcs.c.
 * The influence of parameters such as gravity and dry_rho are incorporated into
 * the coefficients.
 *
 * partial_rho: The mass concentration of the cloud particles [kg/m^3], ie.,
 *                       dry_rho*Q(is,ip). This is typically much less than the particle
 *                       density itself, which is not needed explicitly.
 *
 * NOTE: The return value should be positive downwards.
 *
 * NOTE: In the layer-value case (kk%2 == 0), input partial_rho on the interface
 *       above rather than averaging with the interface below, in upstream fashion.
 *       This is especially important for the leading edge (bottom) of the precipition.
 */

double terminal_velocity(int    is, 
                         int    ip,
                         double temperature,
                         double pressure,
                         double partial_rho)
{
  double
    Re_c,
    w,lambda_r,lambda_s,
    N_0S,N_0I,D_I;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms = 0;
  static char
    dbmsname[] = "terminal_velocity";

  w = 0.;
  switch(ip) {
    case SNOW:
      N_0S     = MIN(2.e+8,2.e+6*exp(0.12*(T_triple_pt(is)-temperature)));
      lambda_s = pow(M_PI*RHO_SNOW(is)*N_0S/partial_rho,0.25);
      w        = COEFF_XS(is)*GAMMA_SNOW(is)*pow(lambda_s,-COEFF_YS(is))
                                            *pow((P_REF/pressure),P_EXP_ICE(is))/6.;
    break;
    case RAIN:
      lambda_r = pow(M_PI*RHO_RAIN(is)*N_0R(is)/partial_rho,0.25);
      w        = COEFF_XR(is)*GAMMA_RAIN(is)*pow(lambda_r,-COEFF_YR(is))
                                            *pow((P_REF/pressure),P_EXP_LIQ(is))/6.;
    break;
    case LIQUID:
      /*
       * NOTE: Need to implement cloud-liquid case.
       *       Currently dropping it into the cloud-ice case.
       */
    case ICE:
      /*  RS 02/13/2019 Cloud ice case. */

      /* 
       * Mean number concentration [1/m^3]
       * COEFF_C and COEFF_D correspond to parameters c and d in Hong et al. (2004),
       * eqns (4) and (5c), for the "single bullet" habit. 
       */
      N_0I = COEFF_C(is)*pow(partial_rho,COEFF_D(is));

      /* 
       * Mean particle diameter [m]
       * COEFF_M and COEFF_N correspond to parameters in Hong et al. (2004),
       * eqn (5b) and M_I = partial_rho/N_OI. 
       */
      D_I = pow(partial_rho/(COEFF_M(is)*N_0I),1./COEFF_N(is));

      /*
       * Terminal velocity [m/s, positive down]
       * COEFF_XI, COEFF_YI, and P_EXP_ICE(is) correspond to parameters x, y, and gamma
       * in Table 1 of Sankar, Klare, and Palotai (2021, Icarus 368, 114589). 
       */
      w = COEFF_XI(is)*pow(D_I,COEFF_YI(is))*pow((P_REF/pressure),P_EXP_ICE(is));
    break;
    case VAPOR:
      w = 0.;
    break;
    default:
      sprintf(Message,"phase %s not yet implemented",var.species[is].phase[ip].info[0].name);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * Ensure that the result is nonnegative.
   */
  return MAX(w,0.);
}

/*======================== end of terminal_velocity() ============================*/

/*========================== moist_convection() ==================================*/

/*
 * Ramanakumar Sankar
 * Relaxed Arakawa Schubert moist convection scheme for EPIC.
 */

#define LVAP(is,T) ((fcmp(T,T_triple_pt(is)) >= 0) ? Lc(is,VAPOR,LIQUID,T) : Ls(is))

#define KBASE(is,j,i)     kbase[is][               i+(j)*Iadim            -Shift2d]
#define QOLD(is,ip,k,j,i) old_spec[is].phase[ip].q[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

/* Maximum iterations for moist convection */
#define MAX_MC_IT grid.max_mc_it

void moist_convection(void)
{
  register int
    K,J,I,
    is,ip,js;
  register double
    dq;
  int
    verbose = FALSE;
  static int 
    *kbase[LAST_SPECIES+1],
    initialized = FALSE;
  static double 
    *a,*b,*c,*r,*u;
  static species_variable
    old_spec[LAST_SPECIES+1];
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int 
    idbms = 0;
  static char
    dbmsname[] = "moist_convection";

  if (grid.moist_convection == OFF) {
    sprintf(Message,"grid.moist_convection is OFF\n");
    epic_error(dbmsname,Message);
  }

  if (grid.coord_type != COORD_ISOBARIC && 
      grid.coord_type != COORD_HYBRID     ) {
    sprintf(Message,"Not implemented for grid.coord_type = %d\n",grid.coord_type);
    epic_error(dbmsname,Message);
  }

  if (!initialized) {
    for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
      /* 
       * Allocate memory.
       */
      if (var.species[is].on) {
        kbase[is] = ivector(0,Nelem2d-1,dbmsname);
        if (grid.coord_type != COORD_ISOBARIC) {
          for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
            old_spec[is].phase[ip].q = dvector(0,Nelem3d-1,dbmsname);
          }
        }
      }
      a = dvector(0,grid.nk-1,dbmsname);
      b = dvector(0,grid.nk-1,dbmsname);
      c = dvector(0,grid.nk-1,dbmsname);
      r = dvector(0,grid.nk-1,dbmsname);
      u = dvector(0,grid.nk-1,dbmsname);
    }

    initialized = TRUE;
  }

  /* 
   * Calculate the cloud base pressure for each species as a function of horizontal position.
   * Also fill in CLOUD_BASE(K,J,I) with the cloud base for each species marked by its index (as a double).
   */
  calculate_base_pressure(kbase);

  /* Exit if only calculating diagnostic values. */
  if (grid.moist_convection == PASSIVE) {
    return;
  }

  /* 
   * Zero global heating array for moist convection (MC).
   */
  memset(var.heat_mc.value,0.,(Nelem3d-1)*sizeof(double));

  /*
   * We treat moist convection for each active species in a time-splitting fashion,
   * where each is done in turn, as if it were the only active species. The other species
   * are updated as if they were passive. 
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      /* 
       * Zero arrays.
       */
      memset(var.species[is].cwf.value,      0.,(Nelem3d-1)*sizeof(double));
      memset(var.species[is].mb_mc.value,    0.,(Nelem3d-1)*sizeof(double));
      memset(var.species[is].dAdt.value,     0.,(Nelem3d-1)*sizeof(double));
      memset(var.species[is].lambda_mc.value,0.,(Nelem3d-1)*sizeof(double));

      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        if (var.species[is].phase[ip].on) {
          if (grid.coord_type != COORD_ISOBARIC) {
            /*
             * Copy previous species value.
             */
            memcpy(old_spec[is].phase[ip].q,var.species[is].phase[ip].q,(Nelem3d-1)*sizeof(double));
          }
          /*
           * Zero tendency.
           */
          memset(var.species[is].phase[ip].dqdt,0.,(Nelem3d-1)*sizeof(double));
        }
      }

      /* 
       * RAS_flux() is the main relaxed Arakawa-Schubert (RAS) driver. It calculates the 
       * entrainment parameter, cloud base mass flux, vapor and potential temperature tendencies.
       */
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          RAS_flux(is,KBASE(is,J,I),J,I,var.species[is].cwf.value);
        }
      }

      /* 
       * The cloud work function (cwf) does not have bc_lateral() applied, 
       * since there are no horizontal derivatives of it.
       */

      /*
       * Update all species, since all are affected when one convects.
       */
      for (js = FIRST_SPECIES; js <= LAST_SPECIES; js++) {
        if (var.species[js].on) {
          for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
            if (var.species[js].phase[ip].on) {
              for (K = KLO; K < KHI; K++) {
                for (J = JLO; J <= JHI; J++) {
                  for (I = ILO; I <= IHI; I++) {
                    Q(js,ip,K,J,I) += DQDT_MC(js,ip,K,J,I)*DT;
                  }
                }
              }
              /* Need to apply bc_lateral() */
              bc_lateral(var.species[js].phase[ip].q,THREEDIM);
              restore_mass(js,ip);
            }
          } /*end of ip loop */
        }
      } /* end of js loop */

      if (grid.coord_type != COORD_ISOBARIC) {
        /* 
         * Update the hybrid density.
         *   h     = hdry + sum(hdry*q)
         *   new h = h    + hdry*sum(dq)
         */
        for (K = KLO; K < KHI; K++) {
          for (J = JLO; J <= JHI; J++) {
            for (I = ILO; I <= IHI; I++) {
              dq = 0.;
              for (js = FIRST_SPECIES; js <= LAST_SPECIES; js++) {
                if (var.species[js].on) {
                  for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
                    if (var.species[js].phase[ip].on) {
                      dq += Q(js,ip,K,J,I)-QOLD(js,ip,K,J,I);
                    }
                  }
                }
              }
              if (fabs(dq) > Q_MIN) {
                H3(K,J,I) += HDRY3(K,J,I)*dq;
              }
            }
          }
        }
        /* Need to call bc_lateral() here. */
        bc_lateral(var.h3.value,THREEDIM);
      }

    }
  } /* end of is loop */

  /* Need to apply bc_lateral() to the moist-convection heating array. */
  bc_lateral(var.heat_mc.value,THREEDIM);


  if (grid.coord_type != COORD_ISOBARIC) {
    /*
     * Solve tridiagonal problem in log(H) to find H on the layers given H3 on the interfaces.
     * Assume the same relations as in set_p2_etc():
     *    H3(K)   = sqrt(H(K)*H(K+1))
     *    H3(KHI) = H(KHI)*H(KHI)/H3(KHI-1)
     */
    for (K = KLO; K < KHI; K++) {
      b[K-KLO] = .5;
      c[K-KLO] = .5;
    }
    K        = KHI;
    b[K-KLO] = 2.;
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        for (K = KLO; K < KHI; K++) {
          r[K-KLO] = log(H3(K,J,I));
        }
        K        = KHI;
        r[K-KLO] = log(H3(K-1,J,I))+log(H3(K,J,I));

        tridiag(grid.nk,a,b,c,r,u,WITH_PIVOTING);

        for (K = KLO; K <= KHI; K++) {
          H(K,J,I) = exp(u[K-KLO]);
        }
      }
    }
    /* No need to call bc_lateral() here. */
  }

  set_p2_etc(DONT_UPDATE_THETA);

  if (verbose) {
    fprintf(stdout,"\n");
  }

  /* RAS has been updated so set the status to FALSE. */
  grid.first_RAS_upd = FALSE;

  return;
}

/*========================= end of moist_convection() ============================*/

/* 
 * NOTE: In many textbooks, mass mixing ratio (dry denominator) and 
 *       specific humidity (moist denominator) are denoted w and q, respectively.
 *       However, in EPIC, Q is mass mixing ratio (dry denominator).  Carrying
 *       mass mixing ratio as the prognostic variable is common in many GCMs 
 *       because the denominator is not changed by precipitation, etc. 
 *       The issue with using W for mass mixing ratio is W also refers to vertical
 *       velocity. Nevertheless, in the following macros, the standard textbook
 *       convention is used rather than the EPIC convention.
 *
 * Macros to convert mass mixing ratio, w = rho_i/rho_dry, to
 * specific humidity, q = rho_i/rho_total, and vice versa. 
 */

#define GET_Q_FROM_W(kk,w) ((w)/(ras[kk].qvj))
#define GET_W_FROM_Q(kk,q) (q*ras[kk].qvj)
#define GET_DQDW(    kk,w) (1./ras[kk].qvj-(w)/SQR(ras[kk].qvj))

#define QTice(t) (MAX(0,MIN(1,(TL-t)/(TL-TLF))))

/*========================= RAS_flux() ===========================================*/

/* 
 * Ramanakumar Sankar
 *
 * Relaxed Arakawa-Schubert (RAS) driver for subgrid-scale moist convection in EPIC.
 * Calculates the entrainment parameter, cloud base mass flux, and the vapor and
 * potential temperature tendencies for a given horizontal position, cloud base,
 * and species.
 * 
 * See documentation of RAS v2 for the algorithm and more details:
 *   https://repository.library.noaa.gov/view/noaa/11400
 */

#undef CWF_IS
#define CWF_IS(k,j,i) cwf_is[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

void RAS_flux(int     is,
              int     kbase_ji,
              int     J, 
              int     I,
              double *cwf_is) 
{
  register int 
    K,kk,k,ktop,jj,
    mc_it,ip,js;
  boolean
    too_hot;
  static RAS_kk_spec
    *ras;
  static RAS_2_spec
    *ras2;
  static RAS_3_spec
    *ras3;
  register double
    Cp,fpara,Lvap,
    t,t1,t2,tmp,
    qvap,qsat,qice,qliq,qt,
    Tv,pdry,psat,
    Pkp,Pk,Pkm,
    lambdai,beta,eps,
    exkk,tem,lc,lp,
    Ft,Gt,Ht,aa,bb,cc,
    zetatop,xitop,hutop,qutop,sutop,ztop,etatop,
    Mb,Ai,Fll,dAdt,zbase,
    dqdt,dwdt,dthetadt,dqsdT,
    Kii,Kii1,
    vTsnow,vTrain;
  double
    rho,
    fgibb,fpe,uoup;
  const double
    alphaf = DT/grid.tau_relax,   /*  relaxation parameter         */
    TL     = -263.16,             /*  RAS v2, for calculating Q(T) */
    TLF    = -30.0;               /*  RAS vs, TL-TF                */
  int
    verbose  = FALSE;
  static int
    initialized = FALSE;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int 
    idbms = 0;
  static char
    dbmsname[] = "RAS_flux";

  if (!initialized) {
    /* 
     * Allocate memory
     */

    /* Array of structures of RAS variables needed on both layer and interface. */
    ras = (RAS_kk_spec *)calloc(2*grid.nk+2,sizeof(RAS_kk_spec));
    if (!ras) {
      sprintf(Message,"Failed to allocate ras structure");
      epic_error(dbmsname,Message);
    }

    /* Likewise, RAS variables needed on layers. */
    ras2 = (RAS_2_spec *)calloc(grid.nk+1,sizeof(RAS_2_spec));
    if (!ras2) {
      sprintf(Message,"Failed to allocate ras2 structure");
      epic_error(dbmsname,Message);
    }

    /* Likewise, RAS variables needed on interfaces. */
    ras3 = (RAS_3_spec *)calloc(grid.nk+1,sizeof(RAS_3_spec));
    if (!ras3) {
      sprintf(Message,"Failed to allocate ras3 structure");
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;
  }

  /* 
   * Don't run moist convection if the base of the cloud is
   * above the tropopause, which we take to be 100 hPa.
   */
  if (PBASE_MC(is,J,I) < 100.e2) {
    return;
  }

  zbase = Z2(kbase_ji,J,I);
  eps   = planet->rgas/GAS_R(is);

  /*
   * NOTE: The current closure scheme does not need to calculate the convective inhibition (CIN), so the K loops
   *       do not go to altitudes below kbase_ji. Should CIN be desired, then the corresponding vertical loops
   *       would go to altitudes below kbase_ji. Going below kbase_ji can reach temperatures that are higher 
   *       than a species' critical temperature, where the enthalpy data tables typically stop (this happens
   *       with ammonia on Jupiter at about 20 bar).
   */

  /* 
   * Calculate the mass factor and virtual temperature contributions from all species.
   */
  for (K = kbase_ji; K > KLO; K--) {
    /* On the interface */
    kk          = 2*K+1;
    ras[kk].nu  = 0.;
    ras[kk].qvj = 1.;
    for (js = FIRST_SPECIES; js <= LAST_SPECIES; js++) {
      if ((var.species[js].on)) {
        for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
          if (var.species[js].phase[ip].on) {
            qvap         = Q(js,ip,K,J,I);
            ras[kk].qvj += qvap;
          }
        }
        ras[kk].nu += (GAS_R(js)/planet->rgas-1.)*(Q(js,VAPOR,K,J,I)+EPSILON)
                                                 /(Q(is,VAPOR,K,J,I)+EPSILON);
      }
    }

    /* On the layer */
    kk          = 2*K;
    ras[kk].nu  = 0.;
    ras[kk].qvj = 1.;
    for (js = FIRST_SPECIES; js <= LAST_SPECIES; js++) {
      if ((var.species[js].on)) {
        for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
          if (var.species[js].phase[ip].on) {
            qvap         = get_var(js,ip,grid.it_h,kk,J,I);
            ras[kk].qvj += qvap;
          }
        }
        ras[kk].nu += (GAS_R(js)/planet->rgas-1.)*(get_var(js,VAPOR,grid.it_h,kk,J,I)+EPSILON)
                                                 /(get_var(is,VAPOR,grid.it_h,kk,J,I)+EPSILON);
      }
    }
  }

  /*
   * Main iteration over cloud types, as distinguished by their detrainment level (cloud top).
   */
  for (mc_it = 0; mc_it < MAX_MC_IT; mc_it++) {
    /* 
     * Pick the cloud top from the iteration number.
     * Skip the first layer above.
     */
    ktop    = kbase_ji-mc_it-2;
    ztop    = Z3(ktop,J,I);
    zetatop = ztop-zbase;
    xitop   = 0.5*zetatop*zetatop;

    if (ktop < KLO+1) {
      continue;
    }
    /* 
     * RAS moist convection has not been implemented in regions serviced by an isentropic coordinate.
     */
    if (ktop < grid.k_sigma) {
      continue;
    }

    /* 
     * Zeroing the ras2[K] structure sets the normalized mass flux to zero at ktop (eta).
     * Zeroing the ras3[K] structure sets the tendencies to zero at kbase_ji (gammas, gammah, gammai, gammal).
     *
     * NOTE: Do not zero the ras[kk] structure. This is not needed and will in fact cause a divide by zero error
     *       in the GET_Q_FROM_W(kk,w) macro.
     */
    memset(ras2,0,(grid.nk+1)*sizeof(RAS_2_spec));
    memset(ras3,0,(grid.nk+1)*sizeof(RAS_3_spec));

    /* 
     * Fill in henv, hsat and satmo, which are static stabilities on both layers and interfaces.
     */
    too_hot = FALSE;
    for (K = kbase_ji; K > KLO; K--) {
      /* * * * * * * * * * * * * *
       * On the lower interface  *
       * * * * * * * * * * * * * */
      kk = 2*K+1;
      t  = T3(K,J,I);

      if (fcmp(t,T_critical_pt(is)) >= 0) {
        /*
         * The RAS scheme is not implemented for T >= Tcrit, so screen for this.
         */
        too_hot = TRUE;

        break;
      }

      fpara = get_var(FPARA_INDEX,NO_PHASE,grid.it_h,kk,J,I);
      Cp    = return_cp(fpara,P3(K,J,I),t);

      /*
       * Calculate satmo = Cp*T+gz => H+PHI
       */
      ras[kk].henv = ras[kk].hsat = ras[kk].satmo = PHI3(K,J,I)
                                                   +return_enthalpy(fpara,P3(K,J,I),t,&fgibb,&fpe,&uoup);

      /* Calculate virtual temperature. */
      qvap = Q(is,VAPOR,K,J,I);
      Tv   = (1.+ras[kk].nu*qvap)*t;

      /*
       * NOTE: Confusingly, EPIC's Q is mixing ratio (dry denominator), which in standard textbooks is denoted w,
       *       such that "GET_Q_FROM_W" is in fact converting mixing ratio into specific humidity (moist denominator).
       */ 
      qvap = GET_Q_FROM_W(kk,Q(is,VAPOR,K,J,I));

      psat = var.species[is].sat_vapor_p(t);
      pdry = PDRY3(K,J,I);
      qsat = eps*psat/pdry;
      qsat = GET_Q_FROM_W(kk,qsat);

      /*
       * Find de/dT, divide by de/dw to get dw/dT, and multiply by dq/dw to get dq/dT.
       * ANT_BL and ANT_BS are the Antoine B coefficients for liquid and solid, respectively, from
       *    log10(esat) = A-B/T .
       */
      dqsdT  = -((fcmp(t,T_triple_pt(is)) >= 0) ? ANT_BL(is) : ANT_BS(is))
               *log(10.)*psat*(1./(t*t));
      dqsdT /= pdry*eps;
      dqsdT *= GET_DQDW(kk,qsat);

      ras[kk].hsat += LVAP(is,t)*qsat;
      ras[kk].henv += LVAP(is,t)*qvap;

      /* 
       * gamma is (Cp/L)*d(qsat)/dT
       * Lt is L-tilde, Eq (28)
       * hst is h**
       */
      ras[kk].gamma  = (Cp/LVAP(is,t))*dqsdT;
      ras[kk].gammaf = ras[kk].gamma/(LVAP(is,t)*(1.+ras[kk].gamma));
      ras[kk].Lt     = Cp*Tv*(1.+ras[kk].gamma)/(1.+ras[kk].gamma*ras[kk].nu*Cp*t/LVAP(is,t));

      /* hst is h** */
      ras[kk].hst = ras[kk].hsat-(ras[kk].nu*ras[kk].Lt)/(1.+ras[kk].nu*qsat)*(qsat-qvap);


      /* * * * * * * * *
       * On the layer  *
       * * * * * * * * */
      kk = 2*K;
      t  = T2(K,J,I);

      if (fcmp(t,T_critical_pt(is)) >= 0) {
        /*
         * The RAS scheme is not implemented for T >= Tcrit, so screen for this.
         */
        too_hot = TRUE;
        break;
      }

      fpara = get_var(FPARA_INDEX,NO_PHASE,grid.it_h,kk,J,I);
      Cp    = return_cp(fpara,P2(K,J,I),t);

      /*
       * Calculate satmo = Cp*T+gz => H+PHI
       */
      ras[kk].henv = ras[kk].hsat = ras[kk].satmo = PHI2(K,J,I)
                                                   +return_enthalpy(fpara,P2(K,J,I),t,&fgibb,&fpe,&uoup);

      /* Calculate virtual temperature. */
      qvap = get_var(is,VAPOR,grid.it_h,kk,J,I);
      Tv   = (1.+ras[kk].nu*qvap)*t;

      /*
       * NOTE: Confusingly, EPIC's Q is mixing ratio (dry denominator), which in standard textbooks is denoted w,
       *       such that "GET_Q_FROM_W" is in fact converting mixing ratio into specific humidity (moist denominator).
       */ 
      qvap = GET_Q_FROM_W(kk,qvap);

      psat = var.species[is].sat_vapor_p(t);
      pdry = onto_kk(P2_INDEX,PDRY3(K-1,J,I),PDRY3(K,J,I),kk);
      qsat = eps*psat/pdry;
      qsat = GET_Q_FROM_W(kk,qsat);

      /*
       * Find de/dT, divide by de/dw to get dw/dT, and multiply by dq/dw to get dq/dT.
       * ANT_BL and ANT_BS are the Antoine B coefficients for liquid and solid, respectively, from
       *    log10(esat) = A-B/T .
       */
      dqsdT  = -((fcmp(t,T_triple_pt(is)) >= 0) ? ANT_BL(is) : ANT_BS(is))
               *log(10.)*psat*(1./(t*t));
      dqsdT /= pdry*eps;
      dqsdT *= GET_DQDW(kk,qsat);

      ras[kk].hsat += LVAP(is,t)*qsat;
      ras[kk].henv += LVAP(is,t)*qvap;

      /* 
       * gamma is (Cp/L)*d(qsat)/dT
       * Lt is L-tilde, Eq (28)
       */
      ras[kk].gamma  = (Cp/LVAP(is,t))*dqsdT;
      ras[kk].gammaf = ras[kk].gamma/(LVAP(is,t)*(1.+ras[kk].gamma));
      ras[kk].Lt     = Cp*Tv*(1.+ras[kk].gamma)/(1.+ras[kk].gamma*ras[kk].nu*Cp*t/LVAP(is,t));

      /* hst is h** */
      ras[kk].hst  = ras[kk].hsat-(ras[kk].nu*ras[kk].Lt)/(1.+ras[kk].nu*qsat)*(qsat-qvap);

      ras2[K].zeta = Z2(K,J,I)-zbase;
      ras2[K].xi   = 0.5*SQR(ras2[K].zeta);
    }

    if (too_hot) {
      /*
       * This mc_it encountered T >= Tcrit, which is not implemented, so move on.
       */
      Mb = 0.;

      continue;
    }

    /* 
     * Entrainment terms are defined on the layer.
     */
    t1 = 0.;
    t2 = 0.;
    for (K = kbase_ji; K > ktop; K--) {
      kk = 2*K;
      t  = T2(K,J,I);

      psat = var.species[is].sat_vapor_p(t);
      pdry = onto_kk(P2_INDEX,PDRY3(K-1,J,I),PDRY3(K,J,I),kk);
      qsat = eps*psat/pdry;
      qsat = GET_Q_FROM_W(kk,qsat);

      /* Eqns. 94a,b,c */
      tmp       = qsat-ras[kk].hsat*ras[kk].gammaf;
      ras2[K].C = tmp+ras[kk].gammaf*ras[2*kbase_ji+1].henv;
      ras2[K].D = ras2[K].zeta*tmp;
      ras2[K].E = ras2[K].xi  *tmp;

      if (K < kbase_ji) {
        t1        += ras[kk+1].henv*(ras2[K].zeta-ras2[K+1].zeta);
        ras2[K].D += ras[kk].gammaf*t1;

        t2        += ras[kk+1].henv*(ras2[K].xi  -ras2[K+1].xi  );
        ras2[K].E += ras[kk].gammaf*t2;
      }
    }
    /* Special case at ktop */
    K  = ktop;
    kk = 2*K+1;
    t  = T3(K,J,I);

    psat = var.species[is].sat_vapor_p(t);
    pdry = PDRY3(K,J,I);
    qsat = eps*psat/pdry;
    qsat = GET_Q_FROM_W(kk,qsat);

    /* Eqns. 94a,b,c */
    tmp       = qsat-ras[kk].hsat*ras[kk].gammaf;
    ras2[K].C = tmp+ras[kk].gammaf*ras[2*kbase_ji+1].henv;
    ras2[K].D = zetatop*tmp;
    ras2[K].E = xitop  *tmp; 

    t1        += ras[kk].henv*(zetatop-ras2[K+1].zeta);
    ras2[K].D += ras[kk].gammaf*t1;

    t2        += ras[kk].henv*(xitop  -ras2[K+1].xi  );
    ras2[K].E += ras[kk].gammaf*t2;

    /*
     * Calculate F, G and H with the cloud values.
     */
    for (K = kbase_ji-1; K > ktop; K--) {
      qt = Q(is,VAPOR,K,J,I)+Q(is,LIQUID,K,J,I)+Q(is,SOLID,K,J,I);
      qt = GET_Q_FROM_W(2*K+1,qt);

      ras2[K].F = ras2[K+1].C-ras2[K].C;
      ras2[K].G = ras2[K+1].D-ras2[K].D+qt*(ras2[K].zeta-ras2[K+1].zeta);
      ras2[K].H = ras2[K+1].E-ras2[K].E+qt*(ras2[K].xi  -ras2[K+1].xi  );
    }
    /* Special case at ktop */
    K  = ktop;
    kk = 2*K+1;

    qt = Q(is,VAPOR,K,J,I)+Q(is,LIQUID,K,J,I)+Q(is,SOLID,K,J,I);
    qt = GET_Q_FROM_W(kk,qt);

    ras2[K].F = ras2[K+1].C-ras2[K].C;
    ras2[K].G = ras2[K+1].D-ras2[K].D+qt*(zetatop-ras2[K+1].zeta);
    ras2[K].H = ras2[K+1].E-ras2[K].E+qt*(xitop  -ras2[K+1].xi  );

    /*
     * Calculate Ftilde, Gtilde and Htilde for getting lambda.
     */
    Ft = 0.;
    Gt = 0.;
    Ht = 0.;
    for (K = kbase_ji-1; K >= ktop; K--) {
      Ft += ras2[K].F;
      Gt += ras2[K].G;
      Ht += ras2[K].H;
    }

    /*
     * Finally, calculate a, b, and c in the quadratic equation for lambda.
     */
    aa = 0.;
    bb = 0.;
    for (K = kbase_ji-1; K > ktop; K--) {
      kk  = 2*K+1;
      aa += (ras2[K].xi  -ras2[K+1].xi  )*ras[kk].henv;
      bb += (ras2[K].zeta-ras2[K+1].zeta)*ras[kk].henv;
    }

    /* Add the cloud contribution at the top. */
    kk   = 2*ktop+1;
    qvap = Q(is,SOLID,ktop,J,I)+Q(is,LIQUID,ktop,J,I);
    qvap = GET_Q_FROM_W(kk,qvap);

    /* Add the constant contributions. */
    aa += (xitop-ras2[ktop+1].xi)*ras[kk].henv
          -xitop*ras[kk].hst
          -ras[kk].Lt*(Ht-xitop*qvap);

    bb += (zetatop-ras2[ktop+1].zeta)*ras[kk].henv
         -zetatop*ras[kk].hst
         -ras[kk].Lt*(Gt-zetatop*qvap);

    cc = ras[2*kbase_ji+1].henv-ras[kk].hst
        -ras[kk].Lt*(Ft-qvap);

    /*
     * Calculate lambda by solving the quadratic equation: 
     *   aa*lam^2 + bb*lam + cc = 0 .
     */
    if (fcmp(aa,0.) == 0) {
      /* If aa = 0., lambdai = -cc/bb */
      lambdai = -cc/bb;
    } 
    else {
      tem = bb*bb-4.*aa*cc;
      if (fcmp(tem,0.) >= 0) {
        tem = sqrt(tem);
        t1  = (-bb+tem)/(2.*aa);
        t2  = (-bb-tem)/(2.*aa);

        /* Limit lambda to at most 100% / km */
        t1 = MIN(t1,1.e-3);
        t2 = MIN(t2,1.e-3);

        lambdai = MAX(t1,t2);
      } 
      else {
        /*
         * No solution -- break out.
         */
        Mb = 0.;

        continue;
      }
    }

    LAMBDA_MC(is,ktop,J,I) = lambdai;

    if (verbose) {
      K = ktop;
      fprintf(stdout,"%d %.3f\n",ktop,lambdai);
    }

    if (fcmp(lambdai,0.) <= 0) {
      Mb = 0.;

      continue;
    }

    /*
     * With lambdai now available, get the updraft properties.
     */
    K  = kbase_ji;
    kk = 2*K+1;
    t  = T3(K,J,I);

    ras2[K].eta = 1.;

    /* Use upwind values. */
    ras3[K].quL = 0.;
    ras3[K].quI = 0.;
    ras3[K].quC = 0.;

    ras3[K].hu = ras[kk].henv;
    ras3[K].su = ras[kk].satmo;

    psat  = var.species[is].sat_vapor_p(t);
    pdry  = PDRY3(K,J,I);
    qsat  = eps*psat/pdry;
    qsat  = GET_Q_FROM_W(kk,qsat);

    for (K = kbase_ji-1; K > ktop; K--) {
      kk = 2*K+1;

      ras2[K].eta = 1.+lambdai*(ras2[K].zeta+lambdai*ras2[K].xi);

      /* 
       * Advect h and s.
       */
      ras3[K].hu = (ras2[K+1].eta*ras3[K+1].hu+ras[kk].henv *(ras2[K].eta-ras2[K+1].eta))/ras2[K].eta;
      ras3[K].su = (ras2[K+1].eta*ras3[K+1].su+ras[kk].satmo*(ras2[K].eta-ras2[K+1].eta))/ras2[K].eta;

      kk = 2*K;
      t  = T2(K,J,I);

      pdry = onto_kk(P2_INDEX,PDRY3(K-1,J,I),PDRY3(K,J,I),kk);
      psat = var.species[is].sat_vapor_p(t);
      qsat = eps*psat/pdry;
      qsat = GET_Q_FROM_W(kk,qsat);

      /*
       * Updraft condensate profile
       */
      ras3[K].quC = (ras2[K+1].eta*ras3[K+1].quC+ras2[K].F+lambdai*(ras2[K].G+lambdai*ras2[K].H))/ras2[K].eta;
      ras3[K].quI = QTice(t)*ras3[K].quC;
      ras3[K].quL = ras3[K].quC-ras3[K].quI;
    }

    K  = ktop;
    kk = 2*K+1;
    t  = T3(K,J,I);

    pdry = PDRY3(K,J,I);
    psat = var.species[is].sat_vapor_p(t);
    qsat = eps*psat/pdry;
    qsat = GET_Q_FROM_W(kk,qsat);

    etatop = 1.+lambdai*(zetatop+lambdai*xitop);

    hutop = (ras2[K+1].eta*ras3[K+1].hu+(etatop-ras2[K+1].eta)*ras[kk].henv )/etatop;
    sutop = (ras2[K+1].eta*ras3[K+1].su+(etatop-ras2[K+1].eta)*ras[kk].satmo)/etatop;

    qutop = qsat+ras[kk].gammaf*(hutop-ras[kk].hsat);

    ras3[K].quC = (Ft+lambdai*(Gt+lambdai*Ht))/etatop;
    ras3[K].quI = QTice(t)*ras3[K].quC;
    ras3[K].quL = ras3[K].quC-ras3[K].quI;

    /*
     * Check top buoyancy condition.
     */
    K  = ktop;
    kk = 2*K+1;

    qliq = GET_Q_FROM_W(kk,Q(is,LIQUID,K,J,I));
    qice = GET_Q_FROM_W(kk,Q(is,SOLID, K,J,I));

    /* 
     * Given lambda, calculate the rates of change of static and moist energies, gammas and gammah.
     */
    for (K = kbase_ji-1; K >= ktop; K--) {
      kk = 2*K+1;

      if (K == ktop) {
        t1 = grid.g[kk][2*J+1]/(P2(K+1,J,I)-P3(K,J,I));
      } 
      else {
        t1 = grid.g[kk][2*J+1]/(P2(K+1,J,I)-P2(K,J,I));
      }

      ras3[K].gammah = t1*(ras2[K  ].eta*(ras[kk-1].henv -ras[kk  ].henv )
                          +ras2[K+1].eta*(ras[kk  ].henv -ras[kk+1].henv ));

      ras3[K].gammas = t1*(ras2[K  ].eta*(ras[kk-1].satmo-ras[kk  ].satmo)
                          +ras2[K+1].eta*(ras[kk  ].satmo-ras[kk+1].satmo));

      qliq = GET_Q_FROM_W(kk,Q(is,LIQUID,K,J,I));
      qice = GET_Q_FROM_W(kk,Q(is,SOLID, K,J,I));

      tem = get_var(is,SOLID,grid.it_h,kk-1,J,I);
      tem = GET_Q_FROM_W(kk-1,tem);
      ras3[K].gammai  = ras2[K  ].eta*(tem-qice);

      tem = get_var(is,SOLID,grid.it_h,kk+1,J,I);
      tem = GET_Q_FROM_W(kk+1,tem);
      ras3[K].gammai -= ras2[K+1].eta*(tem-qice);

      tem = get_var(is,LIQUID,grid.it_h,kk-1,J,I);
      tem = GET_Q_FROM_W(kk-1,tem);
      ras3[K].gammal  = ras2[K  ].eta*(tem-qliq);

      tem = get_var(is,LIQUID,grid.it_h,kk+1,J,I);
      tem = GET_Q_FROM_W(kk+1,tem);
      ras3[K].gammal -= ras2[K+1].eta*(tem-qliq);

      if (K == ktop) {
        /*
         * Add the detraining at the top
         */
        ras3[K].gammah += t1*etatop*(hutop      -ras[kk].henv );
        ras3[K].gammas += t1*etatop*(sutop      -ras[kk].satmo);
        ras3[K].gammai +=    etatop*(ras3[K].quI-qice         );
        ras3[K].gammal +=    etatop*(ras3[K].quL-qliq         );
      }

      ras3[K].gammai *= t1;
      ras3[K].gammal *= t1;
    }

    /*
     * Find the CAPE (convective available potential energy) using the trapezoidal rule.
     */
    Ai = 0.;
    for (K = kbase_ji-1; K >= ktop+1; K--) {
      kk = 2*K+1;

      qice = get_var(is,SOLID, grid.it_h,kk-1,J,I);
      qice = GET_Q_FROM_W(kk-1,qice);
      qliq = get_var(is,LIQUID,grid.it_h,kk-1,J,I);
      qliq = GET_Q_FROM_W(kk-1,qliq);
      t1   = ras3[K  ].hu-ras[kk-1].hst-ras[kk-1].Lt*(ras3[K  ].quC-qice-qliq);

      qice = get_var(is,SOLID, grid.it_h,kk+1,J,I);
      qice = GET_Q_FROM_W(kk+1,qice);
      qliq = get_var(is,LIQUID,grid.it_h,kk+1,J,I);
      qliq = GET_Q_FROM_W(kk+1,qliq);
      t2   = ras3[K+1].hu-ras[kk+1].hst-ras[kk+1].Lt*(ras3[K+1].quC-qice-qliq);

      Ai += 0.5*(PHI2(K,J,I)-PHI2(K+1,J,I))*(t1/ras[kk-1].Lt+t2/ras[kk+1].Lt);
    }
    K    = ktop;
    kk   = 2*K+1;

    qice = get_var(is,SOLID,grid.it_h,kk, J,I);
    qice = GET_Q_FROM_W(kk,qice);
    qliq = get_var(is,LIQUID,grid.it_h,kk, J,I);
    qliq = GET_Q_FROM_W(kk,qliq);
    t1   = hutop-ras[kk].hst-ras[kk].Lt*(ras3[K].quC-qice-qliq);

    qice = get_var(is,SOLID, grid.it_h,kk+1,J,I);
    qice = GET_Q_FROM_W(kk+1,qice);
    qliq = get_var(is,LIQUID,grid.it_h,kk+1,J,I);
    qliq = GET_Q_FROM_W(kk+1,qliq);
    t2   = ras3[K+1].hu-ras[kk+1].hst-ras[kk+1].Lt*(ras3[K+1].quC-qice-qliq);

    Ai += 0.5*(PHI3(K,J,I)-PHI2(K+1,J,I))*(t1/ras[kk].Lt+t2/ras[kk+1].Lt);

    /* 
     * Ignore all non-positive CAPE.
     */
    if (fcmp(Ai,0.) <= 0) {
      Mb = 0.;

      continue;
    }

    /* 
     * Find (1./Mb) d(CAPE)/dt.
     */
    Fll = 0.;
    for (K = kbase_ji-1; K >= ktop+1; K--) {
      kk = 2*K;

      Fll += 0.5*(PHI3(K-1,J,I)-PHI3(K,J,I))
                  *(ras2[K-1].eta*((1.+ras[kk-1].gamma)*ras3[K-1].gammas/ras[kk-1].Lt+ras3[K-1].gammai+ras3[K-1].gammal)
                   +ras2[K  ].eta*((1.+ras[kk+1].gamma)*ras3[K  ].gammas/ras[kk+1].Lt+ras3[K  ].gammai+ras3[K  ].gammal));
    }

    dAdt = (Ai-CWF_IS(ktop,J,I))/DT;
    Mb   = -dAdt/Fll;

    /*
     * Set the cloud work function to the new value.
     */
    CWF_IS(ktop,J,I) = Ai;

    /*
     * Exit if this is the first time this function is called.
     */
    if (grid.first_RAS_upd) {
      Mb = 0.;

      continue;
    }

    /*
     * Ignore non-increasing CAPE.
     */
    if (fcmp(dAdt,0.) <= 0) {
      Mb = 0.;

      continue;
    }

    /*
     * Ignore non-positive mass flux.
     */
    if (fcmp(Mb,0.) <= 0) {
      Mb = 0.;

      continue;
    }

    if (verbose) {
      K = ktop;
      fprintf(stdout,"%d Ai: %.3e kernel: %.3e dA/dt: %.3e Mb: %.3e \n",
                      mc_it,Ai,Fll,dAdt,Mb);
    }

    /*
     * Calculate rhodry(kbase).
     */
    rho = PDRY3(kbase_ji,J,I)/(planet->rgas*T3(kbase_ji,J,I));

    /*
     * Limit the mass flux to at most density times 100 m/s large scale motion.
     */
    Mb = MIN(Mb,rho*100.);

    DADT(is, ktop,J,I) = dAdt;
    MB_MC(is,ktop,J,I) = Mb;

    /*
     * Relax the scheme.
     */
    Mb *= alphaf;

    /*
     * Update the specific humidities, but not the thermal profile.
     */
    for (K = kbase_ji; K >= ktop; K--) {
      kk = 2*K+1;
      t  = T3(K,J,I);

      /* Calculate the moist convection heating rate, Cp*(dT/dt). */
      HEAT_MC(K,J,I) += Mb*ras3[K].gammas;

      /*
       * Calculate dqdt from moist convection for the different phases,
       * for diagnostic purposes. Convert into mass-mixing-ratio rates.
       *
       * Calculate the rate of conversion for vapor.
       */
      dqdt  = (ras3[K].gammah-ras3[K].gammas)*Mb/LVAP(is,t);
      dqdt /= GET_DQDW(kk,Q(is,VAPOR,K,J,I));

      DQDT_MC(is,VAPOR,K,J,I) += dqdt;

      /*
       * Calculate the rate of conversion for ice.
       */
      dqdt = ras3[K].gammai*Mb/GET_DQDW(kk,Q(is,SOLID,K,J,I));

      DQDT_MC(is,SOLID,K,J,I) += dqdt;

      /*
       * Calculate the rate of conversion for liquid.
       */
      dqdt = ras3[K].gammal*Mb/GET_DQDW(kk,Q(is,LIQUID,K,J,I));

      DQDT_MC(is,LIQUID,K,J,I) += dqdt;

      if (fabs(DT*HEAT_MC(K,J,I)/EXNER3(K,J,I)) > 10.) {
        sprintf(Message,
                "HEAT: %.3e > 10 [K]; ptop: %.3f Mb: %.3e dAdt: %.3e Ki: %.3e Ai: "
                "%.3e\n",
                DT*HEAT_MC(K,J,I)/EXNER3(K,J,I),P3(ktop,J,I),Mb,dAdt,Fll,Ai);
        epic_warning(dbmsname,Message);
      }
    }

    /* 
     * Advect the vapor from other species from this upwelling. This is done
     * similarly to how h is advected, by reusing the algorithms for quC for 
     * the updraft q profile and gammai for the updraft tendency.
     */

    ip = VAPOR;
    for (js = FIRST_SPECIES; js <= LAST_SPECIES; js++) {
      if ((js != is) && var.species[js].on && var.species[js].phase[ip].on) {
        /* Zero ras3. */
        memset(ras3,0,(grid.nk+1)*sizeof(RAS_3_spec));

        /*
         * First, build the convective profile.
         */
        K  = kbase_ji;
        kk = 2*K+1;

        ras3[K].quC = GET_Q_FROM_W(kk,Q(js,ip,K,J,I));

        for (K = kbase_ji-1; K > ktop; K--) {
          kk = 2*K+1;

          qvap = Q(js,ip,K,J,I);
          qvap = GET_Q_FROM_W(kk,qvap);

          /* Updraft condensate profile. */
          ras3[K].quC = (ras2[K+1].eta*ras3[K+1].quC+qvap*(ras2[K].eta-ras2[K+1].eta))/ras2[K].eta;
        }
        K  = ktop;
        kk = 2*K+1;

        ras3[K].quC = (ras2[K+1].eta*ras3[K+1].quC+GET_Q_FROM_W(kk,Q(js,ip,K,J,I))*(etatop-ras2[K+1].eta))/etatop;

        for (K = kbase_ji - 1; K >= ktop; K--) {
          kk = 2*K+1;

          if (K == ktop) {
            t1 = grid.g[kk][2*J+1]/(P2(K+1,J,I)-P3(K,J,I));
          } 
          else {
            t1 = grid.g[kk][2*J+1]/(P2(K+1,J,I)-P2(K,J,I));
          }

          qvap = GET_Q_FROM_W(kk,Q(js,ip,K,J,I));

          tem = get_var(js,ip,grid.it_h,kk-1,J,I);
          tem = GET_Q_FROM_W(kk-1,tem);
          ras3[K].gammai  = ras2[K  ].eta*(tem-qvap);

          tem = get_var(js,ip,grid.it_h,kk+1,J,I);
          tem = GET_Q_FROM_W(kk+1,tem);
          ras3[K].gammai -= ras2[K+1].eta*(tem-qvap);

          if (K == ktop) {
            ras3[K].gammai += etatop*(ras3[K].quC-GET_Q_FROM_W(kk,Q(js,ip,K,J,I)));
          }
          ras3[K].gammai *= t1;
        }

        for (K = kbase_ji; K >= ktop; K--) {
          kk = 2*K+1;

          dqdt = ras3[K].gammai*Mb/GET_DQDW(kk,Q(js,ip,K,J,I));

          DQDT_MC(js,ip,K,J,I) += dqdt;
        }
      }
    }
  } /* end of mc_it loop */

  return;
}

/*========================== end of RAS_flux() ===================================*/

/*========================= calculate_base_pressure() ============================*/

/*
 * Finds the index and cloud base pressure for all species.
 * This is defined as the layer where the relative humidity is greater than 90%.
 */

void calculate_base_pressure(int **kbase) 
{
  register int
    K,J,I,kk,is;
  register double
    pp,psat,rh;

  /*
   * Zero CLOUD_BASE(K,J,I).
   */
  memset(var.cloud_base.value,0,Nelem3d*sizeof(double));

  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if ((var.species[is].on)) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          /*
           * The base of the cloud is determined as
           * the first interface where RH > 90%.
           *
           * For the case of no clouds in the column, initialize
           * kbase to the second interface below the top.
           */
          KBASE(   is,J,I) = KLO+1;
          PBASE_MC(is,J,I) = P3(KLO+1,J,I);

          for (K = KHI; K > KLO; K--) {
            kk = 2*K+1;

            /*
             * Start with the "old school" definition of RH (relative humidity)
             * as the ratio of vapor partial pressure to saturation pressure.
             */
            pp   = get_p(is,kk,J,I);
            psat = var.species[is].sat_vapor_p(T3(K,J,I));
            rh   = pp/psat;

            /*
             * Factor in the term that converts this to the ratio of mixing
             * ratio to saturation mixing ratio.
             */

            rh *= 1.-(psat-pp)/PDRY3(K,J,I);

            if (rh > 0.9) {
              KBASE(   is, J,I) = K;
              PBASE_MC(is, J,I) = P3(K,J,I);
              CLOUD_BASE(K,J,I) = (double)is;

              break;
            }
          }
        }
      }
      /* No need to apply bc_lateral() */
    }
  }

  return;
}

/*========================= end of calculate_base_pressure() =====================*/

/* * * * * * * * * * * * * end of epic_microphysics.c * * * * * * * * * * * * * * */
