library(pomp)

source("code/common/Csnippet_service.R")



covid_constants = Csnippet(paste0(
  mat2cpp( phi, "phi"),
  mat2cpp( YHR, "YHR"),
#//  vec2cpp( totalpop, "N"),
  vec2cpp( omega_P_overall, "omega_P"),
  vec2cpp( HFR, "HFR")
))

# -------------------

covid_indexed_states = c("S", "E", "PY", "PA", "IY", "IA",
                         "H", "R",
                         "new_H", "leaving_H",
                         "new_E", "leaving_E")


covid_statenames = c(
  expand_state( covid_indexed_states, c(2,5) ),
  "NI",
  "NH",
  "LH",
  "NIY",
  "Beta",
  "Z",
  "leave_H_rate",
    "total_new_H",
    "total_leave_H",
  "total_H",
  "Z_mu",
  "Z_R",
  "b1",  ## DLM
  "b2",  ## DLM
  "b3",  ## DLM
  "bH"  ## DLM
)

covid_accum_vars = c(
   "total_new_H",
    "total_leave_H"
)

covid_init_pars =
  c(
    rNH = 15.19,
    rH  = 15,
    psi_Z = 0.97,    #// AR(1) coefficient 
    psi_bH = 0.97,   #// AR(1) for bH, chance to leave hospital
    sig_Z = 0.05,   #// for the transmission prob
#//    log_mu_0 = log(1/12),
#//    log_gamma_H_0 = log(1/8.7),
    log_beta_0 = log(0.06),
    Z_0=0.0,
 #//   b1_0 = 0.25,         
    sig_delta_b1 = 0.02, 
    sig_delta_b2 = 0.02, 
    sig_delta_b3 = 0.02, 
    bH_0 = logit(0.1432) , # Set from Austin daily release rate. logit(0.1432)         
    sig_delta_bH = 0.03, ### <- change for not admits was 0.5 ???
    sigma = sigma, #// Incubation period (no transmission)
    tau= tau, #// symp_prop
    eta = eta,  #// t_symp_onset_rate
    omega_Y = omega_Y,
    omega_A = omega_A,
    rho_Y = rho_Y,  #// pre-symptomatic rate
    rho_A = rho_A,
    gamma_Y = gamma_Y,  #// recovery rate
    gamma_A = gamma_Y,
    init_scale_0 = 10.0 ,
    hospitalized_0 = 1.0 
  )

covid_paramnames = names(covid_init_pars)

covid_partrans = parameter_trans(log=c("sig_Z", "rNH","rH",  "sig_delta_b1","sig_delta_b2","sig_delta_b3","sig_delta_bH"),
  logit=c("init_scale_0","hospitalized_0")
                                       
                                 )
#// logit=c( "psi_Z")
#//rw.sd=rw.sd(a = 0.03, b1=0.03, b12=0.02, bH=0.03, b3=0.03, b4=0.02, b5=0.02,
covid_rw.sd =rw.sd( 
                 log_beta_0 = 0.03, 
                 rNH= 0.02,
                 rH = 0.02, #// if optimize r, uncomment the prior in dmeas
                 sig_Z=0.01, 
                 init_scale_0 = 0.05, 
                 hospitalized_0 = 0.005 ,
                 sig_delta_b1 = 0.01, 
                 sig_delta_b2 = 0.01, 
                 sig_delta_b3 = 0.01, 
                 sig_delta_bH = 0.00,  ### <- change for not admits was 0.01
                 bH_0 = 0.00 # DLM
                 #// b1_0 = 0.00, 
                 #// psi_Z=0.0, 
                 #// log_mu_0=0.0,
                 #// log_gamma_H_0=0.0
                ) 
    



# c("<5", "5-17", "18-49", "50-64", "65+"
init_infected = matrix( dimnames = list( risk=1:2,
  age=c(
    "<5", "5-17", "18-49", "50-64", "65+")),
  c(   0,      0,       0,       0,     0,
       0,      0,       1,       0,     0), 
       2, 5, byrow = TRUE)
zero_mat = matrix(0, 2, 5)



rinit = Csnippet(paste0(
  init_mat_var( pop     ,      name="S"),
  #// The following initialize the infected pop, and have to be converted to actual
  #// individuals. They are initialized to a value that sums to one altogether.
  init_mat_var( init_E_ratio,  name= "E"),
  init_mat_var( init_PA_ratio, name= "PA"),
  init_mat_var( init_PY_ratio, name= "PY"),
  init_mat_var( init_IA_ratio, name= "IA"),
  init_mat_var( init_IY_ratio, name= "IY"),
  init_mat_var( init_H_ratio,  name= "H"),
  init_mat_var( zero_mat,      name= "new_H"),
  init_mat_var( zero_mat,      name= "leaving_H"),
  init_mat_var( zero_mat,      name= "R"),
  enable_indexes(covid_indexed_states, dims=c(2, 5)),
  #'
  '
  NI=0;
  NIY=0;
  NH=0;
  LH=0;
  Z=Z_0;
  Z_mu=0;
  Z_R=0;
  total_new_H=0.0;
  total_leave_H=0.0;
  total_H=0.0;
  b1 = 0.0; // no need for b1_0, because this is added to log_beta_0. 
  b2 = 0.0; // no need for b2_0, because this is added to log_beta_0. 
  b3 = 0.0; // no need for b3_0, because this is added to log_beta_0. 
  bH = bH_0; 
//  Rprintf("%g, \\n",init_scale_0) ;

  // Calculate sum of H[g][i] to normalize
  double sH = 0.0 ;
  for(int g=0; g<2; g++) {
    for(int i=0; i<5; i++) {
      sH += H[g][i] ;
    }
  }
  
  // hospitalized_0 is total number of in H.
  // normalize H[g][i] and then multiply by hospitalized_0
  // Actually, we do not need to normalize... fitting hospitalized_0 should be able
  // to handle this.
  for(int g=0; g<2; g++) {
    for(int i=0; i<5; i++) {
      H[g][i] = rpois( H[g][i] / sH * hospitalized_0 * S[g][i]) ;
    }
  }
//  Rprintf("%g\\n",init_scale_0);
     
  for(int g=0; g<2; g++) {
    for(int i=0; i<5; i++) {
      E[g][i]  = rpois( E[g][i]  * init_scale_0* S[g][i]) ;
      PY[g][i] = rpois( PY[g][i] * init_scale_0* S[g][i]) ;
      PA[g][i] = rpois( PA[g][i] * init_scale_0* S[g][i]) ;
      IY[g][i] = rpois( IY[g][i] * init_scale_0* S[g][i]) ;
      IA[g][i] = rpois( IA[g][i] * init_scale_0* S[g][i]) ;
  // Rprintf("E[%d,%d]=%g\\n",g,i,E[g][i]);
  //    H[g][i] = 0 ;
  //    H[g][i] = rpois( H[g][i] * init_scale_0) ;
      total_H += H[g][i];
      NI      += E[g][i]+PY[g][i]+PA[g][i]+IY[g][i]+IA[g][i] ;
      NH      += H[g][i] ;
      NIY     += IY[g][i];
      S[g][i] -= E[g][i]+PY[g][i]+PA[g][i]+IY[g][i]+IA[g][i] + H[g][i] ;
      S[g][i] = round( S[g][i]) ;
 // Rprintf("S[%d,%d]=%g\\n",g,i,S[g][i]);
    }
  }
  ', #'
    register_indexed_changes(covid_indexed_states, dims=c(2, 5))
))









# time t=1,..,T
# data Y_t   
# latent Z_t
# reconstruct Z_t from the likelihood



# Z_t = (S_t, E_t, I_t, ...)
# filtered states:  p(Z_T | Y_{1:T}) latent states today given the history
# smoothed states:  p(Z_1:T | Y_{1:T}) latent state at any point t given all history
# forecasted states:

# ---------------------

covid_rprocess = Csnippet(paste0(
  enable_indexes(covid_indexed_states, dims=c(2, 5)),
  covid_constants, #'
  '
  double new_PA;
  double new_PY;
  double new_IA;
  double new_IY;
  double recovering_IA;
  double recovering_IY;
  double leaving_IY ;


//#define BIN_dt( N, rate) ( ( (N) > 0.5) ? rbinom( (N) , 1.0 - exp(         -(rate)  * dt) ) : 0.0 )
// The code below gives 1-(1-rate)^dt instead of 1-exp(-rate*dt) above.
  #define BIN_dt( N, rate) ( ( (N) > 0.5) ? rbinom( (N) , 1.0 - exp(  log(1.0-(rate)) * dt) ) : 0.0 )
  #define BIN(    N, rate) ( ( (N) > 0.5) ? rbinom( (N) , (rate)                   ) : 0.0 )

  b1       = b1    +  rnorm( 0.0, sig_delta_b1 * sqrt(dt)) *(1-future) ; //DLM
  b2       = b2    +  rnorm( 0.0, sig_delta_b2 * sqrt(dt)) *(1-future) ; //DLM
  b3       = b3    +  rnorm( 0.0, sig_delta_b3 * sqrt(dt)) *(1-future) ; //DLM
  bH       = bH    + (rnorm( 0.0, sig_delta_bH * sqrt(dt)) + bH  * (psi_bH-1.0)  ) *(1-future) ; // leaving
  
  Z = future==1? Z: psi_Z * Z + rnorm( 0.0, sig_Z        * sqrt(dt)) ;

  Beta = exp(   log_beta_0  
              + b1 * PC1 + b2 * PC2 + b3 * PC3
              + Z  
              );
// leave_H_rate = expit( bH + bH_0 ) ;
 leave_H_rate = expit(  bH_0 ) ;
// Calc total pop. Doesnt really need to happen every step, but better than relying on a single calc from pop.
  double N[5] ;

  for( int j=0; j<5; j++) {
    N[j] = 0 ;
    for( int g=0; g<2; g++)
      N[j] += S[g][j]  + E[g][j] +
              PA[g][j] + PY[g][j] + 
              IA[g][j] + IY[g][j] + 
              //D[g][j]  +
              R[g][j] ;
  }

  double inf[5] ;
  double leaving_S_rate[5] ;

  for (int j=0; j<5; j++) {
    inf[j] = 
      Beta *              // j is infecting i
            (   omega_A * omega_P[j] * (PA[0][j] + PA[1][j]) +
                omega_Y * omega_P[j] * (PY[0][j] + PY[1][j]) +
                omega_A              * (IA[0][j] + IA[1][j]) +
                omega_Y              * (IY[0][j] + IY[1][j])
            ) / N[j]; // should subtract dead?
}


  // rate for being infected
  // calculate leaving_S_rate[g][i]
  for (int i=0; i<5; i++) {
    // rate of infectiousness
    leaving_S_rate[i] = 0.0;
    for (int j=0; j<5; j++) {
      leaving_S_rate[i] +=  phi[i][j] *  inf[j] ;
    }
  }

//  double mu      = exp( log_mu_0      + Z_mu);
//  double gamma_H = exp( log_gamma_H_0 + Z_R );

  total_H = 0.0;


  for (int g=0; g<2; g++) {
    for (int i=0; i < 5; i++) {
      // Rprintf("in: S[%d,%d]=%g\\n",g,i,S[g][i]);
      double pi =  gamma_Y * YHR[g][i] / (eta + (gamma_Y - eta) * YHR[g][i]);
      double leaving_IY_rate = ((1.0 - pi) * gamma_Y + pi * eta);

      // newly infected
      new_E[g][i]     = BIN_dt( S[g][i], leaving_S_rate[i] );    
      leaving_E[g][i] = BIN_dt( E[g][i], sigma ) ; // end of incubation period

      new_PA    = BIN   ( leaving_E[g][i],   1.0 - tau ) ;
      new_PY    = leaving_E[g][i] - new_PA;
      
      new_IA    = BIN_dt( PA[g][i], rho_A ) ;
      new_IY    = BIN_dt( PY[g][i], rho_Y );

      recovering_IA = BIN_dt(IA[g][i], gamma_A );  // Recovering asymptomatics 
      // Hospitalization and recovering symptomatics
      leaving_IY = BIN_dt( IY[g][i], leaving_IY_rate );
      new_H[g][i]      = BIN( leaving_IY, pi * eta / leaving_IY_rate );
      recovering_IY = leaving_IY - new_H[g][i];


      leaving_H[g][i]  = BIN_dt( H[g][i]+new_H[g][i], leave_H_rate );

      //////////////////////////////////////////////////////////////////////////////////////        
      S[g][i]  +=              - new_E[g][i];
      E[g][i]  += new_E[g][i]  - leaving_E[g][i];
      PA[g][i] += new_PA - new_IA;
      PY[g][i] += new_PY - new_IY;
      IA[g][i] += new_IA - recovering_IA;
      IY[g][i] += new_IY - leaving_IY;
      H[g][i]  += new_H[g][i]  - leaving_H[g][i];
      R[g][i] += recovering_IA + recovering_IY + leaving_H[g][i] ;
//      D[g][i] += dying_H[g][i];
      NI  += leaving_E[g][i];
      NH  += new_H[g][i];
      LH  += leaving_H[g][i];
      NIY += new_IY ;
      
      // Totals for the likelihood
      total_H            += H[g][i];
      total_new_H        += new_H[g][i];
 //     Rprintf("%g,%g:%g ",new_H[g][i],leaving_IY[g][i],total_new_H) ;
      total_leave_H      += leaving_H[g][i];

    }
  }
  
  // Update states
  
  
//  Rprintf("Here %g") ;
 
 // Rprintf("\\n");
  ', #'
  register_indexed_changes(covid_indexed_states, dims=c(2, 5))
))


#// dmeas is f(y_n|x_n;par)
covid_dmeas_admitdischarge_nb = Csnippet(paste0( #'
  '
  double loglike = 0.0;  
  // clamp r for stability
  double like=0 ;

// EQ:  hosp_diff = new_H - leaving_H
//  new_H = leaving_H + hosp_diff ;
// We have leaving_H and new_H in the model.

// We can do two things: 
// 1. What is chance to observe diff?
// leaving_H is binomial on H with rate Z_leave
// new_H is a binomial on IY, but we do it as a negative binomial
// So, we could take all possible leaving_H, for each calculate new_H, and then do 
// dnbinom(new_H, r, IY*) * binom( leaving_H, H, Z_leave)
// 
// We can also draw leave_H from binomial, then calculate new_H, and then do dnbinom. That seems easiest.
  double eps = 1e-6 ;
  /// ********* disabling use of admits
//**  if( !R_IsNA(adm_total) )
//**    loglike += dnbinom_mu( adm_total, rNH +eps,  total_new_H, 1 ) ;
//  if( !R_IsNA(leave_total) )
//    loglike += dnbinom_mu( leave_total, r_,  total_leave_H, 1 ) ;
  // Since hosp_dif = total_new_H - total_leave_H
  //  Rprintf(\"dnbinom_mu %g %g %g = %g\\n\", 
   //           total_leave_H, hosp_dif, total_new_H, 
    //          dnbinom_mu( total_leave_H + hosp_dif, r_,  total_new_H, 1 )); //"



 
  if( !R_IsNA( hospitalized) )
      loglike +=   dnbinom_mu( hospitalized, rH +eps, total_H,            1); // annealed cause not necessary // Michael: stronger anneal, was 0.5*
  // prior
  loglike += dexp(rH , 1.0, 1) / maxT; 
  loglike += dexp(rNH, 1.0, 1) / maxT; 
  loglike += dgamma(sig_Z, 1.1, 1.1, 1) / maxT; //
  // loglike += dgamma(sig_mu, 1.1, 1.1, 1) / maxT; //
 //  loglike += dgamma(sig_R, 1.1, 1.1, 1) / maxT; //
  loglike += dgamma(sig_delta_b1, 1.1, 1/1.1, 1) / maxT; //## DLM
  loglike += dgamma(sig_delta_b2, 1.1, 1/1.1, 1) / maxT; //## DLM
  loglike += dgamma(sig_delta_b3, 1.1, 1/1.1, 1) / maxT; //## DLM
  loglike += dgamma(sig_delta_bH, 1.1, 1/1.1, 1) / maxT; //## DLM
  
  //loglike += dnorm(b1_0, 0.0, 1.0, 1) / maxT; //
  loglike += dnorm(bH_0, 0.0, 1.0, 1) / maxT; //
//  loglike += dnorm(b2, 0.0, 1.0, 1) / maxT; //
//  loglike += dnorm(b3, 0.0, 1.0, 1) / maxT; //
//  loglike += dnorm(b5, 0.0, 1.0, 1) / maxT; //

 // Rprintf("%g %g %g %g \\n", loglike,hospitalized, total_H, total_new_H) ;
  

  lik = (give_log) ? loglike : exp(loglike);
  ' #'
))

covid_rmeas = Csnippet( #'
  '
      adm_total = rnbinom_mu(total_new_H, rNH + 1e-6);       
  '
  #'
)