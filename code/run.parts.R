

#// Quite ugly, but currently city read after model, because some parameters need both.

#source("code/models/austin_data.R")
#source("code/models/agegroups+mu+pres+risk+vacc.R.cpp")
source("code/models/agegroups+mu+pres+risk.R")


# This function fits a mif2 object to the data
# data contains the hospitalization data, which needs a column called 
# hospitalized 
    covid_init_pars["init_scale_0"] = (hospitalization_data$hospitalized[1] * 25 +25)/sum(pop)
    covid_init_pars["hospitalized_0"] = (hospitalization_data$hospitalized[1] +1)/sum(pop)
    


find.mif = function( data, covars, N.steps=200, N.particles=2000, clust.n=14, dt = 1) {
  my.cl = makeCluster( clust.n, outfile="" )
  registerDoParallel( my.cl )
  
  covid_mod =
        pomp(data=data,
            times="day", 
            t0 = 0,
            rinit=rinit,
            rmeasure=covid_rmeas,
            dmeasure=covid_dmeas_admitdischarge_nb, 
            covar=covariate_table(select(covars, -date), times="day"),
            rprocess = pomp::euler(covid_rprocess,dt),
            statenames= covid_statenames,
            paramnames= covid_paramnames,
            params = covid_init_pars,
            accumvars = covid_accum_vars
        )


    # Simulate from the prior pomp model and plot ------------------------------------
  covid_mod %>%
       pomp::simulate(nsim = 1, format="data.frame") ->a
    #   View

  
    #save.session(file=paste0("before.mf2.",run.name,".Rda"))

  mif2.N.steps= N.steps
  mif2.iter = min( 10, N.steps)
  # browser() ;

  mf2s = foreach( cpu=1:clust.n, 
                        .packages = c("magrittr", "pomp"), 
                        .export=c("covid_partrans","covid_init_pars","covid_paramnames","covid_rw.sd"),
                        .verbose = T,
                        .errorhandling = "remove"
  ) %dopar% {
    mf2 <- covid_mod %>%
    mif2(
        Np= N.particles,
        Nmif= mif2.iter,
        paramnames= covid_paramnames,
        params = covid_init_pars,
        partrans=covid_partrans,
        rw.sd=covid_rw.sd,
        cooling.fraction.50=0.5
    )


    for( i in (attributes(mf2)$Nmif/mif2.iter):(mif2.N.steps/mif2.iter) ) {
        print( c(cpu=cpu, mif2.iter=i, time=system.time( { mf2 %>% continue(Nmif= mif2.iter) ->mf2 } )[1] ) )
    }
    mf2
  }
  ll = foreach( i=1:length(mf2s) , 
                .packages = c("magrittr", "pomp")
                ) %dopar% { 
                  replicate( 10, mf2s[[i]] %>% pfilter() %>% logLik ) %>% mean
                }
  stopCluster(my.cl)
  mf2s[[ which.max(ll)]]
}


make_smoothed_states = function( mf2, covars,  N.draws=500, N.particles=2500, clust.n = 1, dt = 1) {
    if( clust.n > 1) {
      my.cl = makeCluster(clust.n )
      registerDoParallel( my.cl )
    }
    ## Estimate a smoothed distribution over beta --------------------------------
    N_smooth_draws = N.draws * 3
    N_particles_per_draw = N.particles

    # Filtered states  p(Z_t | Y_1:t)
    # Smoothed states  p(Z_t | Y_1:T)  t <= T
    # Forecasted states p(Z_t | Y_1:T) = int p(Z_t | Z_T)p(Z_T | Y_{1:T}) dZ_T t > T

    # construct the raw filtered states for the smoothed distribution
    smooth_dist = foreach(i=1:N_smooth_draws, 
                        .packages = c("magrittr", "pomp")
    ) %dopar% {
        print( c( smooth.iter=as.character(i), time=as.character( Sys.time())  ) )
        mf_filter = mf2 %>%
            pfilter(Np = N_particles_per_draw, save.states=TRUE, filter.mean=TRUE, filter.traj=TRUE)
        st=drop(filter.traj(mf_filter))
        cur.i = length(covars$date)
        st = st[,1:cur.i]  # This is just because sometimes st has one more column than needed ?
        
        list(ll = logLik(mf_filter), states=st)
    }

   if( clust.n > 1) {
      stopCluster( my.cl )
    }

    n_states = nrow(smooth_dist[[1]]$states)

    # extract the weights for each filtered trajectory
    smoothing_weights = smooth_dist %>%
    lapply(function(x) x$ll) %>%
    unlist %>%
    (function(x) {exp(x - max(x))})

    # bind up the filtered trajectories for each state in an array
    filtered_states = smooth_dist %>%
      lapply(function(x) x$states) %>%
      unlist %>%
    #  array(dim=c(n_states, maxT, length(smooth_dist)))
      array(dim=c(n_states, dim(smooth_dist[[1]]$states)[2], length(smooth_dist)))
    dimnames(filtered_states)[[1]] = dimnames( smooth_dist[[1]]$states )[[1]]
    i = sample( pomp::systematic_resample(smoothing_weights), N.draws, rep=F )
    filtered_states[,,i]
}

add_future_states = function( coefs, data, covars, filtered_states, dt = 1 ) {
  future.pomp = pomp(
    data=data,
    times= "day", 
    t0 = covars$day[1],
    rinit=rinit,
    rmeasure=covid_rmeas,
    dmeasure=covid_dmeas_admitdischarge_nb, 
    covar=covariate_table(select(covars, -date), times="day"),
    rprocess = pomp::euler(covid_rprocess,dt),
    statenames= covid_statenames,
    paramnames= covid_paramnames,
    params = coefs,
    accumvars = covid_accum_vars
  ) 

  sim_future_states = sapply( 1:dim(filtered_states)[3], simplify="array",function(i) {
    future.pomp %>% 
    pomp:::simulate(nsim=1, 
                    format="arrays",
                    rinit=function(...) { setNames( filtered_states[ , dim(filtered_states)[2], i],
                                                dimnames(filtered_states)[[1]] ) },
                    statenames=covid_statenames,
                    times=future_covars$day,
                    t0 = future_covars$day[1] ) %>% 
    .[["states"]] %>%
    .[,1,]
  })

  filtered_states_new = abind::abind(filtered_states, sim_future_states[,-(1:2),], along=2)


}


make_future_states = function() {
    # ------------------------------------------------------------------------
# Augment with forecasts counterfactual using sim
max_dat_date = max(hospitalization_data$date)
min_dat_date = min(hospitalization_data$date)

# Modify for forecasting
new_dates = seq(max_dat_date + 1, max_dat_date + num_forecast_days, by=1)


# makes a counterfactual covariate matrix for the future

type = "stay"

if (type == "stay") {
  if (days_average == 7) { # if we already averaged, we don't need to average, just take last row
    last_row = tail(covars, 1)  
  } else {
    last_row = tail(covars %>% filter(date <= covars_last_date), 7)  %>% 
      summarise_all(mean)
  }
  
  forecast_color = "orange"
} else if (type == "decrease") {
  last_row = covars %>% filter(date == "2020-02-24")
  forecast_color = "red"
} else if (type == "increase") {
  last_row = covars %>% filter(date == "2020-04-14")
  forecast_color = "blue"
}

new_covars = covars %>% 
  bind_rows(purrr::map(new_dates, ~ {
    D = last_row
    D$date = .x
    D}) %>% 
      bind_rows()
  ) %>% 
  mutate(day = 0:(n() - 1)) %>%
  mutate(weekend = as.numeric(chron::is.weekend(date))) %>%
  mutate( future = ifelse( date > covars_last_date, 1, 0))



sim_future_states =  pomp(
  data=filter(hospitalization_data, day==maxT - 1),
  times="day", 
  t0 = maxT - 1,
  rinit=rinit,
  rmeasure=covid_rmeas,
  dmeasure=covid_dmeas_admitdischarge_nb, 
  covar=covariate_table(select(new_covars, -date), times="day"),
  rprocess = pomp::euler(covid_rprocess,dt),
  statenames= covid_statenames,
  paramnames= covid_paramnames,
  params = coef(mf2),
  accumvars = covid_accum_vars
) %>% 
  pomp:::simulate(nsim=N_smooth_draws, 
                  format="arrays",
                  rinit=function(...) setNames(filtered_states[ , maxT, sample(N_smooth_draws, 1)],
                                               dimnames(smooth_dist[[1]]$states)$variable),
                  statenames=covid_statenames,
                  times=maxT:(maxT + num_forecast_days - 1)) %>% 
  .[["states"]] %>% 
  aperm(c(1, 3, 2)) 


filtered_states_new = abind::abind(filtered_states, sim_future_states, along=2)
my_dim_names = dimnames(smooth_dist[[1]]$states)
my_dim_names$sample = 1:length(smooth_dist)

my_dim_names$time = c(my_dim_names$time, maxT:(maxT + num_forecast_days - 1))
dimnames(filtered_states_new) = my_dim_names

# resample indices and get smoothed states
sample_ind = pomp::systematic_resample(smoothing_weights)
smoothed_states = filtered_states_new[,,sample_ind]

}


library(abind)
add_r0_rt = function( smoothed_states, clust.n = 1 ) {
   if( clust.n > 1) {
      my.cl = makeCluster(clust.n, outfile="" )
      registerDoParallel( my.cl )
    }
# compute r0
source("code/common/nextgen.R")
source("code/common/NGM_R_with_S.R")

pop_totals = colSums(pop)
pop_ratio_mat = get_pop_ratio(pop_totals)

myfun =  function(beta)         get_r0(beta, (phi), pop_totals, pop_ratio_mat, sigma, tau, rho_A, rho_A, 
                                      gamma_A, gamma_Y, omega_A, omega_Y, omega_P_overall)
myfun2 = function(beta,S_prop) get_rt(beta, (phi), pop, pop_ratio_mat, sigma, tau, rho_A, rho_Y, 
                                      gamma_A, gamma_Y, omega_A, omega_Y, omega_P, S_prop) 


  cat("r0\n")
  r0 = foreach(i=1:dim(smoothed_states)[3], 
    .combine = cbind,
    .export = c("get_r0", "sigma", "tau", "phi", "rho_A", "rho_Y",
                "gamma_A", "gamma_Y", "omega_A", "omega_Y", "omega_P_overall")
  ) %dopar% {
    r0_i = purrr::map_dbl(smoothed_states["Beta", , i], myfun)
  }
  
  state.names = gsub("_1_1","",dimnames(smoothed_states)[[1]][grep("^[A-Z]+_1_1",dimnames(smoothed_states)[[1]])])
  all.states = dimnames(smoothed_states)[[1]][grep("^[A-Z]+_[0-9]+_[0-9]+",dimnames(smoothed_states)[[1]])]
  S.states = matrix( paste( "S",matrix(1:2,2,5), matrix(1:5,2,5,byrow = T), sep="_"), 2, 5)
  cat("rt\n")

  rt  = foreach(i=1:dim(smoothed_states)[3], 
      .combine = cbind,
      .export = c("get_rt", "sigma", "tau", "phi", "rho_A", "rho_Y",
                "gamma_A", "gamma_Y", "omega_A", "omega_Y", "omega_P"),
      .packages = c("magrittr")
          ) %dopar% {
    Age.Pop.N = sapply(1:5,function(g) smoothed_states[ all.states[grep(paste0("_",g,"$"),all.states)]  ,,i] %>% colSums ) 
    Age.S.N = sapply(1:5,function(g) smoothed_states[ S.states[,g]  ,,i] %>% colSums ) 
    Age.S.prop = Age.S.N / Age.Pop.N
    rt_i = sapply( seq_len(nrow(Age.S.N)), function(j) myfun2(beta=smoothed_states["Beta",j , i], S_prop=Age.S.prop[j,])   )
  }
  if( clust.n > 1) {
    stopCluster(my.cl)
  }
  r0 = array( r0, dim = c(1,dim(smoothed_states)[2:3]))
  rt = array( rt, dim = c(1,dim(smoothed_states)[2:3]))
  dimnames(r0)=c(list("R0"), dimnames(smoothed_states)[2:3])
  dimnames(rt)=c(list("Rt"), dimnames(smoothed_states)[2:3])
  abind( smoothed_states, r0, rt, along=1)
}

