run.name="7-9-Matt.data"
#date.range=c( 5.11, 8.20)
date.range=c( 4.1)

num_forecast_days = 60
global.dt=1.0

# Timing - record how long each part takes
timing=list()
timing$start=Sys.time()


library(session)

library(foreach)
library(doParallel)

library(tidyverse)
library(readxl)
library(pomp)
library(lubridate)
library(magrittr)

#png()

#################################
#read command line arguments
args = commandArgs(trailingOnly=TRUE)

print( "args:")
print(args[])


if( grepl("frontera",system("hostname",intern=T)) ) {
  clust.n = 30
} else if( grepl("stampede",system("hostname",intern=T)) ) {
  clust.n = 90
} else {
  clust.n = 14
}

if( length(args) == 0) {
  indexes.to.run = 3
} else  {
  indexes.to.run = args[1] %>% strsplit(split=",") %>% unlist %>% as.numeric
}  
#pdf()

# calculating r0 and rt
source("code/common/nextgen.R")
source("code/common/NGM_R_with_S.R")
source("code/common/small_func.R")

# Fill in the data --------------------------------------------------------
# Read in covariates

source("code/common/tsa_mobility.R")
source("code/common/data-processing.R")





# read hosp data
#TSA = "Amarillo"
#TSA = "Lubbock"
# TSA = "Witchita Falls"
#TSA = "Abilene"
#TSA = "Dallas/Ft. Worth"
#TSA = "Paris"
#TSA = "Longview/Tyler"
#TSA = "Lufkin"
#TSA = "El Paso"
#TSA = "Midland/Odessa"
#TSA = "San Angelo"
#TSA = "Belton/Killeen"
#TSA = "Waco"
#TSA = "Bryan/College Station"
#TSA = "Austin"
#TSA = "San Antonio"
#TSA = "Houston"
TSA = "Galveston"
#TSA = "Victoria"
#TSA = "Laredo"
#TSA = "Corpus Christi"
#TSA = "Lower Rio Grande Valley"

mob = get_mobility_data_tsa( min.date = da(4.10), rerun=c() ) 
# %>%
#    filter( date < ymd("2020-12-18"))

a=read.csv("data/bed_and_vents.csv",as.is=T)
TSAs = unique(a$location) %>% sort
TSAs.let = sapply(TSAs, function(x) a$tsa[a$location==x][1] )

################# for loop
for( TSA in TSAs[indexes.to.run])  {
  cat("TSA = ", TSA, "\n") ;
  TSA.let = TSAs.let[TSA]
  source("code/models/model_parameters.R") 
  source("code/TSA_data.R")
#  hospitalization_data[ hospitalization_data$date == da(7.27),"adm_total"] = NA 




# if we want to read Austin data
if(F) {
      table.file = "../seir_regression/data/UT Template for COVID Hospital Admissions + Discharges 7.22.2020.xlsx"

      hospitalization_data = 
          read_xlsx(table.file) %>% 
        rename(date=`Date (yyyy-mm-dd)`,
              adm_total=`New Admits`,
              recov_total=`New Alive Discharges`,
              deaths_total=`New Dead Discharges`) %>% 
        mutate(date=ymd(date))
        
      hospitalization_data$hospitalized = (hospitalization_data %$% cumsum( adm_total - deaths_total - recov_total))
        
      hospitalization_data %<>% 
        filter( date >= da(4.10) ) 
}

  hospitalization_data %<>% 
  #  left_join( dat, by="date")%>%
    mutate( hosp_dif = diff(c(0,hospitalized))) %>%
    mutate( leave_total = adm_total - hosp_dif) %>%
        filter( date > da( date.range[1]))


  hospitalization_data  %<>% 
  #  filter( date > da(6.01)) %>%
    mutate( day=1:(n() )  )

#  if( TSA == "Lower Rio Grande Valley") {
#      hospitalization_data$adm_total[ hospitalization_data$date > da(7.21) ] = NA
#  }
#  if( TSA == "Lower Rio Grande Valley") {
#      hospitalization_data$adm_total[ hospitalization_data$date > da(7.27) ] = NA
#  }

  mob.tsa = mob %>% 
    filter( TSA == TSA.let) %>%
#    mutate(PC1=1,PC2=1,PC3=1) %>%
    mutate( PC1=rolling(PC1),
            PC2=rolling(PC2),
            PC3=rolling(PC3) ) %>%
    select(date, PC1, PC2, PC3, weekend=weekend, cases=case_count) %>%
    fill_dates(dates_from=hospitalization_data$date)

  covars = hospitalization_data %>% 
              { data.frame( date=daS( min(.$date)-1, max( .$date)), 
                            day =seq_along ( daS( min(.$date)-1, max( .$date)) )-1) 
              } %>%
      mutate( maxT   = n() ) %>%
      mutate( future = 0  ) %>% 
      left_join( mob.tsa, by="date")




  source("code/common/Csnippet_service.R")

  # read ratio of states for initialization at first time point
  load(file="initial_state_ratio.Rda")
  state_ratio %<>% {. / sum(.)}

  init_E_ratio  = state_ratio[expand_state( "E" ,c(2,5))] %>% matrix( ., 2,5)
  init_PA_ratio = state_ratio[expand_state( "PA",c(2,5))] %>% matrix( ., 2,5)
  init_PY_ratio = state_ratio[expand_state( "PY",c(2,5))] %>% matrix( ., 2,5)
  init_IA_ratio = state_ratio[expand_state( "IA",c(2,5))] %>% matrix( ., 2,5)
  init_IY_ratio = state_ratio[expand_state( "IY",c(2,5))] %>% matrix( ., 2,5)
  init_H_ratio  = state_ratio[expand_state( "H" ,c(2,5))] %>% matrix( ., 2,5)



  source("code/run.parts.R")

  print(TSA)
  print(tail(hospitalization_data))

  mf2 = find.mif( data=hospitalization_data, covars=covars, N.steps = 300, N.particles = 3500, clust.n = clust.n, dt=global.dt)
  #mf2 = find.mif( data=hospitalization_data[1:10,], covars=covars[1:10,], N.steps = 1, N.particles = 2, clust.n = 1, dt=global.dt)
  cat("mf2 done\n")

  save.session(file=paste0("after.mf2.Rda"))
  #library(session); restore.session(file=paste0("after.mf2.Rda"))

  ## Estimate a smoothed distribution over beta --------------------------------

  filtered_states = make_smoothed_states( mf2, covars=covars , N.draws = 500 , clust.n=clust.n, dt=global.dt)
  cat("smooth done\n")

  save.session(file=paste0("after.filter.Rda"))
  #run.name="6-26-simple-4";library(session); restore.session(file=paste0("after.filter.",run.name,".Rda"))




  # ------------------------------------------------------------------------
  # Augment with forecasts counterfactual using sim
  cov.range = range( covars$date)

  # Modify for forecasting
  new_dates = seq( from= cov.range[1], to = cov.range[2]+ num_forecast_days,  by=1)



  mob.tsa.future = mob %>% 
    filter( TSA == TSA.let) %>%
    select(date, PC1, PC2, PC3, weekend=weekend) %>%
    fill_dates(dates_from=new_dates) 


  
  new_covars = hospitalization_data %>% 
              { data.frame( date=daS( min(.$date)-1, max( .$date)+num_forecast_days), 
                            day =seq_along ( daS( min(.$date)-1, max( .$date)+num_forecast_days) )-1) 
              } %>%
      mutate( maxT   = n() ) %>%
      mutate( future = date > max(hospitalization_data$date)  ) %>% 
      left_join( mob.tsa.future, by="date") %>%
      filter( date >= da(4.10) )

  future_covars = tail(new_covars, num_forecast_days+2)


  smoothed_states = add_future_states(  coefs= coef(mf2), dt = global.dt,
        data = tail( hospitalization_data,1 ),
        covars = future_covars,
        filtered_states =  filtered_states
  )

  smoothed_states = add_r0_rt( smoothed_states, clust.n = clust.n  )






  TSA.4f = gsub(pat="/",rep="-",TSA)
  save( TSA, mf2, hospitalization_data, covars, new_covars, smoothed_states, file=paste0("RDA/res-10.8-",TSA.4f,".Rda"))
 
 
}
