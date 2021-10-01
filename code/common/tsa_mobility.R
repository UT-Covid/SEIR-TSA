library(magrittr)
library(tidyverse)
library(lubridate)

source("code/common/conditional_run.R")

get_mobility_data_tsa = function (average_mode=c("trailing", "centered"),
                              days_average=1, rerun=c() , min.date = da(4.11) ) 
{
  ## condition by TACC





SG.dir=c( "/work2/projects/utprojections/safegraph_data/final-metrics/", "data/safegraph_data/final_metrics/" )
visits.f = "visits_by_category_bycounty_baselined.csv"
sdmetric.f = "sd-metrics_bycounty_baselined.csv"

i=1
repeat {
  if( file.exists( paste0( SG.dir[i],visits.f)) & file.exists( paste0( SG.dir[i], sdmetric.f))  ) {
    SG.dir=SG.dir[i]
    break
  }
  if( i == length(SG.dir)) {
    stop( paste("couldn't find",visits.f,"or",sdmetric.f,"in",SG.dir,sep=" ",collapse="\n"  ))
  }
  i = i+1
}

SG.dir=SG.dir[i]
visits.f %<>% paste0(SG.dir,.)
sdmetric.f %<>% paste0(SG.dir,.)




  county_population.url = "https://raw.githubusercontent.com/JieYingWu/COVID-19_US_County-level_Summaries/master/data/counties.csv" 
  rural_urban.f = "data/counties_basicdata.csv"
  TSA_county_table.f = "data/TSA_counties_csv.csv"
  res.f = "data/tsa_mobility_data.Rda"

   #### Check if data in files is newer than our latest data
   if( length(rerun)==0) {
      need.run =  check.files.for.rerun.needed( files=c( visits.f, sdmetric.f), data.file = res.f)
   }


  ## if not, just load it from before
  
  if( (length(rerun)==0 && !need.run) || (length(rerun)>0 && !rerun) ) {
      load( res.f)
      return( tsa_data)
  } 


  
  cases_data = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv" %>% 
    read_csv()    %>%
    pivot_longer(cols = ends_with("20") | ends_with("21") ,
                 names_to = "date_mdy",
                 values_to = "cumulative_cases") %>%
    mutate(date = mdy(date_mdy) %>% ymd()) %>%
    rename( fips=FIPS) %>%
    mutate( fips = as.character(fips))
  
  
  

  average_mode=c("trailing")
  days_average = 1
    
  top_categories <- c(
      "Colleges, Universities, and Professional Schools"   ="colleges",   
      "Drinking Places (Alcoholic Beverages)"              ="drinking",
      "Elementary and Secondary Schools"                   ="schools",
      "General Medical and Surgical Hospitals"             ="medical",
      "Grocery Stores"                                     ="grocery",
      "Museums, Historical Sites, and Similar Institutions"="museums/parks",
      "Restaurants and Other Eating Places"                ="restaurants"
  )
 
  
  county_population = county_population.url  %>% 
    data.table::fread() %>% 
    as_tibble() %>% 
    select(fips=FIPS, state_short=State, county=Area_Name, pop=POP_ESTIMATE_2018) %>% 
    mutate(fips=sprintf("%05d", fips))
  
  rural_urban = rural_urban.f %>%
    read_csv(.) %>% tail(n=-1) %>% 
    mutate(fips = sprintf("%05d", fips)) %>% 
    mutate(rural_urban = Rural.urban_Continuum.Code_2013 / 9) %>% 
    select(fips, rural_urban)


  visits =    read_csv( visits.f   )[ -1]
  sdmetrics = read_csv( sdmetric.f )[,-1]
  


  cp_TX = county_population %>%
      filter( state_short == "TX")    %>% 
      mutate(county = county %>%
              str_remove(" County") )

  TSA_county_table =read.csv(TSA_county_table.f,header = F,as.is=T)

  tsa = lapply(1:dim( TSA_county_table )[1],function(i) {
      x= strsplit(TSA_county_table[i,3],split=",")[[1]] %>%
          gsub(pat="^ ",rep="",perl=T)
      data.frame( TSA=TSA_county_table[i,2],county=x)
  } )  %>% 
      Reduce( f=rbind) %>%
      mutate( TSA = gsub(TSA,pat=" ",rep="") )
      

  tsa   %<>% 
      left_join( cp_TX, by="county" )

  tsa %<>%     left_join( group_by(., TSA) %>% 
                summarize(
                  tsa_pop=sum(pop)) %>% 
                ungroup() %>% 
                mutate(pop_rank=rank(-tsa_pop)),
        by="TSA") %>%
      mutate(county_weight = pop / tsa_pop)
  

  county_poi = visits %>% 
    # A bit of conversion and renaming to make things nicer
    mutate(date = ymd(date)) %>%
    rename(fips = fips5d,
          state_short = region,
          `county_full` = County) %>% 
    mutate(top_category = top_categories[top_category]) %>% 
    # Add the TSA label to each row, using fips
    filter(fips %in% unique( tsa$fips)) %>% 
    left_join(select( tsa, fips, TSA), by="fips") %>%

    group_by(date, top_category, TSA) %>% 
    summarize(pop = sum(pop),
              count = sum(count),
              count_baseline = sum(count_baseline)) %>% 
    mutate(count_pc = count / pop,
          count_relative = count / count_baseline) %>% 
    ungroup() %>%
    mutate( base.date = date > da(1.15) & date < da(2.15) )  %>%
    left_join( 
      group_by(., base.date, top_category, TSA, .add=T) %>% 
      summarize( count_base = mean(count)) %>%
      ungroup( ) 
    ) %>%
    filter( date >= min.date) %>%
    mutate( count_rel = count / count_base)
  
    cases_tsa = cases_data %>%
        filter( fips %in% unique(tsa$fips) ) %>%
      left_join(select( tsa, fips, TSA), by="fips") %>%
      group_by(date, TSA) %>% 
      summarize( case_count = sum(cumulative_cases)) %>%
      ungroup()
  


  tsa_totals = tsa %>% 
      select( TSA, state_short, tsa_pop ) %>% 
      distinct()    
  
  sd_cols = c("median_home_dwell_time_relative")

  sd = sdmetrics %>% 
      select(date, fips=fips5d, all_of(sd_cols)) %>% 
      left_join(county_population, by="fips") %>% 
      filter(fips %in% unique(tsa$fips)) %>% 
      left_join(select(tsa, fips, TSA, county_weight), by="fips") %>% 
      group_by(date, TSA) %>% 
      summarize(median_home_dwell_time_relative = sum(median_home_dwell_time_relative * county_weight))
      # summarize(median_home_dwell_time = sum(median_home_dwell_time * county_weight),
      #           full_time_work_behavior_devices = sum(full_time_work_behavior_devices))

# browser()
  tsa_data = county_poi %>%
    select(date, TSA, top_category, count_pc, count_rel) %>% 
    pivot_wider(names_from = top_category, values_from = c(count_pc, count_rel)) %>%
#    left_join(sd, by=c("date", "TSA")) %>% 
    right_join( cases_tsa, by=c("date", "TSA")) %>%
    # na.omit() %>% 
    # group_by(TSA) %>% 
    # group_modify(~ { # add 7-week average
    #   df = arrange(.x, date)
    #   # add T-day running averages
    #   cur_cols = c(paste0("count_relative_", unname(top_categories)), sd_cols)
    #   new_cols = paste0(cur_cols, "_av")
    #   for (c in new_cols)
    #     df[[c]] = NA_real_  # df[[c]] = NA gave errors - could not add doubles later.
    #   offset = ifelse(average_mode == "trailing", 0, days_average %/% 2 + 1)
    #   for (i in days_average:nrow(df)) {
    #     start_index = i - days_average
    #     end_index = i
    #     df[i - offset, new_cols] = t(apply(df[start_index:end_index, cur_cols, drop=FALSE], 2, mean)) # added t() to avoid errors.
    #   }
    #   # output of group modify
    #   df
    # }) %>% 
    # ungroup() %>% 
    left_join(tsa_totals, by="TSA") %>%
  #   left_join(death_data, by=c("msa", "date")) %>% 
  #   right_join(cases_data, by=c("msa", "date")) %>% 
  #   mutate(deaths_per_cap = deaths / msa_pop,
  #          cumulative_deaths_per_cap = cumulative_deaths / msa_pop) %>% 
#    left_join(
#      group_by(., TSA) %>%
#      summarize(threshold_day = date[min(which(cumulative_deaths_per_cap >= 3 / 1e7))]),
#      by = "msa") %>% 
#    mutate(days_since_thresh = as.integer(date - threshold_day)) %>% 
    mutate(weekend = chron::is.weekend(date)) # %>%


    # add pca
#    cur_cols = c(paste0("count_rel_", unname(top_categories)), sd_cols)
    cur_cols = c(paste0("count_rel_", unname(top_categories)))
    new_cols = paste0(cur_cols)
    browser()
    pca = prcomp(data.matrix(na.omit(log(tsa_data[ ,new_cols]))), scale=FALSE, rank=3)   # <-- log of relative components
#    pca = prcomp(data.matrix(na.omit(tsa_data[ ,new_cols])), center = TRUE, scale=TRUE, rank=3)
    # comps = data.matrix(msa_data[ ,row.names(pca$rotation)]) %*% pca$rotation
    comps = (predict(pca, newdata = (tsa_data)))  #<--- not exp back 
    tsa_data[ , paste0("PC", 1:3)] = comps

    save( tsa_data, file=res.f)

    tsa_data
                                                           
}















#death_data = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv' %>% 
#    read_csv() 





   
