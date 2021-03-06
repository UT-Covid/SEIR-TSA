###########################################
## Processing Hospitalization Data
###########################################
library(purrr)
library(zoo)


get_hospitalization_data <- function(nyc = FALSE){
  require(googlesheets4)
  require(tidyverse)
  require(lubridate)
  # browser()
  if(!nyc){
    # df <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/16QcLWs5kytxSjiZr8Ze6RtILvReLRh1FGysSXAJRncc/edit?usp=sharing",
    #                                 sheet=1) %>% 
    #   mutate(zt = c(0, diff(hospitalized))) %>% 
    #   select(date, hospitalized, zt)
    # df <- df %>% 
    #   filter(ifelse(seq_along(hospitalized) %in% 1:(rle(df$hospitalized)$lengths[1]-1),
    #                 FALSE, TRUE))  
    df <- data.frame(
      new_admits = c(1L,0L,1L,0L,0L,1L,1L,1L,
                     2L,3L,1L,3L,4L,3L,5L,12L,8L,4L,9L,9L,13L,9L,
                     10L,8L,10L,7L,7L,5L,14L,11L,4L,4L,14L,10L,
                     10L,7L,7L,11L,8L,5L,13L,9L,6L,10L,9L,3L,8L,7L,
                     9L),
      new_discharges = c(0L,0L,0L,0L,0L,0L,0L,0L,
                         0L,1L,0L,1L,1L,0L,1L,6L,6L,1L,3L,4L,5L,3L,
                         7L,12L,3L,8L,4L,5L,10L,11L,7L,6L,7L,8L,9L,8L,
                         6L,12L,5L,4L,12L,5L,9L,13L,16L,0L,3L,5L,7L),
      hospitalized = c(1L,1L,2L,2L,2L,3L,4L,5L,
                       7L,9L,10L,12L,15L,18L,22L,28L,30L,33L,39L,44L,
                       52L,58L,61L,57L,64L,63L,66L,66L,70L,70L,67L,
                       65L,72L,74L,75L,74L,75L,74L,77L,78L,79L,83L,
                       80L,77L,70L,73L,78L,80L,82L),
      date = c("2020-03-10","2020-03-11",
               "2020-03-12","2020-03-13","2020-03-14","2020-03-15",
               "2020-03-16","2020-03-17","2020-03-18","2020-03-19",
               "2020-03-20","2020-03-21","2020-03-22","2020-03-23",
               "2020-03-24","2020-03-25","2020-03-26","2020-03-27",
               "2020-03-28","2020-03-29","2020-03-30","2020-03-31",
               "2020-04-01","2020-04-02","2020-04-03","2020-04-04",
               "2020-04-05","2020-04-06","2020-04-07","2020-04-08",
               "2020-04-09","2020-04-10","2020-04-11","2020-04-12",
               "2020-04-13","2020-04-14","2020-04-15","2020-04-16","2020-04-17",
               "2020-04-18","2020-04-19","2020-04-20","2020-04-21",
               "2020-04-22","2020-04-23","2020-04-24","2020-04-25",
               "2020-04-26","2020-04-27")
    ) %>% 
      as_tibble() %>% 
      mutate(date = lubridate::ymd(date),
             zt = c(0, diff(hospitalized)))
  } else {
    df <- read_csv("https://raw.githubusercontent.com/nychealth/coronavirus-data/master/case-hosp-death.csv") %>%
      mutate(date = lubridate::mdy(DATE_OF_INTEREST),
             new_admits = HOSPITALIZED_COUNT) %>% 
      select(date, new_admits) %>% 
      filter(!is.na(new_admits))
  }
  return(df)
}


get_data_and_parms <- function(days_to_fill_in, fitting_data, parms){

    data_to_fit <- tibble(date = rev(min(fitting_data$date) - days(1:days_to_fill_in)),
           hospitalized = 0, 
           new_admits = 0,
           new_discharges = 0,
           zt = 0) %>% bind_rows(fitting_data) %>% 
      mutate(day = seq_along(date)-1)
    
    return(list(data = data_to_fit,
                parms = parms))
}



get_hospitalization_data_old <- function(nyc = FALSE){
  require(googlesheets4)
  require(tidyverse)
  require(lubridate)
  # browser()
  if(!nyc){
    df <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/16QcLWs5kytxSjiZr8Ze6RtILvReLRh1FGysSXAJRncc/edit?usp=sharing",
                                    sheet=1) %>% 
      mutate(zt = c(0, diff(hospitalized))) %>% 
      select(date, hospitalized, zt)
    df <- df %>% 
      filter(ifelse(seq_along(hospitalized) %in% 1:(rle(df$hospitalized)$lengths[1]-1),
                    FALSE, TRUE))  
  } else {
    df <- read_csv("https://raw.githubusercontent.com/nychealth/coronavirus-data/master/case-hosp-death.csv") %>%
      mutate(date = lubridate::mdy(DATE_OF_INTEREST),
             hospitalized = HOSPITALIZED_CASE_COUNT) %>% 
      select(date, hospitalized) 
    df <- tibble(date = min(df$date)-lubridate::days(1),
           hospitalized = 0) %>% 
      bind_rows(df) %>% 
      mutate(zt = c(0, diff(hospitalized)))
  }
  return(df)
}


get_data_and_parms_reg <- function(days_to_fill_in, fitting_data, parms){
  ## Maybe want to remove the added data eventually (make NAs)
  data_to_fit <- tibble(date = rev(min(fitting_data$date) - days(1:days_to_fill_in)),
         hospitalized = 0, 
         zt = 0) %>% bind_rows(fitting_data) %>% 
    mutate(day = seq_along(date)-1) %>% 
    mutate(date = lubridate::ymd(date))

  return(list(data = data_to_fit,
              parms = parms))
}

replace_parms <- function(new_vals, parms){
  val_names <- names(new_vals)
  for(i in 1:length(new_vals)){
    parms[val_names[i]] <- new_vals[i]  
  }
  return(parms)
}




fill_dates = function(df, dates_from, days_average = 7) {
  df %>% 
    arrange(date) %>% 
    as.list %>% 
    lapply( FUN=function(x){ 
      x%>%na.trim%>% tail(n=days_average) %>% mean ->m.e
      x%>%na.trim%>% head(n=1) ->m.s 
      x %>% na.fill(fill=c(m.s,"extend",m.e))  }) %>% 
    as.data.frame   -> df
  missing_dates = dates_from[!(dates_from %in% unique(df$date))]
  # fill prev dates
  oldest = df %>% 
    arrange(date) %>% 
    head(1)
  missing_dates_prev = missing_dates[missing_dates < min(df$date)]
  prev = missing_dates_prev %>% 
    purrr::map(~ bind_cols(tibble(date=.x), select(oldest, -date))) %>% 
    bind_rows()
  # fill newest dates
  newest_av = df %>% 
    arrange(date) %>% 
    na.trim() %>%
    tail( days_average ) %>%
    summarise_all(mean)
  missing_dates_new = missing_dates[missing_dates > max(df$date)]
  new = missing_dates_new %>% 
    purrr::map(~ bind_cols(tibble(date=.x), select(newest_av, -date))) %>% 
    bind_rows()
  bind_rows(prev, df, new)
}
