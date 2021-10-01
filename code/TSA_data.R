library(openxlsx)
source("code/common/small_func.R")
source("code/common/Csnippet_service.R")



#TSA = "Amarillo"
#TSA = "Lubbock"
#TSA = "Witchita Falls"
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
# TSA = "Austin"
#TSA = "San Antonio"
#TSA = "Houston"
#TSA = "Galveston"
#TSA = "Victoria"
#TSA = "Laredo"
#TSA = "Corpus Christi"
#TSA = "Lower Rio Grande Valley"


if( T) {
  a=read.csv("data/bed_and_vents.csv",as.is=T)
  a$date %<>% gsub(pat="T00:00:00Z",rep="") %>% ymd

  i = a$location == TSA
  b=a[i,]
  tsa.let = b$tsa[1]

  hospitalization_data=data.frame( 
      date=b$date, 
      adm_total = b$susp_covid_gen_admitted_24, 
      hospitalized = b$total_laboratory_confirmed,
      icu          = b$lab_con_icu )

  if( F ) {
  b=read.csv("data/new_hosp_data_dshs.csv")
  b$date %<>% ymd
  i = b$tsa == tsa.let
  b=b[i,]
  b$tot_confirmed_covid_admits_last_24hrs[ is.na(b$tot_confirmed_covid_admits_last_24hrs)] = 0
  hospitalization_data2=data.frame( 
      date=b$date, 
#      adm_total = b$tot_confirmed_covid_admits_last_24hrs + b$tot_suspected_covid_admits_last_24hrs, 
      adm_total = b$tot_confirmed_covid_admits_last_24hrs, 
      hospitalized = b$tot_confirmed_covid_patients_in_hosp,
      icu          = NA )
  }


  sheets = getSheetNames("data/TexasCOVID-19HospitalizationsOverTimebyTSA.xlsx")
  # Read DSHS data in xcel sheet
  if( "COVID-19 Hospitalizations" %in% sheets) {
    x=read.xlsx("data/TexasCOVID-19HospitalizationsOverTimebyTSA.xlsx",startRow=3, sheet="COVID-19 Hospitalizations")
  } else {
    stop( "Couldn't find sheet 'COVID-19 Hospitalizations'")
  }  
  x=x[1:22,]
  d=daS( from=4.12, len=dim(x)[2]-2)
  x.let = x[,1] %>% as.character %>% gsub(pat="[.]",rep="")
  i= x.let == tsa.let 
  h = x[i,-(1:2)] %>% unlist
  if( "COVID-19 ICU" %in% sheets ) {
    x=read.xlsx("data/TexasCOVID-19HospitalizationsOverTimebyTSA.xlsx",startRow=3,sheet="COVID-19 ICU")
  } else if( "Adult COVID-19 ICU" %in% sheets) {
    x=read.xlsx("data/TexasCOVID-19HospitalizationsOverTimebyTSA.xlsx",startRow=3,sheet="Adult COVID-19 ICU")
  } else {
    stop( "Couldn't find 'COVID-19 ICU' sheet or 'Adult COVID-19 ICU' sheet" )
  }  
  x=x[1:22,]
  x.let = x[,1] %>% as.character %>% gsub(pat="[.]",rep="")
  i= x.let == tsa.let 
  icu = x[i,-(1:2)] %>% unlist

  hospitalization_data3 = data.frame(
    date = d,
    adm_total = NA,
    hospitalized = h,
    icu=icu
  )

  i = hospitalization_data$date < da(7.22)
  hospitalization_data = hospitalization_data[i,]
  j = !(hospitalization_data3$date %in% hospitalization_data$date) # Those of hosp3 that are not in hosp
  hospitalization_data = rbind( hospitalization_data, hospitalization_data3[j,])
  hospitalization_data %<>% { .[order(.$date),] }

} else {
   
   
   library(openxlsx)
   library(magrittr)
   source("code/common/small_func.R")
   a=read.xlsx("data/Texas COVID-19 Hospitalizations by TSA.xlsx")
   
   i=which(a[,2]==TSA)
   tsa.let = a[i,1] %>% as.character %>% gsub(pat="[.]",rep="")
   
   x = a[i,-(1:2)] %>% t
   hospitalization_data=data.frame(date = daS(4.08,6.25), hospitalized=c(x))
   
}

a = read.csv("data/TSAPopulation_RiskGroups.csv")

pop = a[a$TSA==tsa.let,-(1:2)] %>% as.matrix

totalpop = colSums(pop)

YHR_overall = (YHR[1, ] * pop[1, ] + YHR[2, ] * pop[2, ]) / colSums(pop)

#// one row per risk group
# high_risk_ratio = c(0.082825, 0.141121, 0.165298, 0.329912, 0.470568)
high_risk = round(high_risk_ratio * totalpop)

relpop = totalpop / sum(totalpop)

 omega_P = P / (1.0 - P) *
   (tau * omega_Y * (YHR / eta + (1 - YHR) / gamma_Y) + (1 - tau) * omega_A / gamma_A) *
   (1 / (tau * omega_Y / rho_Y + (1 - tau)  * omega_A / rho_A))
 omega_P_overall = P / (1.0 - P) *
   (tau * omega_Y * (YHR_overall / eta + (1 - YHR_overall) / gamma_Y) + (1 - tau) * omega_A / gamma_A) *
   (1 / (tau * omega_Y / rho_Y + (1 - tau)  * omega_A / rho_A))



vacc = local({
  # diff starting from 0 - keeps vector length
  diff0 = function(x) {
    x=c(0,x)
    diff(x) %>%  ifelse( .<0, 0, .)
  }

  d = hospitalization_data$date

  vacc.one  = read.csv("data/vaccination/TSA_vaccine_one_dose.csv")[,-1] %>% filter(tsa==tsa.let) %>%
            mutate( across( !c("date","tsa"), diff0) ) 
            
  vacc.one$date = ymd( vacc.one$date)

  vacc.full = read.csv("data/vaccination/TSA_vaccine_fully_vaccinated.csv")[,-1] %>% filter( tsa==tsa.let) 
  vacc.full$date = ymd( vacc.full$date)

  min.d = min(d)

  vacc.one %<>% filter( date > min.d ) 
  vacc.full %<>% filter( date > min.d ) 

  colnames(vacc.one) %<>% gsub(pat="^X",rep="V1_")
  colnames(vacc.full) %<>% gsub(pat="^X",rep="V2_")


  vac = data.frame( vacc.one[,-2], vacc.full[,-(1:2)])
    # Add empty rows before first vaccination
  min.vac.date = as.Date(min( vac$date))

  if( min.d < min.vac.date) {
    empty_vac = data.frame( date = seq( from=min.d, to=min.vac.date-1, by=1))
    empty_vac[,colnames(vac)[-1]] = 0
    vac = rbind(empty_vac,vac)
  }
  expand_state("nV1",c(5,2)) %>% order %>% order ->o

  colnames(vac) = c("date", expand_state("nV1",c(2,5))[o], expand_state("nV2",c(2,5))[o] )
  vac
})
