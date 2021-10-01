
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)

## County populations (some cleaning required...)
countypop <- read.csv("../data/county-pop.csv", stringsAsFactors = FALSE) %>%
  mutate(county = str_sub(county, 2, -1),
         state = str_split_fixed(county, ",", 2)[, 2] %>% #State name
           str_sub(2, -1),
         county = gsub("(.*),.*", "\\1", county), #Remove state name
         county = str_remove(county, " County"),
         county = str_remove(county, " Parish"),
         county = str_remove(county, " city"),
         pop2019 = pop2019 %>% #Convert string to integer
           str_remove_all(",") %>%
           as.integer()) %>%
  select(county, state, pop2019)

glimpse(countypop)

## List of MSAs

## Note: CBSA includes all MSAs and \mu-SAs;
## This .csv file contains an indicator for which type the CBSA is
## Large MSAs are sometimes divided into Metropolitan Divisions
## CSAs are combinations of MSAs or \mu-SAs

cbsa <- read.csv("../data/csa.csv", stringsAsFactors = FALSE)
  
glimpse(cbsa)

## Counties within Houston MSA
cbsa %>% filter(Metropolitan.Micropolitan.Statistical.Area == "Metropolitan Statistical Area",
               CBSA.Title == "Houston-The Woodlands-Sugar Land, TX")

## Narrow down to counties of MSAs, without Puerto Rico
msa <- cbsa %>%
  filter(Metropolitan.Micropolitan.Statistical.Area == "Metropolitan Statistical Area",
         State.Name != "Puerto Rico")

## Create 5-digit fips codes for counties, select columns, simplify
## county/parish names, and add population counts
msa <- msa %>%
  mutate(fips = FIPS.State.Code * 1000 + FIPS.County.Code) %>%
  select(msa = CBSA.Title,
         county = County.County.Equivalent,
         state = State.Name,
         central.outlying = Central.Outlying.County,
         fips) %>%
  mutate(county = county %>%
           str_remove(" County") %>%
           str_remove(" city") %>%
         str_remove(" Parish")) 

## Add population counties
msa <- msa %>% left_join(countypop)

## All counties in MSA dataset have a population
sum(is.na(msa$pop2019 ))

glimpse(msa)

## Read in Counties death count data
us_counties = read.csv('https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv', stringsAsFactors = FALSE)

us_counties$state %>% unique()

## Check overlap in county fips
Match <- (msa$fips %in% us_counties$fips)

mean(Match)

## Besides the five counties of NYC, most counties in MSAs not present
## in the NYTimes county case data are pretty small, so I'm assuming
## these are counties that have not had a case yet
msa[!Match, ] %>%
  select("msa", "county", "state", "pop2019", "central.outlying") %>%
  arrange(desc(pop2019)) %>%
  head(31)


## Note that the NYTimes data combines the five counties.  This is a
## unique case, because NYC is the only city to be split into separate
## counties [citation needed]
us_counties %>%
filter(state == "New York",
       str_detect(county, "New York"))

## Create population weights
msa <- msa %>%
  group_by(msa) %>% 
  mutate(county.weight = pop2019 / sum(pop2019)) %>%
  ungroup()

###############################################################################
                                        #          County basic data          #
###############################################################################


## county_income <- read.csv("../data/counties/income_unemployment.csv", stringsAsFactors = FALSE)

## glimpse(county_income)

## county_education <- read.csv("../data/counties/education.csv", stringsAsFactors = FALSE)

## glimpse(county_education)



## msa2 <- msa %>%
##   left_join(county_income)

## msa3 <- msa %>%
##   left_join(county_education)

## ## Check for completeness of mapping
## msa2[is.na(msa2$State), c("county", "state")]

## msa3[is.na(msa2$State), c("county", "state")]
