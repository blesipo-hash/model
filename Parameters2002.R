library(readxl)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(viridis)
library(broom)
library(purrr)
library(sf)
library(sp)
library(raster)
library(ggmap)
library(rgdal)
library(tmap)
library(haven)
library(foreign)
library(rdhs)
library(Hmisc)

## Parameters 

#### Mean worm burden, calories and cost by commondity ####
setwd("/Users/Lexi/Documents/Data_Cost-Benefit_MDA/world-administrative-boundaries/") 
country_borders <- st_read("world-administrative-boundaries.shp") # shapefile with country borders
setwd("/Users/Lexi/Documents/Data_Cost-Benefit_MDA/shapes/") 
region_borders <- st_read("region.shp")%>% # shapefile with region borders
  st_as_sf()%>%
  mutate(REGION = ifelse(REGION == "Southeastern Asia", "South-Eastern Asia", REGION))%>% # changing region name to match all other files
  rename(Region=REGION)

country_borders <- country_borders %>% 
  st_as_sf() %>%
  mutate(name = ifelse(name =="Côte d'Ivoire", "Ivory Coast", #Changing country names to match other files
                       ifelse(name=="Bolivia", "Bolivia (Plurinational State of)", 
                              ifelse(name=="Democratic People's Republic of Korea", "Dem. People's Republic of Korea", 
                                     ifelse(name=="Guinea-Bissau", "Guinea Bissau",
                                            ifelse(name=="Vietnam", "Viet Nam",
                                                   ifelse(name=="The former Yugoslav Republic of Macedonia", "North Macedonia",
                                                          ifelse(name=="Lao People's Democratic Republic","Lao Peoples Democratic Republic",
                                                                 ifelse(name=="Libyan Arab Jamahiriya","Libya",
                                                                        ifelse(name=="Venezuela ", "Venezuela (Bolivarian Republic of)", name))))))))))%>%
           rename(country_title = name, Region=region)

setwd("/Users/Lexi/Documents/Data_Cost-Benefit_MDA/") 

intensity_infection<-read_excel("Parameters2.xlsx", sheet = "Intensity_Infection") #This sheet contains data on worm burden of each parasite.
parasite_parameters<-read_excel("Parameters2.xlsx", sheet = "Parameters2") # multiple studies for weight differences and blood loss (anthropometric studies) - all intestinal schisto treated the same
diet <- read_excel("Parameters2.xlsx", sheet = "Diet")
efficacy <- read_excel("Parameters2.xlsx", sheet = "drug_efficacy")
efficacy <- efficacy %>%
  dplyr::select(Parasite, Drug_Efficacy)
kcal_kg_blood = 1300 # 2* 650 - Nureye and Tekalign 2020
kcal_kg_tissue = 5000	#http://www.fao.org/3/M2885E/M2885E00.htm # number of calories needed to produce a kg of tissue. 
                                                                # Assuming these are the same number of calories needed to produce a liter of blood

food_balance <- read.csv("FoodBalanceSheets.csv") %>%  # Food balance sheets from FAO -> contain info on quantities of main commodities consumed in a country in a given year
  filter(Element == "Food supply (kcal/capita/day)")%>% 
  filter(Item != "Grand Total")
food_balance2 <-subset(food_balance,food_balance$Avg>0.00 & food_balance$Item.Code<2900)

  # take the mean of each food item across years 2014-2018 

TotValue <- food_balance2 %>% # the purpose of this pipe is to sum the values of all commodities by country and use this value as the denominator in the pipe below  
  ungroup()%>%
  group_by(Area)%>%
  summarise(TotValue = sum(Avg))

food_balance3 <- food_balance2%>% 
  left_join(TotValue, by = "Area")%>%
  group_by(Area, Item)%>%
  summarise(food_comp = Avg / TotValue)

food_balance4 <-food_balance3 %>%  
  subset(food_balance3$food_comp>0.005 )%>%
  filter(Item != "Miscellaneous")

CheckPercent <- food_balance4%>% 
  group_by(Area)%>%
  summarise(sum(food_comp))

#%>% # obtain the proportion of each food item by country


Food <- food_balance4%>% 
  group_by(Item)%>%
  summarise(sum(food_comp))
write.csv(Food,"Food2.csv")

setwd("/Users/Lexi/Documents/Data_Cost-Benefit_MDA/")
food_cals <- read.csv("Food_Cal.csv") 

food_price <- read.csv("food_price.csv") %>%
  group_by(country_title, Commodity)%>%
  summarise(USD_Value = mean(USD_Value))%>% # get mean value of commodities by country. Note these values are in tonnes and combine values from retail and wholesale
  rename(Item = Commodity) %>% # get mean value of commodities by country. Note these values are in tonnes and combine values from retail and wholesale
  rename(Area = country_title) # rename commodity to match food balance dataframe

# These lines estimate the mean worm burden of each parasite. 
burden_by_parasite <- intensity_infection %>%    
  group_by(Parasite) %>%
  summarise(mean_burden = mean(MeanBurden),
            sd_burden = sd(MeanBurden),
            burden_error = qnorm(0.975)*sd_burden/sqrt(n()), # this is assuming worm burden is normally distributed - it isn't and errors are one  of the major things needing fixing
            burden_high95 = mean_burden + burden_error,
            burden_low95 = mean_burden - burden_error)%>%
  dplyr::select(Parasite, mean_burden)


# These lines are for estimating the mean, SD and error of the anthropometric data
params_morb <- parasite_parameters %>%
  mutate(WeightDiff = as.numeric(WeightDiff),
         BloodLoss = as.numeric(BloodLoss),
         WeightDiff_SD = as.numeric(WeightDiff_SD),
         BloodLoss_SD = as.numeric(BloodLoss_SD),
         WD_error = qnorm(0.975)*WeightDiff_SD/sqrt(n_weightdiff),
         BL_error = qnorm(0.975)*BloodLoss_SD/sqrt(n_bloodloss),
         WeightDiff_low95 = WeightDiff - WD_error,
         WeightDiff_high95 = WeightDiff + WD_error,
         BloodLoss_low95 = BloodLoss - BL_error,
         BloodLoss_high95 = BloodLoss + BL_error)%>%
  replace_na(list(BloodLoss=0, WeightDiff=0,WD_error=0,BL_error=0))%>%
  group_by(Parasite)%>%
  mutate( kcal_year_wastedblood = BloodLoss * kcal_kg_blood,
          kcal_year_wastedtissue = WeightDiff * kcal_kg_tissue,
          kcal_year_wastedHost = kcal_year_wastedblood + kcal_year_wastedtissue,# values are summed for parasites that cause weight loss and blood loss 
         Tot_tissue_loss_error = rowSums(cbind(WeightDiff_SD, BloodLoss_SD), na.rm=TRUE),
         kcal_year_wastedHost_error = Tot_tissue_loss_error * kcal_kg_tissue)%>%
  left_join(efficacy, by = "Parasite")%>%
  dplyr::select(Parasite, kcal_year_wastedHost , kcal_year_wastedHost_error , Drug_Efficacy)


#### Epi data ####

pop <- read.csv("Pop_2000_2020.csv") # mean pop from 2000-2015
pop_by_age_year <- read.csv("annual_pop_by_age.csv") # https://population.un.org/wpp/Download/Standard/Population/ # pop by country and age classes for years 1990-2020
high_income <- read_excel("Parameters2.xlsx", sheet ="high_income") #World Economic Situation and Prospects 2022 United Nations
high_income <- high_income %>% # to identify high-income countries 
  mutate(income = "high")

pop <- pop %>%
  dplyr::select(country_title,Index, country_code,MeanPop, Region)%>%
  left_join(high_income, by="country_title")%>%
  replace_na(list(income="low-med"))%>%
  filter(!income == "high")%>%
  filter(!country_title == "Oman")
  

regions <- pop %>%
  dplyr::select(country_title, Region)%>%
  filter(!duplicated(country_title))


mean_pop <- pop_by_age_year %>%
  group_by(country_title)%>%
  summarise(Pop0_4 = mean(Pop0_4), # Mean pop estimates for 0-4 age class
            Pop5_14 = mean(Pop5_14), # Mean pop estimates for 5-14 age class
            Pop15plus = mean(Pop15plus))%>% # Mean pop estimates for 15 plus age class
  gather("AgeClass", "Pop", Pop0_4:Pop15plus, na.rm = TRUE)%>%
  mutate(AgeClass=gsub("Pop0_4","0_4", AgeClass), # rename age classes
         AgeClass=gsub("Pop5_14","5_14", AgeClass),
         AgeClass=gsub("Pop15plus","15plus", AgeClass))%>%
  left_join(high_income, by="country_title")%>%
  replace_na(list(income="low-med"))%>%
  filter(!income == "high")%>%
  filter(!country_title == "Oman")%>%# remove high-income countries
  left_join(regions, by = "country_title")

epi_helminth <- read.csv("epi_helminth2.csv") # data from Global Atlas of Helminth Infections (GAHI)
SCH_new<- read.csv("SCH_sitelevel.csv")
STH_new<- read.csv("STH_sitelevel.csv")

# given that prevalence data is not available for each age class/country combination, 
# I'm using a multiplicator from Lo et al 2016 
# that used an age class as a reference - typically the 5-14 age class 
# and provides a multiplicator to estimate the prevalence of the other age classes 
prev_multiplicator <- read_excel("Parameters2.xlsx", sheet = "Age-specific_multiplicator")%>% 
  dplyr::select(Parasite, AgeClass, Multiplicator)

# these lines are used to clean the GAHI database
# GAHI collects prevalence data from multiple studies which do not follow a standard age class classification
# so I made my own classification based on the age_end column. 
# If a study collected prev data on individuals not older than 4 these were classified as "Prev0_4"
# if age-end between 5 and 14, then "Prev5_14". If age_end = or greater than 15, the "Prev15plus"

SCH_2000 <- SCH_new %>%
  mutate(Year = SurveyYear, age_start=Age_start, age_end=Age_end,country_title=Country) %>%
  filter(age_start > 0, age_end > 0, Year >= 2000,Source!="GAHI") %>%
  mutate(AgeClass = ifelse(age_end <= 4, "Prev0_4",
                           ifelse(age_end >= 5 & age_end <= 14, "Prev5_14",
                                  ifelse(age_end >= 15, "Prev15plus", age_end))))%>%
  rename(Parasite = Species)%>%
  left_join(pop, by = "country_title")%>% # Note that "country_title" is the name I use to name countries in all data frames
  replace_na(list(income="high"))%>%
  filter(!income == "high")%>%
  filter(!country_title == "Oman")

STH_2000 <- STH_new %>%
  mutate(age_start=Age_start, age_end=Age_end,country_title=Country) %>%
  filter(age_start > 0, age_end > 0, Year >= 2000,Source!="GAHI") %>%
  mutate(AgeClass = ifelse(age_end <= 4, "Prev0_4",
                           ifelse(age_end >= 5 & age_end <= 14, "Prev5_14",
                                  ifelse(age_end >= 15, "Prev15plus", age_end))))%>%
  left_join(pop, by = "country_title")%>% # Note that "country_title" is the name I use to name countries in all data frames
  replace_na(list(income="high"))%>%
  filter(!income == "high")%>%
  filter(!country_title == "Oman")

Asc<-STH_2000 %>%
  mutate(Examined=Asc_examined, Positive=Asc_positive,Prevalence=Asc_prevalence)
Asc$Parasite<-"Ascaris"
HK<-STH_2000 %>%
  mutate(Examined=HK_examined, Positive=HK_positive,Prevalence=HK_prevalence)
HK$Parasite<-"Hookworm"
TT<-STH_2000 %>%
  mutate(Examined=TT_examined, Positive=TT_positive,Prevalence=TT_prevalence)
TT$Parasite<-"Trichuris"

STH_2000Long<-rbind(Asc,TT)
STH_2000Long<-rbind(STH_2000Long, HK)

STH_2000Long<-STH_2000Long[,c("Country","Year","AgeClass","Examined","Positive","Prevalence","Parasite")]
SCH_2000<-SCH_2000[,c("Country","Year","AgeClass","Examined","Positive","Prevalence","Parasite")]

New2000<-rbind(STH_2000Long,SCH_2000)
New2000 <- New2000 %>%
  rename(country_title = Country)

 table(epi_helminth$worm_type_title) 
epi_helminth <- epi_helminth %>%
  mutate(Year = year_end) %>%
  filter(age_start >= 0, age_end > 0, Year >= 2000) %>%
  # LOOK INTO CROSSING AGES
  mutate(AgeClass = ifelse(age_end <= 4, "Prev0_4",
                           ifelse(age_end >= 5 & age_end <= 14, "Prev5_14",
                                  ifelse(age_end >= 15, "Prev15plus", age_end))))%>%
  dplyr::select(Year, AgeClass,number_examined,
                number_positive, prevalence, worm_type_title,country_title)%>%
  rename(Parasite = worm_type_title,Positive=number_positive,Prevalence=prevalence,Examined=number_examined )
epi_helminth2<-rbind(epi_helminth,New2000)
table(epi_helminth$Parasite)
epi_helminth2 <- epi_helminth2 %>%
  left_join(pop, by = "country_title")%>% # Note that "country_title" is the name I use to name countries in all data frames
  replace_na(list(income="high"))%>%
  filter(!income == "high")%>%
  filter(!country_title == "Oman")

table(epi_helminth2$Parasite)
epi_helminth2 <-subset(epi_helminth2,epi_helminth2$Parasite!="not specified")
epi_helminth2[epi_helminth2 == "haematobium"] <- "S_haematobium"
epi_helminth2[epi_helminth2 == "mansoni"] <- "S_mansoni"
epi_helminth2[epi_helminth2 == "S. intercalatum"] <- "S_intercalatum"
epi_helminth2[epi_helminth2 == "S. japonicum"] <- "S_japonicum"
epi_helminth2[epi_helminth2 == "S. mekongi"] <- "S_mekongi"
epi_helminth2<-subset(epi_helminth2, is.na(epi_helminth2$Prevalence)=="FALSE")
epi_helminth2<-subset(epi_helminth2, is.na(epi_helminth2$Examined)=="FALSE")
table(is.na(epi_helminth2$Examined))
table(epi_helminth2$Parasite)
epi_helminth2$Examined<-as.numeric(epi_helminth2$Examined)
epi_helminth2$Prevalence<-as.numeric(epi_helminth2$Prevalence)
prev <- epi_helminth2 %>%
  group_by(country_title, Parasite, AgeClass) %>%
  summarise(mean_prevalence = mean(Prevalence), # obtain mean prevalence for the age classes created above
            sd_prevalence = sd(Prevalence, na.rm = TRUE),
            n=sum(Examined))%>%
  group_by(country_title, Parasite, AgeClass) %>%
  gather(Metric, Epi_value, mean_prevalence:sd_prevalence)%>%
  spread(AgeClass, Epi_value)%>% 
  replace_na(list(Prev0_4 = 1, Prev5_14 = 1, Prev15plus = 1)) %>% # for age classes for which there is no info, replace with 1 so that below we can use the multiplicator to estimate prevalence based on a reference value  
  gather(AgeClass,Epi_value, Prev0_4:Prev5_14)%>%
  left_join(prev_multiplicator, by = c("Parasite", "AgeClass"))# join the multiplicator dataframe 

#Combine to figure out which foods need prices
Country_List<-prev %>%
  group_by(country_title) %>%
  summarise(n=n())

table(food_balance5$Area)

food_balance5<-food_balance4%>%
  mutate(Area = ifelse(Area =="Slovakia", "Slovak Republic", 
                       ifelse(Area =="United Kingdom of Great Britain and Northern Ireland", "United Kingdom",
                              ifelse(Area =="Saint Kitts and Nevis", "St. Kitts and Nevis",
                                     ifelse(Area =="C\xf4te d'Ivoire", "Ivory Coast",
                                            ifelse(Area =="China, mainland", "China",
                                                   ifelse(Area =="Democratic People's Republic of Korea", "Republic of Korea",
                       ifelse(Area=="United States of America", "United States", Area))))))))
high_income2<-unique(high_income)
food_balance_cal_LMI <- food_balance5 %>%
  left_join(food_cals, by = c("Item")) %>%
  rename(country_title = Area) %>%
  left_join(high_income2, by = c("country_title")) %>%
  filter(is.na(income)==TRUE)
high_income2<-unique(high_income)
table(high_income$country_title)

food_balance_cal_LMI <-subset(food_balance_cal_LMI, food_balance_cal_LMI$country_title!="New Zealand")
food_balance_cal_LMI <-subset(food_balance_cal_LMI, food_balance_cal_LMI$country_title!="Oman")
food_balance_cal_LMI <-subset(food_balance_cal_LMI, food_balance_cal_LMI$country_title!="New Caledonia")
  
table(food_balance_cal_LMI$country_title)  # check to make sure are all LMI countries
  
# make a dataframe for each age class and rename variables 
prevalence_0_4 <- prev%>%
  filter(AgeClass == "Prev0_4")%>%
  rename(M0_4 = Multiplicator,
         Class0_4 = AgeClass,
         P0_4=Epi_value)

prevalence_5_14 <- prev %>%
  filter(AgeClass == "Prev5_14")%>%
  rename(M5_14 = Multiplicator,
         Class5_14 = AgeClass,
         P5_14 = Epi_value)

prevalence_15p <- prev %>%
  filter(AgeClass == "Prev15plus")%>%
  rename(M15p = Multiplicator,
         Class15plus = AgeClass,
         P15p= Epi_value)

# for age specific prevalence,
# join the dataframes
# estimate prevalence for age classes with missing info
# use first the 5-14 age class as reference if available 
prev_age_country <- prevalence_0_4%>%
  left_join(prevalence_5_14, by = c("country_title", "Parasite", "Metric", "n"))%>%
  left_join(prevalence_15p, by = c("country_title", "Parasite", "Metric", "n"))%>%
  mutate(P5_14 = ifelse(P5_14 == 1, (M5_14+(M5_14-M15p))*P15p, P5_14),# if no data for 5-14, use prevalence of 15plus times multiplicator of 5-14 plus multiplicator of 15p. same for lines below
         P5_14 = ifelse(P5_14 == 1, (M5_14+(M5_14-M0_4))*P0_4, P5_14),
         P0_4 = ifelse(P0_4 == 1, M0_4*P5_14, P0_4),
         P0_4 = ifelse(P0_4 == 1, M0_4*P15p, P0_4),
         P15p = ifelse(P15p == 1, M15p*P5_14, P15p),
         P15p  = ifelse(P15p  == 1, M15p*P0_4, P15p))%>%
  gather(AgeClass,Epi_value, c(P0_4,P5_14,P15p))%>%
  dplyr::select(country_title, Parasite, AgeClass, Metric, Epi_value, n)%>%
  spread(Metric, Epi_value)%>%
  mutate(AgeClass=gsub("P0_4","0_4", AgeClass),
         AgeClass=gsub("P5_14","5_14", AgeClass),
         AgeClass=gsub("P15p","15plus", AgeClass),
         error_prevalence = qnorm(0.975)*sd_prevalence/sqrt(n()))%>%
  group_by(country_title, Parasite, AgeClass)%>%
  summarise(mean_prevalence = mean(mean_prevalence),
            sd_prevalence = mean(sd_prevalence),
            error_prevalence = mean(error_prevalence))%>%
  left_join(mean_pop, by =c("country_title","AgeClass"))%>%
  mutate(mean_prevalence = ifelse(mean_prevalence >100, 99, mean_prevalence), # for observations that ended up with prevalence greater than 100, change to 99. Thiis is an artifact of the multiplicator method
         TotInfected = Pop * (mean_prevalence/100),
         TotInfected_error = Pop * (error_prevalence/100))%>%
  na.omit

# country prevalence, not taking into account age prevalence
prevalence_country <- epi_helminth2 %>%
  group_by(country_title, Parasite) %>%
  summarise(mean_prevalence = mean(Prevalence),
            sd_prevalence = sd(Prevalence, na.rm=TRUE),
            error_prevalence = qnorm(0.975)*sd_prevalence/sqrt(n()),
            low95_prev = mean_prevalence - error_prevalence,
            high95_prev = mean_prevalence + error_prevalence)%>%
  left_join(high_income, by="country_title")%>%
  replace_na(list(income="low-med"))%>%
  filter(!income == "high")%>%
  filter(!country_title == "Oman")


#### MDA costs  ####

# estimate costs of MDA 
# because drugs target all parasites, prevalence must be combined to determine whether prevalence meets the threshold for MDA  
mda_age_country <- prev_age_country %>%
  dplyr::select(country_title, Region, Parasite, AgeClass, TotInfected, Pop, mean_prevalence, error_prevalence, TotInfected_error)%>%
  group_by(country_title, Parasite, Region)%>%
  summarise(prev = sum(mean_prevalence), 
            Infected = sum(TotInfected),
            Population = sum(Pop))%>%
  mutate(Helminth = ifelse(Parasite %in% c("Ascaris", "Trichuris", "Hookworm"), "STH",
                           ifelse(Parasite %in% c("S_haematobium", "S_mansoni","S_intercalatum", "S_japonicum", "S_mekongi"), "WTH", Parasite)))%>%
  group_by(country_title, Helminth, Region) %>%
  summarise(prevalence_comb = mean(prev), # note, the summed prevalence results in prevalence values much greater than 100 because combining prevalence of multiple diseases. I'm not absolutely certain this is the best approach for this
            Infected = sum(Infected),
            Population = sum(Population),
            true_prev = (Infected/Population)*100)%>%
  rename(Parasite = Helminth) 

MDA <- read_excel("Parameters2.xlsx", sheet = "prev_thresholds_MDA") # WHO thresholds fro MDA and Lo et al 2016 thresholds

# Lo et al 2016 and WHO comparison of thresholds
WHO_LO <- MDA %>%
  dplyr::select(Parasite, prev_thresh_community_Lo, prev_thresh_community_WHO)%>%
  rename(Prevalence_Threshold_Lo_etal_2016 = prev_thresh_community_Lo,
         Prevalence_Threshold_WHO = prev_thresh_community_WHO)

# select only WHO's
MDA_WHO_Comm <- MDA %>%
  dplyr::select(Parasite, prev_thresh_community_WHO, coverage, cost_community, cost_community_low, cost_community_high)

# find countries that have STH and WTH and select those that fall at or above the threshold for co-endemic countries
Tx_coendemic_WHO <- mda_age_country %>%
  left_join(MDA_WHO_Comm, by = "Parasite")%>%
  filter(prevalence_comb  >= prev_thresh_community_WHO)%>%
  group_by(country_title)%>%
  filter(duplicated(country_title))%>%
  mutate(Parasite = "Co-endemic")

# find countries that have STH and select those that fall at or above the threshold prevalence for STH
Tx_STH_WHO <- mda_age_country %>%
  filter(Parasite == "STH" )%>%
  left_join(MDA_WHO_Comm, by = "Parasite")%>%
  filter(prevalence_comb  >= prev_thresh_community_WHO)

# find countries that have WTH and select those that fall at or above the threshold prevalence for WTH
Tx_WTH_WHO <- mda_age_country %>%
  filter(Parasite == "WTH")%>%
  left_join(MDA_WHO_Comm, by = "Parasite")%>%
  filter(prevalence_comb  >= prev_thresh_community_WHO)

#same process as above but with Lo et al 2016 proposed guideliines
MDA_Lo_Comm <- MDA %>%
  dplyr::select(Parasite, prev_thresh_community_Lo, coverage, cost_community, cost_community_low, cost_community_high)

Tx_coendemic_Lo <- mda_age_country %>%
  filter(prevalence_comb >= 4)%>%
  group_by(country_title)%>%
  filter(duplicated(country_title))%>%
  mutate(Parasite = "Co-endemic")

Tx_STH_Lo <-mda_age_country %>%
  filter(Parasite == "STH")%>%
  left_join(MDA_Lo_Comm, by = "Parasite")%>%
  filter(prevalence_comb  >= prev_thresh_community_Lo)

Tx_WTH_Lo <- mda_age_country %>%#prevalence_combined %>%
  filter(Parasite == "WTH")%>%
  left_join(MDA_Lo_Comm, by = "Parasite")%>%
  filter(prevalence_comb  >= prev_thresh_community_Lo)

# combine by region and estimate cost by multipliying number of people by cost of community MDA times coverage (75% or 100%)
# Regional Costs for Lo et al. thresholds
Country_MDA_Lo <- rbind(Tx_coendemic_Lo, Tx_STH_Lo, Tx_WTH_Lo)%>%
  group_by(country_title)%>%
  filter(!duplicated(country_title))%>%
  dplyr::select(country_title, Parasite, prevalence_comb)%>%
  left_join(MDA_Lo_Comm, by = "Parasite")%>%
  left_join(pop, by="country_title")%>%
  mutate(Average_Cost75= cost_community * (MeanPop * coverage),
         Low_Cost75 = cost_community_low * (MeanPop * coverage),
         High_Cost75 = cost_community_high * (MeanPop * coverage),
         Average_Cost100= cost_community * MeanPop ,
         Low_Cost100 = cost_community_low * MeanPop,
         High_Cost100 = cost_community_high * MeanPop)

Region_MDA_Lo <- Country_MDA_Lo%>%
  group_by(Region)%>%
  summarise(Average_Cost75 = sum(Average_Cost75),
            Low_Cost75 = sum(Low_Cost75),
            High_Cost75 = sum(High_Cost75),
            Average_Cost100 = sum(Average_Cost100),
            Low_Cost100 = sum(Low_Cost100),
            High_Cost100 = sum(High_Cost100))%>%
  mutate(Source = "Lo_2016")

# Same as above but for WHO thresholds
Country_MDA_WHO <- rbind(Tx_coendemic_WHO, Tx_STH_WHO, Tx_WTH_WHO)%>%
  filter(!duplicated(country_title))%>%
  dplyr::select(country_title, Parasite, prevalence_comb)%>%
  left_join(MDA_WHO_Comm, by = "Parasite")%>%
  left_join(pop, by="country_title")%>%
  mutate(Average_Cost75= cost_community * (MeanPop * coverage),
         Low_Cost75 = cost_community_low * (MeanPop * coverage),
         High_Cost75 = cost_community_high * (MeanPop * coverage),
         Average_Cost100= cost_community * MeanPop ,
         Low_Cost100 = cost_community_low * MeanPop,
         High_Cost100 = cost_community_high * MeanPop)

Region_MDA_WHO <- Country_MDA_WHO%>%
  group_by(Region)%>%
  summarise(Average_Cost75 = sum(Average_Cost75),
            Low_Cost75 = sum(Low_Cost75),
            High_Cost75 = sum(High_Cost75),
            Average_Cost100 = sum(Average_Cost100),
            Low_Cost100 = sum(Low_Cost100),
            High_Cost100 = sum(High_Cost100))%>%
  mutate(Source = "WHO")

# Global costs WHO
MDA_WHO_integrated <- Region_MDA_WHO%>%
  mutate(Source = "WHO")%>%
  group_by(Source)%>%
  summarise(Average_Cost75 = sum(Average_Cost75),
            Low_Cost75 = sum(Low_Cost75),
            High_Cost75 = sum(High_Cost75),
            Average_Cost100 = sum(Average_Cost100),
            Low_Cost100 = sum(Low_Cost100),
            High_Cost100 = sum(High_Cost100))

# Global costs Lo et al 2016
MDA_Lo_integrated <- Region_MDA_Lo%>%
  mutate(Source = "Lo_2016")%>%
  group_by(Source)%>%
  summarise(Average_Cost75 = sum(Average_Cost75),
            Low_Cost75 = sum(Low_Cost75),
            High_Cost75 = sum(High_Cost75),
            Average_Cost100 = sum(Average_Cost100),
            Low_Cost100 = sum(Low_Cost100),
            High_Cost100 = sum(High_Cost100))

# Global costs WHO and Lo et al
Global_MDA_Costs <- rbind(MDA_Lo_integrated,MDA_WHO_integrated)%>%
  mutate(Region = "Global")%>% # add column Region = Global to integrate with regional costs (see below)
  print()

# combine regional costs with global costs
Region_MDA_Costs <- rbind(Region_MDA_WHO,Region_MDA_Lo, Global_MDA_Costs)%>%
  print()

#### Food Savings   ####

# to estimate food savings need to estimate number of parasites ppl are infected with. Might want to talk with Jason, Alex and Sean to see how best to approach this.

# worm morbidity thresholds from Chan et al 1993
# Schistosoma data from WHO's definition of heavy infection (50 epg) for Schistosoma
a1 <- seq(0,9) # ascaris 0-4 year olds
a2 <- seq(0,15) # ascaris 5-14 year olds
a3 <- seq(0,19) # ascaris 15 plus
t1 <- seq(0,89) # trichuris 0-4 year olds
t2 <- seq(0,129) #trichuris 5-14 year olds
t3 <- seq(0,170) # trichuris 15 plus
h1 <- seq(0,19) # hookworm 0-4 year olds
h2 <- seq(0,29) # hookworm 5-14 year olds
h3 <- seq(0,39) # hookworm 15 plus
s1 <- seq(0,170) # Schistosoma 0-4 yo. using the highest worm burden that does not result in NaN. Morbidity threshold burdens could be in the 200-400 based on WHO's definition of heavy infection (50 epg) and Guriiare et al worm to egg ratios  
s2 <- seq(0,170) # Schistosoma 5-14 yo. using the highest worm burden that does not result in NaN. Morbidity threshold burdens could be in the 200-400 based on WHO's definition of heavy infection (50 epg) and Guriiare et al worm to egg ratios
s3 <- seq(0,170) # Schistosoma 15 plus. using the highest worm burden that does not result in NaN. Morbidity threshold burdens could be in the 200-400 based on WHO's definition of heavy infection (50 epg) and Guriiare et al worm to egg ratios

# these lines are for estimating mean burden based on prevalence and dispersion parameter. Formula is obtained from Chan et al 1993
burden <- prev_age_country%>%
  mutate(prevalence = mean_prevalence/100,
         k = ifelse(Parasite == "Ascaris", 0.54,
                    ifelse(Parasite == "Trichuris", 0.23,
                           ifelse(Parasite == "Hookworm", 0.34,
                                  ifelse(Parasite %in% c("S_haematobium", "S_mansoni","S_intercalatum", "S_japonicum", "S_mekongi"), 0.23, Parasite)))),
         k = as.numeric(k),
         mean_burden = k*(((1-prevalence)^(-1/k)) -1)) #NEED TO CHECK

# dispersion parameters for parasites. Obtained from Chan et al 1993
a_k=0.54
t_k=0.23
h_k=0.34
s_k=0.23

# this function is for estimating the proportion of individuals that are infected with number of parasites below the morbidity threshold. Based on Chan et al 1993 
p_calc<-function(x,k,bcnts){
  sum(((1+x/k)^(-k))*
        ((gamma(k+bcnts))/(factorial(bcnts)*gamma(k)))*
        (x/(x+k))^bcnts)  
}

# use function and worm thresholds by parasite and age class to estimate number of individuals below and above threshold
# ascaris 0-4 yo
b_a1 <- burden %>%
  filter(Parasite == "Ascaris" & AgeClass == "0_4")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=a_k,bcnts=a1),
         aboveT = 1-belowT, 
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 10) 

# ascaris 5-14 yo
b_a2 <- burden %>%
  filter(Parasite == "Ascaris" & AgeClass == "5_14")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=a_k,bcnts=a2),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 20) 

# ascaris 15 plus 
b_a3 <- burden %>%
  filter(Parasite == "Ascaris" & AgeClass == "15plus")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=a_k,bcnts=a3),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 20)

# combine data for ascaris
burden_ascaris <- rbind(b_a1,b_a2,b_a3)%>%
  dplyr::select(country_title, Parasite, AgeClass, belowT, aboveT, TotInfected, Infected_aboveT, Infected_belowT,burden)

# trichuris 0-4
b_t1 <- burden %>%
  filter(Parasite == "Trichuris" & AgeClass == "0_4")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=t_k,bcnts=t1),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 90)

#trichuris 5-14
b_t2 <- burden %>%
  filter(Parasite == "Trichuris" & AgeClass == "5_14")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=t_k,bcnts=t2),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 130)

#trichuris 15 plus
b_t3 <- burden %>%
  filter(Parasite == "Trichuris" & AgeClass == "15plus")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=t_k,bcnts=t3),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 171)

# combine burden for trichuris
burden_trichuris <- rbind(b_t1,b_t2,b_t3)%>%
  dplyr::select(country_title, Parasite, AgeClass, belowT, aboveT, TotInfected, Infected_aboveT, Infected_belowT,burden)

# hookworm 0-4
b_h1 <- burden %>%
  filter(Parasite == "Hookworm" & AgeClass == "0_4")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=h_k,bcnts=h1),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 20)

# hookworm 5-14
b_h2 <- burden %>%
  filter(Parasite == "Hookworm" & AgeClass == "5_14")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=h_k,bcnts=h2),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 30)

# hookworm 15 plus
b_h3 <- burden %>%
  filter(Parasite == "Hookworm" & AgeClass == "15plus")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=h_k,bcnts=h3),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 40)

# combine hookworm data
burden_hook <- rbind(b_h1,b_h2,b_h3)%>%
  dplyr::select(country_title, Parasite, AgeClass, belowT, aboveT,TotInfected, Infected_aboveT, Infected_belowT,burden)

# schisto 0-4
b_s1 <- burden %>%
  filter(Parasite %in% c("S_mansoni","S_haematobium","S_intercalatum", "S_japonicum", "S_mekongi") & AgeClass == "0_4")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=s_k,bcnts=s1),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 171)
# schisto 5-14
b_s2 <- burden %>%
  filter(Parasite %in% c("S_mansoni","S_haematobium","S_intercalatum", "S_japonicum", "S_mekongi") & AgeClass == "5_14")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=s_k,bcnts=s2),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 171)
# schisto 15 plus
b_s3 <- burden %>%
  filter(Parasite %in% c("S_mansoni","S_haematobium","S_intercalatum", "S_japonicum", "S_mekongi") & AgeClass == "15plus")%>%
  mutate(belowT = sapply(mean_burden,p_calc,k=s_k,bcnts=s3),
         aboveT = 1-belowT,
         Infected_aboveT = aboveT * TotInfected,
         Infected_belowT = belowT * TotInfected,
         burden = 171)

# combine schisto data
burden_schisto <- rbind(b_s1,b_s2,b_s3)%>%
  dplyr::select(country_title, Parasite, AgeClass, belowT, aboveT, TotInfected, Infected_aboveT, Infected_belowT, burden)

# combine all parasite burden data - Look at distribution instead of just cut off
burden_morb <- burden_ascaris %>%
  rbind(burden_trichuris)%>%
  rbind(burden_hook)%>%
  rbind(burden_schisto)%>%
  left_join(params_morb, by = "Parasite")%>%
  left_join(burden_by_parasite, by="Parasite")%>% 
  mutate(kcals_consumed_year_aboveT = burden * kcal_year_wastedHost, # number of kcals for individuals infected with burden above morbidity threshold
         kcals_consumed_year_belowT = mean_burden * kcal_year_wastedHost) # number of kcals for individuals infected with  mean burden --> below morbidity threshold
table(epi_helminth2$country_title)
# food savings by country
# use the prevalence by age data and combine with burden_morb database to estimate the total number of calories consumed by parasites
# for people infected with burden above and below threshold

fs_country <- Country_MDA_Lo %>%
  dplyr::select(country_title, Region)%>%
  left_join(burden_morb,by="country_title")

#Read in csv with cost per cal
Cost_cal<-read.csv("FoodPriceLMIFull2.csv") 
Cost_cal$prop<-(Cost_cal$food_comp*Cost_cal$PriceCal)

Cost_cal <- Cost_cal%>%
  group_by(country_title)%>%
  summarise(TotalCalProp = sum(prop),TotalFoodcomp = sum(food_comp))%>%
  mutate(CostCal=TotalCalProp/TotalFoodcomp)
range(Cost_cal$CostCal)

fs_country <- fs_country %>%
  left_join(Cost_cal,by="country_title") %>%
  mutate(Tot_kcals_consumed_year = (((Infected_aboveT * kcals_consumed_year_aboveT) + (Infected_belowT * kcals_consumed_year_belowT)) * Drug_Efficacy),
         Food_Savings =  (Tot_kcals_consumed_year*CostCal))
table(is.na(fs_country$Food_Savings),fs_country$Parasite)
range(fs_country$Food_Savings)

# food savings by region and age class
fs_region <- fs_country %>%
  group_by(Region, AgeClass)%>%
  summarise(Food_Savings = sum(Food_Savings), Tot_kcals_consumed_year=sum(Tot_kcals_consumed_year))

fs_species <- fs_country %>%
  group_by(Region, Parasite)%>%
  summarise(Food_Savings = sum(Food_Savings), Tot_kcals_consumed_year=sum(Tot_kcals_consumed_year))

# global food savings by age class
fs_global <- fs_region %>%
  group_by(AgeClass)%>%
  summarise(Food_Savings = sum(Food_Savings), Tot_kcals_consumed_year=sum(Tot_kcals_consumed_year))%>%
  mutate(Region = "Global")

fs_global2 <- fs_global %>%
  group_by(Region)%>%
  summarise(Food_Savings = sum(Food_Savings), Tot_kcals_consumed_year=sum(Tot_kcals_consumed_year))

fs_global2$PeopleFed<-fs_global2$Tot_kcals_consumed_year/(2300*365)


Food_Savings <- rbind(fs_global,fs_region)%>%
  group_by(Region)%>%
  summarise(Cost_Food_Averted = sum(Food_Savings),
            Cost_Food_Averted_75 = sum(Food_Savings* 0.75))%>%
  print()



#### DALYs ####

DALY <- read.csv("DALYS_Helminth_IHME.csv") %>% # use DALYs from the Global Burden of Disease database - years 2015-2020
  filter(metric == "Number")%>%
  group_by(country_title, Parasite, age)%>%
  summarise(DALYs = mean(DALYs)) # take mean dalys across years 

GNIpc <- read.csv("GNIpc.csv") %>% # Mean Gross National Income per capita of 2000-2020
  dplyr::select(country_title, GNIpc_2000_2020)

# estimate DALYs for countries that meet the MDA thresholds based on Lo et al 2016 guidelines (most conservative)
DALYs_USD <- Country_MDA_Lo %>%
  dplyr::select(country_title, Region)%>%
  left_join(DALY, by = "country_title")%>%
  group_by(country_title, Parasite, Region)%>%
  summarise(DALYs = sum(DALYs))%>%
  left_join(efficacy, by = "Parasite")%>%
  left_join(GNIpc, by = "country_title")%>%
  mutate(DALY_USD = (GNIpc_2000_2020 * DALYs) * Drug_Efficacy)

# estimate total DALYs by combining parasite specific dalys
DALYs_Combined <- DALYs_USD %>%
  group_by(country_title, Region)%>%
  summarise(DALY_USD = sum(DALY_USD))

# estimate dalys by region
DALYs_averted_Region <- DALYs_Combined  %>%
  group_by(Region)%>%
  summarise(DALY_USD = sum(DALY_USD, na.rm=TRUE),
            DALY_USD75 = DALY_USD *0.75)%>%
  print()

# estimate dalys at global scale
DALYs_averted_Global <- DALYs_averted_Region %>%
  ungroup()%>%
  summarise(DALY_USD = sum(DALY_USD),
            DALY_USD75 = sum(DALY_USD75))%>%
  mutate(Region= "Global")%>%
  print()

DALYs_AVERTED_USD <- rbind(DALYs_averted_Region, DALYs_averted_Global)


#### Net Benefits (DALYs - MDA) ####

# estimate Net benefits considering ONLY DALYs and MDA

NetBenefits_country <- Country_MDA_Lo %>%
  dplyr::select(country_title, Average_Cost100:High_Cost100)%>%
  left_join(DALYs_Combined, by = "country_title")%>%
  mutate(NB_Low = DALY_USD - Low_Cost100,
         NB_Avg = DALY_USD - Average_Cost100,
         NB_High = DALY_USD - High_Cost100)

NetBenefits_global <- NetBenefits_country %>%
  ungroup()%>%
  summarise(NB_Low = sum(NB_Low),
            NB_Avg = sum(NB_Avg),
            NB_High = sum(NB_High))%>%
  print()

# DALYs benefits at global scale outweigh MDA costs when MDA cost per treatment is $0.5 or $0.75



#### Plots MDA vs (Food Savings + DALYs) ####

#  benefit to cost ratios by region comparing MDA guidelines proposed by WHO and LOo et al 2016

BC_Ratio <- Region_MDA_Costs %>%
  left_join(Food_Savings, by = "Region")%>%
  left_join(DALYs_AVERTED_USD, by = "Region")%>%
  #left_join(ES_ha_saved_region, by = "Region")%>%
  filter(!Region == "Micronesia")%>%
  group_by(Region, Source)%>%
  summarise(Cost_Food_Averted = sum(Cost_Food_Averted),
            DALY_USD = mean(DALY_USD),
            #ES = mean(ES),
            Low_Cost100= mean(Low_Cost100),
            Average_Cost100=mean(Average_Cost100),
            High_Cost100 = mean(High_Cost100))%>%
  mutate(Ratio_Av = (Cost_Food_Averted + DALY_USD)/ Average_Cost100,
         Ratio_Low = (Cost_Food_Averted+DALY_USD ) / Low_Cost100,
         Ratio_High =(Cost_Food_Averted+DALY_USD) /High_Cost100)%>%
  gather("MDA_Type", "CB_Ratio", Ratio_Av:Ratio_High)%>%
  mutate(MDA_Type=gsub("Ratio_Av","Average Cost (USD$ 0.75) @ 100%", MDA_Type),
         MDA_Type=gsub("Ratio_Low","Low Cost (USD$ 0.50) @ 100%", MDA_Type),
         MDA_Type=gsub("Ratio_High","High Cost (USD$ 3) @ 100%", MDA_Type))%>%
  replace_na(list(Average_Cost=0, Low_Cost=0, High_Cost=0))%>%
  na.omit

BC_Ratio$Region<- factor(BC_Ratio$Region,levels = c("Melanesia", "Caribbean", "Central America","South America",
                                                                          "Northern Africa","Middle Africa","Western Africa","Eastern Africa",
                                                                          "Southern Africa","Southern Asia", "South-Eastern Asia","Western Asia","Eastern Asia","Global"))
BC_Ratio$MDA_Type<- factor(BC_Ratio$MDA_Type,levels = c("Low Cost (USD$ 0.50) @ 100%","Average Cost (USD$ 0.75) @ 100%", "High Cost (USD$ 3) @ 100%"))


ggplot(BC_Ratio, aes(y= CB_Ratio, x = Region, shape = Source))+
  geom_segment( aes(x=Region, xend=Region, y=0, yend=CB_Ratio), color="grey") +
  geom_point( size=2) +
  geom_hline(yintercept = 1,colour="black", linetype = "longdash")+
  scale_shape_manual(values=c(1,4), name = "Guidelines")+
  facet_wrap(~MDA_Type, scales="free")+
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=90, size =12),
    axis.text.y = element_text(size = 12)
  ) +
  xlab("Region") +
  ylab("Expected Annual Benefit-Cost Ratio")

ggsave(
  "Figure1-2000.tiff",
  width = 12,
  height = 9,
  dpi = 300)

#### Code for Maps #####
mda_country <- Country_MDA_Lo%>%
  dplyr::select(country_title, Region, Average_Cost100, Low_Cost100, High_Cost100, MeanPop)%>%
  rename(MDA_Low = Low_Cost100,
         MDA_Avg = Average_Cost100,
         MDA_High = High_Cost100)

dalys_country <- DALYs_USD %>%
  group_by(country_title)%>%
  summarise(dalys_avert = sum(DALY_USD))

dalys_USDworm <- DALYs_USD %>%
  group_by(Parasite)%>%
  summarise(dalys_avert = sum(DALY_USD))

dalys_worm <- DALY %>%
  group_by(Parasite)%>%
  summarise(DALYs = sum(DALYs))

fs_country_sum <- fs_country %>%
  group_by(country_title)%>%
  summarise(FS = sum(Food_Savings))

fs_mda_map <- fs_country_sum %>%
  left_join(mda_country, by="country_title")%>%
  left_join(dalys_country, by="country_title")%>%
  mutate(NetBenefits_low = (FS+dalys_avert) - MDA_Low,
         NetBenefits_avg = (FS+dalys_avert) - MDA_Avg,
         NetBenefits_high = (FS+dalys_avert) - MDA_High)
  
fs_mda_map_region <- fs_mda_map %>%
  group_by(Region)%>%
  summarise(FS = sum(FS, na.rm = TRUE),
            dalys_avert = sum(dalys_avert, na.rm = TRUE),
            MDA_Low = sum(MDA_Low, na.rm = TRUE),
            MDA_Avg = sum(MDA_Avg, na.rm = TRUE),
            MDA_High = sum(MDA_High, na.rm = TRUE))%>%
  mutate(NetBenefits_low = (FS+dalys_avert) - MDA_Low,
         NetBenefits_avg = (FS+dalys_avert) - MDA_Avg,
         NetBenefits_high = (FS+dalys_avert) - MDA_High)

fs_mda_map_region_Globaal<-colSums(fs_mda_map_region[,-1])
fs_mda_map_region_Globaal
fs_mda_map_region <- fs_mda_map %>%
  group_by(Region)%>%
  summarise(FS = sum(FS, na.rm = TRUE),
            dalys_avert = sum(dalys_avert, na.rm = TRUE),
            MDA_Low = sum(MDA_Low, na.rm = TRUE),
            MDA_Avg = sum(MDA_Avg, na.rm = TRUE),
            MDA_High = sum(MDA_High, na.rm = TRUE))%>%
  mutate(NetBenefits_low = (FS+dalys_avert) - MDA_Low,
         NetBenefits_avg = (FS+dalys_avert) - MDA_Avg,
         NetBenefits_high = (FS+dalys_avert) - MDA_High)



country_map <- country_borders %>%
  st_as_sf()%>%
  left_join(fs_mda_map, by = "country_title")%>%
  mutate(NetBenefits_low = NetBenefits_low,
         NetBenefits_avg = NetBenefits_avg,
         NetBenefits_high = NetBenefits_high,
         FS = FS,
         dalys_avert = dalys_avert,
         MDA_Low = MDA_Low,
         MDA_Avg = MDA_Avg,
         MDA_High = MDA_High,
         MeanPop=MeanPop,
         NetBenefits_lowPC = NetBenefits_low/MeanPop,
         NetBenefits_avgPC = NetBenefits_avg/MeanPop,
         NetBenefits_highPC = NetBenefits_high/MeanPop)

region_map <- region_borders %>%
  st_as_sf()%>%
  left_join(fs_mda_map_region, by = "Region")%>%
  mutate(NetBenefits_low = NetBenefits_low/1000000,
         NetBenefits_avg = NetBenefits_avg/1000000,
         NetBenefits_high = NetBenefits_high/1000000,
         FS = FS/1000000,
         dalys_avert = dalys_avert/1000000,
         MDA_Low = MDA_Low/1000000,
         MDA_Avg = MDA_Avg/1000000,
         MDA_High = MDA_High/1000000)


below_thresh <- prev_age_country %>%
  ungroup()%>%
  dplyr::select(country_title)%>%
  filter(!duplicated(country_title))%>%
  anti_join(fs_country, by="country_title")%>%
  mutate(Threshold = "Below threshold for MDA")

country_borders_thresh <- country_borders %>%
  st_as_sf()%>%
  right_join(below_thresh, by = "country_title")%>%
  filter(country_title!="Palestine")

  

FS_map <- tm_shape(country_map) +
  tm_fill(col = "FS", style = "quantile",n=4,
          palette = "Blues", title ="Food Savings", 
          textNA = "Not endemic",  
          colorNA = "gray", midpoint = NA) +
  tm_borders() +
  tm_shape(country_borders_thresh)+
  tm_symbols(size = 0.1, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE)  #+tm_compass() 

FS_map

dalys_map <- tm_shape(country_map)+
  tm_borders()+
  tm_fill(col="dalys_avert", style="quantile",n=4,palette = "Purples", title = "DALYs Averted", textNA = "Not endemic",
          colorNA = "gray")+
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.1, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+ sf::sf_use_s2(FALSE)

dalys_map



savings_map <- tmap_arrange(dalys_map,FS_map, ncol=1)

savings_map
tmap_save(savings_map, "Maps_SavingCat.png", height = 7)

#ES_map <- tm_shape(country_map)+
  tm_borders()+
  tm_fill(col="ES", style="quantile",n=4, 
          palette = "Greens", title = "ES Benefits", textNA = "Not endemic",
          colorNA = "gray")+
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.1, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+ sf::sf_use_s2(FALSE)

#ES_map - don't use yet

# Low cost MDA
# country
MDA_low_map <- tm_shape(country_map) +
  tm_borders() +
  tm_fill(col ="MDA_Low",style="quantile",n=4,
          palette = "Reds", title = "Cost of MDA ($0.5)", textNA = "Not endemic", 
          colorNA = "gray") +
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.1, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE) # +tm_compass() 
MDA_low_map


NB_low_map <- tm_shape(country_map) +
  tm_borders() +
  tm_fill(col ="NetBenefits_low", style = "quantile",n=4,
          palette = "GnBu", title = "Low Cost Net Benefits", textNA = "Not endemic", 
          colorNA = "gray", midpoint =NA) +
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.2, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE) # +tm_compass() 
NB_low_map

#tmap_save(NB_low_dalys_map, "Map_NB_low_dalys.png", height = 9, width = 9)

country_maps_ben <- tmap_arrange(FS_map, dalys_map, ES_map, ncol = 1)

#tmap_save(country_maps_ben, "Maps_Ben.png", height = 9, width = 9)

country_maps_Low <- tmap_arrange(MDA_low_map, NB_low_map, ncol=1)

#tmap_save(country_maps_Low, "Maps_DALYS_Low.png", height = 9, width = 9)

# Average cost MDA
# country
MDA_avg_map <- tm_shape(country_map) +
  tm_borders() +
  tm_fill(col ="MDA_Avg",palette = "Reds", style = "quantile", n=4, 
          title = "Cost of MDA ($1.5)", textNA = "Not endemic", 
          colorNA = "gray", midpoint=NA) +
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.2, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE) 
MDA_avg_map

NB_avg_map <- tm_shape(country_map) +
  tm_borders() +
  tm_fill(col ="NetBenefits_avg",palette = "RdYlGn", style = "quantile",n=4,
          title = "Avg Cost Net Benefits", textNA = "Not endemic", 
          colorNA = "gray", midpoint=NA) +
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.2, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE) # +tm_compass() 
NB_avg_map

country_maps_avg <- tmap_arrange(MDA_avg_map, NB_avg_map, ncol=1)

#tmap_save(country_maps_avg, "Maps_Avg.png", height = 7)

# High cost ($3)
MDA_high_map <- tm_shape(country_map) +
  tm_borders() +
  tm_fill(col ="MDA_High",palette = "Reds", style = "quantile", n=4, 
          title = "Cost of MDA ($3.0)", textNA = "Not endemic", 
          colorNA = "gray", midpoint=NA) +
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.2, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE) 
MDA_high_map

NB_high_map <- tm_shape(country_map) +
  tm_borders() +
  tm_fill(col ="NetBenefits_high",palette = "RdYlGn", style = "quantile",n=4,
          title = "High Cost Net Benefits", textNA = "Not endemic", 
          colorNA = "gray", midpoint=NA) +
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.2, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE) # +tm_compass() 
NB_high_map

country_maps_high <- tmap_arrange(MDA_high_map, NB_high_map, ncol=1)

#tmap_save(country_maps_high, "Maps_Higg.png", height = 7)


country_maps_net <- tmap_arrange(NB_low_map,NB_avg_map, NB_high_map, ncol=1)

country_maps_net
tmap_save(country_maps_net, "Maps_Net.png", height = 7)


global_Sums <- fs_mda_map %>%
  group_by()%>%
  summarise(Food_Savings = sum(FS),dalys=sum(dalys_avert),Low=sum(NetBenefits_low),Avg=sum(NetBenefits_avg),High=sum(NetBenefits_high) )



# Per capita

NB_low_mapPC <- tm_shape(country_map) +
  tm_borders() +
  tm_fill(col ="NetBenefits_lowPC", style = "quantile",n=4,
          palette = "GnBu", title = "Low Cost Net Benefits Per Capita", textNA = "Not endemic", 
          colorNA = "gray", midpoint =NA) +
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.2, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE) # +tm_compass() 
NB_low_mapPC

NB_avg_mapPC <- tm_shape(country_map) +
  tm_borders() +
  tm_fill(col ="NetBenefits_avgPC",palette = "RdYlGn", style = "quantile",n=4,
          title = "Avg Cost Net Benefits Per Capita", textNA = "Not endemic", 
          colorNA = "gray", midpoint=NA) +
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.2, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE) # +tm_compass() 
NB_avg_mapPC

NB_high_mapPC <- tm_shape(country_map) +
  tm_borders() +
  tm_fill(col ="NetBenefits_highPC",palette = "RdYlGn", style = "quantile",n=4,
          title = "High Cost Net Benefits Per Capita", textNA = "Not endemic", 
          colorNA = "gray", midpoint=NA) +
  tm_shape(country_borders_thresh)+
  tm_dots(size = 0.2, col = "Threshold",shape =18,
          pal = "black",
          popup.vars = TRUE)+
  tm_layout(legend.outside = TRUE)+sf::sf_use_s2(FALSE) # +tm_compass() 
NB_high_mapPC

country_maps_netPC <- tmap_arrange(NB_low_mapPC,NB_avg_mapPC, NB_high_mapPC, ncol=1)
country_maps_netPC
