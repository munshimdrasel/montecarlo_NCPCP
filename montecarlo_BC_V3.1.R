#Dynamic montecarlo ammonia limit calculation for Norman Cole Pollution Control Plant, Lorton, VA

rm(list=ls(all=TRUE))

#Loading libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(zoo)
library(lubridate)
library(reshape2)
library(data.table)

#setting up working directory

setwd("C:/Users/mrasel/OneDrive - Fairfax County Government/R-analysis/tss-spin/tss-spin")


#taking last permit cycle data
# Effective Date: November1, 2019
# Expiration Date: October 31, 2024

# ==============================summer analysis===========================================
df <- read_csv(("./data/ph_data/ph_data_NMCPCP.csv")) %>% 
  mutate(date=as.Date(date),   year=year(date), month= month(date),
         day= day(date)) %>%  filter( date >="2019-11-01" & date <= "2023-11-30")

#We want to use all pH sample data 
#combining max and min pH, and temperature of plants data
df_max <- df %>% 
  dplyr:: select(date, max_pH_outfall, Temperature_Out, final_effluent_flow, NH3_eff) %>% 
  na.omit()

names(df_max)[names(df_max) == "max_pH_outfall"]     <- "pH"
names(df_max)[names(df_max) == "Temperature_Out"]     <- "temperature"

df_min <- df %>% 
  dplyr:: select(date, min_pH_outfall, Temperature_Out, final_effluent_flow, NH3_eff) %>% 
  na.omit()


names(df_min)[names(df_min) == "min_pH_outfall"]     <- "pH"
names(df_min)[names(df_min) == "Temperature_Out"]     <- "temperature"

#combining max and min pH, temperature and NH3 values
df<-rbind(df_max,df_min)  %>% 
  mutate(date=as.Date(date),   year=year(date), month= month(date),
         day= day(date)) 


#arranged data by date
df <- df %>% arrange(date)  

#adding seasons
# april-October: Summer
# November-March : Winter

df$season <-  as.character(ifelse(df$month %in% c(11, 12,1,2, 3), "winter", 
                                  ifelse (df$month %in% c(4:10), "summer",  "" )))


#considering constant 67 MGD design flow for Norman Cole plant 
df_combined <- df %>% mutate (final_effluent_flow=67)

# ================= generating 30 years data using actual plants effluent data================

df_mon <- df_combined %>% filter(season=="summer")

#repeating pH, temperature data and matching to have 10960 days data 
# equivalent to 30 yrs data

df_mon_repeat<- df_mon[rep(1:nrow(df_mon),6),]

df_mon_repeat_extra <- df_mon [1:1240,]

df_mon_30yr <- rbind (df_mon_repeat, df_mon_repeat_extra)

# Chronic concentration CCC calculation
# https://www.epa.gov/sites/default/files/2015-08/documents/aquatic-life-ambient-water-quality-criteria-for-ammonia-freshwater-2013.pdf

df_mon_30yr <- df_mon_30yr  %>% dplyr:: select( pH, temperature, NH3_eff, 
                                                final_effluent_flow )

#coefficient of variation of plants effluent NH3 data
original_cv <- sd(df_mon_30yr$NH3_eff) / mean(df_mon_30yr$NH3_eff)
original_sd <- sd(df_mon_30yr$NH3_eff) 

# =========================================================================================
# =============================mixing with stream flow===========================
# =========================================================================================

# stream flow data was collected from USGS website (see below) for following station and 
# Interpolated flow data based on drainage area in Pohick Creek at Lower Potomac (32 sq mile)
# USGS 01654000 ACCOTINK CREEK NEAR ANNANDALE, VA (23 sq mile)

# Raw data
# https://waterdata.usgs.gov/nwis/dv?cb_00010=on&cb_00060=on&cb_00400=on&format=rdb&site_no=01654000&legacy=&referred_module=sw&period=&begin_date=2013-01-01&end_date=2023-11-30


#adding stream data
df_mon_stream <- read_csv("./data/ammonia_limit_monte_carlo/stream_data_all.csv") %>% 
  dplyr::select(-date)

date.vec <- data.frame(date=seq( as.Date("2013-01-01"), by=1, len=3986))


#for stream flow we are using last 10 years stream data for better representation of stream flow
#selecting last permit cycle stream flow data
df_mon_stream <- cbind(date.vec, df_mon_stream) %>%  
  filter( date >="2013-01-01" & date <= "2023-11-30") %>% 
  dplyr::select(date, discharge_cft_sc, temperature_max, 
                temperature_min, pH_max, pH_min) %>% 
  mutate(date=as.Date(date),   year=year(date), month= month(date),
         day= day(date))


#adding seasons
df_mon_stream$season <-  as.factor(ifelse(df_mon_stream$month %in% c(11, 12,1,2, 3), "winter", 
                                          ifelse (df_mon_stream$month %in% c(4:10), "summer",  "" )))


#taking summer stream data
df_mon_stream_summer <- df_mon_stream %>% filter(season== "summer")

days_repeat <- nrow(df_mon_stream_summer)


#stream data does not have ammonia value 
#Assuming low ammonia of stream and generating data based on lognormal distribution

set.seed(10000) 

mean_nh3 <- 0.0305
std_nh3 <-  0.0178
lower_limit <- 0.004
upper_limit <- 0.11
num_days <- 10960

# Generate lognormal distribution of stream ammonia value
lognormal_values <- rlnorm(num_days, meanlog = log(mean_nh3^2 / sqrt(mean_nh3^2 + std_nh3^2)), sdlog = sqrt(log(1 + std_nh3^2 / mean_nh3^2)))

# Adjust values to fit within the specified range
lognormal_values <- ifelse(lognormal_values<lower_limit, lower_limit, lognormal_values)
lognormal_values <- ifelse(lognormal_values>upper_limit, upper_limit, lognormal_values)
nh3_random_stream <- lognormal_values

df_summer_stream <- df_mon_stream_summer

#seperating and then combinging max and min pH, temperature of stream data   
df_summer_stream_max <- df_summer_stream %>% 
  dplyr::select(date, discharge_cft_sc, temperature_max,  pH_max)

#renaming variables
names(df_summer_stream_max)[names(df_summer_stream_max) == "discharge_cft_sc"]     <- "final_effluent_flow_stream"
names(df_summer_stream_max)[names(df_summer_stream_max) == "temperature_max"]     <- "temperature_stream"
names(df_summer_stream_max)[names(df_summer_stream_max) == "pH_max"]     <- "pH_stream"

df_summer_stream_min <- df_summer_stream %>% 
  dplyr::select(date, discharge_cft_sc, 
                temperature_min, pH_min)

names(df_summer_stream_min)[names(df_summer_stream_min) == "discharge_cft_sc"]     <- "final_effluent_flow_stream"
names(df_summer_stream_min)[names(df_summer_stream_min) == "temperature_min"]     <- "temperature_stream"
names(df_summer_stream_min)[names(df_summer_stream_min) == "pH_min"]     <- "pH_stream"

#binding max and min values of stream data
df_summer_stream <- rbind(df_summer_stream_max, df_summer_stream_min)

#arranging data by date
df_summer_stream <- df_summer_stream %>% arrange(date)  #df_summer_stream[sample(nrow(df_summer_stream)), ]


#repeating pH, temperature data and matching to have 10960 equivalent to 30 yrs data
df_stream_mon<- df_summer_stream[rep(1:nrow(df_summer_stream),2),]

df_stream_mon_extra <- df_summer_stream [1:1544,]

df_stream_mon <- rbind (df_stream_mon, df_stream_mon_extra) %>% dplyr::select(-date)

df_stream_mon$nh3_eff_stream <- nh3_random_stream

ggplot(df_stream_mon, aes(x = nh3_eff_stream)) +
  geom_density(alpha = 0.8) +
  labs(title = "Log-normal distributions of stream ammonia data",
       x = "Stream ammonia (mg/L NH3-N)",
       y = "Density") + theme_bw() + scale_x_continuous(trans='log10')  


#combining plant data with stream data
new_df_summer <-cbind (df_mon_30yr, df_stream_mon) 

#conversion from ft3/sec to MGD -->0.6463
#drainage area based flow interpolation
new_df_summer$interpolated_effluent_flow_stream <- new_df_summer$final_effluent_flow_stream * (32/23.5)*0.6463 

# pH= -log[H+]
#quantifying H+ ion concentration

new_df_summer <- new_df_summer %>% mutate(outfall_pH_H = 10^(-pH),
                                          outfall_temp = (temperature),
                                          stream_pH_H = 10^(-pH_stream),
                                          stream_temp = (temperature_stream),
                                          NH3_eff_plant= NH3_eff,
                                          nh3_eff_stream= nh3_eff_stream )

#mixing formula

# H+_mixing = (("H+_plant_eff ∗ Q_plant_flow +H+_steam ∗ Q_stream" _𝑓𝑙𝑜𝑤)/("Q_plant_flow +Q_stream" _𝑓𝑙𝑜𝑤)) = 
# pH=-log(H+)
# Temperature_mixing = (("temp_plant_eff ∗ Q_plant_flow +temp_steam ∗ Q_stream" _𝑓𝑙𝑜𝑤)/("Q_plant_flow +Q_stream" _𝑓𝑙𝑜𝑤))


new_df_summer <- new_df_summer%>% mutate(
  pH_mix= (final_effluent_flow * outfall_pH_H + interpolated_effluent_flow_stream * stream_pH_H )/ (final_effluent_flow + interpolated_effluent_flow_stream ),
  temp_mix= (final_effluent_flow * temperature + interpolated_effluent_flow_stream * temperature_stream)/ (final_effluent_flow + interpolated_effluent_flow_stream ),
  nh3_mix= (final_effluent_flow * NH3_eff_plant + interpolated_effluent_flow_stream * nh3_eff_stream)/ (final_effluent_flow + interpolated_effluent_flow_stream )
) %>%  dplyr::select( pH_mix, temp_mix, nh3_mix)  %>% 
  mutate(pH_mix= -log10(pH_mix)) 

# Replacing any NA values with mean value
new_df_summer$pH_mix <-  ifelse(new_df_summer$pH_mix%in% NA, 
                                rnorm(new_df_summer$pH_mix, mean = mean(new_df_summer$pH_mix, na.rm=T), sd = sd(new_df_summer$pH_mix, na.rm=T)), 
                                new_df_summer$pH_mix )

new_df_summer$temp_mix <-  ifelse(new_df_summer$temp_mix%in% NA, 
                                  rnorm(new_df_summer$temp_mix, mean = mean(new_df_summer$temp_mix, na.rm=T), sd = sd(new_df_summer$temp_mix, na.rm=T)), 
                                  new_df_summer$temp_mix)

new_df_summer$nh3_mix <-  ifelse(new_df_summer$nh3_mix%in% NA, 
                                 rnorm(new_df_summer$nh3_mix, mean = mean(new_df_summer$nh3_mix, na.rm=T), sd = sd(new_df_summer$nh3_mix, na.rm=T)), 
                                 new_df_summer$nh3_mix)

#exceedance criteria check
# The ammonia LTA that is associated with no more than 1-in-3 year exceedance of the chronic ammonia criteria. This is determined by setting up a 30-year mixing simulation with randomly generated effluent ammonia concentration (based on an assumed mean and CV, lognormal) , and iteratively changing the assumed mean effluent ammonia concentration until the Monte Carlo analysis show the 1-in-3 year average exceedance rate of the 30-day CCC. 

#Moving average from daily data
mc_pH_mix <- rollmean(new_df_summer$pH_mix , k = 30)
mc_temp_mix <- rollmean(new_df_summer$temp_mix , k = 30)
mc_nh3_effluent <- rollmean(new_df_summer$nh3_mix , k = 30)
mc_30_mavg <- data.frame(cbind(mc_pH_mix, mc_temp_mix, mc_nh3_effluent))

#chronic concentration value calculate for daily data
mc_30_mavg$ccc <-  (0.8876*(0.0278/(1+10^(7.688-mc_30_mavg$mc_pH_mix ))+1.1994/(1+10^(mc_30_mavg$mc_pH_mix -7.688))))*2.126*(10^(0.028*(20-max(mc_30_mavg$mc_temp_mix ,7))))


#Calculating exceedance for each day
mc_30_mavg$exceedance <-  ifelse(mc_30_mavg$mc_nh3_effluent>mc_30_mavg$ccc,  1, 0)

#need to adjust 30 day period exceedance to show 10 allowable exceedance in 30 year 
# to represent 1 in every three year exceedance

#splitting 30 years data to 30 days segments
n <- 30
nr <- nrow(mc_30_mavg)
x <- split(mc_30_mavg, gl(ceiling(nr/n), n, nr))

datalist=list ()

for (i in 1: length(x)) {
  y <- x[[i]]
  dt <- y
  dt$final.exceedance <- ifelse(cumsum(dt$exceedance) > 0, c(1, rep(0, nrow(dt)-1)), 0)
  datalist[[i]] <- dt
}

dt_fin <- do.call(rbind, datalist)


# mean effluent ammonia nitrogen exceeding once in every three year chronic criterion 

# Set seed for reproducibility
set.seed(10000)

# Target sum for the 'exceedance' column
target_exceedance_sum <-10

# Maximum attempts to regenerate data (to avoid infinite loop)
max_attempts <- 10000

# Initial sum of exceedance in 30 day period
initial_exceedance_sum <- sum(dt_fin$final.exceedance)


#adjusting mean effluent ammonia value to have 10 exceedance in 30 years to represent 1 exceedance in every three years
attempts <- 0

df_mon_30yr_x <- df_mon_30yr

while (initial_exceedance_sum != target_exceedance_sum && attempts < max_attempts) {
  
  df_mon_30yr2<- df_mon_30yr_x #plant data
  df_stream_mon2<- df_stream_mon #stream data
  target_cv <- original_cv 
  
  mean_amm <- mean(df_mon_30yr2$NH3_eff)
  std_amm <- mean(df_mon_30yr2$NH3_eff)  * target_cv
  
  
  if (initial_exceedance_sum < target_exceedance_sum) {
   
    # Increase mean.nh3 and regenerate data

    df_mon_30yr2$NH3_eff <-rlnorm(num_days, meanlog = log(mean_amm^2 / sqrt(mean_amm^2 + std_amm^2)) + 0.1, sdlog = sqrt(log(1 + std_amm^2 / mean_amm^2)))
    
    #replacing negative with 0
    df_mon_30yr2$NH3_eff <-  ifelse(df_mon_30yr2$NH3_eff<0, 0 ,df_mon_30yr2$NH3_eff )
    
    } else {
    # Decrease NH3_eff and regenerate data
    df_mon_30yr2$NH3_eff <-rlnorm(num_days, meanlog = log(mean_amm^2 / sqrt(mean_amm^2 + std_amm^2)) - 0.1, sdlog = sqrt(log(1 + std_amm^2 / mean_amm^2)))
    #replacing negative 
    df_mon_30yr2$NH3_eff <-  ifelse(df_mon_30yr2$NH3_eff<0, 0 ,df_mon_30yr2$NH3_eff )
    
    }
  
  #combining plant data with stream data
  new_df_summer <-cbind (df_mon_30yr2, df_stream_mon2) 
  
  
  new_df_summer$interpolated_effluent_flow_stream <- new_df_summer$final_effluent_flow_stream * (32/23.5)*0.6463 
  
  # pH= -log[H+]
  #quantifying H+ ion concentration
  
  new_df_summer <- new_df_summer %>% mutate(outfall_pH_H = 10^(-pH),
                                            outfall_temp = (temperature),
                                            stream_pH_H = 10^(-pH_stream),
                                            stream_temp = (temperature_stream),
                                            NH3_eff_plant= NH3_eff,
                                            nh3_eff_stream= nh3_eff_stream )
  
  #mixing formula
  
  # H+_mixing = (("H+_plant_eff ∗ Q_plant_flow +H+_steam ∗ Q_stream" _𝑓𝑙𝑜𝑤)/("Q_plant_flow +Q_stream" _𝑓𝑙𝑜𝑤)) = 
  # pH=-log(H+)
  
  # Temperature_mixing = (("temp_plant_eff ∗ Q_plant_flow +temp_steam ∗ Q_stream" _𝑓𝑙𝑜𝑤)/("Q_plant_flow +Q_stream" _𝑓𝑙𝑜𝑤))
  
  
  new_df_summer <- new_df_summer%>% mutate(
    pH_mix= (final_effluent_flow * outfall_pH_H + interpolated_effluent_flow_stream * stream_pH_H )/ (final_effluent_flow + interpolated_effluent_flow_stream ),
    temp_mix= (final_effluent_flow * temperature + interpolated_effluent_flow_stream * temperature_stream)/ (final_effluent_flow + interpolated_effluent_flow_stream ),
    nh3_mix= (final_effluent_flow * NH3_eff_plant + interpolated_effluent_flow_stream * nh3_eff_stream)/ (final_effluent_flow + interpolated_effluent_flow_stream )
  ) %>%  dplyr::select( pH_mix, temp_mix, nh3_mix)  %>% 
    mutate(pH_mix= -log10(pH_mix)) 
  
  # Replacing any NA values with mean value
  new_df_summer$pH_mix <-  ifelse(new_df_summer$pH_mix%in% NA, 
                                  rnorm(new_df_summer$pH_mix, mean = mean(new_df_summer$pH_mix, na.rm=T), sd = sd(new_df_summer$pH_mix, na.rm=T)), 
                                  new_df_summer$pH_mix )
  
  new_df_summer$temp_mix <-  ifelse(new_df_summer$temp_mix%in% NA, 
                                    rnorm(new_df_summer$temp_mix, mean = mean(new_df_summer$temp_mix, na.rm=T), sd = sd(new_df_summer$temp_mix, na.rm=T)), 
                                    new_df_summer$temp_mix)
  
  new_df_summer$nh3_mix <-  ifelse(new_df_summer$nh3_mix%in% NA, 
                                   rnorm(new_df_summer$nh3_mix, mean = mean(new_df_summer$nh3_mix, na.rm=T), sd = sd(new_df_summer$nh3_mix, na.rm=T)), 
                                   new_df_summer$nh3_mix)
  
  
  mc_pH_mix <- rollmean(new_df_summer$pH_mix , k = 30)
  mc_temp_mix <- rollmean(new_df_summer$temp_mix , k = 30)
  mc_nh3_effluent <- rollmean(new_df_summer$nh3_mix , k = 30)
  mc_30_mavg <- data.frame(cbind(mc_pH_mix, mc_temp_mix, mc_nh3_effluent))
  mc_30_mavg$ccc <-  (0.8876*(0.0278/(1+10^(7.688-mc_30_mavg$mc_pH_mix ))+1.1994/(1+10^(mc_30_mavg$mc_pH_mix -7.688))))*2.126*(10^(0.028*(20-max(mc_30_mavg$mc_temp_mix ,7))))
  

  
  mc_30_mavg$exceedance <-  ifelse(mc_30_mavg$mc_nh3_effluent>mc_30_mavg$ccc,  1, 0)
  
  #splitting 30 years data to 30 days segments
  n <- 30
  nr <- nrow(mc_30_mavg)
  x <- split(mc_30_mavg, gl(ceiling(nr/n), n, nr))
  
  datalist=list ()

  for (i in 1: length(x)) {
    y <- x[[i]]
    
    dt <- y
    dt$final.exceedance <- ifelse(cumsum(dt$exceedance) > 0, c(1, rep(0, nrow(dt)-1)), 0)
    datalist[[i]] <- dt
    
  }
  
  dt_fin <- do.call(rbind, datalist)
  
  initial_exceedance_sum <- sum(dt_fin$final.exceedance)
  
  df_mon_30yr_x$NH3_eff <- df_mon_30yr2$NH3_eff #adjusted plant effluent NH3
  # Increment attempts counter
  attempts <- attempts + 1
}

#checking coefficeint of variations of final effluent ammonia data
cv <- sd(df_mon_30yr2$NH3_eff)/ mean (df_mon_30yr2$NH3_eff)

cat(" coefficient of variation effluent ammonia data is", cv, "\n")

cat(" coefficient of variation original effluent ammonia data is", original_cv, "\n")

difference_from_original_cv <- abs(original_cv-cv)*100/original_cv

cat(" Difference of coefficient of variation from original CV is", difference_from_original_cv, "% \n")

# Print the results
cat(" ammonia LTA exceedance in plant effluent ammonia", sum(dt_fin$final.exceedance), "\n")

# cat("running 4 day average total ammonia concentration exceedance", sum(meta_exceedance$final.exceedance), "\n")

cat("Allowable Exceedance Sum:", target_exceedance_sum, "\n")


#check 2
## running four day average ammonia nitrogen not exceeding 2.5 times of chronic criterion within a 30 day period

mc_pH_mix <- rollmean(new_df_summer$pH_mix , k = 4)
mc_temp_mix <- rollmean(new_df_summer$temp_mix , k = 4)
mc_nh3_effluent <- rollmean(new_df_summer$nh3_mix , k = 4)
mc_4_mavg <- data.frame(cbind(mc_pH_mix, mc_temp_mix, mc_nh3_effluent))
mc_4_mavg$ccc <-  (0.8876*(0.0278/(1+10^(7.688-mc_4_mavg$mc_pH_mix ))+1.1994/(1+10^(mc_4_mavg$mc_pH_mix -7.688))))*2.126*(10^(0.028*(20-max(mc_4_mavg$mc_temp_mix ,7))))
mc_4_mavg<- mc_4_mavg %>% mutate(ccc2.5x= 2.5*ccc)


mc_4_mavg$exceedance <-  ifelse(mc_4_mavg$mc_nh3_effluent>mc_4_mavg$ccc2.5x,  1, 0)

#splitting 30 years data to 30 days segments
n <- 30
nr <- nrow(mc_4_mavg)
x <- split(mc_4_mavg, gl(ceiling(nr/n), n, nr))

datalist=list ()

for (i in 1: length(x)) {
  y <- x[[i]]
  dt <- y
  dt$final.exceedance <- ifelse(cumsum(dt$exceedance) > 0, c(1, rep(0, nrow(dt)-1)), 0)
  datalist[[i]] <- dt
  
}

dt_fin <- do.call(rbind, datalist)

cat(" ammonia LTA exceedance in plant effluent ammonia", sum(dt_fin$final.exceedance), "\n")
cat("Allowable Exceedance Sum:", target_exceedance_sum, "\n")

#This means running four day average ammonia nitrogen not exceeding 2.5 times of chronic criterion within a 30 day period condition also satisfied [allowable 10 exceedance in 30 years]

#mean ammonia
mean_ammonia <- mean(df_mon_30yr2$NH3_eff)

plot_a <-ggplot(df_mon_30yr2, aes(x = NH3_eff)) +
  geom_density(alpha = 0.8) +
  labs(title = "Lognormal distributions of effluent ammonia during summer",
       x = "Effluent ammonia (mg/L NH3-N)",
       y = "Density") + theme_bw() +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(limits = c(0,0.8)) +
  theme(axis.text= element_text(size=15),
        axis.title= element_text(size=20)) 


#monthly and weekly limit mixing of stream flow	
monthly<-30
weekly<-7
Z97<-	1.881	#Z-score

CV<-	original_cv #Coefficient of variation

σ2 <- log(CV^2+1)	
σ <- 	sqrt(σ2) 

σ230 <- log((CV^2/monthly)+1) 
σ30 <-  sqrt(σ230 )   

σ27 <- log((CV^2/weekly)+1) 
σ7 <-  sqrt(σ27 ) 


LTAc <-  mean_ammonia   #Chronic Long term average

LTAmin <-  min( LTAc, na.rm=T)	

n <- 	monthly	
n_7 <- weekly

σ2n_sample_30 <- log(CV^2/n+1)	
σn_sample_30 <- 	 sqrt(σ2n_sample_30) 

σ2n_sample_7 <- log(CV^2/n_7+1)	
σn_sample_7 <- 	 sqrt(σ2n_sample_7) 

#maximum daily limit
MDL <- 	LTAmin*exp((Z97*σ)-(0.5*σ2))

#average monthly limit
AML <-  LTAmin*exp(Z97*σn_sample_30-0.5*σ2n_sample_30)

AML
#average weekly limit
AWL <-  LTAmin*exp(Z97*σn_sample_7-0.5*σ2n_sample_7)

AWL

# =============================winter montecarlo calculation ==============================

df <- read_csv(("./data/ph_data/ph_data_NMCPCP.csv")) %>% 
  mutate(date=as.Date(date),   year=year(date), month= month(date),
         day= day(date)) %>%  filter( date >="2019-11-01" & date <= "2023-11-30")

#combining max and min pH, and temperature of plants data
df_max <- df %>% 
  dplyr:: select(date, max_pH_outfall, Temperature_Out, final_effluent_flow, NH3_eff) %>% 
  na.omit()

names(df_max)[names(df_max) == "max_pH_outfall"]     <- "pH"
names(df_max)[names(df_max) == "Temperature_Out"]     <- "temperature"

df_min <- df %>% 
  dplyr:: select(date, min_pH_outfall, Temperature_Out, final_effluent_flow, NH3_eff) %>% 
  na.omit()


names(df_min)[names(df_min) == "min_pH_outfall"]     <- "pH"
names(df_min)[names(df_min) == "Temperature_Out"]     <- "temperature"

#combining max and min pH, temperature and NH3 values
df<-rbind(df_max,df_min)  %>% 
  mutate(date=as.Date(date),   year=year(date), month= month(date),
         day= day(date)) 


df <- df %>% arrange(date)


#adding seasons
# april-October: Summer
# November-March : Winter

df$season <-  as.character(ifelse(df$month %in% c(11, 12,1,2, 3), "winter", 
                                  ifelse (df$month %in% c(4:10), "summer",  "" )))


#considering constant 67 MGD design flow for Norman Cole plant 
df_combined <- df %>% mutate (final_effluent_flow=67)

# ================= generating 30 years data using actual plants data================

df_mon <- df_combined %>% filter(season=="winter")


#repeating pH, temperature data and matching to have 10960 days data 
# equivalent to 30 yrs data

df_mon_repeat<- df_mon[rep(1:nrow(df_mon),9),]

df_mon_repeat_extra <- df_mon [1:484,]

df_mon_30yr <- rbind (df_mon_repeat, df_mon_repeat_extra)

# Chronic concentration CCC calculation
# https://www.epa.gov/sites/default/files/2015-08/documents/aquatic-life-ambient-water-quality-criteria-for-ammonia-freshwater-2013.pdf

df_mon_30yr <- df_mon_30yr  %>% dplyr:: select( pH, temperature, NH3_eff, 
                                                final_effluent_flow )

original_cv <- sd(df_mon_30yr$NH3_eff) / mean(df_mon_30yr$NH3_eff)
original_sd <- sd(df_mon_30yr$NH3_eff) 

# =========================================================================================
# =============================mixing with stream flow===========================
# =========================================================================================

# stream flow data was collected from USGS website (see below) for following station and 
# Interpolated flow data based on drainage area in Pohick Creek at Lower Potomac (32 sq mile)
# USGS 01654000 ACCOTINK CREEK NEAR ANNANDALE, VA (23 sq mile)

# Raw data
# https://waterdata.usgs.gov/nwis/dv?cb_00010=on&cb_00060=on&cb_00400=on&format=rdb&site_no=01654000&legacy=&referred_module=sw&period=&begin_date=2013-01-01&end_date=2023-11-30


#adding stream data
df_mon_stream <- read_csv("./data/ammonia_limit_monte_carlo/stream_data_all.csv") %>% 
  dplyr::select(-date)

date.vec <- data.frame(date=seq( as.Date("2013-01-01"), by=1, len=3986))


#selecting last permit cycle stream flow data

df_mon_stream <- cbind(date.vec, df_mon_stream) %>%  
  filter( date >="2013-01-01" & date <= "2023-11-30") %>% 
  dplyr::select(date, discharge_cft_sc, temperature_max, 
                temperature_min, pH_max, pH_min) %>% 
  mutate(date=as.Date(date),   year=year(date), month= month(date),
         day= day(date))


#adding seasons
df_mon_stream$season <-  as.factor(ifelse(df_mon_stream$month %in% c(11, 12,1,2, 3), "winter", 
                                          ifelse (df_mon_stream$month %in% c(4:10), "summer",  "" )))


#taking summer stream data
df_mon_stream_winter <- df_mon_stream %>% filter(season== "winter")

days_repeat <- nrow(df_mon_stream_winter)

#stream data does not have ammonia value 
#Assuming low ammonia of stream and generating data based on lognormal distribution

set.seed(10000) 

mean_nh3 <- 0.0305
std_nh3 <-  0.0178
lower_limit <- 0.004
upper_limit <- 0.11
num_days <- 10960

# Generate lognormal distribution of stream ammonia value
lognormal_values <- rlnorm(num_days, meanlog = log(mean_nh3^2 / sqrt(mean_nh3^2 + std_nh3^2)), sdlog = sqrt(log(1 + std_nh3^2 / mean_nh3^2)))

# Adjust values to fit within the specified range
lognormal_values <- ifelse(lognormal_values<lower_limit, lower_limit, lognormal_values)
lognormal_values <- ifelse(lognormal_values>upper_limit, upper_limit, lognormal_values)
nh3_random_stream <- lognormal_values

df_winter_stream <- df_mon_stream_winter

#seperating and then combinging max and min pH, temperature stream data   
df_winter_stream_max <- df_winter_stream %>% 
  dplyr::select(date, discharge_cft_sc, temperature_max,  pH_max)

#renaming variables
names(df_winter_stream_max)[names(df_winter_stream_max) == "discharge_cft_sc"]     <- "final_effluent_flow_stream"
names(df_winter_stream_max)[names(df_winter_stream_max) == "temperature_max"]     <- "temperature_stream"
names(df_winter_stream_max)[names(df_winter_stream_max) == "pH_max"]     <- "pH_stream"

df_winter_stream_min <- df_winter_stream %>% 
  dplyr::select(date, discharge_cft_sc, 
                temperature_min, pH_min)

names(df_winter_stream_min)[names(df_winter_stream_min) == "discharge_cft_sc"]     <- "final_effluent_flow_stream"
names(df_winter_stream_min)[names(df_winter_stream_min) == "temperature_min"]     <- "temperature_stream"
names(df_winter_stream_min)[names(df_winter_stream_min) == "pH_min"]     <- "pH_stream"

df_winter_stream <- rbind(df_winter_stream_max, df_winter_stream_min)


#randomly shuffled data
df_winter_stream <- df_winter_stream %>% arrange(date)


#repeating pH, temperature data and matching to have 10960 equivalent to 30 yrs data
df_stream_mon<- df_winter_stream[rep(1:nrow(df_winter_stream),3),]

df_stream_mon_extra <- df_winter_stream [1:1168,]

df_stream_mon <- rbind (df_stream_mon, df_stream_mon_extra) %>% dplyr::select(-date)

df_stream_mon$nh3_eff_stream <- nh3_random_stream

ggplot(df_stream_mon, aes(x = nh3_eff_stream)) +
  geom_density(alpha = 0.8) +
  labs(title = "Lognormal distributions stream ammonia data during summer",
       x = "Stream ammonia (mg/L NH3-N)",
       y = "Density") + theme_bw() +
  scale_x_continuous(trans='log10') +
   theme(axis.text= element_text(size=15),
        axis.title= element_text(size=20))


#combining plant data with stream data
new_df_winter <-cbind (df_mon_30yr, df_stream_mon) 

#average seasonal stream flow
# stream_mean_flow_winter <- mean(new_df_winter$final_effluent_flow_stream)


# pH= -log[H+]
#quantifying H+ ion concentration

new_df_winter <- new_df_winter %>% mutate(outfall_pH_H = 10^(-pH),
                                          outfall_temp = (temperature),
                                          stream_pH_H = 10^(-pH_stream),
                                          stream_temp = (temperature_stream),
                                          NH3_eff_plant= NH3_eff,
                                          nh3_eff_stream= nh3_eff_stream )

#conversion from ft3/sec to MGD
#drainage area based flow interpolation
new_df_winter$interpolated_effluent_flow_stream <- new_df_winter$final_effluent_flow_stream * (32/23.5)*0.6463 


#mixing formula

# H+_mixing = (("H+_plant_eff ∗ Q_plant_flow +H+_steam ∗ Q_stream" _𝑓𝑙𝑜𝑤)/("Q_plant_flow +Q_stream" _𝑓𝑙𝑜𝑤)) = 
# pH=-log(H+)
# Temperature_mixing = (("temp_plant_eff ∗ Q_plant_flow +temp_steam ∗ Q_stream" _𝑓𝑙𝑜𝑤)/("Q_plant_flow +Q_stream" _𝑓𝑙𝑜𝑤))


new_df_winter <- new_df_winter%>% mutate(
  pH_mix= (final_effluent_flow * outfall_pH_H + interpolated_effluent_flow_stream * stream_pH_H )/ (final_effluent_flow + interpolated_effluent_flow_stream ),
  temp_mix= (final_effluent_flow * temperature + interpolated_effluent_flow_stream * temperature_stream)/ (final_effluent_flow + interpolated_effluent_flow_stream ),
  nh3_mix= (final_effluent_flow * NH3_eff_plant + interpolated_effluent_flow_stream * nh3_eff_stream)/ (final_effluent_flow + interpolated_effluent_flow_stream )
) %>%  dplyr::select( pH_mix, temp_mix, nh3_mix)  %>% 
  mutate(pH_mix= -log10(pH_mix)) 



# Replacing any NA value with mean value
new_df_winter$pH_mix <-  ifelse(new_df_winter$pH_mix%in% NA, 
                                rnorm(new_df_winter$pH_mix, mean = mean(new_df_winter$pH_mix, na.rm=T), sd = sd(new_df_winter$pH_mix, na.rm=T)), 
                                new_df_winter$pH_mix )

new_df_winter$temp_mix <-  ifelse(new_df_winter$temp_mix%in% NA, 
                                  rnorm(new_df_winter$temp_mix, mean = mean(new_df_winter$temp_mix, na.rm=T), sd = sd(new_df_winter$temp_mix, na.rm=T)), 
                                  new_df_winter$temp_mix)

new_df_winter$nh3_mix <-  ifelse(new_df_winter$nh3_mix%in% NA, 
                                 rnorm(new_df_winter$nh3_mix, mean = mean(new_df_winter$nh3_mix, na.rm=T), sd = sd(new_df_winter$nh3_mix, na.rm=T)), 
                                 new_df_winter$nh3_mix)

#2 exceedance criteria
# a.	The ammonia LTA that is associated with no more than 1-in-3 year exceedance of the chronic ammonia criteria. This is determined by setting up a 30-year mixing simulation with randomly generated effluent ammonia concentration (based on an assumed mean and CV, lognormal) , and iteratively changing the assumed mean effluent ammonia concentration until the Monte Carlo analysis show the 1-in-3 year average exceedance rate of the 30-day CCC. 
# b.	The running four-day average concentration of total ammonia nitrogen (in mg N/L), for verification that it would not exceed 2.5 times the chronic criterion within a 30-day period more than once every three years on the average.

#Moving average from daily data
mc_pH_mix <- rollmean(new_df_winter$pH_mix , k = 30)
mc_temp_mix <- rollmean(new_df_winter$temp_mix , k = 30)
mc_nh3_effluent <- rollmean(new_df_winter$nh3_mix , k = 30)
mc_30_mavg <- data.frame(cbind(mc_pH_mix, mc_temp_mix, mc_nh3_effluent))

#chronic concentration value calculate for daily data
mc_30_mavg$ccc <-  (0.8876*(0.0278/(1+10^(7.688-mc_30_mavg$mc_pH_mix ))+1.1994/(1+10^(mc_30_mavg$mc_pH_mix -7.688))))*2.126*(10^(0.028*(20-max(mc_30_mavg$mc_temp_mix ,7))))



#Exceedance each day
mc_30_mavg$exceedance <-  ifelse(mc_30_mavg$mc_nh3_effluent>mc_30_mavg$ccc,  1, 0)

#need to adjust 30 day period exceedance to show 10 allowable exceedance in 30 year 
# to represent 1 in every three year exceedance

#splitting 30 years data to 30 days segments
n <- 30
nr <- nrow(mc_30_mavg)
x <- split(mc_30_mavg, gl(ceiling(nr/n), n, nr))

datalist=list ()

for (i in 1: length(x)) {
  y <- x[[i]]
  dt <- y
  dt$final.exceedance <- ifelse(cumsum(dt$exceedance) > 0, c(1, rep(0, nrow(dt)-1)), 0)
  datalist[[i]] <- dt
}

dt_fin <- do.call(rbind, datalist)


# mean effluent ammonia nitrogen exceeding once in every three year chronic criterion 

# Set seed for reproducibility
set.seed(10000)

# Target sum for the 'exceedance' column
target_exceedance_sum <-10

# Maximum attempts to regenerate data (to avoid infinite loop)
max_attempts <- 1000

# Initial sum of exceedance in 30 day period
initial_exceedance_sum <- sum(dt_fin$final.exceedance)

#Check if initial exceedance sum is not equal to the target
attempts <- 0

df_mon_30yr_x <- df_mon_30yr

while (initial_exceedance_sum != target_exceedance_sum && attempts < max_attempts) {
  
  df_mon_30yr2<- df_mon_30yr_x #plant data
  df_stream_mon2<- df_stream_mon #stream data
  target_cv <- original_cv 
  
  mean_amm <- mean(df_mon_30yr2$NH3_eff)
  std_amm <- mean(df_mon_30yr2$NH3_eff)  * target_cv
  
  
  if (initial_exceedance_sum < target_exceedance_sum) {
    
    # Increase mean.nh3 and regenerate data
    
    df_mon_30yr2$NH3_eff <-rlnorm(num_days, meanlog = log(mean_amm^2 / sqrt(mean_amm^2 + std_amm^2)) + 0.1, sdlog = sqrt(log(1 + std_amm^2 / mean_amm^2)))
    
    #replacing negative 
    df_mon_30yr2$NH3_eff <-  ifelse(df_mon_30yr2$NH3_eff<0, 0 ,df_mon_30yr2$NH3_eff )
    
  } else {
    # Decrease NH3_eff and regenerate data
    # df_mon_30yr2$NH3_eff <-  rnorm(nrow(df_mon_30yr2), mean(df_mon_30yr2$NH3_eff) , sd = (mean(df_mon_30yr2$NH3_eff)) * target_cv)
    
    df_mon_30yr2$NH3_eff <-rlnorm(num_days, meanlog = log(mean_amm^2 / sqrt(mean_amm^2 + std_amm^2)) - 0.1, sdlog = sqrt(log(1 + std_amm^2 / mean_amm^2)))
    #replacing negative 
    df_mon_30yr2$NH3_eff <-  ifelse(df_mon_30yr2$NH3_eff<0, 0 ,df_mon_30yr2$NH3_eff )
    
  }
  
  #combining plant data with stream data
  new_df_winter <-cbind (df_mon_30yr2, df_stream_mon2) 
  
  #conversion from ft3/sec to MGD
  #drainage area based flow interpolation
  new_df_winter$interpolated_effluent_flow_stream <- new_df_winter$final_effluent_flow_stream * (32/23.5)*0.6463 
  
  
  #average seasonal stream flow
  # stream_mean_flow_winter <- mean(new_df_winter$final_effluent_flow_stream)
  
  
  # pH= -log[H+]
  #quantifying H+ ion concentration
  
  new_df_winter <- new_df_winter %>% mutate(outfall_pH_H = 10^(-pH),
                                            outfall_temp = (temperature),
                                            stream_pH_H = 10^(-pH_stream),
                                            stream_temp = (temperature_stream),
                                            NH3_eff_plant= NH3_eff,
                                            nh3_eff_stream= nh3_eff_stream )
  
  #mixing formula
  
  # H+_mixing = (("H+_plant_eff ∗ Q_plant_flow +H+_steam ∗ Q_stream" _𝑓𝑙𝑜𝑤)/("Q_plant_flow +Q_stream" _𝑓𝑙𝑜𝑤)) = 
  #   pH=-log(H+)
  # 
  # Temperature_mixing = (("temp_plant_eff ∗ Q_plant_flow +temp_steam ∗ Q_stream" _𝑓𝑙𝑜𝑤)/("Q_plant_flow +Q_stream" _𝑓𝑙𝑜𝑤))
  
  
  new_df_winter <- new_df_winter%>% mutate(
    pH_mix= (final_effluent_flow * outfall_pH_H + interpolated_effluent_flow_stream * stream_pH_H )/ (final_effluent_flow + interpolated_effluent_flow_stream ),
    temp_mix= (final_effluent_flow * temperature + interpolated_effluent_flow_stream * temperature_stream)/ (final_effluent_flow + interpolated_effluent_flow_stream ),
    nh3_mix= (final_effluent_flow * NH3_eff_plant + interpolated_effluent_flow_stream * nh3_eff_stream)/ (final_effluent_flow + interpolated_effluent_flow_stream )
  ) %>%  dplyr::select( pH_mix, temp_mix, nh3_mix)  %>% 
    mutate(pH_mix= -log10(pH_mix)) 
  
  # Replacing any NA value with mean value
  # Replacing any NA value with mean value
  new_df_winter$pH_mix <-  ifelse(new_df_winter$pH_mix%in% NA, 
                                  rnorm(new_df_winter$pH_mix, mean = mean(new_df_winter$pH_mix, na.rm=T), sd = sd(new_df_winter$pH_mix, na.rm=T)), 
                                  new_df_winter$pH_mix )
  
  new_df_winter$temp_mix <-  ifelse(new_df_winter$temp_mix%in% NA, 
                                    rnorm(new_df_winter$temp_mix, mean = mean(new_df_winter$temp_mix, na.rm=T), sd = sd(new_df_winter$temp_mix, na.rm=T)), 
                                    new_df_winter$temp_mix)
  
  new_df_winter$nh3_mix <-  ifelse(new_df_winter$nh3_mix%in% NA, 
                                   rnorm(new_df_winter$nh3_mix, mean = mean(new_df_winter$nh3_mix, na.rm=T), sd = sd(new_df_winter$nh3_mix, na.rm=T)), 
                                   new_df_winter$nh3_mix)
  
  
  mc_pH_mix <- rollmean(new_df_winter$pH_mix , k = 30)
  mc_temp_mix <- rollmean(new_df_winter$temp_mix , k = 30)
  mc_nh3_effluent <- rollmean(new_df_winter$nh3_mix , k = 30)
  mc_30_mavg <- data.frame(cbind(mc_pH_mix, mc_temp_mix, mc_nh3_effluent))
  mc_30_mavg$ccc <-  (0.8876*(0.0278/(1+10^(7.688-mc_30_mavg$mc_pH_mix ))+1.1994/(1+10^(mc_30_mavg$mc_pH_mix -7.688))))*2.126*(10^(0.028*(20-max(mc_30_mavg$mc_temp_mix ,7))))
  

  
  mc_30_mavg$exceedance <-  ifelse(mc_30_mavg$mc_nh3_effluent>mc_30_mavg$ccc,  1, 0)
  
  #splitting 30 years data to 30 days segments
  n <- 30
  nr <- nrow(mc_30_mavg)
  x <- split(mc_30_mavg, gl(ceiling(nr/n), n, nr))
  
  datalist=list ()
  
  for (i in 1: length(x)) {
    y <- x[[i]]
    
    dt <- y
    dt$final.exceedance <- ifelse(cumsum(dt$exceedance) > 0, c(1, rep(0, nrow(dt)-1)), 0)
    datalist[[i]] <- dt
    
  }
  
  dt_fin <- do.call(rbind, datalist)
  
  initial_exceedance_sum <- sum(dt_fin$final.exceedance)
  
  df_mon_30yr_x$NH3_eff <- df_mon_30yr2$NH3_eff #adjusted plant effluent NH3
  # Increment attempts counter
  attempts <- attempts + 1
}

#checking coefficeint of variations of final effluent ammonia data
cv <- sd(df_mon_30yr2$NH3_eff)/ mean (df_mon_30yr2$NH3_eff)

cat(" coefficient of variation effluent ammonia data is", cv, "\n")

cat(" coefficient of variation original effluent ammonia data is", original_cv, "\n")

difference_from_original_cv <- abs(original_cv-cv)*100/original_cv

cat(" Difference of coefficient of variation from original CV is", difference_from_original_cv, "% \n")

# Print the results
cat(" ammonia LTA exceedance in plant effluent ammonia", sum(dt_fin$final.exceedance), "\n")

# cat("running 4 day average total ammonia concentration exceedance", sum(meta_exceedance$final.exceedance), "\n")

cat("Target Exceedance Sum:", target_exceedance_sum, "\n")

#check 2
## running four day average ammonia nitrogen not exceeding 2.5 times of chronic criterion within a 30 day period

mc_pH_mix <- rollmean(new_df_winter$pH_mix , k = 4)
mc_temp_mix <- rollmean(new_df_winter$temp_mix , k = 4)
mc_nh3_effluent <- rollmean(new_df_winter$nh3_mix , k = 4)
mc_4_mavg <- data.frame(cbind(mc_pH_mix, mc_temp_mix, mc_nh3_effluent))
mc_4_mavg$ccc <-  (0.8876*(0.0278/(1+10^(7.688-mc_4_mavg$mc_pH_mix ))+1.1994/(1+10^(mc_4_mavg$mc_pH_mix -7.688))))*2.126*(10^(0.028*(20-max(mc_4_mavg$mc_temp_mix ,7))))
mc_4_mavg<- mc_4_mavg %>% mutate(ccc2.5x= 2.5*ccc)


mc_4_mavg$exceedance <-  ifelse(mc_4_mavg$mc_nh3_effluent>mc_4_mavg$ccc2.5x,  1, 0)

#splitting 30 years data to 30 days segments
n <- 30
nr <- nrow(mc_4_mavg)
x <- split(mc_4_mavg, gl(ceiling(nr/n), n, nr))

datalist=list ()

for (i in 1: length(x)) {
  y <- x[[i]]
  
  dt <- y
  dt$final.exceedance <- ifelse(cumsum(dt$exceedance) > 0, c(1, rep(0, nrow(dt)-1)), 0)
  datalist[[i]] <- dt
  
}

dt_fin <- do.call(rbind, datalist)

sum(dt_fin$final.exceedance)

cat(" ammonia LTA exceedance in plant effluent ammonia", sum(dt_fin$final.exceedance), "\n")
cat("Allowable Exceedance Sum:", target_exceedance_sum, "\n")

#This means running four day average ammonia nitrogen not exceeding 2.5 times of chronic criterion within a 30 day period condition also satisfied

#mean ammonia
mean_ammonia <- mean(df_mon_30yr2$NH3_eff)

plot_b <- ggplot(df_mon_30yr2, aes(x = NH3_eff)) +
  geom_density(alpha = 0.8) +
  labs(title = "Lognormal distributions effluent ammonia during winter",
       x = "Effluent ammonia (mg/L NH3-N)",
       y = "Density") + theme_bw() +
  scale_x_continuous(trans='log10', limits = c(0.01,50)) +
  scale_y_continuous(limits = c(0,0.8)) +
  theme(axis.text= element_text(size=15),
        axis.title= element_text(size=20)) 

cowplot::plot_grid(plot_a, plot_b,
                   # 1 column and two rows - stacked on top of each other
                   ncol = 2,
                   nrow = 1,
                   # top plot is 2/3 as tall as second
                   rel_heights = c(0.5, 0.5))

ggsave("distributions.png", path = "./plots/", width=12, height=5, units="in")



#monthly and weekly limit without mixing of stream flow	
monthly<-30
weekly<-7
Z97<-	1.881	
CV<-	original_cv #Coefficient of variaion
σ2 <- log(CV^2+1)	
σ <- 	sqrt(σ2) 
σ230 <- log((CV^2/monthly)+1) 
σ30 <-  sqrt(σ230 )   
σ27 <- log((CV^2/weekly)+1) 
σ7 <-  sqrt(σ27 ) 


LTAc <-  mean_ammonia #Chronic Long term average

LTAmin <-  LTAc

n <- 	monthly	
n_7 <- weekly

σ2n_sample_30 <- log(CV^2/n+1)	
σn_sample_30 <- 	 sqrt(σ2n_sample_30) 

σ2n_sample_7 <- log(CV^2/n_7+1)	
σn_sample_7 <- 	 sqrt(σ2n_sample_7) 

#maximum daily limit
MDL <- 	LTAmin*exp((Z97*σ)-(0.5*σ2))

#average monthly limit
AML <-  LTAmin*exp(Z97*σn_sample_30-0.5*σ2n_sample_30)

AML
#average weekly limit
AWL <-  LTAmin*exp(Z97*σn_sample_7-0.5*σ2n_sample_7)

AWL

