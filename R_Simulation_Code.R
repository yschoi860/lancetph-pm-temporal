################################################################################
#
# R code for the analysis in:
#  Assessing the Temporal Change in Mortality Risk from Particulate Matter (PM10): An Analysis Across 143 Cities in 26 Countries
# 
#  ## IMPORTANT ##
#  ## THE DATASET INCLUDED IN THIS ANALYSIS IS NOT THE DATA USED FOR THE PAPER
#  ## THIS IS A SIMULATION DATASET MADE FOR ILLUSTRATIVE PURPOSES
#  ## THE RESULTS OBTAINED FROM THIS R CODE WILL DIFFER FROM THE MAIN STUDY
#
################################################################################

# LOAD PACKAGES
library(tidyverse); library(dlnm); library(splines); library(mixmeta); library(readr)



################################################################################
# 1ST STAGE: ESTIMATE LOCATION AND TIME SPECIFIC ESTIMATES
################################################################################

# LOAD DATA
df <- read_csv("1_Data/R_SIM_DATA.csv", col_types = cols(date = col_date(format = "%Y-%m-%d")))

# SCALE PM10 LEVEL BY 10
df <- df %>% mutate(pm10 = pm10/10)


################################################################################
### 1-1. ESTIMAT TIME-SPECIFIC ESTIMAT FOR A SINGLE CITY

# FILTER DATA FOR ONE CITY 
df_example <- df %>% filter(cityname == "City1") %>%
   arrange(date)

# SEPERATE THE DATASET INTO 3-YEAR TIME WINDOWS
years <- first(df_example$year):last(df_example$year)
intervals <- split(years, ceiling(seq_along(years)/3))
df_example_tw <- map2(intervals, intervals, ~{
   df_example %>%
      filter(year %in% .x) %>%
      mutate(mv_period = paste(min(.y), "~", max(.y)))  # ADD A NEW COLUMN CONTAINING INFORMATION ABOUT THE TIME WINDOWS
})


# ITERATE THE ANALYSIS OVER DIFFERENT PERIODS.
# USED THE MAP_DFR FUNCTION, WHICH IS SIMILAR TO A FOR LOOP. HOWEVER, THE DIFFERENCE IS THAT THE MAP FUNCTION ITERATES OVER A USER-DEFINED FUNCTION.
# MORE INFORMATION ABOUT THE MAP FUNCTION CAN BE FOUND HERE: https://www.rdocumentation.org/packages/purrr/versions/0.2.5/topics/map

fin_example <- df_example_tw %>% 
   map_dfr(
      function(x){
         
         df_mf <- as.data.frame(x)
         
         cb.pm <- crossbasis(df_mf$pm10, lag=2, argvar=list(fun="identity"))
         cb.tmean <- crossbasis(df_mf$tmean, lag=4,
                                argvar=list(fun="ns", df = 6),
                                arglag=list(fun ="identity"))
         sp_model <- glm(all ~ cb.pm + cb.tmean + ns(date, df=7*length(unique(df_mf$year))) + as.factor(dow), data=df_mf, family=quasipoisson)
         
         fin_cp <- crossreduce(cb.pm, sp_model, type="overall", cen=0, from=0, to=200, by = 1)
         
         mid_year <- as.numeric(unique(substr(df_mf$mv_period, 1, 4))) + 1
         afin <- as.data.frame(matrix(data = c(unique(df_mf$countryname), unique(df_mf$cityname), unique(df_mf$mv_period), mid_year,
                                               fin_cp$RRfit[2], fin_cp$RRlow[2], fin_cp$RRhigh[2], 
                                               coef(fin_cp), fin_cp$se[2], vcov(fin_cp)), nrow=1)) %>%
            setNames(c("country", "city", "period", "mid_year",
                       "RR", "RRlow", "RRhigh", 
                       "coef", "se",  "vcov")) %>%
            mutate(across(mid_year:vcov, ~as.numeric(.)))
         
         return(afin)
      }
   )

### 1-1 FIN
################################################################################



################################################################################
### 1-2. ESTIMATE LOCATION-AND TIME-SPECIFIC ESTIMATES USING FUNCTION (ESTIMATIONS ARE THE SAME AS IN 1-1)

# CLEAN ENVIRONMENT
rm(list=ls())

# LOAD DATA
df <- read_csv("1_Data/R_SIM_DATA.csv", col_types = cols(date = col_date(format = "%Y-%m-%d")))

# SCALE PM10 LEVEL BY 10
df <- df %>% mutate(pm10 = pm10/10)

# DEFINE FUNCTION
f.temporal_lv1 <- function(m_area){
   
   # DF_INF CONTAINS DATA FOR ONLY ONE CITY
   df_inf <- df %>% filter(cityname == m_area) %>%
      arrange(date)
   
   # SEPERATE THE DATASET INTO 3-YEAR TIME WINDOWS
   years <- first(df_inf$year):last(df_inf$year)
   intervals <- split(years, ceiling(seq_along(years)/3))
   df_inf_tw <- map2(intervals, intervals, ~{
      df_inf %>%
         filter(year %in% .x) %>%
         mutate(mv_period = paste(min(.y), "~", max(.y)))  # ADD A NEW COLUMN CONTAINING INFORMATION ABOUT THE TIME WINDOWS
   })
   
   
   # ITERATE THE RESULT OVER DIFFERENT TIME PERIOD 
   # SAME FUNCTION USED IN 1-1.
   fin_mf <- df_inf_tw %>% 
      map_dfr(
         function(x){
            df_mf <- as.data.frame(x)
            
            cb.pm <- crossbasis(df_mf$pm10, lag=2, argvar=list(fun="identity"))
            cb.tmean <- crossbasis(df_mf$tmean, lag=4,
                                   argvar=list(fun="ns", df = 6),
                                   arglag=list(fun ="identity"))
            sp_model <- glm(all ~ cb.pm + cb.tmean + ns(date, df=7*length(unique(df_mf$year))) + as.factor(dow), data=df_mf, family=quasipoisson)
            
            fin_cp <- crossreduce(cb.pm, sp_model, type="overall", cen=0, from=0, to=200, by = 1)
            
            mid_year <- as.numeric(unique(substr(df_mf$mv_period, 1, 4))) + 1
            afin <- as.data.frame(matrix(data = c(unique(df_mf$countryname), unique(df_mf$cityname), unique(df_mf$mv_period), mid_year,
                                                  fin_cp$RRfit[2], fin_cp$RRlow[2], fin_cp$RRhigh[2], 
                                                  coef(fin_cp), fin_cp$se[2], vcov(fin_cp)), nrow=1)) %>%
               setNames(c("country", "city", "period", "mid_year",
                          "RR", "RRlow", "RRhigh", 
                          "coef", "se",  "vcov")) %>%
               mutate(across(mid_year:vcov, ~as.numeric(.)))
            
            return(afin)
         }
      )
   return(fin_mf)
}


# DO ITERATION
df_lv2 <- unique(df$cityname) %>% map_dfr(f.temporal_lv1)

# ADD SES DATA (ONLY INCLUDED PERCENTAGE AGED 65 AND MORE JUST FOR ILLUSTRATIVE PURPOSE)
perc_aged <- read.csv("1_Data/R_SIM_DATA_SES.csv") 

# MERGE SES DATASET WITH THE MAIN RESULT
df_lv2 <- left_join(df_lv2, perc_aged, by = c("city" = "city", "mid_year" = "year"))



### 1-2 FIN
################################################################################




################################################################################
# 2ND STAGE: ESTIMATE TEMPORAL TREND OF PM10 MORTALITY EFFECTS USING LONGITUDINAL META REGRESSION MODEL
################################################################################

# CLEAR THE ENVIRONMENT, RETAINING ONLY THE RESULTS FROM THE FIRST STAGE
rm(list=setdiff(ls(), ls(pattern = "df_lv2")))


# ESTIMATE TEMPORAL TREND OF PM10 MORTALITY EFFECTS USING LONGITUDINAL META REGRESSION MODEL
m.main  <- mixmeta(coef ~ mid_year,
                   vcov, method ="ml", bscov="diag", data = df_lv2, random = ~1|country/city)


# CALCULATE PERCENT INCREASE OVER YEAR AND ITS 95% CONFIDENCE INTERVAL
sprintf("%.4f", (exp(summary(m.main)$coefficients[2,1])-1)*100)
sprintf("%.4f", (exp(summary(m.main)$coefficients[2,1] - 1.96 * summary(m.main)$coefficients[2,2])-1)*100)
sprintf("%.4f", (exp(summary(m.main)$coefficients[2,1] + 1.96 * summary(m.main)$coefficients[2,2])-1)*100)


# PREDICT RESULTS THROUGHOUT THE THE STUDY PERIOD
pred <- exp(predict(m.main, newdata=list(mid_year=1979:2019), ci=T)) %>%
   as.data.frame() %>% 
   mutate(year=1979:2019) %>%
   mutate(across(c("fit", "ci.lb", "ci.ub"), ~(.x-1)*100))

# DRAW PLOT
ggplot(pred, aes(year, fit)) + 
   geom_ribbon(aes(ymin=ci.lb, ymax=ci.ub), fill="lightblue1") + geom_line(color="black") +
   theme_classic() + ylab(expression(atop("Percent increase in mortality risk", paste("per 10", mu, "g/m"^3, " increase of PM"[10]))))


# EXPLAIN THE TEMPORAL TREND USING A SES VARIABLE (PERCENTAGE OF POPULATION AGED 65 AND OVER)
m2.main  <- mixmeta(coef ~ mid_year + perc_aged,
                    vcov, method ="ml", bscov="diag", data = df_lv2, random = ~1|country/city)
summary.mixmeta(m2.main)

# PERFORM THE LR (LIKELIHOOD RATIO) TEST 
drop1(m2.main, test="Chisq")


