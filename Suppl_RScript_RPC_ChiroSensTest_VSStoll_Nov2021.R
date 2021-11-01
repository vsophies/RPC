############### RPC: The effect of Bacillus thuringiensis subsp. israelensis on Chironomus riparius ##
#################     under different exposure conditions:                                          ##
###############How application time, treatment and larval stage affect toxicity. #####################
############### Supplementary Material: R Script #####################################################
############### Author: Victoria Sophie Stoll ########################################################
############### Date: 11-01-2021 #####################################################################

# Prepare data ----
setwd('C:/Users/Sophie/Desktop/M_Sc_Ecotoxicology/Research_Project_Course/Data/R_Analysis')
data = read.csv('C:/Users/Sophie/Desktop/M_Sc_Ecotoxicology/Research_Project_Course/Data/R_Analysis/Sensitivity_data_kumulativ_DAUER_CORR.csv') # load data

##load packages
require(data.table)
require(ggplot2)
require(dplyr)
require(gridExtra)
require(ggpubr)
require(tidyverse)
library(multcomp)
library(knitr) 
library(rstatix)
library(ggpubr)
library(grid)
library(pBrackets)
library(broom)

data <- setDT(data) # set to data.table

data[ , Date := as.Date(Date, format = "%d/%m/%Y")]

data[ , Day := yday(Date) - 61 ] # extract day after start of emergence
data[ , Percent_Ind := Ind/20] # percent individuals of 20 individuals in total
data[ , Percent_Ind_M := M/20] # male
data[ , Percent_Ind_F := F/20] # female
data$Treatment[data$Treatment == '5%'] <- '05%'
data[ , Rpl := as.character(Rpl)] 
data[ , Bti_application := as.factor(Bti_application)] 
data[ , Treatment := as.factor(Treatment)] 
data[ , Larval_stage := as.factor(Larval_stage)] 
data[ , Date := as.character(Date)] 
data[is.na(data)] <- 0  # convert NA to 0

levels(data$Treatment) # check if levels are in the right order
data <- data %>%
  reorder_levels(Treatment, order = c("Ctrl", "05%", "10%", "25%", "50%")) # reorder

levels(data$Bti_application) # check if levels are in the right order
data <- data %>%
  reorder_levels(Bti_application, order = c("none", "d0", "d05", "d10", "d15")) # reorder

levels(data$Larval_stage) # check if levels are in the right order
data <- data %>%
  reorder_levels(Larval_stage, order = c("L1", "Lx")) # reorder

application <- c("none" = "green", "d0" = "dark blue", "d05" = "blue", "d10" = "light blue", "d15" = "grey" )
treatment <- c("Ctrl" = "green", "05%" = "yellow", "10%" = "orange", "25%" = "red", "50%" = "dark red")
 

# Overview ----

##**Description of the (relevant) variables in the data set.** 
  
#Treatment                  = Bti concentration (5, 10, 25, 50% of field rate, Ctrl)   
# Bti_application/Time      = Day of application before test start (=d0) (d0, d05, d10, d15)   
#Larval_stage               = Larval_stage at test start (L1 (earlier), Lx (later))   
#Day                        = Day afer start of emergence   
#Ind/All                    = Number of individuals emerged     
#M                          = Number of males emerged in one day per test unit   
#F                          = Number of females emerged   
#Emergence_all              = Duration from start of emergence till last day of emergence (per test unit)   
#Emergence_M                = Duration from start of emergence till last day of emergence (per test unit) (only males)    
#Emergence_F                = Duration from start of emergence till last day of emergence (per test unit) (only females)   
#Mean_LD_Ind                = Mean larval development (Duration from hatching till emergence)   
#Mean_LD_M                  = Mean larval development (Duration from hatching till emergence) (only males)    
#Mean_LD_F                  = Mean larval development (Duration from hatching till emergence) (only females)  
#Percentage_M               = Percentage of males emerged  
  
###### The investigated endpoints are: ####### 
  #Number of individuals emerged (all, M, F)
  #Duration of emergence (all, M, F)
  # Mean duration of larval development (all, M, F)
  #Percentage of emerged males (gender ratio, only M)

###### The independent variables are: #######
  #Treatment (Bti amount)
  #Bti application (day of application)
  #Larval stage


##### Hypotheses: ########
  #Bti effect mitigates with increasing application time but is still visible at d15
  #The effect is dependent on the treatment (increasing effect with increasing treatment)
  #Earlier larval stages are more sensitive than later larval stages
  #There are different effects in males and females 







# ANOVAs   ----
#### --> note: plot() was used to check assumptions of each ANOVA model but is out-commanded right now because it can be annoying when running the script. 
#### --> remove "#" before plot() to see the plots.

## No. of individuals ----
 
#Treatment
aov_Treatment <- aov(Percent_Ind ~ Treatment, data = data[Day == 17])
summary(aov_Treatment)
####### Check model assumptions:
#plot(aov_Treatment) 
#Residuals vs. fitted values: approx. around 0, horizontal line. Die Säulen-Anordnung ist normal für ANOVA
#Normal QQ values around dashed line, no other patterns.
#even distribution around approx. horizontal line. Between 2.5 and -2.5.
#All data points inside cooks distance (lines) (outliers are no influential points).
TukeyHSD(aov_Treatment)

aov_Treatment_M <- aov(Percent_Ind_M ~ Treatment, data = data[Day == 17])
summary(aov_Treatment_M)
#plot(aov_Treatment_M) 
TukeyHSD(aov_Treatment_M)

aov_Treatment_F <- aov(Percent_Ind_F ~ Treatment, data = data[Day == 17])
summary(aov_Treatment_F)
#plot(aov_Treatment_F) 
TukeyHSD(aov_Treatment_F)

#Time (Bti_application)
aov_Bti_application <- aov(Percent_Ind ~ Bti_application, data = data[Day == 17])
summary(aov_Bti_application)
#plot(aov_Bti_application)
TukeyHSD(aov_Bti_application)

aov_Bti_application_M <- aov(Percent_Ind_M ~ Bti_application, data = data[Day == 17])
summary(aov_Bti_application_M)
#plot(aov_Bti_application_M)
TukeyHSD(aov_Bti_application_M)

aov_Bti_application_F <- aov(Percent_Ind_F ~ Bti_application, data = data[Day == 17])
summary(aov_Bti_application_F)
#plot(aov_Bti_application_F)

#Larval_stage
aov_Larval_stage <- aov(Percent_Ind ~ Larval_stage, data = data[Day == 17])
summary(aov_Larval_stage)
#plot(aov_Larval_stage)
TukeyHSD(aov_Larval_stage)

aov_Larval_stage_M <- aov(Percent_Ind_M ~ Larval_stage, data = data[Day == 17])
summary(aov_Larval_stage_M)
#plot(aov_Larval_stage_M) 
TukeyHSD(aov_Larval_stage_M)

aov_Larval_stage_F <- aov(Percent_Ind_F ~ Larval_stage, data = data[Day == 17])
summary(aov_Larval_stage_F)
#plot(aov_Larval_stage_F)
TukeyHSD(aov_Larval_stage_F)

#Treatment x Bti_application
aov_Treatment_appl <- aov(Percent_Ind ~ Treatment * Bti_application, data = data[Day == 17])
summary(aov_Treatment_appl)
#plot(aov_Treatment_appl) 
TukeyHSD(aov_Treatment_appl)

aov_Treatment_appl_M <- aov(M ~ Treatment * Bti_application, data = data[Day == 17])
summary(aov_Treatment_appl_M)
#plot(aov_Treatment_appl_M)
shapiro_test(residuals(aov_Treatment_appl_M)) 
TukeyHSD(aov_Treatment_appl_M)

aov_Treatment_appl_F <- aov(Percent_Ind_F ~ Treatment * Bti_application, data = data[Day == 17])
summary(aov_Treatment_appl_F)
TukeyHSD(aov_Treatment_appl_F)
#plot(aov_Treatment_appl_F)

#Larval_stage x Bti_application
aov_Larval_stage_appl <- aov(Ind ~ Larval_stage * Bti_application, data = data[Day == 17])
summary(aov_Larval_stage_appl)
#plot(aov_Larval_stage_appl) 
TukeyHSD(aov_Larval_stage_appl)

aov_Larval_stage_appl_M <- aov(M ~ Larval_stage * Bti_application, data = data[Day == 17])
summary(aov_Larval_stage_appl_M)
#plot(aov_Larval_stage_appl_M)
TukeyHSD(aov_Larval_stage_appl_M)

aov_Larval_stage_appl_F <- aov(F ~ Larval_stage * Bti_application, data = data[Day == 17])
summary(aov_Larval_stage_appl_F)
#plot(aov_Larval_stage_appl_F)
TukeyHSD(aov_Larval_stage_appl_F)

#Larval_stage x Treatment
aov_Larval_stage_Treatment <- aov(Percent_Ind ~ Larval_stage * Treatment, data = data[Day == 17])
summary(aov_Larval_stage_Treatment)
#plot(aov_Larval_stage_Treatment)

aov_Larval_stage_Treatment_M <- aov(Percent_Ind_M ~ Larval_stage * Treatment, data = data[Day == 17])
summary(aov_Larval_stage_Treatment_M)
#plot(aov_Larval_stage_Treatment_M)

aov_Larval_stage_Treatment_F <- aov(Percent_Ind_F ~ Larval_stage * Treatment, data = data[Day == 17])
summary(aov_Larval_stage_Treatment_F)
#plot(aov_Larval_stage_Treatment_F)

#Larval_stage x Treatment x Bti_application
aov_Larval_stage_Treatment_Bti_application <- aov(Percent_Ind ~ Larval_stage * Treatment * Bti_application, data = data[Day == 17])
summary(aov_Larval_stage_Treatment_Bti_application)
#plot(aov_Larval_stage_Treatment_Bti_application)

aov_Larval_stage_Treatment_Bti_application_M <- aov(Percent_Ind_M ~ Larval_stage * Treatment * Bti_application, data = data[Day == 17])
summary(aov_Larval_stage_Treatment_Bti_application_M)
#plot(aov_Larval_stage_Treatment_Bti_application_M)

aov_Larval_stage_Treatment_Bti_application_F <- aov(Percent_Ind_F ~ Larval_stage * Treatment * Bti_application, data = data[Day == 17])
summary(aov_Larval_stage_Treatment_Bti_application_F)
#plot(aov_Larval_stage_Treatment_Bti_application_F)
 

## Duration of emergence ----
 
# Treatment
aov_Treatment_Em <- aov(Emergence_all ~ Treatment, data = data[Day == 17 & Ind > 0])
summary(aov_Treatment_Em)
#plot(aov_Treatment_Em)

aov_Treatment_EmM <- aov(Emergenz_M ~ Treatment, data = data[Day == 17 & M > 0])
summary(aov_Treatment_EmM)
#plot(aov_Treatment_EmM) 

aov_Treatment_EmF <- aov(Emergenz_F ~ Treatment, data = data[Day == 17 & F > 0])
summary(aov_Treatment_EmF)
#plot(aov_Treatment_EmF)

#Time (Bti_application)
aov_Bti_application_Em <- aov(Emergence_all ~ Bti_application, data = data[Day == 17 & Ind > 0])
summary(aov_Bti_application_Em)
#plot(aov_Bti_application_Em)
TukeyHSD(aov_Bti_application_Em)

aov_Bti_application_EmM <- aov(Emergenz_M ~ Bti_application, data = data[Day == 17 & M > 0])
summary(aov_Bti_application_EmM)
#plot(aov_Bti_application_EmM) 
TukeyHSD(aov_Bti_application_EmM)

aov_Bti_application_EmF <- aov(Emergenz_F ~ Bti_application, data = data[Day == 17 & F > 0])
summary(aov_Bti_application_EmF)
#plot(aov_Bti_application_EmF)
TukeyHSD(aov_Bti_application_EmF)

#Larval_stage
aov_Larval_stage_Em <- aov(Emergence_all ~ Larval_stage, data = data[Day == 17 & Ind > 0])
summary(aov_Larval_stage_Em)
#plot(aov_Larval_stage_Em)
TukeyHSD(aov_Larval_stage_Em)

aov_Larval_stage_EmM <- aov(Emergenz_M ~ Larval_stage, data = data[Day == 17 & M > 0])
summary(aov_Larval_stage_EmM)
#plot(aov_Larval_stage_EmM) 
TukeyHSD(aov_Larval_stage_EmM )

aov_Larval_stage_EmF <- aov(Emergenz_F ~ Larval_stage, data = data[Day == 17 & F > 0])
summary(aov_Larval_stage_EmF)
#plot(aov_Larval_stage_EmF)
TukeyHSD(aov_Larval_stage_EmF)

#Treatment x Bti_application
aov_Treatment_appl_Em <- aov(Emergence_all ~ Treatment * Bti_application, data = data[Day == 17 & Ind > 0])
summary(aov_Treatment_appl_Em)
#plot(aov_Treatment_appl_Em)

aov_Treatment_appl_EmM <- aov(Emergenz_M ~ Treatment * Bti_application, data = data[Day == 17 & M > 0])
summary(aov_Treatment_appl_EmM)
#plot(aov_Treatment_appl_EmM)

aov_Treatment_appl_EmF <- aov(Emergenz_F ~ Treatment * Bti_application, data = data[Day == 17 & F > 0])
summary(aov_Treatment_appl_EmF)
#plot(aov_Treatment_appl_EmF)

#Larval_stage x Bti_application
aov_Larval_stage_appl_Em <- aov(Emergence_all ~ Larval_stage * Bti_application, data = data[Day == 17 & Ind > 0])
summary(aov_Larval_stage_appl_Em)
#plot(aov_Larval_stage_appl_Em) 

aov_Larval_stage_appl_EmM <- aov(Emergenz_M ~ Larval_stage * Bti_application, data = data[Day == 17 & M > 0])
summary(aov_Larval_stage_appl_EmM)
#plot(aov_Larval_stage_appl_EmM) 

aov_Larval_stage_appl_EmF <- aov(Emergenz_F ~ Larval_stage * Bti_application, data = data[Day == 17 & F > 0])
summary(aov_Larval_stage_appl_EmF)
#plot(aov_Larval_stage_appl_EmF)

#Larval_stage x Treatment
aov_Larval_stage_Treatment_Em <- aov(Emergence_all ~ Larval_stage * Treatment, data = data[Day == 17 & Ind > 0])
summary(aov_Larval_stage_Treatment_Em)
#plot(aov_Larval_stage_Treatment_Em)

aov_Larval_stage_Treatment_EmM <- aov(Emergenz_M ~ Larval_stage * Treatment, data = data[Day == 17 & M > 0])
summary(aov_Larval_stage_Treatment_EmM)
#plot(aov_Larval_stage_Treatment_EmM) 

aov_Larval_stage_Treatment_EmF <- aov(Emergenz_F ~ Larval_stage * Treatment, data = data[Day == 17 & F > 0])
summary(aov_Larval_stage_Treatment_EmF)
#plot(aov_Larval_stage_Treatment_EmF)

#Larval_stage x Treatment x Bti_application
aov_Larval_stage_Treatment_Bti_application_Em <- aov(Emergence_all ~ Larval_stage * Treatment * Bti_application, data = data[Day == 17 & Ind > 0])
summary(aov_Larval_stage_Treatment_Bti_application_Em)
#plot(aov_Larval_stage_Treatment_Bti_application_Em)

aov_Larval_stage_Treatment_Bti_application_EmM <- aov(Emergenz_M ~ Larval_stage * Treatment * Bti_application, data = data[Day == 17 & M > 0])
summary(aov_Larval_stage_Treatment_Bti_application_EmM)
#plot(aov_Larval_stage_Treatment_Bti_application_EmM) 

aov_Larval_stage_Treatment_Bti_application_EmF <- aov(Emergenz_F ~ Larval_stage * Treatment * Bti_application, data = data[Day == 17 & F > 0])
summary(aov_Larval_stage_Treatment_Bti_application_EmF)
#plot(aov_Larval_stage_Treatment_Bti_application_EmF)
TukeyHSD(aov_Larval_stage_Treatment_Bti_application_EmF)
 

## Mean duration of Larval Development  ----
 
# Treatment
aov_Treatment_LD <- aov(Mean_LD_Ind ~ Treatment, data = data[Day == 17 & Ind > 0])
summary(aov_Treatment_LD)
##plot(aov_Treatment_LD)

aov_Treatment_LDM <- aov(Mean_LD_M ~ Treatment, data = data[Day == 17 & M > 0])
summary(aov_Treatment_LDM)
##plot(aov_Treatment_LDM)

aov_Treatment_LDF <- aov(Mean_LD_F ~ Treatment, data = data[Day == 17 & F > 0])
summary(aov_Treatment_LDF)
##plot(aov_Treatment_LDF)

#Time (Bti_application)
aov_Bti_application_LD <- aov(Mean_LD_Ind ~ Bti_application, data = data[Day == 17 & Ind > 0])
summary(aov_Bti_application_LD)
#plot(aov_Bti_application_LD)

aov_Bti_application_LDM <- aov(Mean_LD_M ~ Bti_application, data = data[Day == 17 & M > 0])
summary(aov_Bti_application_LDM)
#plot(aov_Bti_application_LDM)

aov_Bti_application_LDF <- aov(Mean_LD_F ~ Bti_application, data = data[Day == 17 & F > 0])
summary(aov_Bti_application_LDF)
#plot(aov_Bti_application_LDF)

#Larval_stage
aov_Larval_stage_LD <- aov(Mean_LD_Ind ~ Larval_stage, data = data[Day == 17 & Ind > 0])
summary(aov_Larval_stage_LD)
#plot(aov_Larval_stage_LD)
TukeyHSD(aov_Larval_stage_LD)

aov_Larval_stage_LDM <- aov(Mean_LD_M ~ Larval_stage, data = data[Day == 17 & M > 0])
summary(aov_Larval_stage_LDM)
#plot(aov_Larval_stage_LDM)
TukeyHSD(aov_Larval_stage_LDM)

aov_Larval_stage_LDF <- aov(Mean_LD_F ~ Larval_stage, data = data[Day == 17 & F > 0])
summary(aov_Larval_stage_LDF)
#plot(aov_Larval_stage_LDF)
TukeyHSD(aov_Larval_stage_LDF)

#Treatment x Bti_application
aov_Treatment_appl_LD <- aov(Mean_LD_Ind ~ Treatment * Bti_application, data = data[Day == 17 & Ind > 0])
summary(aov_Treatment_appl_LD)
#plot(aov_Treatment_appl_LD) 

aov_Treatment_appl_LDM <- aov(Mean_LD_M ~ Treatment * Bti_application, data = data[Day == 17 & M > 0])
summary(aov_Treatment_appl_LDM)
#plot(aov_Treatment_appl_LDM) 

aov_Treatment_appl_LDF <- aov(Mean_LD_F ~ Treatment * Bti_application, data = data[Day == 17 & F > 0])
summary(aov_Treatment_appl_LDF)
#plot(aov_Treatment_appl_LDF)


#Larval_stage x Bti_application
aov_Larval_stage_appl_LD <- aov(Mean_LD_Ind ~ Larval_stage * Bti_application, data = data[Day == 17 & Ind > 0])
summary(aov_Larval_stage_appl_LD)
#plot(aov_Larval_stage_appl_LD)

aov_Larval_stage_appl_LDM <- aov(Mean_LD_M ~ Larval_stage * Bti_application, data = data[Day == 17 & M > 0])
summary(aov_Larval_stage_appl_LDM)
#plot(aov_Larval_stage_appl_LDM) 

aov_Larval_stage_appl_LDF <- aov(Mean_LD_F ~ Larval_stage * Bti_application, data = data[Day == 17 & F > 0])
summary(aov_Larval_stage_appl_LDF)
#plot(aov_Larval_stage_appl_LDF)

#Larval_stage x Treatment
aov_Larval_stage_Treatment_LD <- aov(Mean_LD_Ind ~ Larval_stage * Treatment, data = data[Day == 17 & Ind > 0])
summary(aov_Larval_stage_Treatment_LD)
#plot(aov_Larval_stage_Treatment_LD)

aov_Larval_stage_Treatment_LDM <- aov(Mean_LD_M ~ Larval_stage * Treatment, data = data[Day == 17 & M > 0])
summary(aov_Larval_stage_Treatment_LDM)
#plot(aov_Larval_stage_Treatment_LDM)

aov_Larval_stage_Treatment_LDF <- aov(Mean_LD_F ~ Larval_stage * Treatment, data = data[Day == 17 & F > 0])
summary(aov_Larval_stage_Treatment_LDF)
#plot(aov_Larval_stage_Treatment_LDF)

#Larval_stage x Treatment x Bti_application
aov_Larval_stage_Treatment_Bti_application_LD <- aov(Mean_LD_Ind ~ Larval_stage * Treatment * Bti_application, data = data[Day == 17 & Ind > 0])
summary(aov_Larval_stage_Treatment_Bti_application_LD)
#plot(aov_Larval_stage_Treatment_Bti_application_LD)

aov_Larval_stage_Treatment_Bti_application_LDM <- aov(Mean_LD_M ~ Larval_stage * Treatment * Bti_application, data = data[Day == 17 & M > 0])
summary(aov_Larval_stage_Treatment_Bti_application_LDM)
#plot(aov_Larval_stage_Treatment_Bti_application_LDM)

aov_Larval_stage_Treatment_Bti_application_LDF <- aov(Mean_LD_F ~ Larval_stage * Treatment * Bti_application, data = data[Day == 17 & F > 0])
summary(aov_Larval_stage_Treatment_Bti_application_LDF)
#plot(aov_Larval_stage_Treatment_Bti_application_LDF)
 

## Percentage of emerged males (gender ratio) ----

 
selectday = 17 # change to see results over time

# Treatment
aov_Treatment_G <- aov(GR_M ~ Treatment, data = data[Day == selectday & Ind > 0])
summary(aov_Treatment_G)
#plot(aov_Treatment_G)

# Time
aov_Time_G <- aov(GR_M ~ Bti_application, data = data[Day == selectday & Ind > 0])
summary(aov_Time_G)
#plot(aov_Time_G)

# Larval_stage
aov_LS_G <- aov(GR_M ~ Larval_stage, data = data[Day == selectday & Ind > 0])
summary(aov_LS_G)
#plot(aov_LS_G)

# Time x Treatment
aov_TT_G <- aov(GR_M ~ Treatment * Bti_application, data = data[Day == selectday  & Ind > 0])
summary(aov_TT_G)
TukeyHSD(aov_TT_G)
#plot(aov_TT_G)

# Larval_stage x Time
aov_LST_G <- aov(GR_M ~ Larval_stage * Bti_application, data = data[Day == selectday & Ind > 0])
summary(aov_LST_G)
#plot(aov_LST_G)

# Larval_stage x Treatment
aov_LSTr_G <- aov(GR_M ~ Larval_stage * Treatment, data = data[Day == selectday & Ind > 0])
summary(aov_LSTr_G)
#plot(aov_LSTr_G)

# Larval_stage x Treatment x Bti_application
aov_LSTrT_G <- aov(GR_M ~ Larval_stage * Treatment * Bti_application, data = data[Day == selectday & Ind > 0])
summary(aov_LSTrT_G)
#plot(aov_LSTrT_G)

 

# Relationship btw. emerged individuals and duration of emergence ----
###The pattern of duration of emergence grouped by application days is similar to the pattern of number of individuals grouped by application days (see RPC report for data visualizations). Therefore, it was further investigated if the duration of emergence is predicted by the number of individuals and if this a stochastic or biological effect.

## Linear regression ----
lm <- lm(Emergence_all ~ Ind, data[Day==17 & Ind > 0])
summary(lm)

ggscatter(
  data[Day == 17 & Ind > 0], x = "Ind", y = "Emergence_all", add = "reg.line") +
  geom_point(aes(colour = Bti_application), size = 3) +
  scale_colour_manual(values = application) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))
  )
 

## ANCOVA  ----
##- ...to control for Bti_application

### Check Assumptions:

#1) linear relationship between covariate and response variable
##response variable = Emergence duration (Emergence_all)
##covariate = Number of Individuals (Ind)

ggscatter(
  data[Day == 17 & Ind > 0], x = "Ind", y = "Emergence_all",
  color = "Bti_application", add = "reg.line", size = 3,
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Bti_application)
  )

lm <- lm(Emergence_all ~ Ind, data[Day==17])
summary(lm)


#2) Homogeneity of regression slopes
 
data[Day == 17 & Ind > 0] %>% anova_test(Emergence_all ~ Bti_application * Ind)
 
data[Day == 17 & Ind > 0 & Bti_application == "none" | Day == 17 & Ind > 0 & Bti_application != "d0"] %>% anova_test(Emergence_all ~ Bti_application * Ind)
#The slopes are not parallel. However, the interaction term is not significant, which means that the regression slopes are homogeneous.

#Fit the model:
# the covariate goes first!!
model.ancova <- lm(Emergence_all ~ Ind + Bti_application, data = data[Day == 17 & Ind > 0])
summary(model.ancova)
plot(model.ancova) # quickly check assumptions
model.metrics <- augment(model.ancova)

#Check assumptions of the model:
#3) Normality of residuals
hist(model.metrics$.resid)
shapiro_test(model.metrics$.resid)
#Normality of residuals is ok. Shapiro test is significant. ANCOVA is relatively robust against violation of normality.

#4) homogeneous variance
model.metrics %>% levene_test(.resid ~ Bti_application)
library(car)
leveneTest(Emergence_all ~ Bti_application, data = data[Day == 17 & Ind > 0])
#Variance is homogeneous.

#5) No outliers
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()
#There are two outliers in the data set (only slightly above 3 (std. err.)). We take them as they are because we know about them and decided to leave them in the data set. They are inside cooks distance (see above), so this should be ok (no influential point).

#Computation:
res.aov <- data[Day == 17 & Ind > 0] %>% anova_test(Emergence_all ~ Ind + Bti_application)
get_anova_table(res.aov)
summary(aov(Emergence_all ~ Ind + Bti_application, data = data[Day == 17 & Ind > 0])) #another way to double check
car::Anova(aov(Emergence_all ~ Ind + Bti_application, data = data[Day == 17 & Ind > 0]), type = "III") #another way to double check
#After adjustment for the number of individuals, there is no statistically significant difference in the duration of emergence between the different Bti_applications.

## Permutation  ----
#Is this correlation between the number of individuals and duration of emergence due to a stochastic or biological effect? We can test this by resampling the data set by using permuation.
#In this case, the resampling is a bit tricky. The single sampled individuals need to shuffled and then the duration of emergence and total number of individuals needs to be calculated again for each sample. This is done manually in the next steps because there is no function in R doing this.
 
#extract information from data set for a better overview
data.perm <- data %>%
  select(Name, Rpl, Bti_application, Ind_Single, Day, Larval_stage)

data.perm.L1 <- data.perm[Larval_stage == "L1" & Day > 2] # general emergence start for L1 was later than for Lx
data.perm.Lx <- data.perm[Larval_stage == "Lx" & Day < 16] # general emergence end for Lx was later than for Lx

data.perm.L1 <- data.perm.L1 %>% arrange(Name) # order by Name 
data.perm.Lx <- data.perm.Lx %>% arrange(Name) # order by Name

set.seed(1979) # for reproducability of results


n1 <- length(data.perm.L1$Name) # the number of observations to sample
nx <- length(data.perm.Lx$Name) # the number of observations to sample

P <- 100 # the number of permutation samples to take

Ind_L1 <- data.perm.L1$Ind_Single # the variable we will resample from (single sampled (not cumulative) individuals!)
Ind_Lx <- data.perm.Lx$Ind_Single # the variable we will resample from (single sampled (not cumulative) individuals!)

PermSamplesL1 <- matrix(0, nrow=n1, ncol=P) # initialize a matrix to store the permutation data
PermSamplesLx <- matrix(0, nrow=nx, ncol=P) # initialize a matrix to store the permutation data
# each column is a permutation sample of data

# get those permutation samples
for(i in 1:P){
  PermSamplesL1[,i] <- sample(Ind_L1, size= n1, replace=FALSE)
}
for(i in 1:P){
  PermSamplesLx[,i] <- sample(Ind_Lx, size= nx, replace=FALSE)
}

PermSamplesL1 <- as.data.table(PermSamplesL1)
PermSamplesLx <- as.data.table(PermSamplesLx)

# combine data 
perm.L1v1 <- cbind(data.perm.L1, PermSamplesL1$V1)
perm.Lxv1 <- cbind(data.perm.Lx, PermSamplesLx$V1)

perm.L1v2 <- cbind(data.perm.L1, PermSamplesL1$V2)
perm.Lxv2 <- cbind(data.perm.Lx, PermSamplesLx$V2)

perm.L1v3 <- cbind(data.perm.L1, PermSamplesL1$V3)
perm.Lxv3 <- cbind(data.perm.Lx, PermSamplesLx$V3)

perm.L1v4 <- cbind(data.perm.L1, PermSamplesL1$V4)
perm.Lxv4 <- cbind(data.perm.Lx, PermSamplesLx$V4)

perm.L1v5 <- cbind(data.perm.L1, PermSamplesL1$V5)
perm.Lxv5 <- cbind(data.perm.Lx, PermSamplesLx$V5)

perm.L1v6 <- cbind(data.perm.L1, PermSamplesL1$V6)
perm.Lxv6 <- cbind(data.perm.Lx, PermSamplesLx$V6)

perm.L1v7 <- cbind(data.perm.L1, PermSamplesL1$V7)
perm.Lxv7 <- cbind(data.perm.Lx, PermSamplesLx$V7)

perm.L1v8 <- cbind(data.perm.L1, PermSamplesL1$V8)
perm.Lxv8 <- cbind(data.perm.Lx, PermSamplesLx$V8)

perm.L1v9 <- cbind(data.perm.L1, PermSamplesL1$V9)
perm.Lxv9 <- cbind(data.perm.Lx, PermSamplesLx$V9)

perm.L1v10 <- cbind(data.perm.L1, PermSamplesL1$V10)
perm.Lxv10 <- cbind(data.perm.Lx, PermSamplesLx$V10)

perm.L1v11 <- cbind(data.perm.L1, PermSamplesL1$V11)
perm.Lxv11 <- cbind(data.perm.Lx, PermSamplesLx$V11)

perm.L1v12 <- cbind(data.perm.L1, PermSamplesL1$V12)
perm.Lxv12 <- cbind(data.perm.Lx, PermSamplesLx$V12)

perm.L1v13 <- cbind(data.perm.L1, PermSamplesL1$V13)
perm.Lxv13 <- cbind(data.perm.Lx, PermSamplesLx$V13)

perm.L1v14 <- cbind(data.perm.L1, PermSamplesL1$V14)
perm.Lxv14 <- cbind(data.perm.Lx, PermSamplesLx$V14)

perm.L1v15 <- cbind(data.perm.L1, PermSamplesL1$V15)
perm.Lxv15 <- cbind(data.perm.Lx, PermSamplesLx$V15)

perm.L1v16 <- cbind(data.perm.L1, PermSamplesL1$V16)
perm.Lxv16 <- cbind(data.perm.Lx, PermSamplesLx$V16)

perm.L1v17 <- cbind(data.perm.L1, PermSamplesL1$V17)
perm.Lxv17 <- cbind(data.perm.Lx, PermSamplesLx$V17)

perm.L1v18 <- cbind(data.perm.L1, PermSamplesL1$V18)
perm.Lxv18 <- cbind(data.perm.Lx, PermSamplesLx$V18)

perm.L1v19 <- cbind(data.perm.L1, PermSamplesL1$V19)
perm.Lxv19 <- cbind(data.perm.Lx, PermSamplesLx$V19)

perm.L1v20 <- cbind(data.perm.L1, PermSamplesL1$V20)
perm.Lxv20 <- cbind(data.perm.Lx, PermSamplesLx$V20)

perm.L1v21 <- cbind(data.perm.L1, PermSamplesL1$V21)
perm.Lxv21 <- cbind(data.perm.Lx, PermSamplesLx$V21)

perm.L1v22 <- cbind(data.perm.L1, PermSamplesL1$V22)
perm.Lxv22 <- cbind(data.perm.Lx, PermSamplesLx$V22)

perm.L1v23 <- cbind(data.perm.L1, PermSamplesL1$V23)
perm.Lxv23 <- cbind(data.perm.Lx, PermSamplesLx$V23)

perm.L1v24 <- cbind(data.perm.L1, PermSamplesL1$V24)
perm.Lxv24 <- cbind(data.perm.Lx, PermSamplesLx$V24)

perm.L1v25 <- cbind(data.perm.L1, PermSamplesL1$V25)
perm.Lxv25 <- cbind(data.perm.Lx, PermSamplesLx$V25)

perm.L1v26 <- cbind(data.perm.L1, PermSamplesL1$V26)
perm.Lxv26 <- cbind(data.perm.Lx, PermSamplesLx$V26)

perm.L1v27 <- cbind(data.perm.L1, PermSamplesL1$V27)
perm.Lxv27 <- cbind(data.perm.Lx, PermSamplesLx$V27)

perm.L1v28 <- cbind(data.perm.L1, PermSamplesL1$V28)
perm.Lxv28 <- cbind(data.perm.Lx, PermSamplesLx$V28)

perm.L1v29 <- cbind(data.perm.L1, PermSamplesL1$V29)
perm.Lxv29 <- cbind(data.perm.Lx, PermSamplesLx$V29)

perm.L1v30 <- cbind(data.perm.L1, PermSamplesL1$V30)
perm.Lxv30 <- cbind(data.perm.Lx, PermSamplesLx$V30)

perm.L1v31 <- cbind(data.perm.L1, PermSamplesL1$V31)
perm.Lxv31 <- cbind(data.perm.Lx, PermSamplesLx$V31)

perm.L1v32 <- cbind(data.perm.L1, PermSamplesL1$V32)
perm.Lxv32 <- cbind(data.perm.Lx, PermSamplesLx$V32)

perm.L1v33 <- cbind(data.perm.L1, PermSamplesL1$V33)
perm.Lxv33 <- cbind(data.perm.Lx, PermSamplesLx$V33)

perm.L1v34 <- cbind(data.perm.L1, PermSamplesL1$V34)
perm.Lxv34 <- cbind(data.perm.Lx, PermSamplesLx$V34)

perm.L1v35 <- cbind(data.perm.L1, PermSamplesL1$V35)
perm.Lxv35 <- cbind(data.perm.Lx, PermSamplesLx$V15)

perm.L1v36 <- cbind(data.perm.L1, PermSamplesL1$V36)
perm.Lxv36 <- cbind(data.perm.Lx, PermSamplesLx$V16)

perm.L1v37 <- cbind(data.perm.L1, PermSamplesL1$V37)
perm.Lxv37 <- cbind(data.perm.Lx, PermSamplesLx$V37)

perm.L1v38 <- cbind(data.perm.L1, PermSamplesL1$V38)
perm.Lxv38 <- cbind(data.perm.Lx, PermSamplesLx$V38)

perm.L1v39 <- cbind(data.perm.L1, PermSamplesL1$V39)
perm.Lxv39 <- cbind(data.perm.Lx, PermSamplesLx$V39)

perm.L1v40 <- cbind(data.perm.L1, PermSamplesL1$V40)
perm.Lxv40 <- cbind(data.perm.Lx, PermSamplesLx$V40)

perm.L1v41 <- cbind(data.perm.L1, PermSamplesL1$V41)
perm.Lxv41 <- cbind(data.perm.Lx, PermSamplesLx$V41)

perm.L1v42 <- cbind(data.perm.L1, PermSamplesL1$V42)
perm.Lxv42 <- cbind(data.perm.Lx, PermSamplesLx$V42)

perm.L1v43 <- cbind(data.perm.L1, PermSamplesL1$V43)
perm.Lxv43 <- cbind(data.perm.Lx, PermSamplesLx$V43)

perm.L1v44 <- cbind(data.perm.L1, PermSamplesL1$V44)
perm.Lxv44 <- cbind(data.perm.Lx, PermSamplesLx$V44)

perm.L1v45 <- cbind(data.perm.L1, PermSamplesL1$V45)
perm.Lxv45 <- cbind(data.perm.Lx, PermSamplesLx$V45)

perm.L1v46 <- cbind(data.perm.L1, PermSamplesL1$V46)
perm.Lxv46 <- cbind(data.perm.Lx, PermSamplesLx$V46)

perm.L1v47 <- cbind(data.perm.L1, PermSamplesL1$V47)
perm.Lxv47 <- cbind(data.perm.Lx, PermSamplesLx$V47)

perm.L1v48 <- cbind(data.perm.L1, PermSamplesL1$V48)
perm.Lxv48 <- cbind(data.perm.Lx, PermSamplesLx$V48)

perm.L1v49 <- cbind(data.perm.L1, PermSamplesL1$V49)
perm.Lxv49 <- cbind(data.perm.Lx, PermSamplesLx$V49)

perm.L1v50 <- cbind(data.perm.L1, PermSamplesL1$V50)
perm.Lxv50 <- cbind(data.perm.Lx, PermSamplesLx$V50)

perm.L1v51 <- cbind(data.perm.L1, PermSamplesL1$V51)
perm.Lxv51 <- cbind(data.perm.Lx, PermSamplesLx$V51)

perm.L1v52 <- cbind(data.perm.L1, PermSamplesL1$V52)
perm.Lxv52 <- cbind(data.perm.Lx, PermSamplesLx$V52)

perm.L1v53 <- cbind(data.perm.L1, PermSamplesL1$V53)
perm.Lxv53 <- cbind(data.perm.Lx, PermSamplesLx$V53)

perm.L1v54 <- cbind(data.perm.L1, PermSamplesL1$V54)
perm.Lxv54 <- cbind(data.perm.Lx, PermSamplesLx$V54)

perm.L1v55 <- cbind(data.perm.L1, PermSamplesL1$V55)
perm.Lxv55 <- cbind(data.perm.Lx, PermSamplesLx$V55)

perm.L1v56 <- cbind(data.perm.L1, PermSamplesL1$V56)
perm.Lxv56 <- cbind(data.perm.Lx, PermSamplesLx$V56)

perm.L1v57 <- cbind(data.perm.L1, PermSamplesL1$V57)
perm.Lxv57 <- cbind(data.perm.Lx, PermSamplesLx$V57)

perm.L1v58 <- cbind(data.perm.L1, PermSamplesL1$V58)
perm.Lxv58 <- cbind(data.perm.Lx, PermSamplesLx$V58)

perm.L1v59 <- cbind(data.perm.L1, PermSamplesL1$V59)
perm.Lxv59 <- cbind(data.perm.Lx, PermSamplesLx$V59)

perm.L1v60 <- cbind(data.perm.L1, PermSamplesL1$V60)
perm.Lxv60 <- cbind(data.perm.Lx, PermSamplesLx$V60)

perm.L1v61 <- cbind(data.perm.L1, PermSamplesL1$V61)
perm.Lxv61 <- cbind(data.perm.Lx, PermSamplesLx$V61)

perm.L1v62 <- cbind(data.perm.L1, PermSamplesL1$V62)
perm.Lxv62 <- cbind(data.perm.Lx, PermSamplesLx$V62)

perm.L1v63 <- cbind(data.perm.L1, PermSamplesL1$V63)
perm.Lxv63 <- cbind(data.perm.Lx, PermSamplesLx$V63)

perm.L1v64 <- cbind(data.perm.L1, PermSamplesL1$V64)
perm.Lxv64 <- cbind(data.perm.Lx, PermSamplesLx$V64)

perm.L1v65 <- cbind(data.perm.L1, PermSamplesL1$V65)
perm.Lxv65 <- cbind(data.perm.Lx, PermSamplesLx$V65)

perm.L1v66 <- cbind(data.perm.L1, PermSamplesL1$V66)
perm.Lxv66 <- cbind(data.perm.Lx, PermSamplesLx$V66)

perm.L1v67 <- cbind(data.perm.L1, PermSamplesL1$V67)
perm.Lxv67 <- cbind(data.perm.Lx, PermSamplesLx$V67)

perm.L1v68 <- cbind(data.perm.L1, PermSamplesL1$V68)
perm.Lxv68 <- cbind(data.perm.Lx, PermSamplesLx$V68)

perm.L1v69 <- cbind(data.perm.L1, PermSamplesL1$V69)
perm.Lxv69 <- cbind(data.perm.Lx, PermSamplesLx$V69)

perm.L1v70 <- cbind(data.perm.L1, PermSamplesL1$V70)
perm.Lxv70 <- cbind(data.perm.Lx, PermSamplesLx$V70)

perm.L1v71 <- cbind(data.perm.L1, PermSamplesL1$V71)
perm.Lxv71 <- cbind(data.perm.Lx, PermSamplesLx$V71)

perm.L1v72 <- cbind(data.perm.L1, PermSamplesL1$V72)
perm.Lxv72 <- cbind(data.perm.Lx, PermSamplesLx$V72)

perm.L1v73 <- cbind(data.perm.L1, PermSamplesL1$V73)
perm.Lxv73 <- cbind(data.perm.Lx, PermSamplesLx$V73)

perm.L1v74 <- cbind(data.perm.L1, PermSamplesL1$V74)
perm.Lxv74 <- cbind(data.perm.Lx, PermSamplesLx$V74)

perm.L1v75 <- cbind(data.perm.L1, PermSamplesL1$V75)
perm.Lxv75 <- cbind(data.perm.Lx, PermSamplesLx$V75)

perm.L1v76 <- cbind(data.perm.L1, PermSamplesL1$V76)
perm.Lxv76 <- cbind(data.perm.Lx, PermSamplesLx$V76)

perm.L1v77 <- cbind(data.perm.L1, PermSamplesL1$V77)
perm.Lxv77 <- cbind(data.perm.Lx, PermSamplesLx$V77)

perm.L1v78 <- cbind(data.perm.L1, PermSamplesL1$V78)
perm.Lxv78 <- cbind(data.perm.Lx, PermSamplesLx$V78)

perm.L1v79 <- cbind(data.perm.L1, PermSamplesL1$V79)
perm.Lxv79 <- cbind(data.perm.Lx, PermSamplesLx$V79)

perm.L1v80 <- cbind(data.perm.L1, PermSamplesL1$V80)
perm.Lxv80 <- cbind(data.perm.Lx, PermSamplesLx$V80)

perm.L1v81 <- cbind(data.perm.L1, PermSamplesL1$V81)
perm.Lxv81 <- cbind(data.perm.Lx, PermSamplesLx$V81)

perm.L1v82 <- cbind(data.perm.L1, PermSamplesL1$V82)
perm.Lxv82 <- cbind(data.perm.Lx, PermSamplesLx$V82)

perm.L1v83 <- cbind(data.perm.L1, PermSamplesL1$V83)
perm.Lxv83 <- cbind(data.perm.Lx, PermSamplesLx$V83)

perm.L1v84 <- cbind(data.perm.L1, PermSamplesL1$V84)
perm.Lxv84 <- cbind(data.perm.Lx, PermSamplesLx$V84)

perm.L1v85 <- cbind(data.perm.L1, PermSamplesL1$V85)
perm.Lxv85 <- cbind(data.perm.Lx, PermSamplesLx$V85)

perm.L1v86 <- cbind(data.perm.L1, PermSamplesL1$V86)
perm.Lxv86 <- cbind(data.perm.Lx, PermSamplesLx$V86)

perm.L1v87 <- cbind(data.perm.L1, PermSamplesL1$V87)
perm.Lxv87 <- cbind(data.perm.Lx, PermSamplesLx$V87)

perm.L1v88 <- cbind(data.perm.L1, PermSamplesL1$V88)
perm.Lxv88 <- cbind(data.perm.Lx, PermSamplesLx$V88)

perm.L1v89 <- cbind(data.perm.L1, PermSamplesL1$V89)
perm.Lxv89 <- cbind(data.perm.Lx, PermSamplesLx$V89)

perm.L1v90 <- cbind(data.perm.L1, PermSamplesL1$V90)
perm.Lxv90 <- cbind(data.perm.Lx, PermSamplesLx$V90)

perm.L1v91 <- cbind(data.perm.L1, PermSamplesL1$V91)
perm.Lxv91 <- cbind(data.perm.Lx, PermSamplesLx$V91)

perm.L1v92 <- cbind(data.perm.L1, PermSamplesL1$V92)
perm.Lxv92 <- cbind(data.perm.Lx, PermSamplesLx$V92)

perm.L1v93 <- cbind(data.perm.L1, PermSamplesL1$V93)
perm.Lxv93 <- cbind(data.perm.Lx, PermSamplesLx$V93)

perm.L1v94 <- cbind(data.perm.L1, PermSamplesL1$V94)
perm.Lxv94 <- cbind(data.perm.Lx, PermSamplesLx$V94)

perm.L1v95 <- cbind(data.perm.L1, PermSamplesL1$V95)
perm.Lxv95 <- cbind(data.perm.Lx, PermSamplesLx$V95)

perm.L1v96 <- cbind(data.perm.L1, PermSamplesL1$V96)
perm.Lxv96 <- cbind(data.perm.Lx, PermSamplesLx$V96)

perm.L1v97 <- cbind(data.perm.L1, PermSamplesL1$V97)
perm.Lxv97 <- cbind(data.perm.Lx, PermSamplesLx$V97)

perm.L1v98 <- cbind(data.perm.L1, PermSamplesL1$V98)
perm.Lxv98 <- cbind(data.perm.Lx, PermSamplesLx$V98)

perm.L1v99 <- cbind(data.perm.L1, PermSamplesL1$V99)
perm.Lxv99 <- cbind(data.perm.Lx, PermSamplesLx$V99)

perm.L1v100 <- cbind(data.perm.L1, PermSamplesL1$V100)
perm.Lxv100 <- cbind(data.perm.Lx, PermSamplesLx$V100)

allframes <- list(perm.L1v1, perm.Lxv1, 
                  perm.L1v2, perm.Lxv2, 
                  perm.L1v3, perm.Lxv3, 
                  perm.L1v4, perm.Lxv4, 
                  perm.L1v5, perm.Lxv5,
                  perm.L1v6, perm.Lxv6, 
                  perm.L1v7, perm.Lxv7,
                  perm.L1v8, perm.Lxv8,
                  perm.L1v9, perm.Lxv9,
                  perm.L1v10, perm.Lxv10,
                  perm.L1v11, perm.Lxv11,
                  perm.L1v12, perm.Lxv12,
                  perm.L1v13, perm.Lxv13,
                  perm.L1v14, perm.Lxv14,
                  perm.L1v15, perm.Lxv15,
                  perm.L1v16, perm.Lxv16,
                  perm.L1v17, perm.Lxv17,
                  perm.L1v18, perm.Lxv18,
                  perm.L1v19, perm.Lxv19,
                  perm.L1v20, perm.Lxv20,
                  perm.L1v21, perm.Lxv21, 
                  perm.L1v22, perm.Lxv22, 
                  perm.L1v23, perm.Lxv23, 
                  perm.L1v24, perm.Lxv24, 
                  perm.L1v25, perm.Lxv25,
                  perm.L1v26, perm.Lxv26, 
                  perm.L1v27, perm.Lxv27,
                  perm.L1v28, perm.Lxv28,
                  perm.L1v29, perm.Lxv29,
                  perm.L1v30, perm.Lxv30,
                  perm.L1v31, perm.Lxv31,
                  perm.L1v32, perm.Lxv32,
                  perm.L1v33, perm.Lxv33,
                  perm.L1v34, perm.Lxv34,
                  perm.L1v35, perm.Lxv35,
                  perm.L1v36, perm.Lxv36,
                  perm.L1v37, perm.Lxv37,
                  perm.L1v38, perm.Lxv38,
                  perm.L1v39, perm.Lxv39,
                  perm.L1v40, perm.Lxv40,
                  perm.L1v41, perm.Lxv41, 
                  perm.L1v42, perm.Lxv42, 
                  perm.L1v43, perm.Lxv43, 
                  perm.L1v44, perm.Lxv44, 
                  perm.L1v45, perm.Lxv45,
                  perm.L1v46, perm.Lxv46, 
                  perm.L1v47, perm.Lxv47,
                  perm.L1v48, perm.Lxv48,
                  perm.L1v49, perm.Lxv49,
                  perm.L1v50, perm.Lxv50,
                  perm.L1v51, perm.Lxv51,
                  perm.L1v52, perm.Lxv52,
                  perm.L1v53, perm.Lxv53,
                  perm.L1v54, perm.Lxv54,
                  perm.L1v55, perm.Lxv55,
                  perm.L1v56, perm.Lxv56,
                  perm.L1v57, perm.Lxv57,
                  perm.L1v58, perm.Lxv58,
                  perm.L1v59, perm.Lxv59,
                  perm.L1v60, perm.Lxv60,
                  perm.L1v61, perm.Lxv61, 
                  perm.L1v62, perm.Lxv62, 
                  perm.L1v63, perm.Lxv63, 
                  perm.L1v64, perm.Lxv64, 
                  perm.L1v65, perm.Lxv65,
                  perm.L1v66, perm.Lxv66, 
                  perm.L1v67, perm.Lxv67,
                  perm.L1v68, perm.Lxv68,
                  perm.L1v69, perm.Lxv69,
                  perm.L1v70, perm.Lxv70,
                  perm.L1v71, perm.Lxv71,
                  perm.L1v72, perm.Lxv72,
                  perm.L1v73, perm.Lxv73,
                  perm.L1v74, perm.Lxv74,
                  perm.L1v75, perm.Lxv75,
                  perm.L1v76, perm.Lxv76,
                  perm.L1v77, perm.Lxv77,
                  perm.L1v78, perm.Lxv78,
                  perm.L1v79, perm.Lxv79,
                  perm.L1v80, perm.Lxv80,
                  perm.L1v81, perm.Lxv81, 
                  perm.L1v82, perm.Lxv82, 
                  perm.L1v83, perm.Lxv83, 
                  perm.L1v84, perm.Lxv84, 
                  perm.L1v85, perm.Lxv85,
                  perm.L1v86, perm.Lxv86, 
                  perm.L1v87, perm.Lxv87,
                  perm.L1v88, perm.Lxv88,
                  perm.L1v89, perm.Lxv89,
                  perm.L1v90, perm.Lxv90,
                  perm.L1v91, perm.Lxv91,
                  perm.L1v92, perm.Lxv92,
                  perm.L1v93, perm.Lxv93,
                  perm.L1v94, perm.Lxv94,
                  perm.L1v95, perm.Lxv95,
                  perm.L1v96, perm.Lxv96,
                  perm.L1v97, perm.Lxv97,
                  perm.L1v98, perm.Lxv98,
                  perm.L1v99, perm.Lxv99,
                  perm.L1v100, perm.Lxv100)
# edit all frames
framies <- lapply(allframes, function(x) {
  x %>%
    mutate(Name = as.factor(Name)) %>%
    group_by(Name) %>%
    mutate(Ind = sum(V2)) %>%
    filter(V2 > 0)%>%
    mutate(max = max(Day)) %>%
    mutate(min = min(Day)) %>%
    mutate(Emergence_Ind = max - min) %>%
    distinct(Name, .keep_all= TRUE)
  
} )

onebigframe <- bind_rows(framies, .id = "column_label")
onebigframe <- as.data.table(onebigframe)
onebigframe %>% anova_test(Emergence_Ind ~ column_label * Ind)
onebigframe %>% anova_test(Emergence_Ind ~ column_label * Ind)

# linear models per data frame
fits <- lapply(framies, function(x) {
  lm(formula = Emergence_Ind ~ Ind, data = x )
} )


# extract information from all models
r.squared <- sapply(fits, function(x) summary(x)$r.squared )
p.value <- sapply(fits, function(x) summary(x)$coefficients[2,4] )
slope <- sapply(fits, function(x) summary(x)$coefficients[2,1] )


# compare with original 
lm <- lm(Emergence_all ~ Ind, data[Day==17 & Ind > 0])
summary(lm)
# make up original data
datalm <- data[Day == 17] %>%
  select(Name, Rpl, Bti_application, Ind, Day, Larval_stage, Emergence_all) %>%
  rename(Emergence_Ind = Emergence_all)
datalm$group <- "A"
# make up permuted data
onebigframelm <- onebigframe %>%
  select(Name, Rpl, Bti_application, Ind, Day, Larval_stage, Emergence_Ind)
onebigframelm$group <- "B"
# merge original and permuted data
allframeslm <- rbind(datalm, onebigframelm)


## Compare permuted/original data  ----
 
# original
modelA <- lm(Emergence_Ind ~ Ind, allframeslm[group == 'A' & Ind > 0])
summary(modelA)
# permuted
modelB <- lm(Emergence_Ind ~ Ind, allframeslm[group == 'B' & Ind > 0])
summary(modelB)
# are the slopes significantly different?
anova <- aov(Emergence_Ind ~ group * Ind, allframeslm[Ind > 0])
summary(anova) # significant
library(sjstats)
effectsize::eta_squared(anova, partial = FALSE) # always report effect sizes

ggscatter(
  allframeslm[Ind > 0], x = "Ind", y = "Emergence_Ind",
  color = "group", add = "reg.line", size = "group")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = group), 
    label.x = c(25,25), label.y = c(3.5,5), size = 3.5) +
  scale_size_manual(values=c(3,1), guide = "none")+
  scale_color_manual(values = c("red3", "royalblue2"), labels = c("A" = "Original data", 
                                                                  "B" = "Permuted data"))+
  ylab("Emergence duration [days]") +
  xlab("Number of emerged individuals") +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), legend.text = element_text(size = 13), 
        legend.key = element_blank(), 
        plot.title = element_text(size = 13, face = "bold"),  
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(size = 11, face = "bold"), strip.background = element_blank(),
        axis.title = element_text(size = 15), 
        axis.title.x = element_text(margin = margin(10, 10, -10, 10)), 
        axis.title.y = element_text(margin = margin(10, 10, -10, 10)),
        axis.text = element_text(size = 13),
        plot.margin = margin(1, 1, 1, 1)
  ) 



