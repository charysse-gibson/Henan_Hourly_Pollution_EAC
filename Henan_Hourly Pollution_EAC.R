######################################################################
## Title: SLU Capstone Project                                    
## Date created: February 14, 2020
######################################################################
# Inputs: hnHour1719.csv
# outputs: gam_model_results.csv
######################################################################

setwd("C:/Users/chary/OneDrive/Documents/SLU/Sping 2020/PUBH 5960 Capstone in Public Health Biostatistics/Project/")

# Variables:
  # date: the observation days from 2017.1.1 through 2019.12.31
  # hour: 0-23 in one day
  # all: all the causes for the emergency call except accidental causes
  # cvd: cardiovascular causes
  # resp: respiratory causes
  # city: the Chinese name of the two cities in this study
  # cityE: the English name of the two cities in this study
  # pm25, o3, pm10, so2, no2 are the concentrations of the air pollutants
  # tm: hourly temperature
  # rh: relative humidity

#----import data----
henan_data <- 
  read.csv("hnHour1719.csv")

library(data.table)
henan_data <- as.data.table(henan_data)

#----explore data----
head(henan_data)
library(Hmisc)
str(henan_data)

# convert date to posix date
library(lubridate)
henan_data[,obs_date:=as.Date(date,format="%Y/%m/%d")]

# verify unique observation dates for each city from 2017.1.1 to 2019.12.31
henan_data[,unique(cityE)]
henan_data[cityE=='Zhengzhou',list(.N,unique_dates=uniqueN(obs_date))]
1095*24 #26280
henan_data[cityE=='Xuchang',list(.N,unique_dates=uniqueN(obs_date))]
# check dates
henan_data[cityE=='Zhengzhou',list(min_date=min(obs_date),max_date=max(obs_date))]
henan_data[cityE=='Xuchang',list(min_date=min(obs_date),max_date=max(obs_date))]

##----descriptive stats----
summary(henan_data)
# number of EACs by city
henan_data[cityE=='Xuchang',list(all_eac=sum(all),cvd_cases=sum(cvd),resp_cases=sum(resp))]
#    all_eac cvd_cases resp_cases
# 1:   36551     16553       2548
henan_data[cityE=='Zhengzhou',list(all_eac=sum(all),cvd_cases=sum(cvd),resp_cases=sum(resp))]
#    all_eac cvd_cases resp_cases
# 1:  269744     89837      21352

library(tableone)
myvars <- c("all","cvd","resp","pm25", "pm10", "o3", "so2", "no2", "tm", "rh")
tab1 <- CreateTableOne(vars = myvars, strata='cityE', data = henan_data)
tab1
summary(tab1) #includes details and standardized mean differences

# Xuchang quantiles
quantile(henan_data[cityE=='Xuchang',pm25])
quantile(henan_data[cityE=='Xuchang',pm10])
quantile(henan_data[cityE=='Xuchang',o3])
quantile(henan_data[cityE=='Xuchang',so2])
quantile(henan_data[cityE=='Xuchang',no2])

quantile(henan_data[cityE=='Xuchang',tm])
quantile(henan_data[cityE=='Xuchang',rh])

quantile(henan_data[cityE=='Xuchang',all])
quantile(henan_data[cityE=='Xuchang',cvd])
quantile(henan_data[cityE=='Xuchang',resp])

# Zhengzhou quantiles
quantile(henan_data[cityE=='Zhengzhou',pm25])
quantile(henan_data[cityE=='Zhengzhou',pm10])
quantile(henan_data[cityE=='Zhengzhou',o3])
quantile(henan_data[cityE=='Zhengzhou',so2])
quantile(henan_data[cityE=='Zhengzhou',no2])

quantile(henan_data[cityE=='Zhengzhou',tm])
quantile(henan_data[cityE=='Zhengzhou',rh])

quantile(henan_data[cityE=='Zhengzhou',all])
quantile(henan_data[cityE=='Zhengzhou',cvd])
quantile(henan_data[cityE=='Zhengzhou',resp])
                  
##----variable distributions----
hist(henan_data$all)
hist(henan_data$cvd)
hist(henan_data$resp)
hist(henan_data$pm25)
hist(henan_data$pm10)
hist(henan_data$o3)
hist(henan_data$no2)
hist(henan_data$so2)

##----ready GAM model variables and datasets----
library(tsModel) # Time Series Modeling for Air Pollution and Health package
library(mgcv) # Mixed GAM Computation Vehicle package
?gam
#GCV score is minimised generalised cross-validation (un-biased risk estimator)

# create t (daily number variable) by city
setkey(henan_data,cityE,obs_date,hour)
t <- henan_data[,list(obs_date=unique(obs_date)),by='cityE']
t <- t[,list(cityE,obs_date,t=1:.N),by='cityE']

# merge t back to henan_data
setkey(henan_data,cityE,obs_date)
setkey(t,cityE,obs_date)
henan_data[t,t:=t]
head(henan_data)

# create time-series variables
henan_data[,':=' (dow=wday(obs_date),
                  rmtemp=runMean(tm,0:72),
                  rmrh=runMean(rh,0:72))]
henan_data[,':=' (pm25lag0=pm25/10,
                  pm10lag0=pm10/10,
                  so2lag0=so2/10, 
                  no2lag0=no2/10, 
                  o3lag0=o3/10)]
henan_data[,pmc:=pm10lag0-pm25lag0]

# subset by city
xuchang <- henan_data[cityE=='Xuchang']
zhengzhou <- henan_data[cityE=='Zhengzhou']

#----Create table for GAM models----
# 1 - create cartesian join of EAC types, air pollutants, and lag times
# 2 - attach cities to two separate cartesian joined tables
# 3 - bind tables together
# 4 - create model value columns (gcv, coef, se, er, lower, upper)

# create EAC types, air pollutants, and lag times variables
eac <- c('all','cvd','resp')
aam <- c('pm25','pm10','so2','no2','o3') #aam = ambient air pollutant
lag <- c(0,1,3)

# create cartesian joins
xuchang_mods <- as.data.table(CJ(eac,aam,lag))
zhengzhou_mods <- as.data.table(CJ(eac,aam,lag))

# attach city names
xuchang_mods[,city:='Xuchang']
setcolorder(xuchang_mods,c('city','eac','aam','lag'))
zhengzhou_mods[,city:='Zhengzhou']
setcolorder(zhengzhou_mods,c('city','eac','aam','lag'))

# bind city tables
gam_models <- 
  rbind(xuchang_mods,zhengzhou_mods)

# create model value columns (gcv, coef, se, er, lower, upper)
gam_models[,':=' (city=city,
                  eac=eac,
                  lag=lag,
                  gcv=numeric(),
                  coefficient=numeric(),
                  standard_error=numeric(),
                  excess_risk=numeric(),
                  lower_CI=numeric(),
                  upper_CI=numeric())]

# sort by city, lag, eac, aam
setkey(gam_models,city,lag,eac,aam)

#----GAM models, no lag---- 
library(splines)

#----Xuchang, all EACs, no lag----
fit.pm25L0=gam(all~pm25lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm25L0)$p.coeff[2] #0.00343
se=summary(fit.pm25L0)$se[2] #0.00121
er=(exp(coef)-1)*100 #0.344
lower=(exp(coef-1.96*se)-1)*100 #1.06
upper=(exp(coef+1.96*se)-1)*100 #0.582


gam_models[city=='Xuchang' & eac=='all' & aam=='pm25' & lag==0,
           ':=' (gcv=fit.pm25L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L0=gam(all~pm10lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm10L0)$p.coeff[2]
se=summary(fit.pm10L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='pm10' & lag==0,
           ':=' (gcv=fit.pm10L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L0=gam(all~o3lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.o3L0)$p.coeff[2]
se=summary(fit.o3L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='o3' & lag==0,
           ':=' (gcv=fit.o3L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L0=gam(all~no2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.no2L0)$p.coeff[2]
se=summary(fit.no2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='no2' & lag==0,
           ':=' (gcv=fit.no2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L0=gam(all~so2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.so2L0)$p.coeff[2]
se=summary(fit.so2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='so2' & lag==0,
           ':=' (gcv=fit.so2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Zhengzhou, all EACs, no lag----
fit.pm25L0=gam(all~pm25lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm25L0)$p.coeff[2]
se=summary(fit.pm25L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='pm25' & lag==0,
           ':=' (gcv=fit.pm25L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L0=gam(all~pm10lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm10L0)$p.coeff[2]
se=summary(fit.pm10L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='pm10' & lag==0,
           ':=' (gcv=fit.pm10L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L0=gam(all~o3lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=zhengzhou)
coef=summary(fit.o3L0)$p.coeff[2]
se=summary(fit.o3L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='o3' & lag==0,
           ':=' (gcv=fit.o3L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L0=gam(all~no2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.no2L0)$p.coeff[2]
se=summary(fit.no2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='no2' & lag==0,
           ':=' (gcv=fit.no2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L0=gam(all~so2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.so2L0)$p.coeff[2]
se=summary(fit.so2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='so2' & lag==0,
           ':=' (gcv=fit.so2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


#----Xuchang, cvd EACs, no lag----
fit.pm25L0=gam(cvd~pm25lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm25L0)$p.coeff[2] 
se=summary(fit.pm25L0)$se[2]
er=(exp(coef)-1)*100
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='pm25' & lag==0,
           ':=' (gcv=fit.pm25L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L0=gam(cvd~pm10lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm10L0)$p.coeff[2]
se=summary(fit.pm10L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='pm10' & lag==0,
           ':=' (gcv=fit.pm10L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L0=gam(cvd~o3lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=xuchang)
coef=summary(fit.o3L0)$p.coeff[2]
se=summary(fit.o3L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='o3' & lag==0,
           ':=' (gcv=fit.o3L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L0=gam(cvd~no2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.no2L0)$p.coeff[2]
se=summary(fit.no2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='no2' & lag==0,
           ':=' (gcv=fit.no2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L0=gam(cvd~so2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.so2L0)$p.coeff[2]
se=summary(fit.so2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='so2' & lag==0,
           ':=' (gcv=fit.so2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Zhengzhou, cvd EACs, no lag----
fit.pm25L0=gam(cvd~pm25lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm25L0)$p.coeff[2]
se=summary(fit.pm25L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='pm25' & lag==0,
           ':=' (gcv=fit.pm25L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L0=gam(cvd~pm10lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm10L0)$p.coeff[2]
se=summary(fit.pm10L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='pm10' & lag==0,
           ':=' (gcv=fit.pm10L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L0=gam(cvd~o3lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=zhengzhou)
coef=summary(fit.o3L0)$p.coeff[2]
se=summary(fit.o3L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='o3' & lag==0,
           ':=' (gcv=fit.o3L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L0=gam(cvd~no2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.no2L0)$p.coeff[2]
se=summary(fit.no2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='no2' & lag==0,
           ':=' (gcv=fit.no2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L0=gam(cvd~so2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.so2L0)$p.coeff[2]
se=summary(fit.so2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='so2' & lag==0,
           ':=' (gcv=fit.so2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Xuchang, resp EACs, no lag----
fit.pm25L0=gam(resp~pm25lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm25L0)$p.coeff[2] 
se=summary(fit.pm25L0)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='pm25' & lag==0,
           ':=' (gcv=fit.pm25L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L0=gam(resp~pm10lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm10L0)$p.coeff[2]
se=summary(fit.pm10L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='pm10' & lag==0,
           ':=' (gcv=fit.pm10L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L0=gam(resp~o3lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=xuchang)
coef=summary(fit.o3L0)$p.coeff[2]
se=summary(fit.o3L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='o3' & lag==0,
           ':=' (gcv=fit.o3L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L0=gam(resp~no2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.no2L0)$p.coeff[2]
se=summary(fit.no2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='no2' & lag==0,
           ':=' (gcv=fit.no2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L0=gam(resp~so2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.so2L0)$p.coeff[2]
se=summary(fit.so2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='so2' & lag==0,
           ':=' (gcv=fit.so2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Zhengzhou, resp EACs, no lag----
fit.pm25L0=gam(resp~pm25lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm25L0)$p.coeff[2]
se=summary(fit.pm25L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='pm25' & lag==0,
           ':=' (gcv=fit.pm25L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L0=gam(resp~pm10lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm10L0)$p.coeff[2]
se=summary(fit.pm10L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='pm10' & lag==0,
           ':=' (gcv=fit.pm10L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L0=gam(resp~o3lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=zhengzhou)
coef=summary(fit.o3L0)$p.coeff[2]
se=summary(fit.o3L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='o3' & lag==0,
           ':=' (gcv=fit.o3L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L0=gam(resp~no2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.no2L0)$p.coeff[2]
se=summary(fit.no2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='no2' & lag==0,
           ':=' (gcv=fit.no2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L0=gam(resp~so2lag0+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.so2L0)$p.coeff[2]
se=summary(fit.so2L0)$se[2]
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='so2' & lag==0,
           ':=' (gcv=fit.so2L0[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----GAM models, lag 1---- 
#----Xuchang, all EACs, lag 1----
fit.pm25L1=gam(all~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm25L1)$p.coeff[2]
se=summary(fit.pm25L1)$se[2] 
er=(exp(coef)-1)*100
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100

gam_models[city=='Xuchang' & eac=='all' & aam=='pm25' & lag==1,
           ':=' (gcv=fit.pm25L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L1=gam(all~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm10L1)$p.coeff[2] 
se=summary(fit.pm10L1)$se[2] 
er=(exp(coef)-1)*100
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='pm10' & lag==1,
           ':=' (gcv=fit.pm10L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L1=gam(all~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.o3L1)$p.coeff[2] 
se=summary(fit.o3L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='o3' & lag==1,
           ':=' (gcv=fit.o3L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L1=gam(all~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.no2L1)$p.coeff[2] 
se=summary(fit.no2L10)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='no2' & lag==1,
           ':=' (gcv=fit.no2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L1=gam(all~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.so2L1)$p.coeff[2] 
se=summary(fit.so2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='so2' & lag==1,
           ':=' (gcv=fit.so2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Xuchang, cvd EACs, lag 1----
fit.pm25L1=gam(cvd~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm25L1)$p.coeff[2] 
se=summary(fit.pm25L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='pm25' & lag==1,
           ':=' (gcv=fit.pm25L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L1=gam(cvd~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm10L1)$p.coeff[2] 
se=summary(fit.pm10L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='pm10' & lag==1,
           ':=' (gcv=fit.pm10L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L1=gam(cvd~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=xuchang)
coef=summary(fit.o3L1)$p.coeff[2] 
se=summary(fit.o3L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='o3' & lag==1,
           ':=' (gcv=fit.o3L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L1=gam(cvd~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.no2L1)$p.coeff[2] 
se=summary(fit.no2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='no2' & lag==1,
           ':=' (gcv=fit.no2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L1=gam(cvd~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.so2L1)$p.coeff[2] 
se=summary(fit.so2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='so2' & lag==1,
           ':=' (gcv=fit.so2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


#----Xuchang, resp EACs, lag 1----
fit.pm25L1=gam(resp~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm25L1)$p.coeff[2] 
se=summary(fit.pm25L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='pm25' & lag==1,
           ':=' (gcv=fit.pm25L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L1=gam(resp~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm10L1)$p.coeff[2] 
se=summary(fit.pm10L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='pm10' & lag==1,
           ':=' (gcv=fit.pm10L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L1=gam(resp~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=xuchang)
coef=summary(fit.o3L1)$p.coeff[2] 
se=summary(fit.o3L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='o3' & lag==1,
           ':=' (gcv=fit.o3L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L1=gam(resp~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.no2L1)$p.coeff[2] 
se=summary(fit.no2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='no2' & lag==1,
           ':=' (gcv=fit.no2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L1=gam(resp~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.so2L1)$p.coeff[2] 
se=summary(fit.so2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='so2' & lag==1,
           ':=' (gcv=fit.so2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Zhengzhou, all EACs, lag 1----
fit.pm25L1=gam(all~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm25L1)$p.coeff[2] 
se=summary(fit.pm25L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='pm25' & lag==1,
           ':=' (gcv=fit.pm25L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L1=gam(all~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm10L1)$p.coeff[2] 
se=summary(fit.pm10L1)$se[2] 
er=(exp(coef)-1)*100
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='pm10' & lag==1,
           ':=' (gcv=fit.pm10L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L1=gam(all~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=zhengzhou)
coef=summary(fit.o3L1)$p.coeff[2] 
se=summary(fit.o3L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='o3' & lag==1,
           ':=' (gcv=fit.o3L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L1=gam(all~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.no2L1)$p.coeff[2] 
se=summary(fit.no2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='no2' & lag==1,
           ':=' (gcv=fit.no2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L1=gam(all~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.so2L1)$p.coeff[2] 
se=summary(fit.so2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='so2' & lag==1,
           ':=' (gcv=fit.so2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Zhengzhou, cvd EACs, lag 1----
fit.pm25L1=gam(cvd~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm25L1)$p.coeff[2] 
se=summary(fit.pm25L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='pm25' & lag==1,
           ':=' (gcv=fit.pm25L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L1=gam(cvd~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm10L1)$p.coeff[2] 
se=summary(fit.pm10L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='pm10' & lag==1,
           ':=' (gcv=fit.pm10L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L1=gam(cvd~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=zhengzhou)
coef=summary(fit.o3L1)$p.coeff[2] 
se=summary(fit.o3L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='o3' & lag==1,
           ':=' (gcv=fit.o3L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L1=gam(cvd~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.no2L1)$p.coeff[2] 
se=summary(fit.no2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='no2' & lag==1,
           ':=' (gcv=fit.no2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L1=gam(cvd~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.so2L1)$p.coeff[2] 
se=summary(fit.so2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='so2' & lag==1,
           ':=' (gcv=fit.so2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


#----Zhengzhou, resp EACs, lag 1----
fit.pm25L1=gam(resp~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm25L1)$p.coeff[2] 
se=summary(fit.pm25L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='pm25' & lag==1,
           ':=' (gcv=fit.pm25L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L1=gam(resp~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm10L1)$p.coeff[2] 
se=summary(fit.pm10L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='pm10' & lag==1,
           ':=' (gcv=fit.pm10L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L1=gam(resp~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=zhengzhou)
coef=summary(fit.o3L1)$p.coeff[2] 
se=summary(fit.o3L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='o3' & lag==1,
           ':=' (gcv=fit.o3L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L1=gam(resp~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.no2L1)$p.coeff[2] 
se=summary(fit.no2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='no2' & lag==1,
           ':=' (gcv=fit.no2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L1=gam(resp~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.so2L1)$p.coeff[2] 
se=summary(fit.so2L1)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='so2' & lag==1,
           ':=' (gcv=fit.so2L1[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----GAM models, lag 3---- 
#----Xuchang, all EACs, lag 3----
fit.pm25L3=gam(all~Lag(pm25lag0,3)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm25L3)$p.coeff[2] 
se=summary(fit.pm25L0)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='pm25' & lag==3,
           ':=' (gcv=fit.pm25L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L3=gam(all~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm10L3)$p.coeff[2] 
se=summary(fit.pm10L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='pm10' & lag==3,
           ':=' (gcv=fit.pm10L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L3=gam(all~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=xuchang)
coef=summary(fit.o3L3)$p.coeff[2] 
se=summary(fit.o3L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='o3' & lag==3,
           ':=' (gcv=fit.o3L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L3=gam(all~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.no2L3)$p.coeff[2] 
se=summary(fit.no2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='all' & aam=='no2' & lag==3,
           ':=' (gcv=fit.no2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L3=gam(all~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.so2L3)$p.coeff[2] 
se=summary(fit.so2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100

gam_models[city=='Xuchang' & eac=='all' & aam=='so2' & lag==3,
           ':=' (gcv=fit.so2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Xuchang, cvd EACs, lag 3----
fit.pm25L3=gam(cvd~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm25L3)$p.coeff[2] #0.00343
se=summary(fit.pm25L3)$se[2] #0.00121
er=(exp(coef)-1)*100 #0.344
lower=(exp(coef-1.96*se)-1)*100 #1.06
upper=(exp(coef+1.96*se)-1)*100 #0.582

gam_models[city=='Xuchang' & eac=='cvd' & aam=='pm25' & lag==3,
           ':=' (gcv=fit.pm25L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L3=gam(cvd~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm10L3)$p.coeff[2] 
se=summary(fit.pm10L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='pm10' & lag==3,
           ':=' (gcv=fit.pm10L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L3=gam(cvd~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=xuchang)
coef=summary(fit.o3L3)$p.coeff[2] 
se=summary(fit.o3L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='o3' & lag==3,
           ':=' (gcv=fit.o3L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L3=gam(cvd~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.no2L3)$p.coeff[2] 
se=summary(fit.no2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='no2' & lag==3,
           ':=' (gcv=fit.no2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L3=gam(cvd~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.so2L3)$p.coeff[2] 
se=summary(fit.so2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='cvd' & aam=='so2' & lag==3,
           ':=' (gcv=fit.so2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


#----Xuchang, resp EACs, lag 3----
fit.pm25L3=gam(resp~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm25L3)$p.coeff[2] 
se=summary(fit.pm25L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='pm25' & lag==3,
           ':=' (gcv=fit.pm25L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L3=gam(resp~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=xuchang)
coef=summary(fit.pm10L3)$p.coeff[2] 
se=summary(fit.pm10L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='pm10' & lag==3,
           ':=' (gcv=fit.pm10L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L3=gam(resp~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=xuchang)
coef=summary(fit.o3L3)$p.coeff[2] 
se=summary(fit.o3L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='o3' & lag==3,
           ':=' (gcv=fit.o3L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L3=gam(resp~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.no2L3)$p.coeff[2] 
se=summary(fit.no2L3)$se[2] 
er=(exp(coef)-1)*100 #0.344
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='no2' & lag==3,
           ':=' (gcv=fit.no2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L3=gam(resp~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=xuchang)
coef=summary(fit.so2L3)$p.coeff[2] 
se=summary(fit.so2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Xuchang' & eac=='resp' & aam=='so2' & lag==3,
           ':=' (gcv=fit.so2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Zhengzhou, all EACs, lag 3----
fit.pm25L3=gam(all~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm25L3)$p.coeff[2] 
se=summary(fit.pm25L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='pm25' & lag==3,
           ':=' (gcv=fit.pm25L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L3=gam(all~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm10L3)$p.coeff[2] 
se=summary(fit.pm10L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='pm10' & lag==3,
           ':=' (gcv=fit.pm10L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L3=gam(all~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=zhengzhou)
coef=summary(fit.o3L3)$p.coeff[2] 
se=summary(fit.o3L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='o3' & lag==3,
           ':=' (gcv=fit.o3L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L3=gam(all~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.no2L3)$p.coeff[2] 
se=summary(fit.no2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='no2' & lag==3,
           ':=' (gcv=fit.no2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L3=gam(all~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.so2L3)$p.coeff[2] 
se=summary(fit.so2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='all' & aam=='so2' & lag==3,
           ':=' (gcv=fit.so2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

#----Zhengzhou, cvd EACs, lag 3----
fit.pm25L3=gam(cvd~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm25L3)$p.coeff[2] 
se=summary(fit.pm25L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='pm25' & lag==3,
           ':=' (gcv=fit.pm25L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L3=gam(cvd~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm10L3)$p.coeff[2] 
se=summary(fit.pm10L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='pm10' & lag==3,
           ':=' (gcv=fit.pm10L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L3=gam(cvd~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=zhengzhou)
coef=summary(fit.o3L3)$p.coeff[2] 
se=summary(fit.o3L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='o3' & lag==3,
           ':=' (gcv=fit.o3L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L3=gam(cvd~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.no2L3)$p.coeff[2] 
se=summary(fit.no2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='no2' & lag==3,
           ':=' (gcv=fit.no2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L3=gam(cvd~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.so2L3)$p.coeff[2] 
se=summary(fit.so2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='cvd' & aam=='so2' & lag==3,
           ':=' (gcv=fit.so2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


#----Zhengzhou, resp EACs, lag 3----
fit.pm25L3=gam(resp~Lag(pm25lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm25L3)$p.coeff[2] 
se=summary(fit.pm25L3)$se[2] 
er=(exp(coef)-1)*100
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='pm25' & lag==3,
           ':=' (gcv=fit.pm25L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.pm10L3=gam(resp~Lag(pm10lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                 ns(rmtemp,6+1)+ns(rmrh,3+1),
               family=quasipoisson(link="log"), 
               scale=-1,
               control=gam.control(epsilon=0.0000001,maxit=1000), 
               na.action=na.omit, 
               data=zhengzhou)
coef=summary(fit.pm10L3)$p.coeff[2] 
se=summary(fit.pm10L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='pm10' & lag==3,
           ':=' (gcv=fit.pm10L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.o3L3=gam(resp~Lag(o3lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
               ns(rmtemp,6+1)+ns(rmrh,3+1),
             family=quasipoisson(link="log"), 
             scale=-1,
             control=gam.control(epsilon=0.0000001,maxit=1000), 
             na.action=na.omit, 
             data=zhengzhou)
coef=summary(fit.o3L3)$p.coeff[2]
se=summary(fit.o3L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='o3' & lag==3,
           ':=' (gcv=fit.o3L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.no2L3=gam(resp~Lag(no2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.no2L3)$p.coeff[2] 
se=summary(fit.no2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='no2' & lag==3,
           ':=' (gcv=fit.no2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]


fit.so2L3=gam(resp~Lag(so2lag0,1)+ns(t,3*6+1)+ns(hour,3)+as.factor(dow)+
                ns(rmtemp,6+1)+ns(rmrh,3+1),
              family=quasipoisson(link="log"), 
              scale=-1,
              control=gam.control(epsilon=0.0000001,maxit=1000), 
              na.action=na.omit, 
              data=zhengzhou)
coef=summary(fit.so2L3)$p.coeff[2] 
se=summary(fit.so2L3)$se[2] 
er=(exp(coef)-1)*100 
lower=(exp(coef-1.96*se)-1)*100 
upper=(exp(coef+1.96*se)-1)*100 

gam_models[city=='Zhengzhou' & eac=='resp' & aam=='so2' & lag==3,
           ':=' (gcv=fit.so2L3[["gcv.ubre"]][["GCV.Cp"]],
                 coefficient=coef,
                 standard_error=se,
                 excess_risk=er,
                 lower_CI=lower,
                 upper_CI=upper)]

##----Outputs----

# Check GAM model results
gam_models

# create GAM models csv 
write.csv(gam_models, file='gam_model_results.csv')
