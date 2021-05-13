packages=c('plyr','dplyr','tidyr','ggplot2','tsModel',
           'chron','lubridate','reshape2',
           'mgcv','stringr','epiDisplay',
           'splines','zoo')
lapply(packages, library, character.only=T)

setwd("C:/Users/chary/OneDrive/Documents/SLU/Sping 2020/PUBH 5960 Capstone in Public Health Biostatistics/Project/")
GAM <- read.csv("gam_model_results.csv")
library(data.table)
GAM <- as.data.table(GAM)
head(GAM)
Xuchang <- GAM[city=='Xuchang',list(eac,aap,lag,excess_risk,lower_CI,upper_CI)][1:15]
Zhengzhou <- GAM[city=='Zhengzhou',list(eac,aap,lag,excess_risk,lower_CI,upper_CI)][1:15]

##----Xuchang Plot----

aap = c('pm[2.5]','pm[2.5]','pm[2.5]',
        'pm[10]','pm[10]','pm[10]',
        'o[3]','o[3]','o[3]',
        'so[2]','so[2]','so[2]',
        'no[2]','no[2]','no[2]')

strata = c('All-Cause','Cardiovascular','Respiratory')

paste(paste(as.character(Xuchang[aap=='pm25',excess_risk]), collapse=", "),paste(','),
paste(as.character(Xuchang[aap=='pm10',excess_risk]), collapse=", "),paste(','),
paste(as.character(Xuchang[aap=='o3',excess_risk]), collapse=", "),paste(','),
paste(as.character(Xuchang[aap=='so2',excess_risk]), collapse=", "),paste(','),
paste(as.character(Xuchang[aap=='no2',excess_risk]), collapse=", "))

er    = c(0.344076975, 0.312704103, 1.256516132,
          0.312909444, 0.341026851, 0.729280482,
          -0.320658203, -0.832002535, -2.64986884,
          0.826331112, 0.549135212, 2.625020735,
          1.259445651, 1.822966656, 3.840776339)

paste(paste(as.character(Xuchang[aap=='pm25',lower_CI]), collapse=", "),paste(','),
paste(as.character(Xuchang[aap=='pm10',lower_CI]), collapse=", "),paste(','),
paste(as.character(Xuchang[aap=='o3',lower_CI]), collapse=", "),paste(','),
paste(as.character(Xuchang[aap=='so2',lower_CI]), collapse=", "),paste(','),
paste(as.character(Xuchang[aap=='no2',lower_CI]), collapse=", "))

lower = c(0.106469531, -0.029911072, 0.43524732, 
          0.146571723, 0.101302168, 0.134722807, 
          -0.66054354, -1.327731905, -3.908727323, 
          -0.080216309, -0.770410924, -0.575474921, 
          0.646146402, 0.910572261, 1.612612064)

paste(paste(as.character(Xuchang[aap=='pm25',upper_CI]), collapse=", "),paste(','),
      paste(as.character(Xuchang[aap=='pm10',upper_CI]), collapse=", "),paste(','),
      paste(as.character(Xuchang[aap=='o3',upper_CI]), collapse=", "),paste(','),
      paste(as.character(Xuchang[aap=='so2',upper_CI]), collapse=", "),paste(','),
      paste(as.character(Xuchang[aap=='no2',upper_CI]), collapse=", "))

upper = c(0.582248392, 0.656493482, 2.08450054, 
          0.479523442, 0.581325631, 1.32736839, 
          0.020390035, -0.333782622, -1.374518488, 
          1.741103413, 1.886228554, 5.928540997, 
          1.876482112, 2.74361057, 6.117799865)

data = data.frame(aap,strata,er,lower,upper)

data$aap_f = factor(data$aap, levels=c('pm[2.5]','pm[10]','o[3]','so[2]','no[2]'))
data$strata_f = factor(data$strata, levels=c('Respiratory','Cardiovascular','All-Cause'))

windows(8,3)
ggplot(data,aes(er,strata_f,color=aap_f))+
  geom_point(size=1)+
  facet_grid(.~aap_f,scales = "free",
             labeller=label_parsed)+
  geom_errorbarh(aes(xmin=lower,xmax=upper),
                 size=0.75,height=0)+
  geom_vline(xintercept=0, size=0.5)+
  theme(legend.position="none")+ 
  xlab("ER (%, 95% CI)")+ylab(NULL)+
  theme(axis.title.x =element_text(size=11),
        axis.title.y=element_text(size=16))


##----Zhengzhou Plot----

aap = c('pm[2.5]','pm[2.5]','pm[2.5]',
        'pm[10]','pm[10]','pm[10]',
        'o[3]','o[3]','o[3]',
        'so[2]','so[2]','so[2]',
        'no[2]','no[2]','no[2]')

strata = c('All-Cause','Cardiovascular','Respiratory')

paste(paste(as.character(Zhengzhou[aap=='pm25',excess_risk]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='pm10',excess_risk]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='o3',excess_risk]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='so2',excess_risk]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='no2',excess_risk]), collapse=", "))

er    = c(0.218752962, 0.179284611, 0.458604809, 
          0.287139071, 0.206010393, 0.308842016, 
          -0.629577518, -0.703068894, -0.244857785, 
          4.751000597, 4.021562458, 4.664626905, 
          1.807657917, 2.052323641, 1.580784775)

paste(paste(as.character(Zhengzhou[aap=='pm25',lower_CI]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='pm10',lower_CI]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='o3',lower_CI]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='so2',lower_CI]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='no2',lower_CI]), collapse=", "))

lower = c(0.09794538, 0.019515807, 0.1670061, 
          0.208222355, 0.099618081, 0.106506187, 
          -0.789753872, -0.923520396, -0.669374466, 
          3.962339891, 3.006282867, 2.618855116, 
          1.549438837, 1.698977953, 0.913482014)

paste(paste(as.character(Zhengzhou[aap=='pm25',upper_CI]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='pm10',upper_CI]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='o3',upper_CI]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='so2',upper_CI]), collapse=", "),paste(','),
      paste(as.character(Zhengzhou[aap=='no2',upper_CI]), collapse=", "))

upper = c(0.339706346, 0.339308625, 0.751052398, 
          0.366117936, 0.312515787, 0.511586807, 
          -0.469142558, -0.482126874, 0.181473184, 
          5.545644102, 5.046849135, 6.751182449, 
          2.066533594, 2.406897003, 2.252500157)

data2 = data.frame(aap,strata,er,lower,upper)

data2$aap_f = factor(data$aap, levels=c('pm[2.5]','pm[10]','o[3]','so[2]','no[2]'))
data2$strata_f = factor(data$strata, levels=c('Respiratory','Cardiovascular','All-Cause'))

windows(8,3)
ggplot(data2,aes(er,strata_f,color=aap_f))+
  geom_point(size=1)+
  facet_grid(.~aap_f,scales = "free",
             labeller=label_parsed)+
  geom_errorbarh(aes(xmin=lower,xmax=upper),
                 size=0.75,height=0)+
  geom_vline(xintercept=0, size=0.5)+
  theme(legend.position="none")+ 
  xlab("ER (%, 95% CI)")+ylab(NULL)+
  theme(axis.title.x =element_text(size=11),
        axis.title.y=element_text(size=16))

