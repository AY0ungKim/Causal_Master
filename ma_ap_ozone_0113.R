library(dplyr);library(sf);library(data.table);library(dlnm);library(survival)
library(ggplot2);library(splines)
library(tsModel);library(lubridate)

#############################1. data preprocess############################
rm(list=ls());gc()
###file import
setwd("C:/Users/user1/Dropbox (개인용)/발표/MCC/data")


cco1521<-fread("230910_death0021_dup.csv")%>%filter(matched_death_dt<=as.Date("2021-12-31")&matched_death_dt>=as.Date("2015-01-01"))%>%arrange(person_id,matched_death_dt)

pm1521<-fread("PM25_1521_행정동_groupby_gam보정.csv")%>%select(c(date,sgg_h,gam))
oz24h1521<-fread("O3_24H_1521_행정동_groupby.csv")%>%select(c(date,sgg_h,gam))
oz8h1521<-fread("O3_8H_1521_행정동_groupby.csv")%>%select(c(date,sgg_h,gam))
names(pm1521)[3]="PM2_5";names(oz24h1521)[3]="OZONE_24H";names(oz8h1521)[3]="OZONE_8H"
holiday<-readxl::read_xlsx("C:/Users/user1/Dropbox (개인용)/취약계층연구/glm+ed/231009_holiday_0021.xlsx")
holiday$date<-as.Date(paste0(holiday$YEAR,"-",holiday$MONTH,"-",holiday$DAY))
holiday<-holiday%>%select(c(holiday,date))

pm1521<-pm1521%>%left_join(oz24h1521,by=c("date","sgg_h"))
pm1521$date<-as.Date(pm1521$date);oz8h1521$date<-as.Date(oz8h1521$date)
pm1521<-pm1521%>%left_join(oz8h1521,by=c("date","sgg_h"))

met1521<-fread("0621_변수_행정동별_평균_229.csv")%>%filter(date>=as.Date("2015-01-01") & date<=as.Date("2021-12-31"))

met1521$temperature_2m_era5_land=met1521$temperature_2m_era5_land-273.15
tlag <- Lag(met1521$temperature_2m_era5_land,0:21,group=list(met1521$sgg_h, met1521$year))
ttlag <- cbind(met1521[,c(3,2,1)], tlag) %>% data.table; rm(tlag)
colnames(ttlag) <- c("year", "sgg_h", "date", paste0("temp_lag", 0:21))
ttlag$date <- ymd(ttlag$date)



###select variables 
cco1521<-cco1521%>%select(person_id,matched_death_dt,icd_cd,age_g,sgg_h,sex,CASE_YN)%>%filter(year(matched_death_dt)%in%c(2015:2021))
substrRight <- function(x, n) {
  substr(x, nchar(x)-n+1, nchar(x))
}
###fill the NA
non<- unique(cco1521$sgg_h)[!unique(cco1521$sgg_h) %in% unique(pm1521$sgg_h)] #check
cco1521$sgg_h<-ifelse(substrRight(cco1521$sgg_h,1)!=0,paste0(substr(cco1521$sgg_h,1,4),0)%>%as.numeric(),cco1521$sgg_h) #integrate
pm1521$sgg_h<-ifelse(substrRight(pm1521$sgg_h,1)!=0,paste0(substr(pm1521$sgg_h,1,4),0)%>%as.numeric(),pm1521$sgg_h) #integrate
pm1521<-pm1521%>%filter(!(pm1521$sgg_h%in%non))

pm1521=distinct(pm1521,date,sgg_h,.keep_all=T)
pm1521<-pm1521[pm1521$sgg_h!=33010,]
cco1521<-cco1521[cco1521$sgg_h!=33010,]
cco1521<-cco1521%>%filter(!(cco1521$sgg_h%in%non))
cco1521$sgg_h[cco1521$sgg_h==23030]<-23090 #integrate namgu+michuhol (same)

#31320 33010 33310 34320 34390 non-matched


###left join
#with cco and temp
cco1521$matched_death_dt<-as.Date(cco1521$matched_death_dt)
cco1<-cco1521%>%left_join(pm1521,by=c('matched_death_dt'='date','sgg_h'))


#with meteorological data
cco1<-cco1%>%filter(!sgg_h%in%c(37430,39010,39020,23090)) #except jeju

met1521$date<-as.Date(met1521$date)
cco1<-cco1%>%left_join(met1521,by=c('matched_death_dt'='date','sgg_h'='sgg_h'))


#temp and met

tm1521<-pm1521%>%left_join(met1521,by=c('date','sgg_h'))
tm1521<-tm1521%>%filter(!sgg_h%in%c(37430,39010,39020,23090,29010,31280,33040)) #except jeju

tm1521$humidity<-100-5*((tm1521$temperature_2m_era5_land-32)*(5/9)
                        -(tm1521$dewpoint_temperature_2m_era5_land-32)*(5/9))
tm1521$wind_speed<-sqrt(tm1521$u_component_of_wind_10m_era5_land^2+tm1521$v_component_of_wind_10m_era5_land^2)
tm1521$prep<-ifelse(tm1521$total_precipitation_era5_land<quantile(tm1521$total_precipitation_era5_land)[3],"normal",
                    ifelse(tm1521$total_precipitation_era5_land<quantile(tm1521$total_precipitation_era5_land)[4],"hard","extreme"))

#OZONE PERCENTILE
sgg<-unique(tm1521$sgg_h)%>%sort
temp<-sgg%>%lapply(function(x)tm1521[tm1521$sgg_h==x,])
temp<-temp%>%lapply(function(x){
  x<-x%>%data.table
  x<-x%>%mutate(OZONE_8H_RANK=round(rank(OZONE_8H)/length(OZONE_8H)*100,1),
                OZONE_24H_RANK=round(rank(OZONE_24H)/length(OZONE_24H)*100,1))
})

tm1521<-bind_rows(temp)%>%arrange(sgg_h,date)


lag01<-tm1521%>%mutate(pm_01=runMean(PM2_5,0:1,group=sgg_h),
                       OZONE_8H_01=runMean(OZONE_8H_RANK,0:1,group=sgg_h),
                       OZONE_24H_01=runMean(OZONE_24H_RANK,0:1,group=sgg_h),
                       temp_01=runMean(temperature_2m_era5_land,0:1,group=sgg_h),
                       prep_01=runMean(total_precipitation_era5_land,0:1,group=sgg_h),
                       wind_01=runMean(wind_speed,0:1,group=sgg_h),
                       hum_01=runMean(humidity,0:1,group=sgg_h),
                       aod_01=runMean(optical_depth_047,0:1,group=sgg_h))

#cco1 data
cco1<-cco1%>%filter(!sgg_h%in%c(37430,39010,39020,23090,29010,31280,33040)) #except jeju

cco1$humidity<-100-5*((cco1$temperature_2m_era5_land-32)*(5/9)
                      -(cco1$dewpoint_temperature_2m_era5_land-32)*(5/9))
cco1$wind_speed<-sqrt(cco1$u_component_of_wind_10m_era5_land^2+cco1$v_component_of_wind_10m_era5_land^2)
cco1$prep<-ifelse(cco1$total_precipitation_era5_land<quantile(cco1$total_precipitation_era5_land)[3],"normal",
                  ifelse(cco1$total_precipitation_era5_land<quantile(cco1$total_precipitation_era5_land)[4],"hard","extreme"))

sum(is.na(cco1)) #1663915

#################################dlnm start################################
###########MoVING Average (0,1)
#make moving average

library(tsModel)

tm1521<-tm1521%>%arrange(sgg_h,date)
tm1521$PM2_5<-tm1521$PM2_5/10
#tm1521$OZONE_8H<-tm1521$OZONE_8H/10
#tm1521$OZONE_24H<-tm1521$OZONE_24H/10


#cco1521<-cco1521%>%filter(month(matched_death_dt)%in%c(6:9))
cco1<-cco1%>%select(c(person_id,matched_death_dt,icd_cd,age_g,sgg_h,sex,CASE_YN))

cco2<-cco1%>%left_join(lag01,by=c("matched_death_dt"="date","sgg_h"))#data.frame

cco2<-cco2%>%left_join(holiday,by=c("matched_death_dt"="date"))
cco2$holiday[is.na(cco2$holiday)]=0
####################Modeling start##################3
##############Model 1: only CCO
rtable<-data.frame()

#model fitting
varfun<-"ns"
vardf<-4
lagknots=c(1,3)
z<-qnorm(0.975)
cco4<-na.omit(cco2)
cco5 <- left_join(cco4, ttlag, by=c("matched_death_dt"="date", "sgg_h")) %>% as_tibble()
x<-cco5[,c("temp_lag0","temp_lag1","temp_lag2","temp_lag3","temp_lag4","temp_lag5","temp_lag6",
           "temp_lag7","temp_lag8","temp_lag9","temp_lag10","temp_lag11","temp_lag12","temp_lag13",
           "temp_lag14","temp_lag15","temp_lag16","temp_lag17","temp_lag18","temp_lag19","temp_lag20","temp_lag21"
)]
#ldat$income1<-ifelse(ldat$income=="middle",1,0)
#  ldat$year1<-ifelse(ldat$year==j,0,1)

cb<-crossbasis(x,argvar=list(fun=varfun,df=vardf),arglag=list(fun=varfun,knots=lagknots))



model1<-survival::clogit(CASE_YN~cco5$pm_01+cb+holiday+strata(person_id),cco5) #only with air pollution

gc()
coef<-summary(model1)$coef
vcov<-vcov(model1)

#get relative risk

coef <- summary(model1)$coef[1]
se <- vcov(model1)[1,1] %>% sqrt()
rr <- c(exp(coef), exp(coef-z*se), exp(coef+z*se), se)


rtable<-rbind(rtable,rr)
#rm(model1,dat,predvar);gc()

##############Model 2: CCO+COVARIATES ADJUSTMENT
model2<-survival::clogit(CASE_YN~pm_01+cb+holiday+OZONE_8H_01+prep_01+wind_01+hum_01+aod_01+strata(person_id),cco5) #only with air pollution
gc()
coef<-summary(model2)$coef[1]
se <- vcov(model2)[1,1] %>% sqrt()

#get relative risk
rr <- c(exp(coef), exp(coef-z*se), exp(coef+z*se), se)
rtable<-rbind(rtable,rr)




##############Model 3: CCO+with GPS(ref=https://www.youtube.com/watch?v=mSueFNm_a9g)
#Adjust total_precipitation_era5_land+surface_pressure_era5_land+humidity+sur_refl_b01
#+wind_speed+optical_depth_047+evi+pop_density+covid
#make GPS

##Make weight using the glm

x_fit<- glm (pm_01~ OZONE_8H_01+cb+prep_01+hum_01+ wind_01+aod_01, data = cco4)

cco4$GPS<- dnorm (cco4$pm_01, mean=x_fit$fitted.values, sd (x_fit$residuals))
ipw1 <- (1/cco4$GPS)
ipw1[ipw1>=quantile(ipw1,0.95)]=quantile(ipw1,0.95)
ipw1[ipw1<=quantile(ipw1,0.05)]=quantile(ipw1,0.05)



##Using WeightIt pacakage
#library(WeightIt)
#library(splines)
#ipw1<-weightit(pm_01~temp_01+prep_01+hum_01+ wind_01+aod_01, data = cco4,method="gbm",stabilize = TRUE,criterion="p.mean")


#######Makeing model#######
model3<-clogit(CASE_YN~pm_01+cb+holiday+strata(person_id),cco5,weights=ipw1,method="efron") 
#There can't be applied to exact partial likelihood
coef<-summary(model3)$coef[1]
se <- vcov(model3)[1,1] %>% sqrt()

#get relative risk
rr <- c(exp(coef), exp(coef-z*se), exp(coef+z*se), se)
rtable<-rbind(rtable,rr)













