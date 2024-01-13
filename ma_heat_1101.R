library(dplyr);library(sf);library(data.table);library(dlnm);library(survival);library(wCorr);
library(ggplot2)

#############################1. data preprocess############################

###file import
cco0021<-fread("C:/data/230910_death0021_dup.csv")%>%filter(matched_death_dt<=as.Date("2021-12-31"))%>%arrange(person_id,matched_death_dt)
temp1120<-fread("C:/data/230826_lagged_tmean_250ver.csv")
met1120<-fread("C:/data/0621_변수_행정동별_평균_229.csv")%>%filter(date>=as.Date("2011-01-01") & date<=as.Date("2020-12-31"))
###select variables 
cco1120<-cco0021%>%select(person_id,matched_death_dt,icd_cd,age_g,sgg_h,sex,CASE_YN)%>%filter(year(matched_death_dt)%in%c(2011:2020))
substrRight <- function(x, n) {
  substr(x, nchar(x)-n+1, nchar(x))
}
###fill the NA
non<- unique(cco1120$sgg_h)[!unique(cco1120$sgg_h) %in% unique(temp1120$sgg_h)] #check
cco1120$sgg_h<-ifelse(substrRight(cco1120$sgg_h,1)!=0,paste0(substr(cco1120$sgg_h,1,4),0)%>%as.numeric(),cco1120$sgg_h) #integrate
temp1120$sgg_h<-ifelse(substrRight(temp1120$sgg_h,1)!=0,paste0(substr(temp1120$sgg_h,1,4),0)%>%as.numeric(),temp1120$sgg_h) #integrate
temp1120<-temp1120%>%filter(!(temp1120$sgg_h%in%non))
temp1120<-temp1120%>%select(-sgg_b)
temp1120=distinct(temp1120,date,sgg_h,.keep_all=T)
temp1120<-temp1120[temp1120$sgg_h!=33010,]
cco1120<-cco1120[cco1120$sgg_h!=33010,]
cco1120<-cco1120%>%filter(!(cco1120$sgg_h%in%non))
cco1120$sgg_h[cco1120$sgg_h==23030]<-23090 #integrate namgu+michuhol (same)

#31320 33010 33310 34320 34390 non-matched


###left join
#with cco and temp
cco1<-cco1120%>%left_join(temp1120,by=c('matched_death_dt'='date','sgg_h'))


#with meteorological data
cco1<-cco1%>%filter(!sgg_h%in%c(37430,39010,39020,23090)) #except jeju
cco1<-cco1%>%left_join(met1120,by=c('matched_death_dt'='date','sgg_h'='sgg_h'))


#temp and met
temp1120_tm<-temp1120%>%select(c(date,sgg_h,tmean))
tm1120<-temp1120_tm%>%left_join(met1120,by=c('date','sgg_h'))
tm1120<-tm1120%>%filter(!sgg_h%in%c(37430,39010,39020,23090,29010,31280,33040)) #except jeju

tm1120$humidity<-100-5*((tm1120$temperature_2m_era5_land-32)*(5/9)
                        -(tm1120$dewpoint_temperature_2m_era5_land-32)*(5/9))
tm1120$wind_speed<-sqrt(tm1120$u_component_of_wind_10m_era5_land^2+tm1120$v_component_of_wind_10m_era5_land^2)
tm1120$prep<-ifelse(tm1120$total_precipitation_era5_land<quantile(tm1120$total_precipitation_era5_land)[3],"normal",
                    ifelse(tm1120$total_precipitation_era5_land<quantile(tm1120$total_precipitation_era5_land)[4],"hard","extreme"))
#cco1 data
cco1$humidity<-100-5*((cco1$temperature_2m_era5_land-32)*(5/9)
                      -(cco1$dewpoint_temperature_2m_era5_land-32)*(5/9))
cco1$wind_speed<-sqrt(cco1$u_component_of_wind_10m_era5_land^2+cco1$v_component_of_wind_10m_era5_land^2)
cco1$prep<-ifelse(cco1$total_precipitation_era5_land<quantile(cco1$total_precipitation_era5_land)[3],"normal",
                  ifelse(cco1$total_precipitation_era5_land<quantile(cco1$total_precipitation_era5_land)[4],"hard","extreme"))

#sum(is.na(cco1)) #0

#only for summertime
cco1<-cco1%>%filter(month(matched_death_dt)%in%c(6:9))
#################################dlnm start################################


######1) only CCO 
rr<-data.frame()
#set parameters
par(mfrow=c(1,4))
varfun<-"ns"
vardf<-4
lagknots<-c(1,3) #among lag days

x<-cco1[,c("tmean","tmean_lag_1","tmean_lag_2","tmean_lag_3",
           "tmean_lag_4","tmean_lag_5")]
cb<-crossbasis(x,argvar=list(fun=varfun,df=vardf),arglag=list(knots=lagknots))

model<-clogit(CASE_YN~cb+strata(person_id),cco1) #only with temp
coef<-summary(model)$coef[1:(vardf*(length(lagknots)+1)),1]
vcov<-vcov(model)[1:(vardf*(length(lagknots)+1)),1:(vardf*(length(lagknots)+1))]

red<-crossreduce(cb,model)

predvar<-quantile(cco1$tmean,c(0:100)/100,na.rm=T)
b1<-onebasis(predvar,fun=varfun,knots=predvar[c(25,50,75)+1])
pred<-crosspred(b1,coef=coef(red),vcov=vcov(red),model.link="logit",at=predvar,cen=predvar[75])
plot(pred)
title(main="without adjustment")
rr<-rbind(rr,cbind(pred$allRRfit,pred$allRRlow,pred$allRRhigh)[99,])
#75 vs. 99 1.1108058 1.0535854 1.1711338


######2) CCO + clogit 보정
#temp+precipiation+surface_pressure+humidity+wind_speed+optical_depth_047+evi+pop_density+covid+sur_refl_b01


model2<-clogit(CASE_YN~cb+total_precipitation_era5_land+surface_pressure_era5_land+humidity+sur_refl_b01
               +wind_speed+optical_depth_047+evi+pop_density+covid+
                 strata(person_id),cco1) 
coef<-summary(model2)$coef[1:(vardf*(length(lagknots)+1)),1]
vcov<-vcov(model2)[1:(vardf*(length(lagknots)+1)),1:(vardf*(length(lagknots)+1))]

red<-crossreduce(cb,model2)

predvar<-quantile(cco1$tmean,c(0:100)/100,na.rm=T)
b1<-onebasis(predvar,fun=varfun,knots=predvar[c(25,50,75)+1])
pred<-crosspred(b1,coef=coef(red),vcov=vcov(red),model.link="logit",at=predvar,cen=predvar[75])
plot(pred)
title(main="with adjustment")

rr<-rbind(rr,cbind(pred$allRRfit,pred$allRRlow,pred$allRRhigh)[99,])

######3) CCO + GPS(ref=https://www.youtube.com/watch?v=mSueFNm_a9g)
#Adjust total_precipitation_era5_land+surface_pressure_era5_land+humidity+sur_refl_b01
#+wind_speed+optical_depth_047+evi+pop_density+covid
#+
#+
#temp and met joining
temp1120_tm<-temp1120%>%select(c(date,sgg_h,tmean))
tm1120<-temp1120_tm%>%left_join(met1120,by=c('date','sgg_h'))
tm1120<-tm1120%>%filter(!sgg_h%in%c(37430,39010,39020,23090,29010,31280,33040)) #except jeju

tm1120$humidity<-100-5*((tm1120$temperature_2m_era5_land-32)*(5/9)
                        -(tm1120$dewpoint_temperature_2m_era5_land-32)*(5/9))
tm1120$wind_speed<-sqrt(tm1120$u_component_of_wind_10m_era5_land^2+tm1120$v_component_of_wind_10m_era5_land^2)
tm1120$prep<-ifelse(tm1120$total_precipitation_era5_land<quantile(tm1120$total_precipitation_era5_land)[3],"normal",
                    ifelse(tm1120$total_precipitation_era5_land<quantile(tm1120$total_precipitation_era5_land)[4],"hard","extreme"))


#make moving average

library(tsModel)

#percentile temp
sgg<-unique(tm1120$sgg_h)%>%sort
temp<-sgg%>%lapply(function(x)tm1120 [tm1120 $sgg_h==x,])
temp<-temp%>%lapply(function(x){
  x<-x%>%data.table
  x<-x%>%mutate(temp_rank=round(rank(tmean)/length(tmean)*100,1))
})

tm1120<-bind_rows(temp)%>%arrange(sgg_h,date)
lag01<-tm1120%>%mutate(temp0_1=runMean(temp_rank,0:1,group=sgg_h),
                       prep0_1=runMean(total_precipitation_era5_land,0:1,group=sgg_h),
                       wind0_1=runMean(wind_speed,0:1,group=sgg_h),
                       hum0_1=runMean(humidity,0:1,group=sgg_h),
                       aod0_1=runMean(optical_depth_047,0:1,group=sgg_h))

cco1120<-cco1120%>%filter(month(matched_death_dt)%in%c(6:9))

cco2<-cco1120%>%left_join(lag01,by=c("matched_death_dt"="date","sgg_h"))

#make GPS
library(WeightIt)
library(splines)
# 
# x_fit<- glm (temp0_1 ~ prep0_1+hum0_1+
#              wind0_1+aod0_1, data = lag01 )
# 
# lag01$GPS<- dnorm (lag01$temp0_1, mean=x_fit$fitted.values, sd (x_fit$residuals))
# 
# ipw1 <- (1/lag01$GPS)
# ipw1[ipw1>=quantile(ipw1,0.95)]=quantile(ipw1,0.95)
# ipw1[ipw1<=quantile(ipw1,0.05)]=quantile(ipw1,0.05)
lag01<-lag01%>%filter(month%in%c(6:9))
#lag01_1<-na.omit(lag01)
ipw1<-weightit(temp0_1~prep0_1+hum0_1+ wind0_1+aod0_1, data = lag01 ,method="glm",stabilize = TRUE,criterion="p.mean")

#dlnm-moving average
model3<-clogit(CASE_YN~ns(temp0_1%>%unlist,4)+strata(person_id),cco2,weights=ipw1$weights) #only with temp


coef<-summary(model3)$coef[1:4]
vcov<-vcov(model3)[(1:4),(1:4)]


predvar<-quantile(cco2$tmean,c(0:100)/100,na.rm=T)
b1<-onebasis(predvar,fun=varfun,knots=predvar[c(25,50,75)+1])
pred<-crosspred(b1,coef=coef(red),vcov=vcov(red),model.link="logit",at=predvar,cen=predvar[75])
plot(pred)
lines(pred,col="red")
title(main="with GPS")
rr<-rbind(rr,cbind(pred$allRRfit,pred$allRRlow,pred$allRRhigh)[99,])
#75 vs. 99 1.1108058 1.0535854 1.1711338


######2) CCO + clogit 보정
#temp+precipiation+surface_pressure+humidity+wind_speed+optical_depth_047+evi+pop_density+covid+sur_refl_b01


model2<-clogit(CASE_YN~cb+total_precipitation_era5_land+surface_pressure_era5_land+humidity+sur_refl_b01
               +wind_speed+optical_depth_047+evi+pop_density+covid+
                 strata(person_id),cco2) 
coef<-summary(model2)$coef[1:(vardf*(length(lagknots)+1)),1]
vcov<-vcov(model2)[1:(vardf*(length(lagknots)+1)),1:(vardf*(length(lagknots)+1))]

red<-crossreduce(cb,model2)

predvar<-quantile(cco2$tmean,c(0:100)/100,na.rm=T)
b1<-onebasis(predvar,fun=varfun,knots=predvar[c(25,50,75)+1])
pred<-crosspred(b1,coef=coef(red),vcov=vcov(red),model.link="logit",at=predvar,cen=predvar[75])
plot(pred)
title(main="with adjustment")

rr<-rbind(rr,cbind(pred$allRRfit,pred$allRRlow,pred$allRRhigh)[99,])

######3) CCO + GPS(ref=https://www.youtube.com/watch?v=mSueFNm_a9g)
#Adjust total_precipitation_era5_land+surface_pressure_era5_land+humidity+wind_speed
#+
#+

model4<- glm (tmean ~ total_precipitation_era5_land+surface_pressure_era5_land+humidity
              +wind_speed, data = cco1 )
#strate적용 how~?


cco1$GPS<- dnorm (cco1$tmean, mean=model4$fitted.values, sd (model4$residuals))

ipw2 <- (1/cco1$GPS)






#################################Balance Check!##################################

##########love plot function
my.love.plot <- function (data, variables, ipw2){
  
  
  test <- data%>%select(variables)  
  test$ipw <- ipw1
  
  
  correlation1 <- c()
  correlation.w1 <- c()
  test<-test%>%as_tibble()
  for (i in 2:(length(colnames(test)))){
    
    correlation   <- wCorr::weightedCorr(x =  test [, "tmean"] %>%as.matrix, 
                                         y = test [, i]%>%as.matrix, 
                                         method = "Spearman")
    coreelation.w <- wCorr::weightedCorr(x = test [, "tmean"]%>%as.matrix, 
                                         y = test [, i]%>%as.matrix, 
                                         method = "Spearman", 
                                         weights = test [,"ipw"]%>%as.matrix)
    
    correlation1   <- c(correlation1,   abs (correlation))
    correlation.w1 <- c(correlation.w1, abs (coreelation.w))
    
  }
  
  A <- colnames (test[2:(dim(test)[2])])
  B <- correlation1   # unweighted
  C <- correlation.w1 # weighted
  
  corr.data <- data.frame( rbind (cbind (A, B, "Unadjusted"), cbind (A,C, "adjusted")))
  
  colnames (corr.data) <- c ("Variables", "abs.corr", "Adjusted.yn") 
  
  my.love.plot <- corr.data %>% ggplot (aes (x = Variables, y = as.numeric (abs.corr), color = Adjusted.yn)) +
    geom_point (position = position_dodge(0.4), size = 4) +
    coord_flip() +
    scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) +
    theme_classic () +
    labs (color  = "Status", y = "Abolute correlation", x = "Variables") +
    geom_hline (aes (yintercept = 0.1), linetype="dashed") 
  
  return (my.love.plot)
  
}


my.love.plot (cco1, variables = c('tmean','total_precipitation_era5_land','surface_pressure_era5_land','humidity','wind_speed','optical_depth_047','evi','pop_density','covid','sur_refl_b01'), 
              ipw2 = ipw1)


my.love.plot (cco1, variables = c('tmean','total_precipitation_era5_land','surface_pressure_era5_land','humidity','wind_speed'), 
              ipw2 = ipw2)
