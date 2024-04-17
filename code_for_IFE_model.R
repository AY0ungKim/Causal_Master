#230414
#Referece:https://jamanetwork.com/journals/jamanetworkopen/article-abstract/2815655
#Use Korean mortality data

library(phtt)
library(splines)
library(dplyr)
library(data.table)

################################ data preparation ###############################

#Load the Korean mortality data
setwd("C:/Users/82105/Dropbox (개인용)/발표/MCC/data")
cco1521<-fread("230910_death0021_dup.csv")%>%filter(matched_death_dt<=as.Date("2020-12-31")&matched_death_dt>=as.Date("2020-01-01")&CASE_YN==1)%>%arrange(person_id,matched_death_dt)
temp<-fread("1521_변수들_행정동별_평균.csv")%>%filter(date<=as.Date("2020-12-31")&date>=as.Date("2020-01-01"))%>%select(c("date","sgg_h","temperature_2m_era5_land","PM2_5"))

#Join with pm data and temperature
sgg=unique(cco1521$sgg_h)
key=data.frame(date=rep(seq(as.Date("2020-01-01"),as.Date("2020-12-31"),by="day"),each=length(unique(sgg))),
               sgg_h=rep(sgg,366))

data=cco1521%>%group_by(sgg_h,death_dt)%>%summarise(Mortality=n())
data=data%>%left_join(temp,by=c("death_dt"="date","sgg_h"))
data$death_dt=as.Date(data$death_dt)


data_key=key%>%left_join(data,by=c("date"="death_dt","sgg_h"))


#NA (mean for neighbors)

replace_na_with_mean_of_neighbors <- function(x) {
  x[1] <- x[2]
  x[length(x)] <- x[length(x) - 1]
  for (i in 2:(length(x) - 1)) {
    if (is.na(x[i])) {
      x[i] <- mean(c(x[i - 1], x[i + 1]), na.rm = TRUE)
    }
  }
  return(x)
}

# 변수 x의 결측치를 앞행과 뒤행의 평균으로 대체
data_key$temperature_2m_era5_land <- replace_na_with_mean_of_neighbors(data_key$temperature_2m_era5_land)
data_key$PM2_5 <- replace_na_with_mean_of_neighbors(data_key$PM2_5)
data_key$Mortality[is.na(data_key$Mortality)] <- 0
data=data_key

rm(cco1521,temp);gc()

###Make the variables

#Total.rate.dif: mortality difference from t-1 to t day
#PM25.lag.dif: PM2.5 difference from t-1 to t day
#TMP.lag.dif: Temperature difference from t-1 to t day
data$Total.rate.dif=0
data$PM25.lag.dif=0
data$TMP.lag.dif=0
for (i in 2:nrow(data)){
  data$Total.rate.dif[i]=data$Mortality[i]-data$Mortality[i-1]
  data$PM25.lag.dif[i]=data$PM2_5[i]-data$PM2_5[i-1]
  data$TMP.lag.dif[i]=data$temperature_2m_era5_land[i]-data$temperature_2m_era5_land[i-1]
}

head(data)


#naming for code
data=data%>%select(c("sgg_h","date","Total.rate.dif","PM25.lag.dif","TMP.lag.dif"))
names(data)[c(1,2)]=c("countycode","time")

save(data,file="Mortality_PM_IFE_21.RData")

gc()





#################################Fit IFE model###############################
load("Mortality_PM_IFE_21.RData")

head(data)


n.date <- length(unique(data$time)) # number of total days -1
n.county <- length(unique(data$countycode)) # number of spatial units in each study region
Total.rate <- matrix(data$Total.rate.dif, nrow=n.date, ncol=n.county) # mortality rate after first-order differencing

PM25 <- matrix(data$PM25.lag.dif/10, n.date, n.county)
Tmean.spline <- ns(data$TMP.lag.dif, df=5)
Tmean.ns1 <- matrix(Tmean.spline[,1], n.date, n.county)
Tmean.ns2 <- matrix(Tmean.spline[,2], n.date, n.county)
Tmean.ns3 <- matrix(Tmean.spline[,3], n.date, n.county)
Tmean.ns4 <- matrix(Tmean.spline[,4], n.date, n.county)
Tmean.ns5 <- matrix(Tmean.spline[,5], n.date, n.county)


## main model for PM2.5
model.PM25 <- Eup(Total.rate ~ PM25 + Tmean.ns1 + Tmean.ns2 + Tmean.ns3 + Tmean.ns4 +
                    Tmean.ns5, additive.effects = "individual") 


summary(model.PM25)
#BAD MODEL...
