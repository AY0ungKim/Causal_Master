evi<-fread("C:/Users/82105/Dropbox (개인용)/대기오염 모델링 자료이용 연구/DATA/NDVI_EVI/EVI_SGG.csv")
evi<-evi%>%filter(year%in%c(2002:2020))


yearly_dat<-evi%>%group_by(sgg_cd)%>%summarise(value=mean(mean))

#data import shp
geo<-st_read("C:/Users/82105/AiMS-CREATE Dropbox/지도/data_badn_sgg, sido_2021/bnd_sigungu_00_2021_2021_2Q.shp")
#join

geo$SIGUNGU_CD=as.numeric(paste0(substr(geo$SIGUNGU_CD,1,4),0))


geo_dat<-geo%>%left_join(yearly_dat,by=c("SIGUNGU_CD"="sgg_cd"))

geo_dat<-na.omit(geo_dat)

sum(is.na(geo_dat))


#plot-NDVI
g1=ggplot(geo_dat) + 
  geom_sf(aes(fill = value))+
  scale_fill_gradient(name = "",low = '#E3F2FD', high ='navy')+
  theme_bw()+
  ggtitle("EVI")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))

#plot-EVI
g2=ggplot(geo_dat) + 
  geom_sf(aes(fill = mean_evi))+
  scale_fill_gradient(name = "", low = '#E8F5E9', high ='#388E3C')+
  theme_bw()+
  ggtitle("EVI(standarized)")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))

#plot-BARE
g3=ggplot(geo_dat) + 
  geom_sf(aes(fill = bare))+
  scale_fill_gradient(name = "", low = '#E3F2FD', high ='navy')+
  theme_bw()+
  ggtitle("Bare Coverfraction(standarized)")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))

#plot-crops
g4=ggplot(geo_dat) + 
  geom_sf(aes(fill =crops))+
  scale_fill_gradient(name = "", low = '#E3F2FD', high ='navy')+
  theme_bw()+
  ggtitle("Crops Coverfraction(standarized)")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))

#plot-crops
g5=ggplot(geo_dat) + 
  geom_sf(aes(fill =tree))+
  scale_fill_gradient(name = "", low = '#E3F2FD', high ='navy')+
  theme_bw()+
  ggtitle("Tree Coverfraction(standarized)")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))

#plot-crops
g6=ggplot(geo_dat) + 
  geom_sf(aes(fill =water))+
  scale_fill_gradient(name = "", low = '#E3F2FD', high ='navy')+
  theme_bw()+
  ggtitle("Water Permanent Coverfraction(standarized)")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))


#plot-sr
g7=ggplot(geo_dat) + 
  geom_sf(aes(fill =mean_surface))+
  scale_fill_gradient(name = "", low = '#E3F2FD', high ='navy')+
  theme_bw()+
  ggtitle("Surface Pressure(standarized)")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))



#plot-aod47
g8=ggplot(geo_dat) + 
  geom_sf(aes(fill =aod_47))+
  scale_fill_gradient(name = "", low = '#E3F2FD', high ='navy')+
  theme_bw()+
  ggtitle("Aerosol Optical Depth in the MODIS Blue band(Standarized)")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))


#plot-aod47
g9=ggplot(geo_dat) + 
  geom_sf(aes(fill =aod_55))+
  scale_fill_gradient(name = "", low = '#E3F2FD', high ='navy')+
  theme_bw()+
  ggtitle("Aerosol Optical Depth in the MODIS Green band(Standarized)")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))


grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,g9, ncol=3,nrow=3)
