# Code to conduct paper analyses.

# note: remove first and/or last point of each burst?


library(ggplot2)
library(dplyr)
library(rayshader)
library(patchwork)
library(lubridate)

setwd("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights")

# read data

dat<-read.csv("analyses/tripdat_4_analyses.csv", h=T)
burst_summary<-read.csv("analyses/burst_summary_dat.csv", h=T)

#attrib burst summary
dat<-left_join(dat, burst_summary[c("burstID", "class")], by="burstID")
dat$DateTime_AEDT<-ymd_hms(dat$DateTime_AEDT)
dat<-dat%>%filter(deployed_ID!="predeployment")     


table(dat$burstID, dat$class)

## calc flight heights, do all but only ware about T and A burst classes

# barometric formula (Berberan Santos et al. 1997)
#h=((k*T)/(m*g))*ln(p/p0)
k=8.31432
m=0.0289644
g=9.80665

# for p0 per burst we want to identify pressure level at sea level, we can assume waves
# so finding upper 95th? quantile could represent burst mean sea level, negative values 
# in resultant flight height would be wave troughs. A more precise method would find the 
# max pressure point of each osscillation (after a short window smooth to remove error),
# and then take the mean of these points, you could even run a gam through these if we
# expect local pressure changes to change over the burst - possible on transiting bursts.
# Question we need answered if if birds always meet the surface during their soaring oscillation
# Anchoring each oscillaition to the sea surface removes negative values but probably 
# overestimates.. hmm also look at a small moving window to clean up 

dat$index<-1:nrow(dat)
dat$p0<-0

dat$pres_alt<-NA
for (i in unique(dat$burstID))
{
  
  #split burst in 5 (= 1 minute sub-bursts) then 
  # get 95% quantile of sub-burst pressure to use as p0 for sub burst
  
  splitter<-data.frame(indy=dat[dat$burstID==i,]$index, split=cut(dat[dat$burstID==i,]$index, 5))
  
  for(j in unique(splitter$split))
  {
    dat[dat$index%in%splitter[splitter$split==j,]$indy,]$p0<-
      quantile(dat[dat$index%in%splitter[splitter$split==j,]$indy,]$pres_pa, probs=0.75) # select quantile
  }
  # temp model
  mod_dat<-dat[dat$burstID==i,]
  mod_dat$pres_above<-ifelse(mod_dat$pres_pa>mod_dat$p0, mod_dat$pres_pa, NA)
  m1<-loess(pres_above~index, data=mod_dat, span=0.75, control = loess.control(surface = "direct")) 
  mod_dat$pred<-predict(m1, mod_dat)
  m2<-loess(pres_above~index, data=mod_dat, span=0.5, control = loess.control(surface = "direct")) 
  mod_dat$pred2<-predict(m2, mod_dat)
  
  dat[dat$burstID==i,]$pres_alt<-(-1*  # *-1 flips negative/positive values
                                    ((k*(mod_dat$temp+273.15))/(m*g))*log(mod_dat$pres_pa/mod_dat$pred)) # using pred as p0

  p2<-ggplot(data=mod_dat)+geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1))+
    geom_point(aes(x=DateTime_AEDT, y=pres_pa), size=1)+geom_line(aes(x=DateTime_AEDT, y=p0), col='red')+
    geom_line(aes(x=DateTime_AEDT, y=pred), col='green')+
    geom_line(aes(x=DateTime_AEDT, y=pred2), col='blue')+scale_y_reverse()+labs(title =unique(mod_dat$class))
  
  p3<-ggplot(data=dat[dat$burstID==i,])+geom_path(aes(x=Longitude, y=Latitude, group=1))+
    geom_point(aes(x=Longitude, y=Latitude, colour=pres_alt, shape=sit_fly))
            
  print(p2/p3)
  print(i)
  print(burst_summary%>%filter(burstID==i))
  readline("")
   
}

  dat[dat$burstID==i,]$pres_alt<-(-1*
    ((k*(dat[dat$burstID==i,]$temp+273.15))/(m*g))*log(dat[dat$burstID==i,]$pres_pa/dat[dat$burstID==i,]$p0)) # *-1 flips negative/positive values




ggplot(data=dat[dat$burstID=="08611854_02_71",])+geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa), size=1)+scale_x_date(date_breaks="minutes")
# do overall distribution for flying birds


ggplot(data=dat%>%filter(class %in% c("T", "L", "A")& sit_fly=="fly"))+geom_histogram(aes(x=pres_alt), binwidth = 1)

ggplot(data=dat%>%filter(class %in% c("T", "L", "A")& sit_fly=="fly"))+geom_density(aes(x=pres_alt, fill=class))




### example plots for talk

# pressure + dyn soaring
ggplot(data=dat[dat$burstID=="08611854_02_71",])+geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa), size=1)+scale_y_continuous(breaks=c(100800, 100850, 100900, 100950, 101000, 101050))+
  labs(y="Pressure (mb)", x="Time")+theme(axis.text=element_text(size=12),
                                                                                                                                                                           axis.title=element_text(size=14,face="bold"))

# 3d plot with rayshader - example track
p1<-ggplot(data=dat[dat$burstID=="08611854_02_74",])+
  geom_point(aes(x=X, y=Y, colour=pres_alt))+scale_color_viridis_b()+
  labs(x="Longitude", y="Latitude", colour="   Altitude (m)")

plot_gg(p1, height=3, width=8, pointcontract = 0.7)

#comparison with GPS
p1<-ggplot(data=dat[dat$burstID=="08611854_02_74",])+
  geom_point(aes(x=X, y=Y, colour=alt))+scale_color_viridis_b()+
  labs(x="Longitude", y="Latitude", colour="   Altitude (m)")

plot_gg(p1, height=3, width=8, pointcontract = 0.7)

# 3d plot to google earth

temp1<-dat
temp1$DESCRIPTION=temp1$class
temp1$TIME<-ymd_hms(temp1$DateTime_AEDT)

sf4d<-temp1[c("DESCRIPTION","Longitude", "Latitude", "pres_alt", "TIME")]%>%
  st_as_sf(coords = c("Longitude", "Latitude", "pres_alt", "TIME"), crs = 4326, dim = "XYZM")

sf3d<-temp1[c("DESCRIPTION","Longitude", "Latitude", "pres_alt", "TIME")]%>%
  st_as_sf(coords = c("Longitude", "Latitude", "pres_alt"), crs = 4326, dim = "XYZ")

#st_write(sf3d, "analyses/GIS/google_earth_3d_vis.kml", driver='kml')

#st_write(sf3d%>%filter(DESCRIPTION=="T"), "analyses/GIS/google_earth_3d_vis_transit.kml", driver='kml')
#st_write(sf3d%>%filter(DESCRIPTION=="L"), "analyses/GIS/google_earth_3d_vis.kml_loop.kml", driver='kml')
#st_write(sf3d%>%filter(DESCRIPTION=="A"), "analyses/GIS/google_earth_3d_vis.kml_alight.kml", driver='kml')

# find good bursts

T_b<-dat%>%filter(class=='T')%>%pull(burstID)%>%unique()

for(j in T_b){
p2<-ggplot(data=dat[dat$burstID==j,])+geom_line(aes(x=tim_UTC, y=pres_pa, group=1))+
  geom_point(aes(x=tim_UTC, y=pres_pa), size=1)+scale_y_reverse()
p3<-ggplot(data=dat[dat$burstID==j,])+geom_line(aes(x=X, y=Y, group=1))+
  geom_point(aes(x=X, y=Y, colour=pres_alt))
p2/p3             
print(p2/p3)
print(j)
readline("")
}

#"08611854_01_64"
"08611854_02_68"
"08611854_02_69"
"08611854_02_70"
"08611854_02_71"

p2<-ggplot(data=dat[dat$burstID==j,])+geom_line(aes(x=tim_UTC, y=pres_pa, group=1))+
  geom_point(aes(x=tim_UTC, y=pres_pa), size=1)+geom_hline(yintercept=upper_95, col='red')+scale_y_reverse()

# test export as kmz

temp1<-dat
temp1$tim_UTC<-ymd_hms(temp1$tim_UTC)

#sf4d<-temp1%>%st_as_sf(coords = c("X", "Y", "pres_alt", "tim_UTC"), crs = 4326, dim = "XYZM")

sf3d<-temp1%>%st_as_sf(coords = c("X", "Y", "pres_alt"), crs = 4326, dim = "XYZ")

st_write(sf3d, "C:/Users/mmil0049/Downloads/temp.kml")