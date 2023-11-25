# Code to conduct paper analyses.

# note: remove first and/or last point of each burst?


library(ggplot2)
library(dplyr)
library(rayshader)
library(patchwork)
library(lubridate)
library(viridis)
library(sf)

setwd("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights")

# read data

dat<-read.csv("analyses/tripdat_4_analyses_all.csv", h=T)
burst_summary<-read.csv("analyses/burst_summary_dat_all.csv", h=T)

dat$DateTime_AEDT<-ymd_hms(dat$DateTime_AEDT)
dat<-dat%>%filter(deployed_ID!="predeployment")     

table(dat$burstID, dat$class)

bad_ids<-c(
  "-1_1",  "08611854_01_30","08611854_01_31", "08611854_01_33", "08611854_01_34",
  "08611854_01_35",  "08611854_01_41", "08611854_01_42", "08611854_01_46",  "08611854_01_47",
  "08611854_01_48",  "08611854_01_49",  "08611854_01_50",  "08611854_01_51",
  "08611854_01_52", "08611854_01_55", "08611854_01_65",  "-1_66",  "08611854_02_79",
  "08611854_02_82",  "08611854_02_85",  "08611854_02_86",  "08611854_02_87",
  "08611854_02_91",  "08611854_02_92",  "08611854_03_99",  "08611854_03_101",
  "08611854_03_102",  "08611854_03_103",  "08611854_03_110",  "-1_111",
  "08611854_04_123",  "08611854_04_124",  "08611854_04_125", # good example of bobbing, then swamp
  "08611854_04_126",  "08611854_05_128",  "08611854_05_129",  "08611854_05_130",
  "08611854_05_131",  "08611854_05_132",  "-1_133",  "08611854_06_134", # same
  "08611854_06_136",  "08611854_06_139",  "08611854_06_140",  "08611854_06_153",
  "08611854_06_157",  "08611854_06_161",  "08611854_06_162",  "08611854_06_163",
  "08611854_06_164",  "08611854_06_166",  "08611854_06_167",  "08611854_06_168",
  "08611854_06_169",  "08611854_06_172",  "08611854_06_175",  "08611854_06_176",
  "08611854_06_178",  "08611854_06_183",  "08611854_06_185",  "08611854_06_188",
  "08611854_07_191",  "08611854_07_192",  "08611854_07_195",  "08611854_07_196",
  "08611854_08_199",  "08611854_08_200",  "08611854_08_201",  "41490936_01_16",
  "08611649_01_11",  "08611649_01_12",  "08611649_01_16",  "08611649_01_16",
  "08611649_01_18",  "08611649_01_19",  "08611649_01_21",  "08611649_01_25",
  "08611649_01_26",  "08611649_01_27",  "08611649_01_28",  "08611649_01_29",
  "08611649_01_37", # "08611649_01_39" allowed thru but includes a bit
  "08611649_01_40",  "08611649_01_41", "08611649_01_17")
  
dat<-dat%>%filter(!burstID %in% bad_ids) # could do extra check based on min/max burst pressure difference  to see if any missed

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
  # original method - 95% upper qunatile of pressure to set p0 for entire burst
  dat[dat$burstID==i,]$p0<-quantile(dat[dat$burstID==i,]$pres_pa, probs=0.95)
}
  
dat$pres_alt<-(-1*  # *-1 flips negative/positive values
                                    ((k*(dat$temp+273.15))/(m*g))*log(dat$pres_pa/dat$p0))



# do overall distribution for flying birds

mean(dat%>%filter(class %in% c("T", "L") & sit_fly=="fly")%>%pull(pres_alt), na.rm=T)

ggplot()+
  geom_vline(xintercept=0, colour='blue', size=1)+
  geom_vline(xintercept=10, colour='blue', size=0.5)+
  geom_vline(xintercept=20, colour='blue', size=0.5)+
  geom_vline(xintercept=4.209, colour='red',linetype='dotted', size=1)+
  geom_label(aes(x=6.9, y=4200, label="Mean\n4.21 m"),size=5, col='red')+
  geom_rect(aes(xmin=30, xmax=50, ymin=0, ymax=4300), fill='red', size=0.5, alpha=0.3)+
  geom_histogram(data=dat%>%filter(class %in% c("T", "L") & sit_fly=="fly"), aes(x=pres_alt),colour=1, binwidth = 1)+
  scale_x_continuous(breaks=seq(-5, 50, 5),minor_breaks=seq(-5, 50, 1), limits=c(-5, 50))+
  labs(x="Altitude (m)", y="Number 3D datapoints (Lat,Lon,Pressure)")+theme_bw()+
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold"))

tl_dat<-dat%>%filter(class %in% c("T", "L") & sit_fly=="fly" & burstID!="-1_93")

mean(tl_dat[tl_dat$pres_alt>0,]$pres_alt) # exclude negative vals
  
nrow(tl_dat[tl_dat$pres_alt<20,])/nrow(tl_dat)  
nrow(tl_dat[tl_dat$pres_alt<30,])/nrow(tl_dat)
nrow(tl_dat[tl_dat$pres_alt<10,])/nrow(tl_dat)
nrow(tl_dat[tl_dat$pres_alt<5,])/nrow(tl_dat)

# write out kmz of bird that crosses the island
sf3d<-dat%>%filter(burstID=="-1_93")%>%st_as_sf(coords = c("Longitude", "Latitude", "pres_alt"), crs = 4326, dim = "XYZ")
st_write(sf3d, "analyses/GIS/over_the_island_burst.kml")

# test GPS difference

sensor_comp<-rbind(data.frame(Sensor="gps", Altitude=tl_dat[tl_dat$alt<2000,]$alt),
      data.frame(Sensor="pressure", Altitude=tl_dat[tl_dat$alt<2000,]$pres_alt))

sensor_comp%>%group_by(Sensor)%>%summarise(mean_alt=mean(Altitude))

ggplot(data=sensor_comp%>%filter(Altitude<90))+
  geom_hline(yintercept=0, colour='blue', size=1)+
  geom_jitter(aes(x=Sensor, y=Altitude), height=0, width=0.4, alpha=0.1, shape=16, size=1.5)+
  geom_violin(aes(x=Sensor, y=Altitude, colour=Sensor), fill=NA, size=1)+
  scale_y_continuous(breaks=seq(-50, 90, 10),minor_breaks=seq(-50, 90, 5))+
  labs(y="Altitude (m)")+theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),legend.position="none")
  
# test windspeed

ggplot(data=tl_dat, aes(x=wind_speed, y=pres_alt))+
  geom_point(alpha=0.1)+geom_smooth(method='lm')+facet_wrap(~b_w_class)

av_dat<-tl_dat%>%group_by(burstID)%>%summarise(mn_alt=mean(pres_alt), mn_w_speed=mean(wind_speed, na.rm=T))

ggplot()+
  geom_hline(yintercept=0, colour='blue', size=1)+
  geom_point(data=tl_dat, aes(x=wind_speed, y=pres_alt), alpha=0.05)+
  geom_smooth(data=av_dat, aes(x=mn_w_speed, y=mn_alt), method='lm', colour='red', fill='pink')+
  scale_y_continuous(breaks=seq(-4,24,2), limits=c(-4, 24))+
  scale_x_continuous(breaks=seq(0,15,1))+
  labs(y="Altitude (m)", x="Windspeed (m/s)")+theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),legend.position="none")

 mjyuy# compare pressure and GPS altitude

ggplot(data=dat%>%filter(class %in% c("T", "L")& sit_fly=="fly"))+geom_density(aes(y=alt, x=pres_alt))



### example plots for talk

# pressure + dyn soaring
ggplot(data=dat[dat$burstID=="08611854_02_71",])+geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa), size=1)+scale_y_continuous(breaks=c(100800, 100850, 100900, 100950, 101000, 101050))+
  labs(y="Pressure (mb)", x="Time")+theme(axis.text=element_text(size=12),
                                                                                                                                                                           axis.title=element_text(size=14,face="bold"))

# 3d plot with rayshader - example track ALSO check 08611854_06_146 for loopy example and "08611854_04_120" for Alight
p1<-ggplot(data=dat[dat$burstID=="41490936_01_24",])+
  geom_point(aes(x=Longitude, y=Latitude, colour=pres_alt))+scale_color_viridis()+
  labs(x="Longitude", y="Latitude", colour="   Altitude (m)")

plot_gg(p1, height=3, width=8, pointcontract = 0.7)

#comparison with GPS
p1<-ggplot(data=dat[dat$burstID=="41490936_01_24",])+
  geom_point(aes(x=Longitude, y=Latitude, colour=alt))+scale_color_viridis()+
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