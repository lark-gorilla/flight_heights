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

dat$pres_alt<-NA
for (i in unique(dat$burstID))
{
  # get 95% quantile of burst - want upper as high pressure = low altitude 
  upper_95<-quantile(burst_pres<-dat[dat$burstID==i,]$pres_pa,
           probs=0.95)
  
  
  dat[dat$burstID==i,]$pres_alt<-(-1*
    ((k*(dat[dat$burstID==i,]$temp+273.15))/(m*g))*log(dat[dat$burstID==i,]$pres_pa/upper_95)) # *-1 flips negative/positive values
}

# do overall distribution for flying birds




# 3d plot

ggplot(data=dat[dat$burstID==i,])+geom_line(aes(x=tim_UTC, y=pres_pa, group=1))+geom_point(aes(x=tim_UTC, y=pres_pa), size=1)+geom_hline(yintercept=upper_95, col='red')+scale_y_reverse()

p1<-ggplot(data=dat[dat$burstID==i,])+geom_point(aes(x=X, y=Y, colour=pres_alt))

plot_gg(p1)

# 3d plot to google earth

temp1<-dat
temp1$tim_UTC<-ymd_hms(temp1$tim_UTC)
temp1$DESCRIPTION=temp1$class

#sf4d<-temp1%>%st_as_sf(coords = c("X", "Y", "pres_alt", "tim_UTC"), crs = 4326, dim = "XYZM")

sf3d<-temp1[c("DESCRIPTION","X", "Y", "pres_alt")]%>%
  st_as_sf(coords = c("X", "Y", "pres_alt"), crs = 4326, dim = "XYZ")

st_write(sf3d, "analyses/GIS/google_earth_3d_vis.kml", driver='kml')

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