library(dplyr)
library(tidyr)
library(lubridate)
library(oceanwaves)
library(readxl)
library(sf)
library(tmap)


setwd("C:/Users/mmil0049/OneDrive - Monash University/fieldwork/Seadragon atsea deployment")

# read data

dat<-NULL
for(i in list.files("05mar_datapull"))
{
d1<-read.csv(paste0("05mar_datapull/", i), skip=1)
dat<-rbind(dat, data.frame(d1, ID=substr(i, 3, 6)))
}

dat<-dat%>%dplyr::select(-c('Number.of.SVs','Number.of.SVs', 'Altitude..metres.',
                            'Sensor.Temperature', 'Time.Underwater', 'Activity.Count'))
#assign species
dat$sp<-"BUAL"
dat[dat$ID %in%c(1183,1200, 1196, 1316, 1318),]$sp<-"SHAL"

#convert datetime to local
dat$UTC.Timestamp<-ymd_hms(dat$UTC.Timestamp, tz="UTC")
dat$DateTime_AEDT<-with_tz(dat$UTC.Timestamp, tz="Australia/Sydney")
#remove some duplicates
dat<-dat[-which(duplicated(paste(dat$ID, dat$DateTime_AEDT))),]


# remove pre-deployment data - manually specify for each device
dat<-dat%>%filter(!(ID=="1183" & DateTime_AEDT<ymd_hms("2023-12-20 13:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1200" & DateTime_AEDT<ymd_hms("2023-12-20 14:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1196" & DateTime_AEDT<ymd_hms("2023-12-20 15:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1313" & DateTime_AEDT<ymd_hms("2024-02-21 11:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1321" & DateTime_AEDT<ymd_hms("2024-02-21 12:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1314" & DateTime_AEDT<ymd_hms("2024-02-21 14:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1315" & DateTime_AEDT<ymd_hms("2024-02-21 14:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1316" & DateTime_AEDT<ymd_hms("2024-02-22 10:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1318" & DateTime_AEDT<ymd_hms("2024-02-22 11:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1322" & DateTime_AEDT<ymd_hms("2024-02-22 11:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1323" & DateTime_AEDT<ymd_hms("2024-02-22 14:00:00", tz="Australia/Sydney")))

# make burst ID
dat$burstID<-as.numeric(as.factor(paste(dat$ID, year(dat$DateTime_AEDT), month(dat$DateTime_AEDT), day(dat$DateTime_AEDT), hour(dat$DateTime_AEDT),
                                        dat$Latitude, dat$Longitude, dat$Heading, dat$Battery.Voltage)))
table(table(dat$burstID))

# write out compiled coord file
#write.csv(dat%>%group_by(ID, burstID)%>%summarise_all(first), "compiled_05mar_coords.csv", quote=F, row.names=F)

id_l<-data.frame(table(dat$burstID))
dat<-dat[dat$burstID%in% id_l[id_l$Freq>60,]$Var1,] # remove bursts less than 60 secs 


dat_sf<-st_as_sf(dat%>%group_by(ID, burstID)%>%summarise_all(first), coords=c("Longitude", "Latitude"), crs=4326)
dat_sf<-dat_sf%>%arrange(ID, DateTime_AEDT)

linez<-dat_sf %>%group_by(ID) %>%
  dplyr::summarize(do_union=F) %>%  # do_union=FALSE doesn't work as well
  st_cast("LINESTRING") 

tmap_mode("view")
tm_shape(linez)+tm_lines(col=as.character("ID"))


# get speed histogram
sp_df<-dat%>%group_by(sp, burstID)%>%summarise(speed=first(Speed..km.h.))

sp_df%>%filter((speed/3.6)>4)%>%
  summarise(mean_s=mean(speed/3.6), sd_s=sd(speed/3.6))
#sp    mean_s  sd_s
#<chr>  <dbl> <dbl>
#  1 BUAL    13.5  5.05
#  2 SHAL    14.0  4.93
 
ggplot()+
  geom_histogram(data=sp_df,
                 aes(x=speed, fill=ifelse(speed<15, "blue", "red")), colour=1, binwidth=1, boundary=0)+
  scale_x_continuous(breaks=0:150)+theme_bw()+ylab("Count of GPS datapoints")+xlab("speed km/h")+
  theme(legend.position="none")+facet_wrap(~sp, nrow=2)

ggplot()+
  geom_histogram(data=sp_df,
                 aes(x=speed/3.6, fill=ifelse((speed/3.6)<4, "blue", "red")), colour=1, binwidth=1, boundary=0)+
  scale_x_continuous(breaks=0:32)+theme_bw()+ylab("Count of GPS datapoints")+xlab("speed m/s")+
  theme(legend.position="none")+facet_wrap(~sp, nrow=2)


zc_summary<-data.frame(dat%>%group_by(sp, ID, burstID)%>%summarise(speed=first(Speed..km.h.), temp=first(Device.Temperature)),
                       Hsig=NA,Hmean=NA, Tmean=NA, Tsig=NA)

for (i in unique(dat$burstID))
{
  possibleError <-tryCatch(
    waveStatsZC(dat[dat$burstID==i,]$Pressure, 1,),
    error=function(e) e)
  if(!inherits(possibleError, "error")){
  d1<-waveStatsZC(dat[dat$burstID==i,]$Pressure, 1,)
  zc_summary[zc_summary$burstID==i,]$Hsig=d1$Hsig
  zc_summary[zc_summary$burstID==i,]$Hmean=d1$Hmean
  zc_summary[zc_summary$burstID==i,]$Tmean=d1$Tmean
  zc_summary[zc_summary$burstID==i,]$Tsig=d1$Tsig
  }else{next}
}
# 11/03 - only 11 NAs


ggplot(data=zc_summary)+geom_histogram(aes(x=speed))
ggplot(data=zc_summary)+geom_histogram(aes(x=temp))  # a few minus values
ggplot(data=zc_summary)+geom_histogram(aes(x=Hsig))
ggplot(data=zc_summary)+geom_histogram(aes(x=Hmean))
ggplot(data=zc_summary)+geom_histogram(aes(x=Tmean))
ggplot(data=zc_summary)+geom_histogram(aes(x=Tsig))

zc_summary[zc_summary$temp<0,]

#summary
na.omit(zc_summary)%>%filter((speed/3.6)>4)%>%select(-c('burstID', 'ID', 'sp'))%>%summarise_all(list(mean=mean, sd=sd))

ggplot(data=zc_summary)+geom_point(aes(x=speed, y=Hmean))
ggplot(data=zc_summary[zc_summary$Hmean<5,])+geom_point(aes(x=speed, y=Hmean), alpha=0.2, shape=1) # looks good to cut on speed alone

#### Calculation of flight height using dynamic soaring method ####

# barometric formula (Berberan Santos et al. 1997)
#h=((k*T)/(m*g))*ln(p/p0)
k=8.31432
m=0.0289644
g=9.80665

dat$p0<-0
dat$alt_DS<-NA
for (i in unique(dat$burstID))
{
  # original method - 95% upper quantile of pressure to set p0 for entire burst
  dat[dat$burstID==i,]$p0<-quantile(dat[dat$burstID==i,]$Pressure, probs=0.95)
}

dat$alt_DS<-(-1*  # *-1 flips negative/positive values
               ((k*(dat$Device.Temperature+273.15))/(m*g))*log(dat$Pressure/dat$p0))

#### ^^^ ####

dat%>%filter((Speed..km.h./3.6)>4)%>%group_by(sp)%>%summarise(mn_alt=mean(alt_DS), sd_alt=sd(alt_DS), median=median(alt_DS),
                                        min=min(alt_DS), max=max(alt_DS),
                                        q25=quantile(alt_DS, 0.25), q75=quantile(alt_DS, 0.75))

ggplot(data=dat%>%filter((Speed..km.h./3.6)>4))+geom_density(aes(x=alt_DS, colour=sp), fill=NA)+
  theme_bw()+geom_vline(xintercept = 0, linetype='dotted')+scale_x_continuous(limits=c(-20, 100))
