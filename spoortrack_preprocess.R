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
for(i in list.files("18nov_datapull"))
{
d1<-read.csv(paste0("18nov_datapull/", i), skip=1)
dat<-rbind(dat, data.frame(d1, ID=substr(i, 3, 6)))
}

dat<-dat%>%dplyr::select(-c('Number.of.SVs', 'Altitude..metres.',
                            'Sensor.Temperature', 'Time.Underwater', 'Activity.Count'))

# remove loggers that failed to turn on
dat<-filter(dat, !ID%in%c(1202, 1319))

#assign species
dat$sp<-"NA"
dat[dat$ID %in%c(1313,1314,1315,1321,1322,1323),]$sp<-"BUAL"
dat[dat$ID %in%c(1183,1200, 1196, 1316, 1318),]$sp<-"WCAL"
dat[dat$ID %in%c(1428),]$sp<-"SHAL"
dat[dat$ID %in%c(1427, 1437, 1433, 1443),]$sp<-"WAAL"
dat[dat$ID %in%c(1324,1424,1425,1430,1431,1432),]$sp<-"NGPE"
dat[dat$ID %in%c(1426),]$sp<-"SGPE"
dat[dat$ID %in%c(1320),]$sp<-"IYNA"
dat[dat$ID %in%c(1326, 1327),]$sp<-"BBAL"

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
dat<-dat%>%filter(!(ID=="1327" & DateTime_AEDT<ymd_hms("2024-03-12 09:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1326" & DateTime_AEDT<ymd_hms("2024-03-12 15:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1320" & DateTime_AEDT<ymd_hms("2024-04-17 13:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1430" & DateTime_AEDT<ymd_hms("2024-06-20 13:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1324" & DateTime_AEDT<ymd_hms("2024-06-28 11:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1426" & DateTime_AEDT<ymd_hms("2024-09-11 10:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1443" & DateTime_AEDT<ymd_hms("2024-09-11 11:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1431" & DateTime_AEDT<ymd_hms("2024-09-11 13:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1432" & DateTime_AEDT<ymd_hms("2024-09-11 13:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1424" & DateTime_AEDT<ymd_hms("2024-09-11 16:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1428" & DateTime_AEDT<ymd_hms("2024-10-21 14:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1425" & DateTime_AEDT<ymd_hms("2024-10-22 10:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1433" & DateTime_AEDT<ymd_hms("2024-10-22 12:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1437" & DateTime_AEDT<ymd_hms("2024-10-22 12:00:00", tz="Australia/Sydney")))
dat<-dat%>%filter(!(ID=="1427" & DateTime_AEDT<ymd_hms("2024-10-22 12:00:00", tz="Australia/Sydney")))

# make burst ID
dat$burstID<-as.numeric(as.factor(paste(dat$ID, year(dat$DateTime_AEDT), month(dat$DateTime_AEDT), day(dat$DateTime_AEDT), hour(dat$DateTime_AEDT),
                                        dat$Latitude, dat$Longitude, dat$Heading, dat$Battery.Voltage)))
table(table(dat$burstID))

# write out compiled coord file
#write.csv(dat%>%group_by(ID, burstID)%>%summarise_all(first), "compiled_18nov_coords.csv", quote=F, row.names=F)

id_l<-data.frame(table(dat$burstID))
dat<-dat[dat$burstID%in% id_l[id_l$Freq>299,]$Var1,] # remove bursts less than 299 secs 

#remove first row of each burst as first pressure reading sometimes iffy
dat<-dat%>%group_by(burstID)%>%slice(-1)%>%as.data.frame()

dat_sf<-st_as_sf(dat%>%group_by(ID, burstID)%>%summarise_all(first), coords=c("Longitude", "Latitude"), crs=4326)
dat_sf<-dat_sf%>%arrange(ID, DateTime_AEDT)

linez<-dat_sf %>%group_by(ID, sp) %>%
  dplyr::summarize(do_union=F) %>%  # do_union=FALSE doesn't work as well
  st_cast("LINESTRING")%>%st_wrap_dateline() 

tmap_mode("view")
tm_shape(linez)+tm_lines(col=as.character("ID"))

tm_shape(linez)+tm_lines(col=as.character("sp"))


# get speed histogram
sp_df<-dat%>%group_by(sp, burstID)%>%summarise(speed=first(Speed..km.h.))

sp_df%>%filter((speed/3.6)>4)%>%
  summarise(mean_s=mean(speed/3.6), sd_s=sd(speed/3.6))
#sp    mean_s  sd_s  -- remember includes some landing bursts for some
#<chr>  <dbl> <dbl>
#1 BBAL    14.8  5.51
#2 BUAL    13.5  4.94
#3 IYNA    11.5  3.88
#4 NGPE    14.7  5.50
#5 SGPE    16.3  5.73
#6 SHAL    13.7  4.71
#7 WAAL    15.0  4.99
#8 WCAL    14.2  4.85
 

ggplot()+
  geom_histogram(data=sp_df,
                 aes(x=speed/3.6, fill=ifelse((speed/3.6)<4, "blue", "red")), colour=1, binwidth=1, boundary=0)+
theme_bw()+ylab("Count of GPS datapoints")+xlab("speed m/s")+
  theme(legend.position="none")+facet_wrap(~sp, nrow=2, scales='free_y')

ggplot()+
  geom_density(data=sp_df%>%filter((speed/3.6)>4),
                 aes(x=speed/3.6, colour=sp))+
 theme_bw()+ylab("Count of GPS datapoints")+xlab("speed m/s")+labs(colour='species')

# Use zero crossing analyses to summarize pressure values per burst.
# Hoping this combined with speed/temp can classify sitting flying bursts.

#make 1 min(almost) sub-bursts
dat<-dat%>%group_by(burstID)%>%mutate(sub_burst = cut(UTC.Timestamp, breaks = 5, labels=c("a", "b", "c", "d", "e" )))
dat$sub_burst<-paste(dat$burstID, dat$sub_burst, sep="_")

zc_summary<-data.frame(dat%>%group_by(sp, ID, sub_burst)%>%summarise(speed=first(Speed..km.h.), temp=first(Device.Temperature), Pvar=var(Pressure)),
                       Hsig=NA,Hmean=NA, Tmean=NA, Tsig=NA)

for (i in unique(dat$sub_burst))
{
  possibleError <-tryCatch(
    waveStatsZC(dat[dat$sub_burst==i,]$Pressure, 1,),
    error=function(e) e)
  if(!inherits(possibleError, "error")){
  d1<-waveStatsZC(dat[dat$sub_burst==i,]$Pressure, 1,)
  zc_summary[zc_summary$sub_burst==i,]$Hsig=d1$Hsig
  zc_summary[zc_summary$sub_burst==i,]$Hmean=d1$Hmean
  zc_summary[zc_summary$sub_burst==i,]$Tmean=d1$Tmean
  zc_summary[zc_summary$sub_burst==i,]$Tsig=d1$Tsig
  }else{next}
}
# 11/03 - only 11 NAs

#check and remove extreme values (temp, Hsig etc)

ggplot(data=zc_summary)+geom_histogram(aes(x=speed))
ggplot(data=zc_summary)+geom_histogram(aes(x=Pvar), bins=50)
ggplot(data=zc_summary)+geom_density(aes(x=temp, col=ifelse((speed/3.6)<4, "sit", "fly")))  # a few minus values
ggplot(data=zc_summary)+geom_density(aes(x=Hsig, col=ifelse((speed/3.6)<4, "sit", "fly")))
ggplot(data=zc_summary)+geom_density(aes(x=Hmean, col=ifelse((speed/3.6)<4, "sit", "fly")))
ggplot(data=zc_summary)+geom_density(aes(x=Tmean, col=ifelse((speed/3.6)<4, "sit", "fly")))
ggplot(data=zc_summary)+geom_density(aes(x=Tsig, col=ifelse((speed/3.6)<4, "sit", "fly")))
ggplot(data=zc_summary[zc_summary$Pvar<1,])+geom_point(aes(x=Pvar, y=Hmean, col=ifelse((speed/3.6)<4, "sit", "fly")), alpha=0.2, shape=1)

ggplot(data=dat%>%filter(burstID==2242))+geom_line(aes(x=DateTime_AEDT, y=Pressure))

#summary
na.omit(zc_summary)%>%filter((speed/3.6)>4)%>%select(-c('burstID', 'ID'))%>%group_by(sp)%>%summarise_all(list(mean=mean, sd=sd))

ggplot(data=zc_summary[zc_summary$Hmean<5,])+geom_point(aes(x=speed, y=Hmean), alpha=0.2, shape=1) # looks good to cut on speed alone
ggplot(data=zc_summary[zc_summary$Hmean<5 & zc_summary$Tmean< 25,])+geom_point(aes(x=Tmean, y=Hmean, col=ifelse((speed/3.6)<4, "sit", "fly")), alpha=0.2, shape=1)
ggplot(data=zc_summary[zc_summary$Hmean<5,])+geom_point(aes(x=temp, y=Hmean, col=ifelse((speed/3.6)<4, "sit", "fly")), alpha=0.2, shape=1)

#remove erroneous data bursts
bad_bursts<-zc_summary%>%filter(Hmean>5)%>%pull(burstID) # these are bursts with dives, not too many so remove
bad_bursts<-c(bad_bursts, zc_summary%>%filter(is.na(Hmean))%>%pull(burstID)) # weird pressure bursts, remove
bad_bursts<-c(bad_bursts, zc_summary%>%filter(temp<0)%>%pull(burstID)) # negative temps remove
bad_bursts<-c(bad_bursts, zc_summary%>%filter(Pvar>3)%>%pull(burstID)) # large variances

bad_bursts<-unique(bad_bursts)
#remove from burst summary data
zc_summary<-na.omit(zc_summary%>%filter(!burstID%in%bad_bursts))

# try kmeans clustering
km1<-kmeans(zc_summary$Hmean, 2)
pam1<-pam(zc_summary$Hmean, 2)
zc_summary$Hmean_class<-km1$cluster # 1= sit, 2= fly
zc_summary$final_class<-ifelse((zc_summary$speed/3.6)>=4 & zc_summary$Hmean_class==2, "fly",
                               ifelse((zc_summary$speed/3.6)<4 & zc_summary$Hmean_class==1, "sit",
                               ifelse((zc_summary$speed/3.6)>=4 & zc_summary$Hmean_class==1, "takeoff",
                               ifelse((zc_summary$speed/3.6)<4 & zc_summary$Hmean_class==2, "land", "NA"))))
table(zc_summary$final_class)

ggplot(data=zc_summary)+geom_point(aes(x=speed, y=Hmean, col=final_class), alpha=0.2, shape=1)

ggplot(data=zc_summary[zc_summary$final_class=='takeoff',])+geom_point(aes(x=Pvar, y=Hmean), alpha=0.2, shape=1)

for(i in zc_summary[zc_summary$final_class=='takeoff',]$burstID)
{
  print(ggplot(data=dat%>%filter(burstID==i))+geom_line(aes(x=DateTime_AEDT, y=Pressure))+labs(title=i))
  print(zc_summary[zc_summary$burstID==i,])
  readline()
}

zc_summary%>%filter(Tmean>18&Tmean<20)%>%pull(burstID)

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
