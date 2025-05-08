library(dplyr)
library(tidyr)
library(lubridate)
library(oceanwaves)
library(readxl)
library(sf)
library(tmap)
library(suncalc)
library(lutz)
library(adehabitatLT)

setwd("C:/Users/mmil0049/OneDrive - Monash University/fieldwork/Seadragon atsea deployment")

# read data

dat<-NULL
for(i in list.files("04mar2025_datapull"))
{
d1<-read.csv(paste0("04mar2025_datapull/", i), skip=1)
dat<-rbind(dat, data.frame(d1, ID=substr(i, 3, 6)))
}

dat<-dat%>%dplyr::select(-c('Number.of.SVs', 'Altitude..metres.',
                            'Sensor.Temperature', 'Time.Underwater', 'Activity.Count'))

# remove loggers that failed to turn on
dat<-dplyr::filter(dat, !ID%in%c(1202, 1319, 1666))

#assign species
dat$sp<-"NA"
dat[dat$ID %in%c(1313,1314,1315,1321,1322,1323),]$sp<-"BUAL"
dat[dat$ID %in%c(1183,1200, 1196, 1316, 1318),]$sp<-"WCAL"
dat[dat$ID %in%c(1428),]$sp<-"SHAL"
dat[dat$ID %in%c(1427, 1433),]$sp<-"WAAL"
dat[dat$ID %in%c(1437, 1660, 1443),]$sp<-"NZAL"
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
dat<-dat%>%filter(!(ID=="1660" & DateTime_AEDT<ymd_hms("2024-11-28 12:00:00", tz="Australia/Sydney")))

# make burst ID
dat$burstID<-as.numeric(as.factor(paste(dat$ID, year(dat$DateTime_AEDT), month(dat$DateTime_AEDT), day(dat$DateTime_AEDT), hour(dat$DateTime_AEDT),
                                        dat$Latitude, dat$Longitude, dat$Heading, dat$Battery.Voltage)))
table(table(dat$burstID))

# make Deployment duration fig
d1<-dat%>%group_by(sp, ID)%>%summarise(st_dep=first(DateTime_AEDT),
                                       end_dep=last(DateTime_AEDT))%>%as.data.frame()
d1$sp<-base::factor(d1$sp, levels=c('WCAL', 'SHAL', 'BUAL', 'BBAL', "IYNA",'NGPE', "SGPE","WAAL", "NZAL"),
                            labels= c('White-capped Albatross', 'Shy Albatross', "Buller's Albatross", 'Black-browed Albatross',
                                      "Indian Yellow-nosed Albatross",'Northern Giant-Petrel', "Southern Giant-Petrel", "Wandering Albatross", "New Zealand Wandering Albatross"))

d1$sp_ID<-factor(paste(d1$sp, d1$ID))

ggplot(data=d1)+
  geom_rect(aes(xmin=st_dep, xmax=end_dep, ymin=paste(sp,ID), ymax=paste(sp,ID)), col=1, size=3)+
  scale_x_datetime(date_minor_breaks="month")

#Duration lengths
d1$end_dep-d1$st_dep
#> mean(d1$end_dep-d1$st_dep)
#Time difference of 43.43571 days
# min(d1$end_dep-d1$st_dep)
#Time difference of 3.211794 days
# max(d1$end_dep-d1$st_dep)
#Time difference of 132.868 days
# sd(d1$end_dep-d1$st_dep)
#[1] 30.15556

# write out compiled coord file

#write.csv(dat%>%group_by(ID, burstID)%>%summarise_all(first), "compiled_04mar25_coords.csv", quote=F, row.names=F)

id_l<-data.frame(table(dat$burstID))
dat<-dat[dat$burstID%in% id_l[id_l$Freq>299,]$Var1,] # remove bursts less than 299 secs 

#remove first row of each burst as first pressure reading sometimes iffy
dat<-dat%>%group_by(burstID)%>%slice(-1)%>%as.data.frame()

dat_sf<-st_as_sf(dat%>%group_by(ID, burstID)%>%summarise_all(first), coords=c("Longitude", "Latitude"), crs=4326)
dat_sf<-dat_sf%>%arrange(ID, DateTime_AEDT)

linez<-dat_sf %>%group_by(ID, sp) %>%
  dplyr::summarize(do_union=F) %>%  # do_union=FALSE doesn't work as well
  st_cast("LINESTRING")%>%st_wrap_dateline()
#st_write(linez, "compiled_13dec_linez.shp")

linez_laea<-st_transform(linez, crs="+proj=laea +lat_0=-46 +lon_0=158 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

ggplot(data=linez_laea)+geom_sf(aes(colour=sp))

tmap_mode("view")
tm_shape(linez)+tm_lines(col=as.character("ID"))

tm_shape(linez)+tm_lines(col=as.character("sp"))


# get speed histogram
sp_df<-dat%>%group_by(sp, burstID)%>%summarise(speed=first(Speed..km.h.))

sp_df%>%filter((speed/3.6)>4)%>%
  summarise(mean_s=mean(speed/3.6), sd_s=sd(speed/3.6))
#sp    mean_s  sd_s  
#<chr>  <dbl> <dbl>
# sp    mean_s  sd_s
#<chr>  <dbl> <dbl>
#  1 BBAL    14.8  5.51
#2 BUAL    13.5  4.94
#3 IYNA    11.5  3.88
#4 NGPE    14.7  5.43
#5 NZAL    14.5  4.72
#6 SGPE    16.3  5.73
#7 SHAL    13.2  4.30
#8 WAAL    15.5  5.16
#9 WCAL    14.2  4.85

sp_df%>%filter((speed/3.6)>4)%>%
summarise(mean_s=mean(speed/3.6), sd_s=sd(speed/3.6), min=min(speed/3.6), q01=quantile(speed/3.6, 0.01), q05=quantile(speed/3.6, 0.05), q25=quantile(speed/3.6, 0.25),
          med=median(speed/3.6), q75=quantile(speed/3.6, 0.75),q95=quantile(speed/3.6, 0.95), q99=quantile(speed/3.6, 0.99), max=max(speed/3.6))

# sp    mean_s  sd_s   min   q01   q05   q25   med   q75   q95   q99   max
#<chr>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#1 BBAL    14.8  5.51  4.17  6.14  7.86 10.7   13.9  17.8  25.6  30.8  31.7
#2 BUAL    13.5  4.94  4.17  4.72  6.67 10     12.8  16.4  22.9  27.2  33.6
#3 IYNA    11.5  3.88  4.44  5     5.86  9.17  11.1  13.3  19.7  23.6  24.7
#4 NGPE    14.7  5.43  4.17  5.68  7.22 10.6   13.9  18.3  24.9  28.8  36.9
#5 NZAL    14.5  4.72  4.17  5.83  8.06 11.4   13.9  16.9  23.1  27.7  61.9
#6 SGPE    16.3  5.73  4.44  6.83  8.61 12.4   15.3  19.4  27.5  32.0  35.3
#7 SHAL    13.2  4.30  4.44  5.28  6.57 10.3   12.8  15.3  21.3  25.9  27.2
#8 WAAL    15.5  5.16  4.17  5.83  8.33 11.7   14.7  18.6  25    28.6  36.9
#9 WCAL    14.2  4.85  4.17  5.28  7.78 10.8   13.3  16.7  23.6  28.7  33.6  

ggplot()+
  geom_histogram(data=sp_df,
                 aes(x=speed/3.6, fill=ifelse((speed/3.6)<4, "blue", "red")), colour=1, binwidth=1, boundary=0)+
theme_bw()+ylab("Count of GPS datapoints")+xlab("speed m/s")+
  theme(legend.position="none")+facet_wrap(~sp, nrow=2, scales='free_y')

ggplot()+
  geom_density(data=sp_df%>%filter((speed/3.6)>4),
                 aes(x=speed/3.6, colour=sp))+
 theme_bw()+ylab("Count of GPS datapoints")+xlab("speed m/s")+labs(colour='species')

# Nocturnal flight analyses
# Add day night classification.. slow loop
NAF_df<-dat%>%group_by(sp,ID, burstID)%>%summarise(speed=first(Speed..km.h.), UTC.Timestamp =first(UTC.Timestamp),
                                          DateTime_AEDT  =first(DateTime_AEDT), Latitude=first(Latitude), Longitude=first(Longitude))

#interpolate to hourly points
trajectories <- as.ltraj(xy=data.frame(NAF_df$Longitude,
                                       NAF_df$Latitude),
                         date=as.POSIXct(NAF_df$UTC.Timestamp, tz="UTC"),
                         id=NAF_df$ID, infolocs=NAF_df[c("burstID","speed")]%>%as.data.frame(), typeII = TRUE)   

trajectories_int <- redisltraj(trajectories, 3600, type="time")

# to data.frame
trajectories_int<-ld(trajectories_int)
plot(y~x, data=trajectories_int, col=id)
#rm some weird points
bad_pz<-c(trajectories_int[trajectories_int$y> -40 &trajectories_int$x < -100&trajectories_int$x > -150 &trajectories_int$id==1433,]$pkey,
          trajectories_int[trajectories_int$x>50 &trajectories_int$x<100,]$pkey)

trajectories_int<-trajectories_int%>%filter(!pkey%in%bad_pz)

#join original data
NAF_df$pkey<-paste(NAF_df$ID, NAF_df$UTC.Timestamp, sep=".")
NAF_df<-left_join(trajectories_int, NAF_df[,c(1,3,4,9)], by='pkey')

NAF_df$UTC.Timestamp<-ymd_hms(NAF_df$date, tz="UTC")
NAF_df$local_tz<-tz_lookup_coords(lat=NAF_df$y, lon=NAF_df$x, method='accurate')

NAF_df$Timestamp.local<-as.character("none")
for(i in 1:nrow(NAF_df)){NAF_df$Timestamp.local[i]<-as.character(with_tz(NAF_df$date[i], tz=NAF_df$local_tz[i]))}

NAF_df$daynight<-"night"
for(i in 1:nrow(NAF_df)){
  dn1<-getSunlightTimes(date = as.Date(NAF_df$Timestamp.local[i],tz = NAF_df$local_tz[i]), lat = NAF_df$y[i], lon = NAF_df$x[i],tz = NAF_df$local_tz[i])
  if(format(as.POSIXct(NAF_df$Timestamp.local[i]), format = '%H:%M:%S')<format(as.POSIXct(dn1$sunset), format = '%H:%M:%S')& 
     format(as.POSIXct(NAF_df$Timestamp.local[i]), format = '%H:%M:%S')>format(as.POSIXct(dn1$dawn), format = '%H:%M:%S')){
    NAF_df$daynight[i]<-"day"};print(i)}                       
# use civil dawn and dusk as per Bonnet-Lebrun
NAF_df%>%group_by(id, daynight)%>%summarise(n())%>%View()
NAF_df%>%filter(id==1433)%>%View()

#need to remove points when on land: Macca, Snares and Solander

NAF_df_sf<-st_as_sf(NAF_df, coords=c("x", "y"), crs=4326)

tmap_mode("view")
tm_shape(dat_sf)+tm_dots(fill='red')+tm_shape(NAF_df_sf[,1])+tm_dots()



# Use zero crossing analyses to summarize pressure values per burst.
# Hoping this combined with speed/temp can classify sitting flying bursts.

#make 1 min(almost) sub-bursts
dat<-dat%>%group_by(burstID)%>%mutate(sub_burst = cut(UTC.Timestamp, breaks = 5, labels=c("a", "b", "c", "d", "e" )))
dat$sub_burst<-paste(dat$burstID, dat$sub_burst, sep="_")

sub_burst_summary<-dat%>%group_by(burstID, sub_burst)%>%
                         summarise(Pmed=median(Pressure))%>%ungroup%>%
  group_by(burstID)%>%mutate(pres_ch=TRUE %in% (abs((median(Pmed)-Pmed))>0.5))

sub_burst_summary<-sub_burst_summary%>%group_by(burstID)%>%slice(1)

zc_summary<-data.frame(dat%>%group_by(sp, ID, burstID)%>%summarise(speed=first(Speed..km.h.),
                      temp=first(Device.Temperature),Psd=sd(Pressure)),
                       Hsig=NA,Hmean=NA, Tmean=NA, Tsig=NA)
zc_summary<-left_join(zc_summary, sub_burst_summary%>%dplyr::select(burstID, pres_ch), by="burstID")

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

#check and remove extreme values (temp, Hsig etc)

ggplot(data=zc_summary)+geom_histogram(aes(x=speed))
ggplot(data=zc_summary)+geom_histogram(aes(x=Psd), bins=50)
ggplot(data=zc_summary)+geom_density(aes(x=temp, col=ifelse((speed/3.6)<4, "sit", "fly")))  # a few minus values
ggplot(data=zc_summary)+geom_density(aes(x=Hsig, col=ifelse((speed/3.6)<4, "sit", "fly")))
ggplot(data=zc_summary)+geom_density(aes(x=Hmean, col=ifelse((speed/3.6)<4, "sit", "fly")))
ggplot(data=zc_summary)+geom_density(aes(x=Tmean, col=ifelse((speed/3.6)<4, "sit", "fly")))
ggplot(data=zc_summary)+geom_density(aes(x=Tsig, col=ifelse((speed/3.6)<4, "sit", "fly")))
ggplot(data=zc_summary[zc_summary$Psd<1,])+geom_point(aes(x=Psd, y=Hmean, col=ifelse((speed/3.6)<4, "sit", "fly")), alpha=0.2, shape=1)

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
bad_bursts<-c(bad_bursts, zc_summary%>%filter(Psd>1.5)%>%pull(burstID)) # large variances
bad_bursts<-c(bad_bursts, zc_summary%>%filter(pres_ch==TRUE)%>%pull(burstID)) # earlier identified step changes/gradients

bad_bursts<-unique(bad_bursts)
#remove from burst summary data
zc_summary<-na.omit(zc_summary%>%filter(!burstID%in%bad_bursts))

# try kmeans clustering
km1<-kmeans(zc_summary$Hmean, 2)
library(cluster)
pam1<-pam(zc_summary$Hmean, 2)
zc_summary$Hmean_class<-km1$cluster # 1= sit, 2= fly
zc_summary$final_class<-ifelse((zc_summary$speed/3.6)>=4 & zc_summary$Hmean_class==2, "fly",
                               ifelse((zc_summary$speed/3.6)<4 & zc_summary$Hmean_class==1, "sit",
                               ifelse((zc_summary$speed/3.6)>=4 & zc_summary$Hmean_class==1, "takeoff",
                               ifelse((zc_summary$speed/3.6)<4 & zc_summary$Hmean_class==2, "land", "NA"))))
table(zc_summary$final_class)

ggplot(data=zc_summary)+geom_point(aes(x=speed, y=Hmean, col=final_class), alpha=0.2, shape=1)
# speed break, pretty clear
#could check takeoff with loop below?
ggplot(data=zc_summary[zc_summary$final_class=='takeoff',])+geom_point(aes(x=Psd, y=Hmean), alpha=0.2, shape=1)
for(i in zc_summary[zc_summary$final_class=='takeoff',]$burstID)
{
  print(ggplot(data=dat%>%filter(burstID==i))+geom_line(aes(x=DateTime_AEDT, y=Pressure))+labs(title=i))
  print(zc_summary[zc_summary$burstID==i,])
  readline()
}
zc_summary$speed_class<-"sit"
zc_summary[(zc_summary$speed/3.6)>4,]$speed_class<-"fly"
ggplot(data=zc_summary)+geom_point(aes(x=speed, y=Hmean, col=speed_class), alpha=0.2, shape=1)
# use simple sit/fly to identify flying bursts for flight height calculation
fly_bursts<-zc_summary[zc_summary$speed_class=='fly',]$burstID

#### Calculation of flight height using dynamic soaring method ####
#remove bar bursts first
dat<-dat%>%ungroup()%>%filter(!burstID%in% bad_bursts)%>%as.data.frame()

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

#### Summarise and compare altitude between species  ####

dat_flying<-dat%>%filter(burstID%in%fly_bursts) # select flying only data

#remove some erroneous bursts
dat_flying<-filter(dat_flying, burstID!=8788)
dat_flying<-filter(dat_flying, burstID!=11207)
dat_flying<-filter(dat_flying, burstID!=11306)

#summarise
dat_comp<-dat_flying%>%group_by(sp)%>%summarise(n_bird=n_distinct(ID), n_bursts=n_distinct(burstID),mn_alt=mean(alt_DS), sd_alt=sd(alt_DS), median=median(alt_DS),
                                        min=min(alt_DS), max=max(alt_DS),
                                        q25=quantile(alt_DS, 0.25), q75=quantile(alt_DS, 0.75),
                                        q5=quantile(alt_DS, 0.05), q95=quantile(alt_DS, 0.95),
                                        q1=quantile(alt_DS, 0.01), q99=quantile(alt_DS, 0.99))
#write.csv(dat_comp, 'reporting/final_reporting_nov24/height_table_mar2025.csv')
#ignore min values as dives!

dat_flying$sp<-as.factor(dat_flying$sp)
dat_flying$sp<-base::factor(dat_flying$sp, levels=c('WCAL', 'SHAL', 'BUAL', 'BBAL', "IYNA",'NGPE', "SGPE","WAAL", "NZAL"),
                            labels= c('White-capped Albatross', 'Shy Albatross', "Buller's Albatross", 'Black-browed Albatross',
                                      "Indian Yellow-nosed Albatross",'Northern Giant-Petrel', "Southern Giant-Petrel", "Wandering Albatross",
                                      "New Zealand Wandering Albatross"))

#bin into 1m bands
bins_1m<-dat_flying%>%group_by(sp)%>%
  mutate(category = cut(alt_DS, breaks = c(-Inf,seq(-10, 50, 1),Inf), labels = seq(-10, 51, 1))) %>%
                        count(category)

ggplot(bins_1m, aes(x=factor(category), y=n)) +
  geom_bar(stat="identity", width=1.0, col=1)+
  geom_vline(xintercept = "0", linetype='dotted')+
  geom_vline(xintercept = "10", linetype='dotted', col='orange')+
  geom_vline(xintercept = "20", linetype='dotted', col='orange')+
  geom_vline(xintercept = "30", linetype='dotted', col='orange')+
  scale_x_discrete(breaks=seq(-10, 50, 5))+
  facet_wrap(~sp, ncol=2,scales='free_y')+theme_bw()+
  labs(x="Flight height (m)", y="Count")

# make sp plots


cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ffff33','#ff7f00','#a65628','#f781bf')

cols.alpha<-c(grDevices::adjustcolor(cols[1], alpha.f = 0.75),
              grDevices::adjustcolor(cols[2], alpha.f = 0.75),
              grDevices::adjustcolor(cols[3], alpha.f = 0.75),
              grDevices::adjustcolor(cols[4], alpha.f = 0.75),
              grDevices::adjustcolor(cols[5], alpha.f = 0.75),
              grDevices::adjustcolor(cols[6], alpha.f = 0.75),
              grDevices::adjustcolor(cols[7], alpha.f = 0.75),
              grDevices::adjustcolor(cols[8], alpha.f = 0.75))

ggplot(data=dat_flying)+geom_density(aes(x=alt_DS, colour=sp), bw="bcv", trim=T, fill=NA, size=2)+
  geom_vline(xintercept = 0, linetype='dotted')+scale_x_continuous(breaks=seq(-10,50,2))+
  coord_cartesian(xlim=c(-10, 50))+
  theme(legend.position= c(0.8,0.8), axis.text=element_text(size=10),axis.title=element_text(size=12),
        legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  scale_colour_manual("Species", values=cols.alpha,labels=c(
    'White-capped Albatross', 'Shy Albatross', "Buller's Albatross", 'Black-browed Albatross',
    "Indian Yellow-nosed Albatross",'Northern Giant Petrel', "Southern Giant Petrel", "Wandering Albatross"))+
  labs(x="Altitude (m)", y="Density")

ggplot(data=dat_flying)+geom_density(aes(x=alt_DS), bw="bcv", fill=NA, size=1, trim=T)+
  geom_vline(xintercept = 0, linetype='dotted')+
  geom_vline(xintercept = 10, linetype='dotted', col='orange')+
  geom_vline(xintercept = 20, linetype='dotted', col='orange')+
  geom_vline(xintercept = 30, linetype='dotted', col='orange')+
  scale_x_continuous(breaks=seq(-10, 50, 5))+
  coord_cartesian(xlim=c(-10, 50))+
  facet_wrap(~sp, ncol=2)+theme_bw()+
  labs(x="Altitude (m)", y="Density")



dat%>%filter((Speed..km.h./3.6)>4)%>%group_by(sp)%>%summarise(mn_alt=mean(alt_DS), sd_alt=sd(alt_DS), median=median(alt_DS),
                                        min=min(alt_DS), max=max(alt_DS),
                                        q25=quantile(alt_DS, 0.25), q75=quantile(alt_DS, 0.75))

ggplot(data=dat%>%filter((Speed..km.h./3.6)>4))+geom_density(aes(x=alt_DS, colour=sp), fill=NA)+
  theme_bw()+geom_vline(xintercept = 0, linetype='dotted')+scale_x_continuous(limits=c(-20, 100))
