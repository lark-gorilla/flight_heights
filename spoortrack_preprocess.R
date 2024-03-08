library(dplyr)
library(tidyr)
library(patchwork)
library(lubridate)
library(oceanwaves)


setwd("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights")

# read data

dat<-read.csv("C:/Users/mmil0049/OneDrive - Monash University/fieldwork/Seadragon atsea deployment/first_two_weeks/albatross_data_121223_to_020124.csv", h=T, sep=";")

dat$burstID=as.numeric(as.factor(paste(dat$product_sn, dat$latitude, dat$longitude, dat$heading)))

id_l<-data.frame(table(dat$burstID))
table(id_l$Freq) # remove few small IDs

dat<-dat[dat$burstID%in% id_l[id_l$Freq==300,]$Var1,]

dat<-dat%>%arrange(product_sn, timestamp)

dat$timestamp2<-lubridate::ymd_hms("2023-12-13 06:04:01")
for(i in unique(dat$burstID))
{
  dat[dat$burstID==i,]$timestamp2<-lubridate::ymd_hms(
    paste(dat[dat$burstID==i,]$timestamp,rep(c("00","01","02","03", "04", "05", "06", "07", "08", "09", 10:59), 5), sep=":"), tz="UTC")
}

#dat_sf<-st_as_sf(dat, coords=c("longitude", "latitude"), crs=4326
#tm_shape(dat_sf%>%filter(timestamp2> lubridate::ymd_hms("2023-12-19 22:00:00") & timestamp2< lubridate::ymd_hms("2023-12-21 #01:00:00")))+tm_dots(col=as.character("product_sn"))

dat<-dat%>%filter(timestamp2>lubridate::ymd_hms("2023-12-20 01:00:00")) #filter out pre deployment
dat<-dat%>%filter(product_sn!='1202') # filter out failed logger

dat$timestamp2<-lubridate::with_tz(dat$timestamp2, "Australia/Sydney")

# get speed histogram
sp_df<-dat%>%group_by(burstID)%>%summarise(speed=first(speed_kmh))

sp_df%>%filter((speed/3.6)>4)%>%
  summarise(mean_s=mean(speed/3.6), sd_s=sd(speed/3.6))
#mean_s     sd_s
#1 14.1 4.64

 
ggplot()+
  geom_histogram(data=sp_df,
                 aes(x=speed, fill=ifelse(speed<15, "blue", "red")), colour=1, binwidth=1, boundary=0)+
  scale_x_continuous(breaks=0:150)+theme_bw()+ylab("Count of GPS datapoints")+xlab("speed km/h")+theme(legend.position="none")

ggplot()+
  geom_histogram(data=sp_df,
                 aes(x=speed/3.6, fill=ifelse((speed/3.6)<4, "blue", "red")), colour=1, binwidth=1, boundary=0)+
  scale_x_continuous(breaks=0:32)+theme_bw()+ylab("Count of GPS datapoints")+xlab("speed m/s")+theme(legend.position="none")


zc_summary<-NULL
for (i in dat$burstID)
{
  tout<-data.frame(burstID=i, speed=unique(dat[dat$burstID==i,]$speed_kmh),
                   temp=unique(dat[dat$burstID==i,]$temperature),
                   ds_hsig=NA,ds_hmean=NA, ds_tmean=NA, ds_tsig=NA)
  
  d1<-waveStatsZC(dat[dat$burstID==i,]$pressure, 1,)
  tout$ds_hsig=d1$Hsig
  tout$ds_hmean=d1$Hmean
  tout$ds_tmean=d1$Tmean
  tout$ds_tsig=d1$Tsig
  
  zc_summary<-rbind(zc_summary, tout)
}

ggplot(data=zc_summary)+geom_histogram(aes(x=speed))
ggplot(data=zc_summary)+geom_histogram(aes(x=temp))
ggplot(data=zc_summary)+geom_histogram(aes(x=Hsig))
ggplot(data=zc_summary)+geom_histogram(aes(x=Hmean))
ggplot(data=zc_summary)+geom_histogram(aes(x=Tmean))
ggplot(data=zc_summary)+geom_histogram(aes(x=Tsig))
