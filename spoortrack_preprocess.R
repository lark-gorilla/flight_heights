library(dplyr)
library(tidyr)
library(patchwork)
library(lubridate)


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

sp_df%>%filter(tripID!="-1"&is.finite(speed_diff) & ColDist>1 & speed_diff>4)%>%
  summarise(mean_s=mean(speed_diff), sd_s=sd(speed_diff))
#mean_s     sd_s
#1 13.4574 4.1042

#plot non-island points 
ggplot()+
  geom_histogram(data=sp_df,
                 aes(x=speed, fill=ifelse(speed<15, "blue", "red")), colour=1, binwidth=1, boundary=0)+
  geom_label(aes(x=13, y=15000, label="Mean flying speed = 13.46±4.1 m/s"), size=5)+
  scale_x_continuous(breaks=0:150)+theme_bw()+ylab("Count of GPS datapoints")+xlab("speed km/h")+theme(legend.position="none")
