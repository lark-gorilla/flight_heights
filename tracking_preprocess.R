# data explore clean

library(ggplot2)
library(lubridate)
library(dplyr)
library(GGally)
library(track2KBA)
library(tmap)
library(sf)

setwd("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights/data/shy_albatross_island")

# Combine data from two tags that uploaded data while on birds

dat<-rbind(read.csv("gps_89460800120108611854.csv"),read.csv("gps_89460800120141490936.csv"))
table(dat$iccid)

dat$DateTime_UTC<-ymd_hms(dat$time_UTC, tz="UTC")
dat$DateTime_AEDT<-with_tz(dat$DateTime_UTC, "Australia/Sydney")

# positional errors
ggplot(data=dat)+geom_point(aes(x=lon, y=lat))
dat[dat$lon==0,] # 3 of them but don't remove as valid pressure and temp
ggplot(data=filter(dat, lat<0))+geom_point(aes(x=lon, y=lat, colour=factor(iccid)))


# convert to spatial
dat_sf<-st_as_sf(dat, coords=c("lon", "lat"), crs=4326)
#view in tmap
tmap_mode("view")
tm_shape(dat_sf)+tm_dots(col=factor("iccid"))+tm_mouse_coordinates()

#filter time to after deployment on each bird
# 854 out 23-04-01 17:12:00
# 936 out 23-04-01 17:59:00

# frequency
ggplot(data=dat)+
  geom_point(aes(x=DateTime_AEDT, y=temp))+
  geom_vline(xintercept=ymd_hms("23-04-01 17:59:00", tz="Australia/Sydney"))+
  facet_wrap(~iccid)

# lots of variation in temp.. maybe correlated with solar heating?

ggpairs(data=dat%>%filter(DateTime_AEDT>ymd_hms("23-04-01 17:59:00", tz="Australia/Sydney")),
        columns=c("DateTime_AEDT", "lat", "iSolar", "alt", "temp",
                            "pres_pa"), ggplot2::aes(colour=factor(iccid)))

# track2kba workflow to get extra fields

dataGroup <- formatFields(
  dataGroup = dat, 
  fieldID   = "iccid", 
  fieldDateTime = "DateTime_AEDT", 
  fieldLon  = "lon", 
  fieldLat  = "lat",
  formatDT="ymd_HMS")

colony <- dataGroup %>% 
  summarise(Longitude = first(Longitude), Latitude  = first(Latitude))
