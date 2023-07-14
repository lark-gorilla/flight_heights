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

dat<-rbind(data.frame(read.csv("gps_89460800120108611854.csv"), ID="08611854"),
           data.frame(read.csv("gps_89460800120141490936.csv"), ID="41490936"))

table(dat$iccid)
table(dat$ID)

dat$DateTime_UTC<-ymd_hms(dat$time_UTC, tz="UTC")
dat$DateTime_AEDT<-with_tz(dat$DateTime_UTC, "Australia/Sydney")

# remove 3 positional errors
ggplot(data=dat)+geom_point(aes(x=lon, y=lat, colour=ID))
dat<-dat[dat$lat<0,] # 3 of them but don't remove as valid pressure and temp

#identify time pre/post deployment on each bird
# 854 out 23-04-01 17:12:00
# 936 out 23-04-01 17:59:00
dat$deployed_ID<-dat$ID
dat[dat$DateTime_AEDT<ymd_hms("23-04-01 17:12:00", tz="Australia/Sydney") & dat$ID== "08611854",]$deployed_ID<-"predeployment"
dat[dat$DateTime_AEDT<ymd_hms("23-04-01 17:59:00", tz="Australia/Sydney") & dat$ID== "41490936",]$deployed_ID<-"predeployment"

# convert to spatial
dat_sf<-st_as_sf(dat, coords=c("lon", "lat"), crs=4326)
#view in tmap
tmap_mode("view")
tm_shape(dat_sf)+tm_dots(col="deployed_ID")+tm_mouse_coordinates()

# track2kba workflow to get extra fields

dataGroup <- formatFields(
  dataGroup = dat, 
  fieldID   = "ID", 
  fieldDateTime = "DateTime_AEDT", 
  fieldLon  = "lon", 
  fieldLat  = "lat",
  formatDT="ymd_HMS")

colony <- data.frame(Longitude = 144.65619, Latitude  = -40.37442)

trips <- tripSplit(
  dataGroup  = dataGroup,
  colony     = colony,
  innerBuff  = 1,      # kilometers
  returnBuff = 3,
  duration   = 1,      # hours
  rmNonTrip  = FALSE)

mapTrips(trips = trips, colony = colony, colorBy = "trip")

str(trips)

sumTrips <- tripSummary(trips = trips, colony = colony)

sumTrips


# use trips object for further analyses
tripdat<-trips@data
names(tripdat)[names(tripdat)=="DateTime"]<-"DateTime_AEDT" # rest name

# quick check
ggplot(data=tripdat)+
  geom_point(aes(x=DateTime_AEDT, y=ColDist, colour=tripID))+
  geom_vline(xintercept=ymd_hms("23-04-01 17:59:00", tz="Australia/Sydney"))+
  facet_wrap(~ID)

##### now segment data to classify each ~ 5 minute burst, add tdiff and distance between pts ####

#temporarily split out birds

tripdat1<-tripdat%>%filter(ID=="08611854")
tripdat2<-tripdat%>%filter(ID=="41490936")

tripdat1$tdiff<-c(diff(tripdat1$DateTime_AEDT),0)
tripdat1$burstID<-NA

counter=1 # slow loop but works
for(i in 1:nrow(tripdat1)){
  tripdat1$burstID[i]<-paste(tripdat1$tripID[i], counter, sep="_")
  if(tripdat1$tdiff[i]>100|is.na(tripdat1$tdiff[i])){counter=counter+1};print(i)}

tripdat1_sf<-st_as_sf(tripdat1, coords=c("Longitude", "Latitude"), crs=4326)
tripdat1_sf<-cbind(tripdat1_sf, c(st_geometry(tripdat1_sf)[2:nrow(tripdat1_sf)], st_geometry(tripdat1_sf)[1]))%>%
  mutate(dist_diff = st_distance(geometry, geometry.1, by_element = T))
tripdat1$dist_diff<-tripdat1_sf$dist_diff

tripdat2$tdiff<-c(diff(tripdat2$DateTime_AEDT),0)
tripdat2$burstID<-NA

counter=1
for(i in 1:nrow(tripdat2)){
  tripdat2$burstID[i]<-paste(tripdat2$tripID[i], counter, sep="_")
  if(tripdat2$tdiff[i]>100|is.na(tripdat2$tdiff[i])){counter=counter+1};print(i)}

tripdat2_sf<-st_as_sf(tripdat2, coords=c("Longitude", "Latitude"), crs=4326)
tripdat2_sf<-cbind(tripdat2_sf, c(st_geometry(tripdat2_sf)[2:nrow(tripdat2_sf)], st_geometry(tripdat2_sf)[1]))%>%
  mutate(dist_diff = st_distance(geometry, geometry.1, by_element = T))
tripdat2$dist_diff<-tripdat2_sf$dist_diff

tripdat<-rbind(tripdat1, tripdat2) # bind up again
tripdat$geometry.1<-NULL
tripdat$speed_diff<-as.numeric(tripdat$dist_diff/as.numeric(tripdat$tdiff)) # get speed between points (don't trust device speed)
# remove 3 erroneous datapoints over 50 m/s
tripdat<-tripdat%>%filter(speed_diff<50)

# check segments
blen<-length(unique(tripdat$burstID))
cols = rainbow(blen, s=.6, v=.9)[sample(1:blen,blen)]

ggplot(data=tripdat)+
  geom_point(aes(x=DateTime_AEDT, y=ColDist, colour=burstID))+
  facet_wrap(~ID, scales="free") + scale_colour_manual(values=cols)+theme(legend.position = "none")

table(tripdat$burstID)

# classify non flying points using speed

ggplot(data=tripdat%>%filter(is.finite(speed_diff)))+geom_histogram(aes(x=speed_diff))+facet_wrap(~tripID)

#plot non-island points 
ggplot(data=tripdat%>%filter(tripID!="-1"&is.finite(speed_diff &ColDist>1)))+
  geom_histogram(aes(x=speed_diff, fill=ifelse(speed_diff<4, "blue", "red")), binwidth=1)+scale_x_continuous(breaks=0:50)

tripdat$sit_fly<-"sit"
tripdat[tripdat$]


# convert to spatial
dat_sf<-st_as_sf(tripdat, coords=c("Longitude", "Latitude"), crs=4326)
#view in tmap
tmap_mode("view")
tm_shape(dat_sf)+tm_dots(col="speed_diff")+tm_mouse_coordinates()



# lots of variation in temp.. maybe correlated with solar heating?

ggpairs(data=dat%>%filter(DateTime_AEDT>ymd_hms("23-04-01 17:59:00", tz="Australia/Sydney")),
        columns=c("DateTime_AEDT", "lat", "iSolar", "alt", "temp",
                  "pres_pa"), ggplot2::aes(colour=ID))

