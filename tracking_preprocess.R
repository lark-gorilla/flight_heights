# data explore clean

library(ggplot2)
library(lubridate)
library(dplyr)
library(GGally)
library(track2KBA)
library(tmap)
library(sf)
library(EMbC)
library(geosphere)
library(suncalc)
library(rWind)


#### Load data and prep ####

setwd("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights")

# Combine data from two tags that uploaded data while on birds

dat<-rbind(data.frame(read.csv("data/shy_albatross_island/gps_89460800120108611854.csv"), ID="08611854"),
           data.frame(read.csv("data/shy_albatross_island/gps_89460800120141490936.csv"), ID="41490936"),
           data.frame(read.csv("data/shy_albatross_island/gps_89460800120108611649.csv"), ID="08611649"))

table(dat$iccid)
table(dat$ID)

dat<-dat[-80141,] # fix

dat$DateTime_UTC<-ymd_hms(dat$time_UTC, tz="UTC")
dat$DateTime_AEDT<-with_tz(dat$DateTime_UTC, "Australia/Sydney")

# remove a few positional errors
#ggplot(data=dat)+geom_point(aes(x=lon, y=lat, colour=ID))
dat<-dat[dat$lat<0,]

#identify time pre/post deployment on each bird
# 854 out 23-04-01 17:12:00
# 936 out 23-04-01 17:59:00
dat$deployed_ID<-as.character(dat$ID)
dat[dat$DateTime_AEDT<ymd_hms("23-04-01 17:12:00", tz="Australia/Sydney") & dat$ID== "08611854",]$deployed_ID<-"predeployment"
dat[dat$DateTime_AEDT<ymd_hms("23-04-01 17:59:00", tz="Australia/Sydney") & dat$ID== "41490936",]$deployed_ID<-"predeployment"
dat[dat$DateTime_AEDT<ymd_hms("23-09-29 03:10:00", tz="Australia/Sydney") & dat$ID== "08611649",]$deployed_ID<-"predeployment"


# convert to spatial
dat_sf<-st_as_sf(dat, coords=c("lon", "lat"), crs=4326)
#view in tmap
tmap_mode("view")
tm_shape(dat_sf)+tm_dots(col="deployed_ID")+tm_mouse_coordinates()

#### ---- ####

#### track2kba workflow to get extra fields ####

dataGroup <- formatFields(
  dataGroup = dat[dat$deployed_ID!='predeployment',], 
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
names(tripdat)[names(tripdat)=="DateTime"]<-"DateTime_AEDT" # reset name

#### ---- ####

##### now segment data to classify each ~ 5 minute burst, also variables that use diff between locs: time, distance, bearing ####

#temporarily split out birds

tripdat1<-tripdat%>%filter(ID=="08611854")
tripdat2<-tripdat%>%filter(ID=="41490936")
tripdat3<-tripdat%>%filter(ID=="08611649")

tripdat1$tdiff<-c(diff(tripdat1$DateTime_AEDT),0)
tripdat1$burstID<-NA

counter=1 # slow loop but works
for(i in 1:nrow(tripdat1)){
  tripdat1$burstID[i]<-paste(tripdat1$tripID[i], counter, sep="_")
  if(tripdat1$tdiff[i]>100|is.na(tripdat1$tdiff[i])){counter=counter+1};print(i)}

tripdat1_sf<-st_as_sf(tripdat1, coords=c("Longitude", "Latitude"), crs=4326)
tripdat1_sf<-cbind(tripdat1_sf, c(st_geometry(tripdat1_sf)[2:nrow(tripdat1_sf)], st_geometry(tripdat1_sf)[1]))%>%
  mutate(dist_diff = st_distance(geometry, geometry.1, by_element = T))
tripdat1$dist_diff<-c(as.numeric(tripdat1_sf$dist_diff[1:(nrow(tripdat1_sf)-1)]), 0)
tripdat1$bearing_diff<-c(apply(tripdat1_sf[1:(nrow(tripdat1_sf)-1),], 1,
                               function(x){(geosphere::bearing(unlist(x["geometry"]), unlist(x["geometry.1"]))+360)%%360}), 0)

tripdat2$tdiff<-c(diff(tripdat2$DateTime_AEDT),0)
tripdat2$burstID<-NA

counter=1
for(i in 1:nrow(tripdat2)){
  tripdat2$burstID[i]<-paste(tripdat2$tripID[i], counter, sep="_")
  if(tripdat2$tdiff[i]>100|is.na(tripdat2$tdiff[i])){counter=counter+1};print(i)}

tripdat2_sf<-st_as_sf(tripdat2, coords=c("Longitude", "Latitude"), crs=4326)
tripdat2_sf<-cbind(tripdat2_sf, c(st_geometry(tripdat2_sf)[2:nrow(tripdat2_sf)], st_geometry(tripdat2_sf)[1]))%>%
  mutate(dist_diff = st_distance(geometry, geometry.1, by_element = T))
tripdat2$dist_diff<-c(as.numeric(tripdat2_sf$dist_diff[1:(nrow(tripdat2)-1)]), 0)
tripdat2$bearing_diff<-c(apply(tripdat2_sf[1:(nrow(tripdat2_sf)-1),], 1,
                               function(x){(geosphere::bearing(unlist(x["geometry"]), unlist(x["geometry.1"]))+360)%%360}), 0)

tripdat3$tdiff<-c(diff(tripdat3$DateTime_AEDT),0)
tripdat3$burstID<-NA

counter=1
for(i in 1:nrow(tripdat3)){
  tripdat3$burstID[i]<-paste(tripdat3$tripID[i], counter, sep="_")
  if(tripdat3$tdiff[i]>100|is.na(tripdat3$tdiff[i])){counter=counter+1};print(i)}

tripdat3_sf<-st_as_sf(tripdat3, coords=c("Longitude", "Latitude"), crs=4326)
tripdat3_sf<-cbind(tripdat3_sf, c(st_geometry(tripdat3_sf)[2:nrow(tripdat3_sf)], st_geometry(tripdat3_sf)[1]))%>%
  mutate(dist_diff = st_distance(geometry, geometry.1, by_element = T))
tripdat3$dist_diff<-c(as.numeric(tripdat3_sf$dist_diff[1:(nrow(tripdat3)-1)]), 0)
tripdat3$bearing_diff<-c(apply(tripdat3_sf[1:(nrow(tripdat3_sf)-1),], 1,
                               function(x){(geosphere::bearing(unlist(x["geometry"]), unlist(x["geometry.1"]))+360)%%360}), 0)


tripdat<-rbind(tripdat1, tripdat2, tripdat3) # bind up again
tripdat$geometry.1<-NULL
tripdat$speed_diff<-as.numeric(tripdat$dist_diff/as.numeric(tripdat$tdiff)) # get speed between points (don't trust device speed)

#### ---- ####

#### clean up identified bursts ####

# remove 3 erroneous datapoints over 50 m/s
tripdat<-tripdat%>%filter(speed_diff<50)

# check burst segments
blen<-length(unique(tripdat$burstID))
cols = rainbow(blen, s=.6, v=.9)[sample(1:blen,blen)]

ggplot(data=tripdat)+
  geom_point(aes(x=DateTime_AEDT, y=ColDist, colour=burstID))+
  facet_wrap(~ID, scales="free") + scale_colour_manual(values=cols)+theme(legend.position = "none")

# looks good

tripdat%>%group_by(burstID)%>%summarise(c=length(Latitude))%>%pull(c)%>%table() 
# most bursts just over 300 s ( 5 min). A few shorter ~50-270 sec, then a bunch 6 second or less. See where these are to decide whether to remove
short_bursts<-tripdat%>%group_by(burstID)%>%summarise(c=length(Latitude))%>%filter(c<50)%>%pull(burstID)

# convert to spatial, small trips
dat_sf<-st_as_sf(tripdat%>%filter(burstID %in%short_bursts),
coords=c("Longitude", "Latitude"), crs=4326)
#view in tmap
tmap_mode("view")
tm_shape(dat_sf)+tm_dots(col="burstID")+tm_mouse_coordinates()

# yes remove bursts < 6 seconds
tripdat<-tripdat%>%filter(!burstID%in%short_bursts)

#### ---- ####

#### classify flying/sitting points using speed ####

ggplot(data=tripdat%>%filter(is.finite(speed_diff)))+geom_histogram(aes(x=speed_diff))+facet_wrap(~tripID)

#average FLYING speed
tripdat%>%filter(tripID!="-1"&is.finite(speed_diff) & ColDist>1 & speed_diff>4)%>%
  summarise(mean_s=mean(speed_diff), sd_s=sd(speed_diff))
#mean_s     sd_s
#1 13.4574 4.1042

#plot non-island points 
ggplot()+
  geom_histogram(data=tripdat%>%filter(tripID!="-1"&is.finite(speed_diff) & ColDist>1),
                 aes(x=speed_diff, fill=ifelse(speed_diff<4, "blue", "red")), colour=1, binwidth=1, boundary=0)+
  geom_label(aes(x=13, y=15000, label="Mean flying speed = 13.46±4.1 m/s"), size=5)+
  scale_x_continuous(breaks=0:40)+theme_bw()+ylab("Count of GPS datapoints")+xlab("speed m/s")+theme(legend.position="none")

tripdat$sit_fly<-"fly"
tripdat[tripdat$speed_diff<4,]$sit_fly<-"sit"

#### ---- ####

#### classify behaviour using EmBC and add day/night classification ####

forembc <- tripdat[,c('DateTime_AEDT', 'Longitude', 'Latitude', 'tripID')] #time, longitude, latitude, ID ( ! order is important !)

BC <- stbc(forembc)
sctr(BC)
stts(BC)

smoothedBC <- smth(BC, dlta=1)
view(smoothedBC)

tripdat$embc <- smoothedBC@A
tripdat$embc <- gsub("2","foraging",tripdat$embc)
tripdat$embc <- gsub("1","resting",tripdat$embc)
tripdat$embc <- gsub("3","commuting",tripdat$embc)
tripdat$embc <- gsub("4","relocating",tripdat$embc)
tripdat$embc <- gsub("5","DD",tripdat$embc)

# convert to spatial
dat_sf<-st_as_sf(tripdat, coords=c("Longitude", "Latitude"), crs=4326)
#view in tmap
tmap_mode("view")
tm_shape(dat_sf)+tm_dots(col="embc")+tm_mouse_coordinates()

# Add day night classification.. slow loop
tripdat$daynight<-"day"
for(i in 1:nrow(tripdat)){
dn1<-getSunlightTimes(date = as.Date(tripdat$DateTime_AEDT[i]), lat = tripdat$Latitude[i], lon = tripdat$Longitude[i],tz = "Australia/Sydney")
if(tripdat$DateTime_AEDT[i]>dn1$nauticalDusk  &   tripdat$DateTime_AEDT[i]<dn1$nightEnd){tripdat$daynight[i]<-"night"};print(i)}                       

#### ---- ####

#### Manually classify each burst based on behavior ####

# make burst summary table

burst_summary<-tripdat%>%group_by(burstID)%>%summarise(
  burst_dur1=sum(tdiff),
  burst_dur2=last(DateTime_AEDT)-first(DateTime_AEDT),
  ave_tdiff=mean(tdiff),
  max_tdiff=max(tdiff),
  ave_speed=mean(speed_diff),
  ave_colD=mean(ColDist),
  sum_dist=sum(dist_diff))

# classes
# B = bad (remove)
# I = Inspect (needs something fixing)
# C = Colony (colony control point)
# S = Sit (whole burst is sitting/drifting)
# T = Transit (whole burst is directed flight)
# L = Looping (During burst bird displays searching behavior, breaking from transiting flight)
# A = Alight (sitting AND Flying observed in burst)


burst_summary$class<-"NA"
for(i in unique(tripdat$burstID))
{
  print(tm_basemap("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png")+
          tm_shape(dat_sf%>%filter(burstID==i))+tm_dots(col="sit_fly"))
  print(burst_summary%>%filter(burstID==i))
  ans1<-readline("enter burst class:")
  burst_summary[burst_summary$burstID==i,]$class<-ans1
  print(burst_summary[burst_summary$burstID==i,])
}

# explore 4 'I' classifications
burst_summary[burst_summary$class=="I",]$burstID

# Tidying code
tripdat<-tripdat%>%filter(!(burstID=="08611854_01_66" &ColDist<500)) # drop some near col pts from burst
burst_summary[burst_summary$burstID=="08611854_01_66",]$class<-"S"

burst_summary[burst_summary$burstID=="08611854_02_73",]$class<-"T" # just has gap in middle
burst_summary[burst_summary$burstID=="08611854_04_121",]$class<-"S" # just has gap in middle
burst_summary[burst_summary$burstID=="08611854_06_155",]$class<-"L" # just has gap in middle

#write burst summary
#write.csv(burst_summary, "analyses/burst_summary_dat_all.csv", quote=F, row.names = F)

#write tripdat to upload to Movebank for for ENV-data extraction 
#write.csv(tripdat[c("Latitude", "Longitude","DateTime_UTC","DateTime_AEDT", "ID", "tripID")], "analyses/tripdat_4_movebank_all.csv", quote=F, row.names = F)

#attrib to dat_sf and write out
#dat_sf<-left_join(dat_sf, burst_summary[c("burstID", "class")], by="burstID")
#st_write(dat_sf, "analyses/GIS/tripdat_all_behav_burst_summary.shp")

#### ---- ####

#### Import ENV-data attributed from Movebank and join ####

move_env<-read.csv("data/shy_albatross_island/ENV_data Movebank Shy Albatross Bass Strait-22228471852684413.csv")

# also read in spatial data with full data stored, but drop spatial part

tripdat<-st_read("analyses/GIS/tripdat_all_behav_burst_summary.shp")%>%st_drop_geometry()

# rename varibles we care about for analyses
names(move_env)[names(move_env)=="ECMWF.ERA5.SL.Mean.Sea.Level.Pressure"]<-"mean_sea_level_pressure"
names(move_env)[names(move_env)=="ECMWF.ERA5.SL.Surface.Air.Pressure"]<-"surface_air_pressure"
names(move_env)[names(move_env)=="ETOPO1.Elevation"]<-"bathy"
names(move_env)[names(move_env)=="NASA.Distance.to.Coast"]<-"dist2coast"

# convert wind u and v to speed and direction

wind_df<-as.data.frame(uv2ds(move_env$ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.,
               move_env$ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component.))
names(wind_df)[1]<-"wind_dir"
names(wind_df)[2]<-"wind_speed"

# combine with tracking and calc bird-wind bearing difference, then classify into tail, cross and head-wind

tripdat<-cbind(tripdat, wind_df)

tripdat$b_w_diff<-ifelse(abs(tripdat$bearing_diff - tripdat$wind_dir)>180, 
                          180-(abs(tripdat$bearing_diff - tripdat$wind_dir)-180), abs(tripdat$bearing_diff - tripdat$wind_dir))

tripdat$b_w_class<-cut(tripdat$b_w_diff, breaks=c(0, 45, 135, 180), labels=c("tailwind", "crosswind","headwind"))

# add other env data of interest to tripdat

tripdat<-cbind(tripdat,move_env[c("mean_sea_level_pressure" , "surface_air_pressure","bathy","dist2coast")])

#### ---- ####

# select important columns and export tripdat

# remember to also select in burst class!

#write.csv(tripdat%>%dplyr::select(c(4:6, 9:14, 20, 21, 24:28, 32:48)), "analyses/tripdat_4_analyses.csv", quote=F, row.names=F)

