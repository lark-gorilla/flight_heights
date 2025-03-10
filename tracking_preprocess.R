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
#st_write(dat_sf, "analyses/GIS/alldat_almost_raw.shp")
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

#average FLYING speed - NOW USING GPS SPEED (knots conversion)
tripdat%>%filter(tripID!="-1"& ColDist>1 & (speed/1.94384)>4)%>%
  summarise(mean_s=mean(speed/1.94384), sd_s=sd(speed/1.94384))
#mean_s     sd_s
#1 13.37805 3.980755

#plot non-island points 
ggplot()+
  geom_histogram(data=tripdat%>%filter(tripID!="-1" & ColDist>1),
                 aes(x=speed/1.94384, fill=ifelse((speed/1.94384)<4, "blue", "red")), colour=1, binwidth=1, boundary=0)+
  geom_label(aes(x=13, y=15000, label="Mean flying speed = 13.38±3.98 m/s"), size=5)+
  scale_x_continuous(breaks=0:32)+theme_bw()+ylab("Count of GPS datapoints")+xlab("speed m/s")+theme(legend.position="none")

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
dat$daynight<-"night"
for(i in 1:nrow(dat)){
  dn1<-getSunlightTimes(date = as.Date(dat$DateTime_AEDT[i],tz = "Australia/Sydney"), lat = dat$Latitude[i], lon = dat$Longitude[i],tz = "Australia/Sydney")
  if(dat$DateTime_AEDT[i]<dn1$nauticalDusk  &   dat$DateTime_AEDT[i]>dn1$nightEnd){dat$daynight[i]<-"day"};print(i)}                       

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

#attrib burst summary
tripdat<-left_join(tripdat, burst_summary[c("burstID", "class")], by="burstID")

#write tripdat to upload to Movebank for for ENV-data extraction 
#write.csv(tripdat[c("Latitude", "Longitude","DateTime_UTC","DateTime_AEDT", "ID", "tripID")], "analyses/tripdat_4_movebank_all.csv", quote=F, row.names = F)

#attrib to dat_sf and write out
#dat_sf<-left_join(dat_sf, burst_summary[c("burstID", "class")], by="burstID")
#st_write(dat_sf, "analyses/GIS/tripdat_all_behav_burst_summary.shp")

#### ---- ####

#### Import ENV-data attributed from Movebank and join ####

move_env<-read.csv("data/shy_albatross_island/shy albatross 4 v2-4056681552610937919.csv")
tripdat<-read.csv("analyses/tripdat_4_movebank_all.csv", h=T)

# ~300 less rows in movebank than tripdat so do join to work out
tripdat$DateTime_AEDT<-as_datetime(tripdat$DateTime_AEDT, tz="Australia/Sydney")
move_env$DateTime_AEDT<-with_tz(as_datetime(move_env$timestamp, tz="UTC"), "Australia/Sydney")

tripdat_env<-left_join(tripdat, move_env, by=c("ID"="tag.local.identifier", "DateTime_AEDT"))
#~300 rows missed off end during join. but all crap pressure sensor data so no problem

tripdat_env<-tripdat_env%>%select(!c("X", "Y", "time", "iccid", "time_UTC", "event.id","visible" , "timestamp" ,                                            
                                      "location.long","location.lat" ,"sensor.type"  , "individual.taxon.canonical.name",                       
                                     "individual.local.identifier","study.name"))

# rename varibles we care about for analyses
names(tripdat_env)[names(tripdat_env)=="ECMWF.ERA5.SL.Mean.Sea.Level.Pressure"]<-"mean_sea_level_pressure"
names(tripdat_env)[names(tripdat_env)=="ECMWF.ERA5.SL.Surface.Air.Pressure"]<-"surface_air_pressure"
names(tripdat_env)[names(tripdat_env)=="ETOPO1.Elevation"]<-"bathy"
names(tripdat_env)[names(tripdat_env)=="NASA.Distance.to.Coast"]<-"dist2coast"
names(tripdat_env)[names(tripdat_env)=="MODIS.Ocean.Aqua.OceanColor.4km.8d.Daytime.SST"]<-"sst"
names(tripdat_env)[names(tripdat_env)=="ECMWF.ERA5.SL.Mean.Wave.Direction"]<-"wave_direction"
names(tripdat_env)[names(tripdat_env)=="ECMWF.ERA5.SL.10.Metre.Wind.Gust"]<-"windgust_10m"
names(tripdat_env)[names(tripdat_env)=="MODIS.Ocean.Aqua.OceanColor.4km.8d.Chlorophyll.A..OCI."]<-"chla"
names(tripdat_env)[names(tripdat_env)=="ECMWF.ERA5.SL.Significant.Wave.Height"]<-"wave_height"
names(tripdat_env)[names(tripdat_env)=="OSU.Ocean.NPP.0.083deg.8d.NPP"]<-"net_pp"
names(tripdat_env)[names(tripdat_env)=="ECMWF.ERA5.SL.Mean.Wave.Period"]<-"wave_period"


# convert wind u and v to speed and direction

wind_df<-as.data.frame(uv2ds(tripdat_env$ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.,
                             tripdat_env$ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component.))

tripdat_env$wind_dir<-wind_df$dir
tripdat_env$wind_speed<-wind_df$speed

tripdat_env$ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.<-NULL
tripdat_env$ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.<-NULL

# combine with tracking and calc bird-wind bearing difference, then classify into tail, cross and head-wind


tripdat_env$b_w_diff<-ifelse(is.na(tripdat_env$wind_dir), NA, ifelse(abs(tripdat_env$bearing_diff - tripdat_env$wind_dir)>180, 
                          180-(abs(tripdat_env$bearing_diff - tripdat_env$wind_dir)-180), abs(tripdat_env$bearing_diff - tripdat_env$wind_dir)))

tripdat_env$b_w_class<-cut(tripdat_env$b_w_diff, breaks=c(0, 45, 135, 180), labels=c("tailwind", "crosswind","headwind"))

#### ---- ####


#write.csv(tripdat_env, "analyses/tripdat_4_analyses_all.csv", quote=F, row.names=F)

## read back in and add MSLA data
library(terra)
library(ncdf4)
setwd("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights")
dat<-read.csv("analyses/tripdat_4_analyses_all.csv", h=T)
dat$DateTime_AEDT<-ymd_hms(dat$DateTime_AEDT, tz="Australia/Sydney")

nc_open("sourced_data/IMOS_aggregation_20240129T043941Z/IMOS_aggregation_20240129T043941Z.nc") #see what we're dealing with
#units: days since 1985-01-01 00:00:00 UTC

msla<-terra::rast("sourced_data/IMOS_aggregation_20240129T043941Z/IMOS_aggregation_20240129T043941Z.nc")

#get SLA product
msla<-msla["GSLA"]

sat_dat_ref<-with_tz(hours(as.numeric(gsub( "GSLA_TIME=", "", names(msla)))*24) + ymd_hms('1985-01-01 00:00:00', tz="UTC"), "Australia/Sydney")

df_tere_topo <- subset(msla, 1) %>%
  as.data.frame(xy = TRUE) %>%
  rename(msla = 'GSLA_TIME=13969.25')

dat_sf<-st_as_sf(dat, coords=c("Longitude", "Latitude"), crs=4326)

ggplot()+
  geom_raster(data = df_tere_topo, aes(x = x, y = y, fill = `msla`))+
  geom_sf(dat_sf, mapping = aes(), color = 'red', fill = NA)

# do extract

dat$msla=-999
dat$key<-paste(day(dat$DateTime_AEDT), month(dat$DateTime_AEDT))
for( i in unique(dat$key))
{
  ext1<-subset(msla, which(i==paste(day(sat_dat_ref), month(sat_dat_ref))))
  
  dat[dat$key==i,]$msla<-extract(ext1, dat[dat$key==i,c('Longitude', 'Latitude')])[,2]
}
# filling few NA bursts with nearst MSLA value  
table(dat[is.na(dat$msla),]$burstID)
dat[dat$burstID=='08611649_01_31',]$msla<-0.1644973
dat[dat$burstID=='08611649_01_16',]$msla<-mean(dat[dat$burstID=='08611649_01_15',]$msla)
dat[dat$burstID=='08611649_01_39',]$msla<-mean(dat[dat$burstID=='08611649_01_40',]$msla)

dat$key<-NULL

#citation: The citation in a list of references is: "IMOS [year-of-data-download], [Title], [data-access-URL], accessed [date-of-access]."

# could belnd with tidal data from jetty
# http://www.bom.gov.au/oceanography/projects/abslmp/data/data.shtml

#write.csv(dat, "analyses/tripdat_4_analyses_all.csv", quote=F, row.names=F)

## look into GPS bias - use all bursts before we remove some

#combine predeployment data with newer stuff

dat<-rbind(data.frame(read.csv("data/shy_albatross_island/gps_89460800120108611854.csv"), ID="08611854"),
           data.frame(read.csv("data/shy_albatross_island/gps_89460800120141490936.csv"), ID="41490936"),
           data.frame(read.csv("data/shy_albatross_island/gps_89460800120108611649.csv"), ID="08611649"))

dat<-dat[-80141,] # fix

dat$DateTime_UTC<-ymd_hms(dat$time_UTC, tz="UTC")
dat$DateTime_AEDT<-with_tz(dat$DateTime_UTC, "Australia/Sydney")

dat<-dat[dat$lat<0,]

dat$deployed_ID<-as.character(dat$ID)
dat[dat$DateTime_AEDT<ymd_hms("23-04-01 17:12:00", tz="Australia/Sydney") & dat$ID== "08611854",]$deployed_ID<-"predeployment"
dat[dat$DateTime_AEDT<ymd_hms("23-04-01 17:59:00", tz="Australia/Sydney") & dat$ID== "41490936",]$deployed_ID<-"predeployment"
dat[dat$DateTime_AEDT<ymd_hms("23-09-29 03:10:00", tz="Australia/Sydney") & dat$ID== "08611649",]$deployed_ID<-"predeployment"

pd_dat<-dat%>%filter(deployed_ID=='predeployment')

processed_dat<-read.csv("analyses/tripdat_4_analyses_all.csv", h=T)
processed_dat$DateTime_AEDT<-ymd_hms(processed_dat$DateTime_AEDT, tz="Australia/Sydney")
processed_dat<-processed_dat%>%filter(class=='C')

ai_dem<-terra::rast("sourced_data/Albatross_is_DEM/TAS Government/DEM/2 Metre/Tasmania_Statewide_2m_DEM_14-08-2021.tif")

dat_sf<-processed_dat%>%st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(st_crs(ai_dem))
processed_dat$island_dem<-(extract(ai_dem, dat_sf))[,2]

dat_sf<-pd_dat%>%st_as_sf(coords=c("lon", "lat"), crs=4326)%>%
  st_transform(st_crs(ai_dem))
pd_dat$island_dem<-(extract(ai_dem, dat_sf))[,2]

ggplot()+
  geom_raster(data = ai_dem %>% as.data.frame(xy = TRUE)%>%rename( alt='Tasmania_Statewide_2m_DEM_14-08-2021'),
              aes(x = x, y = y, fill = `alt`))+
  geom_point(data=processed_dat%>%filter(!is.na(island_dem)), aes(x=Longitude, y=Latitude), shape=1)+
  geom_point(data=pd_dat%>%filter(!is.na(island_dem)), aes(x=lon, y=lat), shape=1, col='red')

d1<-rbind(pd_dat%>%filter(!is.na(island_dem))%>%dplyr::select(ID, alt, island_dem, vdop, hdop, nSats),
          processed_dat%>%filter(!is.na(island_dem))%>%dplyr::select(ID, alt, island_dem, vdop, hdop, nSats))

ggplot(data=d1)+geom_point(aes(x=island_dem, y=alt, colour=ID))

d1$alt_diff=d1$alt-d1$island_dem

# difference on colony per device
d1[d1$ID=='08611649',]$ID<-'8611649'
d1[d1$ID=='08611854',]$ID<-'8611854'
d1%>%group_by(ID)%>%summarise(n=n(), med_d=median(alt_diff), mn_d=mean(alt_diff))
summary(d1$alt_diff)

summary(d1$alt_diff-9) #brings difference to within 0.3m of DEM values for island

# difference sitting on water per device
processed_all_dat<-read.csv("analyses/tripdat_4_analyses_all.csv", h=T)
processed_all_dat$DateTime_AEDT<-ymd_hms(processed_all_dat$DateTime_AEDT, tz="Australia/Sydney")
processed_all_dat%>%filter(class=='S' & sit_fly=='sit')%>%
  group_by(ID)%>%summarise(n=n(), med_d=median(alt), mn_d=mean(alt))
summary(processed_all_dat%>%filter(class=='S' & sit_fly=='sit')%>%pull(alt))

ggplot(data=processed_all_dat%>%filter(class=='S' & alt<200))+
  geom_histogram(aes(x=alt), binwidth=1)+facet_wrap(~ID, scales='free')+
  geom_vline(xintercept=7, col='red')+
  geom_vline(xintercept=8, col='blue')+
  geom_vline(xintercept=9, col='green')

ggplot(data=processed_all_dat%>%filter(class=='S' & alt<200))+
  geom_histogram(aes(x=alt), binwidth=1)+
  geom_vline(xintercept=7, col='red')+
  geom_vline(xintercept=8, col='blue')+
  geom_vline(xintercept=9, col='green')

ggplot(data=processed_all_dat%>%filter(alt<200))+
  geom_density(aes(x=alt-9, colour=class)) # ok -9 is the offset
