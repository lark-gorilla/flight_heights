
resample <- function(data=data, timeStep = 0.25) # note limit set to 12 for interpolation 
{
  
  out <- NULL
  Temp <- NULL
  
  TimeRemaining <- 0
  
  for(i in 1:nrow(data))
  {	
    
    O <- data[i,]
    B <- data[i+1,]
    
    #print(paste("at point ", i, " ", TimeRemaining, " Time Remains"))
    
    if(TimeRemaining == 0)
    {
      print(O)
      out <- rbind(out, O)
    }	
    
    if(i >= nrow(data)-1) break
    
    OBDist <- disttraj(O$Latitude, O$Longitude, B$Latitude, B$Longitude)
    OBTime <- (B$TrackTime - O$TrackTime)/3600
    StepDist <- ((OBDist/OBTime) * timeStep)
    
    if(OBTime > 12) next
    
    if(OBDist == 0) SplitRatio <- 0 else SplitRatio <- StepDist / OBDist
    
    if(B$Longitude < -90 & O$Longitude > 90) 
    {	
      LongDiff <- ((B$Longitude - O$Longitude) + 360);
      print(paste("Track crossed the dateline W > E at point ", i, sep=""))
    } else
      if(B$Longitude > 90 & O$Longitude < -90)	
      {
        LongDiff <- (360 - (B$Longitude - O$Longitude));
        print(paste("Track crossed the dateline E > W at point ", i, sep=""))
      } else 
        LongDiff <- (B$Longitude - O$Longitude)
    
    StepLat <- (B$Latitude - O$Latitude) * SplitRatio
    StepLong <- LongDiff * SplitRatio
    StepTrackTime <- (B$TrackTime - O$TrackTime) * SplitRatio
    
    Temp <- B
    
    Lat <- O$Latitude
    Long <- O$Longitude
    TrackTime <- O$TrackTime
    
    if(TimeRemaining > 0)
    {
      RemStep <- StepDist * (TimeRemaining/3600)
      if(OBDist == 0) RemRatio <- 0 else RemRatio <- RemStep / OBDist
      RemStepLat <- (B$Latitude - O$Latitude) * RemRatio
      RemStepLong <- LongDiff * RemRatio
      
      Lat <- Lat + RemStepLat
      Long <- Long + RemStepLong
      TrackTime <- TrackTime + TimeRemaining
      
      Temp$Latitude <- Lat
      Temp$Longitude <- Long
      Temp$TrackTime <- TrackTime
      
      if(TrackTime < B$TrackTime)
      {
        out <- rbind(out, Temp)
      }
    }
    
    
    while(TrackTime < B$TrackTime)
    {
      Lat <- Lat + StepLat
      Long <- Long + StepLong
      TrackTime <- TrackTime + (timeStep * 3600)
      
      if(TrackTime >= B$TrackTime) 
      {
        break
      } else
        
        Temp$Latitude <- Lat;
      Temp$Longitude <- Long;
      Temp$TrackTime <- TrackTime;
      
      out <- rbind(out, Temp)
    }
    
    #plot(Latitude~Longitude, data=out, asp=1)
    
    TimeRemaining <- TrackTime - B$TrackTime
  }
  
  return(out)
  
}





#### Great Circle Distance Calculation for Latitude and Longitude values ####


 disttraj <- function (lat1, long1, lat2, long2) 	## Great Circle Distance
 	{
    rad <- pi/180
    a1 <- lat1 * rad
    a2 <- long1 * rad
    b1 <- lat2 * rad
    b2 <- long2 * rad
    dlon <- b2 - a2
    dlat <- b1 - a1
    a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
    c <- 2 * atan2(sqrt(a), sqrt(1 - a))
    R <- 40003/(2 * pi)
    d <- R * c
    return(d)
 	}



 backforVel <- function(trip, point, n=4, alt=FALSE, filter=FALSE)	## Average velocity calculation, can switch Lat&Longs and ignores Comment == "Remove"
 	{
 library(geosphere)
 
 if(!"Latitude" %in% names(trip)) {print("No Latitude field")}
 if(!"Longitude" %in% names(trip)) {print("No Longitude field")}
 if(!"TrackTime" %in% names(trip)) {print("No TrackTime field")}
 
 if("geosphere" %in% installed.packages()) {} else {install.packages("geosphere")}
 library(geosphere)

 mf <- n/2
 mb <- n/2
 
 pnt.cntr <- data.frame(Latitude=0, Longitude=0, TrackTime=0)

 if( alt == TRUE )
	{
	pnt.cntr$Latitude <- trip[point,]$AltLat;
	pnt.cntr$Longitude <- trip[point,]$AltLon;
	pnt.cntr$TrackTime <- trip[point,]$TrackTime;
	} else 
	{pnt.cntr$Latitude <- trip[point,]$Latitude;
	pnt.cntr$Longitude <- trip[point,]$Longitude;
	pnt.cntr$TrackTime <- trip[point,]$TrackTime}
 
 ## forward

 if(point == 1 | point == nrow(trip)) return(0) else if(point == 2 | point == (nrow(trip) - 1)) m <- 1

 fwd.cnt <- 0
 fwd.vel <- 0
 fwd.inc <- 0 

 while(fwd.inc < mf)
	{
	fwd.cnt <- fwd.cnt + 1
	if((point + fwd.cnt) > nrow(trip)) {mf <- mf - 1; next}
	nxt.pnt <- trip[point + fwd.cnt,]
	if(filter == TRUE) {if(nxt.pnt$Comment == "Remove") {next}}
	fwd.dist <- distCosine(data.frame(pnt.cntr$Longitude, pnt.cntr$Latitude), data.frame(nxt.pnt$Longitude, nxt.pnt$Latitude))/1000
	fwd.time <- Mod((nxt.pnt$TrackTime - pnt.cntr$TrackTime)/3600)
	vel <- fwd.dist/fwd.time
	vel.sqr <- vel^2
	fwd.vel <- sum(fwd.vel, vel.sqr)
	fwd.inc <- fwd.inc + 1
	}

  ## backward

 bck.cnt <- 0
 bck.vel <- 0
 bck.inc <- 0 

 while(bck.inc < mb)
	{
	bck.cnt <- bck.cnt + 1
	if((point - bck.cnt) < 1) {mb <- mb - 1; next}
	nxt.pnt <- trip[point - bck.cnt,]
	if(filter == TRUE) {if(nxt.pnt$Comment == "Remove") {next}}
	bck.dist <- distCosine(data.frame(pnt.cntr$Longitude, pnt.cntr$Latitude), data.frame(nxt.pnt$Longitude, nxt.pnt$Latitude))/1000
	bck.time <- Mod((nxt.pnt$TrackTime - pnt.cntr$TrackTime)/3600)
	vel <- bck.dist/bck.time
	vel.sqr <- vel^2
	bck.vel <- sum(bck.vel, vel.sqr)
	bck.inc <- bck.inc + 1
	}
 
 ave.vel <- (sum(fwd.vel, bck.vel)/(mb + mf))
 vel <- sqrt(ave.vel)

 return(vel)

 }


