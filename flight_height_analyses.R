# Code to conduct paper analyses.

# note: remove first and/or last point of each burst?


library(ggplot2)
library(dplyr)
library(rayshader)

setwd("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights")

# read data

dat<-read.csv("analyses/tripdat_4_analyses.csv", h=T)

table(dat$class)

## calc flight heights, do all but only ware about T and A burst classes

# barometric formula (Berberan Santos et al. 1997)
#h=((k*T)/(m*g))*ln(p/p0)
k=8.31432
m=0.0289644
g=9.80665

dat$pres_alt<-NA
for (i in unique(dat$burstID))
{
  # get 95% CI - want upper as high pressure = low altitude 
  burst_pres<-dat[dat$burstID==i,]$pres_pa
  a <- mean(burst_pres)
  s <- sd(burst_pres)
  n <- length(burst_pres)
  error <- qnorm(0.975)*s/sqrt(n)
  lwr_CI<- a-error
  
  #cant get ci to work using min for the moment
  lwr_CI<-max(burst_pres)
  
  dat[dat$burstID==i,]$pres_alt<-abs(
    ((k*(dat[dat$burstID==i,]$temp+273.15))/(m*g))*log(dat[dat$burstID==i,]$pres_pa/lwr_CI))
}


p0 = 1013.25# not a constant but rough estimate of sea level pressure # 1016 for melbourne at time?
#p and p0 in mbar
#T in kelvin
dat$p0_melb<-1016
dat[ grep(" 23:", dat$Time_orig), ]$p0_melb<-1015
dat[ grep(" 00:", dat$Time_orig), ]$p0_melb<-1015

dat$temp_K=dat$temp+ 273.15

dat$height=((k*dat$temp_K)/(m*g))*log(dat$pressure/p0)

# 3d plot

p1<-ggplot(data=dat[dat$burstID==i,])+geom_point(aes(x=Longitude, y=Latitude, colour=pres_alt))

plot_gg(p1)
