# Code to conduct analyses for manuscript: Altimeters on albatrosses quantifying flight heights for dynamic soaring seabirds

#### Housekeeping and libraries ####

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(rayshader)
library(patchwork)
library(lubridate)
library(viridis)
library(performance)
library(oceanwaves)
library(figpatch)
library(e1071)
library(zoo)
library(pracma)
library(mgcv)

setwd("C:/myloc")

#### ^^ ####

#### Read in data ####

# download from public repository: https://doi.org/10.26180/33829756

dat<-read.csv("Miller_BCI_Albatross.csv", h=T)

dat$datetime_UTC<-ymd_hms(dat$datetime_UTC, tz="UTC")
dat$DateTime_AEDT<-with_tz(dat$datetime_UTC, "Australia/Sydney")

#Add -9m correction to GPS elevation due to geoid model mismatch
dat$alt_gps<-dat$GPSaltitude-9

#remove first GPS fix of each burst as higher error
dat<-dat %>% group_by(burstID) %>%
  filter(row_number()!=1)%>%ungroup()%>%as.data.frame()

dat_flying<-dat%>%filter(burst_class=="flying")

#### ^^ ####

#### dynamic soaring cycle segmentation ####
dat_flying$ds_seg_pressure<-NA
dat_flying$ds_seg_gps<-NA
for(i in unique(dat_flying$burstID))
{
  i_burst<-dat_flying[dat_flying$burstID==i,]
  
  pres_smth <- rollmean(i_burst$pres_pa, k = 3, fill = NA) # apply moving window 3-pt
  gps_smth <- rollmean(i_burst$alt_gps, k = 3, fill = NA)
  
  # negative sign finds valleys (for pressure only)
  pres_valz<-findpeaks(-pres_smth, minpeakdistance = 5, nups=2, ndowns=2, zero='+') 
  # 5 sec = min peak-peak distance min DS duration from Schoombie et al 2023.
  gps_valz<-findpeaks(gps_smth, minpeakdistance = 5, nups=2, ndowns=2, zero='+')
  
  i_burst[pres_valz[,4]%>%sort(),]$ds_seg_pressure<-seq(1:length(pres_valz[,4]))
  i_burst[gps_valz[,4]%>%sort(),]$ds_seg_gps<-seq(1:length(gps_valz[,4]))
  
  i_burst$ds_seg_pressure[1]<-0
  i_burst<-i_burst%>%fill(ds_seg_pressure, .direction='down')
  i_burst$ds_seg_gps[1]<-0
  i_burst<-i_burst%>%fill(ds_seg_gps, .direction='down')
  
  dat_flying[dat_flying$burstID==i,]$ds_seg_pressure<-i_burst$ds_seg_pressure
  dat_flying[dat_flying$burstID==i,]$ds_seg_gps<-i_burst$ds_seg_gps
  print(i)
  
  #hashed code to view each burst output
  #cols = rainbow(nrow(pres_valz)+1, s=.6, v=.9)[sample(1:nrow(pres_valz)+1,nrow(pres_valz)+1,replace=T)]
  #p1<-ggplot()+geom_line(aes(x=i_burst$DateTime_AEDT, y=i_burst$pres_pa))+
  #  geom_point(aes(x=i_burst$DateTime_AEDT, y=i_burst$pres_pa, col=factor(i_burst$ds_seg_pressure)))+
  #  scale_colour_manual(values=cols)+scale_y_reverse()
  
  #cols = rainbow(nrow(gps_valz)+1, s=.6, v=.9)[sample(1:nrow(gps_valz)+1,nrow(gps_valz)+1,replace=T)]
  #p2<-ggplot()+geom_line(aes(x=i_burst$DateTime_AEDT, y=i_burst$alt_gps))+
  #  geom_point(aes(x=i_burst$DateTime_AEDT, y=i_burst$alt_gps, col=factor(i_burst$ds_seg_gps)))+
  #  scale_colour_manual(values=cols)
  #print(p1/p2)
  #readline("")
}

#summarise DS durations
dat_flying%>%group_by(burstID, ds_seg_pressure)%>%summarise(seg_dur=n())%>%
  ungroup()%>%summarise(mn_seg_dur=mean(seg_dur), sd_seg_dur=sd(seg_dur))

dat_flying%>%group_by(burstID, ds_seg_gps)%>%summarise(seg_dur=n())%>%
  ungroup()%>%summarise(mn_seg_dur=mean(seg_dur), sd_seg_dur=sd(seg_dur))

wilcox.test(x= dat_flying%>%group_by(burstID, ds_seg_gps)%>%summarise(seg_dur=n())%>%pull(seg_dur),
            y=dat_flying%>%group_by(burstID, ds_seg_pressure)%>%summarise(seg_dur=n())%>%
              group_by(burstID)%>%summarise(mn_seg_dur=mean(seg_dur))%>%pull(mn_seg_dur))

#Wilcoxon rank sum test with continuity correction
#W = 76295, p-value = 1.031e-06

#### ^^ ####

#### Set p0 for lower, upper and central scenarios ####

# use absolute maximum for p0 (upper scenario)
dat_flying<-dat_flying%>%group_by(burstID)%>%mutate(p0_mx=max(pres_pa)) 

# reset p0 for each dynamic soaring segment (above loop) (lower scenario)
temp_d<-dat_flying %>% group_by(burstID) %>% slice(c(1,1:(n()-1)))%>%ungroup()%>%arrange(birdID, burstID,ds_seg_pressure)
dat_flying<-dat_flying%>%ungroup()%>%arrange(birdID, burstID,ds_seg_pressure) %>%mutate(pres_pa_1=temp_d$pres_pa)
dat_flying<-dat_flying%>%group_by(burstID, ds_seg_pressure)%>%mutate(p0_ds_segSRT=max(pres_pa_1), p0_ds_segEND=max(pres_pa),
                                                                     p0_ds_seg=max(p0_ds_segSRT, p0_ds_segEND) ) 
#fit gam between lower and upper scenario p0s (central scenario)
dat_flying$p0_diff<-dat_flying$p0_mx-dat_flying$p0_ds_seg
dat_flying<-dat_flying%>%group_by(burstID)%>%mutate(p0_gam=
                                                      fitted(gam((p0_mx-(p0_diff/2))~s(DateTime_AEDT%>%as.numeric(), k=7))))%>%ungroup
#make sure no gam predictions higher than upper scenario
dat_flying$p0_gam<-ifelse(dat_flying$p0_gam>dat_flying$p0_mx, dat_flying$p0_mx, dat_flying$p0_gam)

#make sure no gam predictions lower than than lower scenario - removes negative values from gam prediction
dat_flying$p0_gam<-ifelse(dat_flying$p0_gam<dat_flying$p0_ds_seg, dat_flying$p0_ds_seg, dat_flying$p0_gam)

# sanity check
for(i in unique(dat_flying$burstID))
{
  i_burst<-dat_flying[dat_flying$burstID==i,]
  cols = rainbow(length(unique(i_burst$ds_seg_pressure)), s=.6, v=.9)%>%sample(length(unique(i_burst$ds_seg_pressure)))
  p1<-ggplot(data=i_burst)+geom_line(aes(x=DateTime_AEDT, y=pres_pa))+
    geom_point(aes(x=DateTime_AEDT, y=pres_pa, col=factor(ds_seg_pressure)))+
    scale_colour_manual(values=cols)+scale_y_reverse()+
    geom_line(aes(x=DateTime_AEDT, y=p0_mx), col=2)+
    geom_line(aes(x=DateTime_AEDT, y=p0_ds_seg), col=4)+
    geom_line(aes(x=DateTime_AEDT, y=p0_gam), col=5)+
    theme(legend.position='none')
  print(p1)
  readline("")
}

## Now calculate flight heights using barometric formula (Berberan Santos et al. 1997)
#h=((k*T)/(m*g))*ln(p/p0)
k=8.31432
m=0.0289644
g=9.80665
# *-1 flips negative/positive values
# calc altitudes from different p0s: max, min, most likely and also p0 error
dat_flying$alt_p0_mx<-(-1*((k*(dat_flying$temp_degC +273.15))/(m*g))*log(dat_flying$pres_pa/dat_flying$p0_mx))
dat_flying$alt_p0_ds_seg<-(-1*((k*(dat_flying$temp_degC +273.15))/(m*g))*log(dat_flying$pres_pa/dat_flying$p0_ds_seg))
dat_flying$alt_p0_gam<-(-1*((k*(dat_flying$temp_degC +273.15))/(m*g))*log(dat_flying$pres_pa/dat_flying$p0_gam))
dat_flying$alt_p0_error<-(-1*((k*(dat_flying$temp_degC +273.15))/(m*g))*log(dat_flying$p0_ds_seg/dat_flying$p0_mx))

# view error distribution

ggplot(data=dat_flying)+geom_histogram(aes(x=alt_p0_error), binwidth=0.5)+
  geom_vline(data=dat_flying%>%summarise(med_er=median(alt_p0_error)),aes(xintercept=med_er), col=3)+
  scale_x_continuous(breaks=0:11)+theme_bw()+xlab("Plausible range of altimeter flight height estimates (m)")+ylab("Count of datapoints")

ggplot(data=dat_flying)+geom_histogram(aes(x=alt_p0_error), binwidth=0.5)+
  geom_vline(data=dat_flying%>%group_by(birdID)%>%summarise(med_er=median(alt_p0_error)),aes(xintercept=med_er), col=3)+
  scale_x_continuous(breaks=0:11)+ theme_bw()+
  facet_wrap(~birdID, scales='free_y')+xlab("Plausible range of altimeter flight height estimates (m)")+ylab("Count of datapoints")

summary(dat_flying$alt_p0_error)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   1.085   2.087   2.376   3.270  10.842
dat_flying%>%group_by(birdID)%>%summarise(med_er=median(alt_p0_error))

#### ^^ ####

#### Summarise altitude from the three DS methods and GPS, and compare  ####

dat_comp<-rbind(data.frame(method='Altimeters (lower scenario)', Altitude=dat_flying$alt_p0_ds_seg, Logger=as.character(dat_flying$birdID), burstID=dat_flying$burstID),
                data.frame(method='Altimeters (upper scenario)', Altitude=dat_flying$alt_p0_mx, Logger=as.character(dat_flying$birdID), burstID=dat_flying$burstID),
                data.frame(method='Altimeters (central scenario)', Altitude=dat_flying$alt_p0_gam, Logger=as.character(dat_flying$birdID), burstID=dat_flying$burstID),
                data.frame(method='GPS', Altitude=dat_flying$alt_gps, Logger=as.character(dat_flying$birdID), burstID=dat_flying$birdID))

dat_comp$method<-factor(dat_comp$method, levels=c("Altimeters (lower scenario)",
                                                  "Altimeters (central scenario)",
                                                  "Altimeters (upper scenario)", "GPS"))
#summarise
dat_comp%>%group_by(method)%>%summarise(mn_alt=mean(Altitude), sd_alt=sd(Altitude), median=median(Altitude),
                                        min=min(Altitude), max=max(Altitude),
                                        q25=quantile(Altitude, 0.25), q75=quantile(Altitude, 0.75), skew=skewness(Altitude))


#method                        mn_alt sd_alt median   min   max   q25   q75   skew
#<fct>                          <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
#  1 Altimeters (lower scenario)     3.23   2.70   2.67     0  26.1  1.26  4.55 1.46  
#2 Altimeters (central scenario)   4.48   2.89   3.91     0  26.1  2.41  5.97 1.23  
#3 Altimeters (upper scenario)     5.61   3.24   4.97     0  26.1  3.28  7.32 1.09  
#4 GPS                             5.10   9.72   4      -63  89    0     9    0.0607

# make plot
cols <- c('#dc267f','#648fff','#fe6100','#ffb000')

cols.alpha<-c(grDevices::adjustcolor(cols[1], alpha.f = 1),
              grDevices::adjustcolor(cols[2], alpha.f = 1),
              grDevices::adjustcolor(cols[3], alpha.f = 1),
              grDevices::adjustcolor(cols[4], alpha.f = 1))

ggplot(data=dat_comp)+geom_histogram(aes(x=Altitude, fill=method),col=1, binwidth=1)+
  geom_vline(xintercept = 0, linetype='dotted')+
  scale_x_continuous(breaks=seq(-60,60,2))+
  coord_cartesian(xlim=c(-20, 40))+
  scale_fill_manual(values=cols.alpha)+labs(x="Flight height (m)", y="Count")+
  facet_wrap(~method, nrow=4, scales='free_y')+theme_bw()+ theme(legend.position = "none",
                                                                 strip.text = element_text(size = 10, face='bold'),
                                                                 axis.title = element_text(size = 12))
wilcox.test(x= dat_flying$alt_gps,
            y=dat_flying$alt_p0_gam, paired=T)

#Wilcoxon signed rank test with continuity correction

#data:  dat_flying$alt_gps and dat_flying$alt_p0_gam
#V = 167397686, p-value = 2.339e-05

wilcox.test(x= dat_flying$alt_gps,
            y=dat_flying$alt_p0_mx, paired=T)

#data:  dat_flying$alt_gps and dat_flying$alt_p0_mx
#V = 138154406, p-value < 2.2e-16

wilcox.test(x= dat_flying$alt_gps,
            y=dat_flying$alt_p0_ds_seg, paired=T)

#data:  dat_flying$alt_gps and dat_flying$alt_p0_ds_seg
#V = 198704999, p-value < 2.2e-16

#### ^^ ####

#### Zero-crossing analysis  ####

# Using Zero-crossing method in oceanwaves package to extract wave height and period from floating bursts

# Use floating data only
dat_floating<-dat%>%filter(burst_class=="floating")

# Use max p0 to calc wave height from pressure data - zero crossing doesn't care which p0 is used, 
# using barometric equation just converts units from pascals to metres to allow comparison with GPS
dat_floating<-dat_floating%>%group_by(burstID)%>%mutate(p0_mx=max(pres_pa)) 
k=8.31432
m=0.0289644
g=9.80665
dat_floating$alt_p0_mx<-(-1*((k*(dat_floating$temp_degC +273.15))/(m*g))*log(dat_floating$pres_pa/dat_floating$p0_mx))

zc_summary<-NULL
for ( i in unique(dat_floating$burstID))
{
  tout<-data.frame(class=unique(dat_floating[dat_floating$burstID==i,]$burst_class), burstID=i, 
                   ds_hsig=NA,ds_hmean=NA, ds_tmean=NA, ds_tsig=NA,
                   gps_hsig=NA,gps_hmean=NA, gps_tmean=NA, gps_tsig=NA)
  
  d1<-waveStatsZC(dat_floating[dat_floating$burstID==i,]$alt_p0_mx, 1,)
  tout$ds_hsig=d1$Hsig
  tout$ds_hmean=d1$Hmean
  tout$ds_tmean=d1$Tmean
  tout$ds_tsig=d1$Tsig
  
  possibleError <-tryCatch(
    waveStatsZC(dat_floating[dat_floating$burstID==i,]%>%filter(!alt_gps %in% boxplot(alt_gps, plot=F)$out)%>%pull(alt_gps), 1,),
    error=function(e) e)
  
  if(!inherits(possibleError, "error")){
    g1<-waveStatsZC(dat_floating[dat_floating$burstID==i,]%>%filter(!alt_gps %in% boxplot(alt_gps, plot=F)$out)%>%pull(alt_gps), 1,) 
    tout$gps_hsig=g1$Hsig
    tout$gps_hmean=g1$Hmean
    tout$gps_tmean=g1$Tmean
    tout$gps_tsig=g1$Tsig
  }else{}
  
  zc_summary<-rbind(zc_summary, tout)
}

ggplot(data=zc_summary)+geom_point(aes(x=ds_hmean, y=gps_hmean))+coord_cartesian()
ggplot(data=zc_summary)+geom_point(aes(x=ds_tmean, y=gps_tmean))

#### ^^^ ####

#### Measuring wave height w/ altimeters and making figure 5  ####

#summarise first

# Significant wave height = the average wave height of the top one-third highest waves

wave_temp<-dat_floating%>%group_by(burstID)%>%summarise(w_height=mean(ERA5_wave_hgt ),
                                                        w_period=mean(ERA5_wave_prd )) 

wave_sum<-left_join(zc_summary, wave_temp, by="burstID")

# lms tell us which outliers to remove b4 pearsons corr
w1<-lm(w_height~ds_hsig, data=wave_sum)
w2<-lm(w_height~gps_hsig, data=wave_sum)
w3<-lm(w_period~ds_tmean, data=wave_sum)
w4<-lm(w_period~gps_tmean, data=wave_sum)

check_model(w1)
check_model(w2)
check_model(w3)
check_model(w4)

wp1<-ggplot(data=wave_sum[-c(13, 17,18,32),])+geom_point(aes(y=w_height, x=ds_hsig))+theme_bw()+
  labs(x='Wave height from albatross altimeters (m)',y='Wave height from satellite ERA5 (m)', size=5)+
  geom_text(aes(x=7, y=1.4), label=expression(italic(r)*" = "*"0.58, "* italic(p) < 0.001), size=4)+
  theme(axis.text=element_text(size=12))

wp2<-ggplot(data=wave_sum[wave_sum$gps_hsig<6,])+geom_point(aes(y=w_height, x=gps_hsig))+theme_bw()+
  labs(x='Wave height from albatross GPS (m)',y='Wave height from satellite ERA5 (m)', size=5)+
  geom_text(aes(x=4, y=1.4), label=expression(italic(r)*" = "*"0.37, "* italic(p)*" = "*0.03), size=4)+
  theme(axis.text=element_text(size=12))

wp3<-ggplot(data=wave_sum[-c(13, 17,18,32),])+geom_point(aes(y=w_period, x=ds_tmean))+theme_bw()+
  labs(x='Wave period from albatross altimeters (s)',y='Wave period from satellite ERA5 (s)', size=5)+
  geom_text(aes(x=9, y=7), label=expression(italic(r)*" = "*"0.86, "* italic(p) < 0.001), size=4)+
  theme(axis.text=element_text(size=12))

wp4<-ggplot(data=wave_sum[wave_sum$gps_hsig<6,])+geom_point(aes(y=w_period, x=gps_tmean))+theme_bw()+
  labs(x='Wave period from albatross GPS (s)',y='Wave period from satellite ERA5 (s)', size=5)+
  geom_text(aes(x=30, y=7), label=expression(italic(p)*" = NS"), size=4)+
  theme(axis.text=element_text(size=12))

cor.test(x=wave_sum[-c(13, 17,18,32),]$w_height, y=wave_sum[-c(13, 17,18,32),]$ds_hsig, method='pearson', na.action=na.omit) 
cor.test(x=wave_sum[wave_sum$gps_hsig<6,]$w_height, y=wave_sum[wave_sum$gps_hsig<6,]$gps_hsig, method='pearson', na.action=na.omit)

cor.test(x=wave_sum[-c(13, 17,18,32),]$w_period, y=wave_sum[-c(13, 17,18,32),]$ds_tmean, method='pearson', na.action=na.omit)
cor.test(x=wave_sum[wave_sum$gps_hsig<6,]$w_period, y=wave_sum[wave_sum$gps_hsig<6,]$gps_tmean, method='pearson', na.action=na.omit)

# additional request for sitting plot from reviewer

sit_expl<-dat[dat$burstID=='41490936_01_10',]
rescale <- function(x_i){max(sit_expl$pres_pa)-x_i}
sit_expl$pres_pa_rev<-rescale(sit_expl$pres_pa)+max(sit_expl$pres_pa)

p1<-ggplot(sit_expl, aes(x=DateTime_AEDT)) +
  
  geom_line(aes(y=pres_pa_rev), color=1) + 
  geom_line(aes(y=(alt_gps*13.3)+102290), color='#ffb000') +
  
  scale_y_continuous(
    name = "Pressure (Pa) - reversed",
    sec.axis = sec_axis(~(.-102290)/13.3, name="GPS altitude (m)")) +
  theme_bw()+scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)') 


p1/(wp1+wp2)/(wp3+wp4)
#### ^^^ ####

#### Make  Fig 1 ####

p1<-ggplot(data=dat[dat$burstID=="08611854_04_122",])+geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa), size=1)+scale_y_reverse()+
  labs(y="Pressure (Pa) - reversed", x="Time")+
  geom_hline(yintercept=unique(dat[dat$burstID=="08611854_04_122",]$p0), colour='red')+theme_bw()+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)')+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17))

i_burst<-dat_flying[dat_flying$burstID=='08611854_04_122',]

cols=c(rep(c('black', "darkgray"), 12), 'black')
p1<-ggplot(data=i_burst)+geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1,colour=factor(ds_seg_pressure)), size=1)+
  scale_y_reverse()+
  labs(y="Pressure (Pa) - reversed", x="Time")+
  geom_line(aes(x=DateTime_AEDT, y=p0_mx), col='#fe6100', size=1)+
  geom_line(aes(x=DateTime_AEDT, y=p0_ds_seg), col='#dc267f', size=1)+
  geom_line(aes(x=DateTime_AEDT, y=p0_gam), col='#648fff', size=1)+
  geom_line(aes(x=DateTime_AEDT, y=p0_mx), col='#fe6100', linetype='dashed', size=1)+
  scale_colour_manual(values=cols)+
  theme_bw()+ guides(colour="none")+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)')+
  theme(panel.grid.minor=element_blank(),axis.text=element_text(size=15),axis.title=element_text(size=17))

p1.5<-ggplot(data=dat[dat$burstID=="08611854_04_122",])+
  geom_point(aes(x=longitude, y=latitude, colour=pres_pa))+scale_color_viridis(trans="reverse")+
  labs(x="Longitude", y="Latitude", colour="Pressure\n(Pa)\n\n")+theme_bw()
plot_gg(p1.5, height=4, width=8, pointcontract = 0.5, sunangle = 40) # make 3d figure

render_snapshot("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights/writeup/3dplot2.png", dpi=900, clear = T)

p2<-ggplot(data=i_burst)+geom_line(aes(x=DateTime_AEDT, y=alt_p0_gam, group=1), col='#648fff')+
  geom_point(aes(x=DateTime_AEDT, y=alt_p0_gam), size=1, col='#648fff')+
  geom_line(aes(x=DateTime_AEDT, y=alt_gps, group=1), col='#ffb000')+
  geom_point(aes(x=DateTime_AEDT, y=alt_gps), size=1, col='#ffb000')+
  labs(y="Altitude (m)")+
  geom_hline(yintercept=0, linetype='dotted')+theme_bw()+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=17),panel.grid.minor=element_blank())+
  scale_y_continuous(limits=c(-2, 26), breaks=seq(-2,26,2), minor_breaks = NULL)+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)')

p1.6 <- fig("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights/writeup/3dplot2.png")

wrap_plots(p1, p1.6, p2, nrow=3)
p1+p1.6+p2 + plot_layout(nrow=3, heights = c(2,2,1))+ 
  plot_annotation(tag_levels = 'a',tag_suffix = ')')&theme(plot.tag = element_text(size = 26)) # then export png @ 1400 x 1700 
# do manually
p1/p2

p1/p1.6/p2

#### ^^^ ####

#### Check suitability floating-satellite method (Johnston et al 2023) ####

dat<-dat%>%group_by(burstID)%>%mutate(ERA5_mn_sea_lvl_pres =mean(ERA5_mn_sea_lvl_pres ,na.rm = T))%>%
  ungroup()%>%as.data.frame() # average satellite pressure dat per burst

#get sitting pressure diff 
dat$sat_sit_pdiff<-NA
dat[dat$burst_class %in% c('takeoff_landing', 'floating') & dat$speed_ms<4,]$sat_sit_pdiff<-
  (dat[dat$burst_class %in% c('takeoff_landing', 'floating') & dat$speed_ms<4,]$pres_pa-
     dat[dat$burst_class %in% c('takeoff_landing', 'floating') & dat$speed_ms<4,]$ERA5_mn_sea_lvl_pres) 

# summarise per burst
dat<-dat%>%group_by(burstID)%>%
  mutate(burstID_sat_sit_pdiff=mean(sat_sit_pdiff,na.rm = T))%>%ungroup()%>%as.data.frame()

sumr<-NULL
for(i in unique(dat$burstID))
{
  dtemp<-dat%>%filter(burstID==i)
  IDtemp<-dat%>%filter(birdID==unique(dtemp$birdID)) # get nearest from same logger
  sit_burst<-IDtemp[! is.na(IDtemp$burstID_sat_sit_pdiff),] # only bursts with sitting diffs
  
  out1<-data.frame(burstID=i, birdID=unique(dtemp$birdID), burst_class=unique(dtemp$burst_class), 
                   min_t=min(abs((sit_burst$DateTime_AEDT- 
                                    median(dtemp$DateTime_AEDT))), na.rm=T))
  sumr<-rbind(sumr, out1)
}

sumr$min_t<-as.numeric(sumr$min_t)
ggplot(data=sumr%>%filter(burst_class=='flying'))+geom_histogram(aes(x=min_t/3600))+facet_wrap(~birdID)
summary(sumr[sumr$burst_class=='flying',]$min_t/3600)
#Average time gap between nearest flying and floating busts > 6 hrs - too long delay for applying Johnston method 
#### ^^^ ####