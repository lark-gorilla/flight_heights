# Code to conduct paper analyses.

library(ggplot2)
library(dplyr)
library(rayshader)
library(patchwork)
library(lubridate)
library(viridis)
library(sf)
library(nlme)
library(emmeans)
library(performance)
library(oceanwaves)
library(gridExtra)
library(figpatch)
library(e1071)
library(fitdistrplus)

setwd("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights")

# read data

dat<-read.csv("analyses/tripdat_4_analyses_all.csv", h=T)
burst_summary<-read.csv("analyses/burst_summary_dat_all.csv", h=T)

dat$DateTime_AEDT<-ymd_hms(dat$DateTime_AEDT, tz="Australia/Sydney")
dat<-dat%>%filter(deployed_ID!="predeployment")     

table(dat$burstID, dat$class)

bad_ids<-c(
  "-1_1",  "08611854_01_30","08611854_01_31", "08611854_01_33", "08611854_01_34",
  "08611854_01_35",  "08611854_01_41", "08611854_01_42", "08611854_01_46",  "08611854_01_47",
  "08611854_01_48",  "08611854_01_49",  "08611854_01_50",  "08611854_01_51",
  "08611854_01_52", "08611854_01_55", "08611854_01_65",  "-1_66",  "08611854_02_79",
  "08611854_02_82",  "08611854_02_85",  "08611854_02_86",  "08611854_02_87",
  "08611854_02_91",  "08611854_02_92",  "08611854_03_99",  "08611854_03_101",
  "08611854_03_102",  "08611854_03_103",  "08611854_03_110",  "-1_111",
  "08611854_04_123",  "08611854_04_124",  "08611854_04_125", # good example of bobbing, then swamp
  "08611854_04_126",  "08611854_05_128",  "08611854_05_129",  "08611854_05_130",
  "08611854_05_131",  "08611854_05_132",  "-1_133",  "08611854_06_134", # same
  "08611854_06_136",  "08611854_06_139",  "08611854_06_140",  "08611854_06_153",
  "08611854_06_157",  "08611854_06_161",  "08611854_06_162",  "08611854_06_163",
  "08611854_06_164",  "08611854_06_166",  "08611854_06_167",  "08611854_06_168",
  "08611854_06_169",  "08611854_06_172",  "08611854_06_175",  "08611854_06_176",
  "08611854_06_178",  "08611854_06_183",  "08611854_06_185",  "08611854_06_188",
  "08611854_07_191",  "08611854_07_192",  "08611854_07_195",  "08611854_07_196",
  "08611854_08_199",  "08611854_08_200",  "08611854_08_201",  "41490936_01_16",
  "08611649_01_11",  "08611649_01_12",  "08611649_01_16",  "08611649_01_16",
  "08611649_01_18",  "08611649_01_19",  "08611649_01_21",  "08611649_01_25",
  "08611649_01_26",  "08611649_01_27",  "08611649_01_28",  "08611649_01_29",
  "08611649_01_37", # "08611649_01_39" allowed thru but includes a bit
  "08611649_01_40",  "08611649_01_41", "08611649_01_17")

length(unique(dat$burstID));length(bad_ids)
  
dat<-dat%>%filter(!burstID %in% bad_ids) # could do extra check based on min/max burst pressure difference  to see if any missed

## Few edits
dat<-dat%>%filter(burstID!='08611854_01_66') # remove small 80 sec burst
dat[dat$burstID=="08611649_01_36",]$class<-"A"
dat<-dat%>%filter(burstID!="-1_93") # remove over the island bird
dat[dat$burstID=="08611854_06_160",]$class<-"S"
dat[dat$burstID=="08611854_01_59",]$class<-"A"
dat[dat$burstID=="08611854_03_100",]$class<-"A"
dat[dat$burstID=="08611854_06_156",]$class<-"A"
#summary(fitted(lm(alt~vdop, data=dat%>%filter(class=='S'&alt<200))))
dat<-dat%>%filter(!burstID %in%
c('08611854_02_73', '08611854_04_121', '08611854_06_155')) # remove 3 bursts with crazy GPS alt
# removing some small swamp sections
dat<-dat%>%filter(!(burstID=="08611854_06_160" & DateTime_AEDT>ymd_hms("2023-04-14 16:59:35", tz="Australia/Sydney")))
dat<-dat%>%filter(!(burstID=="08611649_01_39" & DateTime_AEDT>ymd_hms("2023-10-02 07:12:25", tz="Australia/Sydney")))
#

#### summary table of env variables ####
tabl1<-dat%>%filter(class %in% c("T", "L", "A"))%>%dplyr::select(ColDist , dist2coast, wind_speed, wind_dir,
                                                          chla, sst, wave_height, wave_period)%>%
                            summarise_all(list(mean=mean, sd=sd, qz=quantile), na.rm=T)

#write_xlsx(tabl1, 'analyses/env_summary_table.xlsx')

dat%>%filter(class %in% c("T", "L", "A"))%>%group_by(b_w_class)%>%summarise(n())
#b_w_class `n()`
#<chr>     <int>
#  1 crosswind 20224
#2 headwind   7540
#3 tailwind   4249

dat%>%filter(class %in% c("T", "L", "A"))%>%group_by(daynight)%>%summarise(n())
# all flights in day!
#### ^^ ####

#### Calculation of flight height using dynamic soaring method and correction to GPS ####

# barometric formula (Berberan Santos et al. 1997)
#h=((k*T)/(m*g))*ln(p/p0)
k=8.31432
m=0.0289644
g=9.80665

#Add -9m correction to GPS elevation 
dat$alt_gps<-dat$alt-9

dat$index<-1:nrow(dat)
dat$p0<-0

dat$alt_DS<-NA
for (i in unique(dat$burstID))
{
  # original method - 95% upper quantile of pressure to set p0 for entire burst
  dat[dat$burstID==i,]$p0<-quantile(dat[dat$burstID==i,]$pres_pa, probs=0.95)
}
  
dat$alt_DS<-(-1*  # *-1 flips negative/positive values
                                    ((k*(dat$temp+273.15))/(m*g))*log(dat$pres_pa/dat$p0))

#### ^^^ ####

#### Sensitivity analysis of dynamic soaring p0 threshold using wave heights ####
dat_sit<-dat%>%filter(class=='S')

for (j in seq(0.7, 1, 0.01))
{
dat_sit$alt_DS<-NA
for (i in unique(dat_sit$burstID))
{
  # original method - 95% upper quantile of pressure to set p0 for entire burst
  dat_sit[dat_sit$burstID==i,]$p0<-quantile(dat_sit[dat_sit$burstID==i,]$pres_pa, probs=j)
}

dat_sit$alt_DS<-(-1*  # *-1 flips negative/positive values
               ((k*(dat_sit$temp+273.15))/(m*g))*log(dat_sit$pres_pa/dat_sit$p0))

names(dat_sit)[which(names(dat_sit)=="alt_DS")]<-paste0("alt_DS", j)
}

# Using Zero-crossing method in oceanwaves package
zc_summary<-NULL

for(i in unique(dat_sit$burstID))
{
  tout<-data.frame(burstID=i, DSp0=seq(0.7, 1, 0.01), Hsig=NA)
   
  for (j in names(dat_sit)[54:84])
  {
  d1<-waveStatsZC(dat_sit[dat_sit$burstID==i, j], 1,)
  tout[which(tout$DSp0==substr(j, 7, nchar(j))),]$Hsig<-d1$Hsig
  }
  zc_summary<-rbind(zc_summary, tout)
}
  
  possibleError <-tryCatch(
    waveStatsZC(dat[dat$burstID==i,]%>%filter(!alt_gps %in% boxplot(alt_gps)$out)%>%pull(alt_gps), 1,),
    error=function(e) e)
  
  if(!inherits(possibleError, "error")){
    g1<-waveStatsZC(dat[dat$burstID==i,]%>%filter(!alt_gps %in% boxplot(alt_gps)$out)%>%pull(alt_gps), 1,) 
    tout$gps_hsig=g1$Hsig
    tout$gps_hmean=g1$Hmean
    tout$gps_tmean=g1$Tmean
    tout$gps_tsig=g1$Tsig
  }else{}
  
  zc_summary<-rbind(zc_summary, tout)
}

#### ^^ ####

#### Calculation of flight height using satellite ocean data (Johnston et al 2023) + 9m GPS offset (final line of code) ####

# Work out difference between pressure of sitting bursts and ECMWF.ERA5.SL.Mean.Sea.Level.Pressure
# then use value to calibrate satellite data to 'true' surface pressure (p0). Apply 'true'
# p0 value sitting bursts to flying bursts within 1 day. Nearest (in time) sitting burst has priority

# analysis run per logger as each will have unique sensor calibration

#### vis helping plots ^^^ ####
#ggplot(data=dat)+geom_point(aes(x=pres_pa, y=mean_sea_level_pressure))+
  #geom_abline()+facet_wrap(~class, scales="free")

#ggplot(data=dat%>%filter(class=="S")%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
#  geom_point(aes(x=index, y=pres_pa), col='red')+
#  facet_wrap(~ID, scales="free")

# check each individually
#p1<-ggplot(data=dat%>%filter(ID==08611649)%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
#  geom_point(aes(x=index, y=pres_pa, colour=embc, shape=class))
#p2<-ggplot(data=dat%>%filter(ID==08611649)%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
#  geom_point(aes(x=index, y=pres_pa, colour=sit_fly, shape=class))

#p1/p2

#p1<-ggplot(data=dat%>%filter(ID==8611854 )%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
#  geom_point(aes(x=index, y=pres_pa, colour=embc, shape=class))
#p2<-ggplot(data=dat%>%filter(ID==8611854 )%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
#  geom_point(aes(x=index, y=pres_pa, colour=sit_fly, shape=class))

#p1/p2

#p1<-ggplot(data=dat%>%filter(ID==41490936   )%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
#  geom_point(aes(x=index, y=pres_pa, colour=embc, shape=class))
#p2<-ggplot(data=dat%>%filter(ID==41490936   )%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
#  geom_point(aes(x=index, y=pres_pa, colour=sit_fly, shape=class))

#p1/p2
#### ^^^ ####

# Ok use sitting points 

#get sitting pressure diff 
dat$sat_sit_pdiff<-NA
dat[dat$class %in% c('A', 'S') & dat$sit_fly=='sit',]$sat_sit_pdiff<-
 (dat[dat$class %in% c('A', 'S') & dat$sit_fly=='sit',]$pres_pa-
      dat[dat$class %in% c('A', 'S') & dat$sit_fly=='sit',]$mean_sea_level_pressure) 

# summarise per burst
dat<-dat%>%group_by(burstID)%>%
  mutate(burstID_sat_sit_pdiff=mean(sat_sit_pdiff,na.rm = T))%>%ungroup()%>%as.data.frame()

dat$nearest_sat_sit_pdiff<-NA
for(i in unique(dat$burstID))
{
  dtemp<-dat%>%filter(burstID==i)
  IDtemp<-dat%>%filter(ID==unique(dtemp$ID)) # get nearest from same logger
  sit_burst<-IDtemp[! is.na(IDtemp$burstID_sat_sit_pdiff),] # only bursts with sitting diffs
  if(min(abs((sit_burst$DateTime_AEDT- 
               median(dtemp$DateTime_AEDT))))>hours(24)){next} #if no sitting within 1 day skip
  
 appl_diff<-sit_burst[which.min(abs((sit_burst$DateTime_AEDT-median(dtemp$DateTime_AEDT)))),]$burstID_sat_sit_pdiff

 dat[dat$burstID==i,]$nearest_sat_sit_pdiff<-appl_diff
}

# check outputs - reproduce Johnston fig 2
ggplot(data=dat%>%filter(ID==08611649)%>%mutate(index=1:nrow(.)))+geom_line(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_line(aes(x=index, y=pres_pa), colour='black')+
  geom_line(aes(x=index, y=ifelse(is.na(sat_sit_pdiff),mean_sea_level_pressure+nearest_sat_sit_pdiff,mean_sea_level_pressure+sat_sit_pdiff)), colour='orange')

ggplot(data=dat%>%filter(ID==8611854)%>%mutate(index=1:nrow(.)))+geom_line(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_line(aes(x=index, y=pres_pa), colour='black')+
  geom_line(aes(x=index, y=ifelse(is.na(sat_sit_pdiff),mean_sea_level_pressure+nearest_sat_sit_pdiff,mean_sea_level_pressure+sat_sit_pdiff)), colour='orange')

ggplot(data=dat%>%filter(ID==41490936)%>%mutate(index=1:nrow(.)))+geom_line(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_line(aes(x=index, y=pres_pa), colour='black')+
  geom_line(aes(x=index, y=ifelse(is.na(sat_sit_pdiff),mean_sea_level_pressure+nearest_sat_sit_pdiff,mean_sea_level_pressure+sat_sit_pdiff)), colour='orange')
# looks good

#calculate p0 and alt for satellite ocean data method
dat$p0_SO<-ifelse(is.na(dat$sat_sit_pdiff),dat$mean_sea_level_pressure+dat$nearest_sat_sit_pdiff,dat$mean_sea_level_pressure+dat$sat_sit_pdiff)
dat$alt_SO<-(-1*  # *-1 flips negative/positive values
               ((k*(dat$temp+273.15))/(m*g))*log(dat$pres_pa/dat$p0_SO))


#### ^^^ ####


#### Zero-crossing analysis  ####

# Using Zero-crossing method in oceanwaves package

zc_summary<-NULL
for ( i in unique(dat[dat$class%in%c("T", "L", "S"),]$burstID))
{
 tout<-data.frame(class=unique(dat[dat$burstID==i,]$class), burstID=i, 
                  ds_hsig=NA,ds_hmean=NA, ds_tmean=NA, ds_tsig=NA,
                  gps_hsig=NA,gps_hmean=NA, gps_tmean=NA, gps_tsig=NA)
  
 d1<-waveStatsZC(dat[dat$burstID==i,]$alt_DS, 1,)
  tout$ds_hsig=d1$Hsig
  tout$ds_hmean=d1$Hmean
  tout$ds_tmean=d1$Tmean
  tout$ds_tsig=d1$Tsig
  
  
  possibleError <-tryCatch(
    waveStatsZC(dat[dat$burstID==i,]%>%filter(!alt_gps %in% boxplot(alt_gps)$out)%>%pull(alt_gps), 1,),
    error=function(e) e)
  
  if(!inherits(possibleError, "error")){
  g1<-waveStatsZC(dat[dat$burstID==i,]%>%filter(!alt_gps %in% boxplot(alt_gps)$out)%>%pull(alt_gps), 1,) 
  tout$gps_hsig=g1$Hsig
  tout$gps_hmean=g1$Hmean
  tout$gps_tmean=g1$Tmean
  tout$gps_tsig=g1$Tsig
  }else{}
  
  zc_summary<-rbind(zc_summary, tout)
}
# two altimeter methods identical SO method not run

ggplot(data=zc_summary)+geom_point(aes(x=ds_hmean, y=gps_hmean))+facet_wrap(~class, scales='free')
ggplot(data=zc_summary)+geom_point(aes(x=ds_tmean, y=gps_tmean))+facet_wrap(~class, scales='free')

na.omit(zc_summary)%>%filter(class!='S')%>%select(-c('burstID', 'class'))%>%summarise_all(mean)
na.omit(zc_summary)%>%filter(class!='S')%>%select(-c('burstID', 'class'))%>%summarise_all(sd)

#ds_hsig ds_hmean ds_tmean  ds_tsig     gps_hsig gps_hmean gps_tmean gps_tsig
#8.018721 5.369986 9.360565 14.20945    9.770262  6.538162  13.35991  22.6568

#ds_hsig ds_hmean ds_tmean  ds_tsig    gps_hsig gps_hmean gps_tmean gps_tsig
#2.861487 1.870761 2.986462 5.957045    4.39308  2.583361  7.725046 20.26673

t.test(zc_summary$ds_hmean, zc_summary$gps_hmean, paired=T)
t.test(zc_summary$ds_tmean, zc_summary$gps_tmean, paired=T)

#### ^^^ ####


#### Make  Fig 1 ####

p1<-ggplot(data=dat[dat$burstID=="08611854_04_122",])+geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa), size=1)+scale_y_reverse()+
  labs(y="Pressure (mb) - reversed", x="Time")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  geom_hline(yintercept=unique(dat[dat$burstID=="08611854_04_122",]$p0), colour='red')+theme_bw()+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)')

p1.5<-ggplot(data=dat[dat$burstID=="08611854_04_122",])+
geom_point(aes(x=Longitude, y=Latitude, colour=pres_pa))+scale_color_viridis(option='plasma',trans="reverse")+
labs(x="Longitude", y="Latitude", colour="Pressure\n(mb)\n-reversed\n\n\n")+theme_bw()
plot_gg(p1.5, height=4, width=8, pointcontract = 0.5, sunangle = 40)

render_snapshot("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights/writeup/3dplot.png", clear = T)

p2<-ggplot(data=dat[dat$burstID=="08611854_04_122",])+geom_line(aes(x=DateTime_AEDT, y=alt_DS, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=alt_DS), size=1)+
geom_line(aes(x=DateTime_AEDT, y=alt_gps, group=1), col='#00A9FF')+
  geom_point(aes(x=DateTime_AEDT, y=alt_gps), size=1, col='#00A9FF')+
  labs(y="Altitude (m)")+
  geom_hline(yintercept=0, linetype='dotted')+theme_bw()+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  scale_y_continuous(limits=c(-2, 26), breaks=seq(-2,26,2), minor_breaks = NULL)+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)')

p1.6 <- fig("C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights/writeup/3dplot_zoom.png")
  
wrap_plots(p1, p1.6, p2, nrow=3)
p1+p1.6+p2 + plot_layout(nrow=3, heights = c(1,2))+ plot_annotation(tag_levels = 'a',tag_suffix = ')')+theme(plot.tag = element_text(size = 14))
 
# do manually
p1/p2

p1/p1.6/p2

#### ^^^ ####

#### Make  plots 1 and 2 of 3 panel figure ####

ggplot(data=dat%>%group_by(ID)%>%mutate(index=1:n())%>%ungroup())+
  geom_line(aes(x=index, y=alt_DS), colour='black')+
  geom_line(aes(x=index, y=alt_gps), colour='#00A9FF', alpha=0.75)+
  geom_line(aes(x=index, y=alt_SO), colour='#E68613', alpha=0.75)+
  facet_wrap(~ID, nrow=3, scales='free')+labs(y='Altitude (m)', x='Index') # ok looks good

# make fig
fig_dat<-dat%>%filter(ID==41490936)%>%mutate(index2=1:n())%>%filter(index2<5109)
p1<-ggplot(data=fig_dat)+
  geom_rect(data=fig_dat%>%filter(class=='S')%>%group_by(burstID)%>%
              summarise(xmin=min(index2), xmax=max(index2)),
            aes(xmin=xmin, xmax=xmax, ymin=-10, ymax=25),fill='darkgrey', alpha=0.5)+
  geom_rect(data=fig_dat%>%filter(class=='A')%>%group_by(burstID)%>%
              summarise(xmin=min(index2), xmax=max(index2)),
            aes(xmin=xmin, xmax=xmax, ymin=-10, ymax=25),fill='lightgrey', alpha=0.5)+
  geom_rect(xmin=4304, xmax=4614 , ymin=-10, ymax=25,colour='purple',fill=NA)+
  geom_line(aes(x=index2, y=alt_DS), colour='black', linewidth=0.1)+
  geom_line(aes(x=index2, y=alt_gps), col='#00A9FF', alpha=0.75, linewidth=0.1)+
  geom_line(aes(x=index2, y=alt_SO), colour='#E68613', alpha=0.75, linewidth=0.1)+
  geom_point(aes(x=index2, y=alt_DS), colour='black', size=0.5)+
  geom_point(aes(x=index2, y=alt_gps), col='#00A9FF', alpha=0.75, size=0.5)+
  geom_point(aes(x=index2, y=alt_SO), colour='#E68613', alpha=0.75, size=0.5)+
  geom_text(aes(x=130, y=24), label="a)", size=8)+
  scale_y_continuous(limits=c(-10, 25), breaks=c(-10,-5,0,5,10,15,20,25),  expand = c(0,0))+labs(y='Altitude (m)', x='5 minute burst index')+
  scale_x_continuous(minor_breaks=NULL,breaks = fig_dat%>%group_by(burstID)%>%summarise(min_i=min(index2))%>%arrange(min_i)%>%pull(min_i), 
                     labels=c("               Sit1", "               Sit2", "               Fly1", "               Fly2",
                              "               Fly3", "                 Land1", "               Fly4", "               Fly5",
                              "               Sit3", "               Sit4", "               Sit5", "                 Land2",
                              "                 Land3","               Sit6", "               Land4", "               Fly7",
                              "          Fly8"),expand = c(0,0))+
  theme_bw()+ theme( axis.text=element_text(size=12),axis.title=element_text(size=14)) 
# may need to tweak labels but OK for now

p2<-ggplot(data=fig_dat[fig_dat$burstID=="41490936_01_17",])+
  geom_hline(aes(yintercept=8), linetype='dotted')+
  geom_hline(aes(yintercept=0), linetype='dotted')+
  geom_point(aes(x=DateTime_AEDT, y=alt_gps), colour='grey', alpha=0.3)+
geom_line(aes(x=DateTime_AEDT, y=alt_DS, group=1))+
geom_point(aes(x=DateTime_AEDT, y=alt_DS, colour=sit_fly))+
  scale_y_continuous(name='Sitting satellite calibrated altitude (m)', breaks=seq(0, 20, 2),
  sec.axis = sec_axis(~.-8, name="Dynamic soaring (flying subset)\ncalibrated altitude (m)",
breaks=seq(0, 10, 2),labels = function(x) {ifelse(x>-1, x, "")}))+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)',  expand = c(0,0))+
theme_bw() +
  geom_text(aes(x=ymd_hms("2023-04-03 09:48:00", tz="Australia/Sydney"), y=18), label="b)", size=8)+
  theme(legend.position= c(0.8,0.8), legend.text=element_text(size=12), axis.text=element_text(size=12),axis.title=element_text(size=14),
        plot.background = element_rect(color = "purple", size = 1), legend.box.background = element_rect(colour = "black"),
        axis.title.y.right=element_text(hjust=0.1))+
  scale_colour_manual("Behaviour", values=c("red", "blue"),labels=c("flying", "sitting"))
#### ^^ ####
  
#### Measuring wave height w/ altimeters and making plot 3 of 3 panel figure  ####

#summarise first
dat%>%filter(class=="S"& wave_height!="NA")%>%
  summarise(mn_gps=mean(alt_gps), sd_gps=sd(alt_gps), mn_ds=mean(alt_DS), sd_ds=sd(alt_DS), 
                mn_wh=mean(wave_height), sd_wh=sd(wave_height)) 

ggplot(data=dat%>%filter(class=="S"& wave_height!="NA"))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa, colour=wave_height))+geom_line(aes(x=DateTime_AEDT, y=pres_pa))+facet_wrap(~burstID, scales="free")+scale_colour_viridis()

ggplot(data=dat%>%filter(class=="S"& wave_height!="NA"))+
  geom_point(aes(x=DateTime_AEDT, y=alt_DS), colour='red')+geom_line(aes(x=DateTime_AEDT, y=alt_DS), colour='red')+
  geom_point(aes(x=DateTime_AEDT, y=alt_gps), colour='green')+geom_line(aes(x=DateTime_AEDT, y=alt_gps), colour='green')+  
  facet_wrap(~burstID, scales="free")

# FYI Significant wave height = the average wave height of the top one-third highest waves

wave_temp<-dat%>%filter(class=="S")%>%group_by(burstID)%>%summarise(w_height=mean(wave_height),
                                                        w_period=mean(wave_period)) 

wave_sum<-left_join(zc_summary%>%filter(class=="S"), wave_temp, by="burstID")

# lms tell us which outliers to remove b4 pearsons corr
w1<-lm(w_height~ds_hsig, data=wave_sum)
w2<-lm(w_height~gps_hsig, data=wave_sum)
w3<-lm(w_period~ds_tmean, data=wave_sum)
w4<-lm(w_period~gps_tmean, data=wave_sum)

check_model(w1)
check_model(w2)
check_model(w3)
check_model(w4)
check_model(w5)

wp1<-ggplot(data=wave_sum[-c(13, 17,18,32),])+geom_point(aes(y=w_height, x=ds_hsig))+theme_bw()+
  labs(x='Wave height from albatross altimeters (m)',y='Wave height from satellite (m)', size=5)+
  geom_text(aes(x=7, y=1.4), label=expression(italic(r)*" = "*"0.58, "* italic(p) < 0.001), size=4)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

wp2<-ggplot(data=wave_sum[wave_sum$gps_hsig<6,])+geom_point(aes(y=w_height, x=gps_hsig))+theme_bw()+
  labs(x='Wave height from albatross GPS (m)',y='Wave height from satellite (m)', size=5)+
  geom_text(aes(x=4, y=1.4), label=expression(italic(r)*" = "*"0.37, "* italic(p)*" = "*0.03), size=4)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

wp3<-ggplot(data=wave_sum[-c(13, 17,18,32),])+geom_point(aes(y=w_period, x=ds_tmean))+theme_bw()+
  labs(x='Wave period from albatross altimeters (s)',y='Wave period from satellite (s)', size=5)+
  geom_text(aes(x=9, y=7), label=expression(italic(r)*" = "*"0.86, "* italic(p) < 0.001), size=4)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

wp4<-ggplot(data=wave_sum[wave_sum$gps_hsig<6,])+geom_point(aes(y=w_period, x=gps_tmean))+theme_bw()+
  labs(x='Wave period from albatross GPS (s)',y='Wave period from satellite (s)', size=5)+
  geom_text(aes(x=30, y=7), label=expression(italic(p)*" = NS"), size=4)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

cor.test(x=wave_sum[-c(13, 17,18,32),]$w_height, y=wave_sum[-c(13, 17,18,32),]$ds_hsig, method='pearson', na.action=na.omit) 
cor.test(x=wave_sum[wave_sum$gps_hsig<6,]$w_height, y=wave_sum[wave_sum$gps_hsig<6,]$gps_hsig, method='pearson', na.acmtion=na.omit)

cor.test(x=wave_sum[-c(13, 17,18,32),]$w_period, y=wave_sum[-c(13, 17,18,32),]$ds_tmean, method='pearson', na.action=na.omit)
cor.test(x=wave_sum[wave_sum$gps_hsig<6,]$w_period, y=wave_sum[wave_sum$gps_hsig<6,]$gps_tmean, method='pearson', na.acmtion=na.omit)

(wp1+wp2)/(wp3+wp4)


# for main fig
w1<-lm(w_height~ds_hsig, data=wave_sum[-c(13, 17,18,32),])
w3<-lm(w_period~ds_tmean, data=wave_sum[-c(13, 17,18,32),])

new_d<-data.frame(ds_hsig=seq(0, 10, 0.2))
new_d<-cbind(new_d, predict(w1, new_d, se.fit = T)[c('fit', 'se.fit')])
new_d$lci<-new_d$fit-(new_d$se.fit*1.96)
new_d$uci<-new_d$fit+(new_d$se.fit*1.96)

p3<-ggplot()+
  geom_point(data=wave_sum, aes(y=w_height, x=ds_hsig))+
  geom_line(data=new_d, aes(x=ds_hsig, y=fit),colour='red', linewidth=1)+
  geom_ribbon(data=new_d, aes(x=ds_hsig, ymin=lci, ymax=uci), alpha=0.5, fill='grey')+
  theme_bw()+
  scale_x_continuous(limits=c(0, 9), breaks=0:9, expand = c(0,0))+scale_y_continuous(limits=c(0.5, 4), expand = c(0,0))+
  labs(x='Wave height from albatross altimeters (m)',y='Wave height from satellite (m)', size=5)+
  geom_text(aes(x=2, y=3.5), label=expression(italic(r)*" = "*"0.58, "* italic(p) < 0.001), size=4)+
  geom_text(aes(x=8, y=1), label="c)", size=8)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

new_d<-data.frame(ds_tmean=seq(4, 13, 0.2))
new_d<-cbind(new_d, predict(w3, new_d, se.fit = T)[c('fit', 'se.fit')])
new_d$lci<-new_d$fit-(new_d$se.fit*1.96)
new_d$uci<-new_d$fit+(new_d$se.fit*1.96)

p4<-ggplot()+
  geom_point(data=wave_sum[-c(13, 17,18,32),], aes(y=w_period, x=ds_tmean))+
  geom_line(data=new_d, aes(x=ds_tmean, y=fit),colour='red', linewidth=1)+
  geom_ribbon(data=new_d, aes(x=ds_tmean, ymin=lci, ymax=uci), alpha=0.5, fill='grey')+
  theme_bw()+
  scale_x_continuous(limits=c(4, 13), breaks=4:13, expand = c(0,0))+scale_y_continuous(limits=c(6, 13), breaks=6:13, expand = c(0,0))+
  labs(x='Wave period from albatross altimeters (s)',y='Wave period from satellite (s)', size=5)+
  geom_text(aes(x=6, y=12), label=expression(italic(r)*" = "*"0.86, "* italic(p) < 0.001), size=4)+
  geom_text(aes(x=12, y=7), label="d)", size=8)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

# Make mega plot!!

areas <- c(patchwork::area(1, 1, 1, 3),patchwork::area(2, 1, 2, 2), patchwork::area(2, 3, 2,3))
p1 + p2 + (p3/p4) + plot_layout(design = areas)

#### ^^ ####

#### Summarise altitude from the three methods and compare  ####

#remove first GPS fix of each burst as higher error
dat<-dat %>% group_by(burstID) %>%
  filter(row_number()!=1)%>%ungroup()%>%as.data.frame()

# compare differences between three methods

# format dataset for comparison # not added 1000 to all alts to make positive for Gamma

dat_flying<-dat%>%filter(class %in% c('T', 'L') & sit_fly=='fly')

cor.test(dat_flying$alt_DS, dat_flying$alt_gps)
ggplot(data=dat_flying)+geom_point(aes(x=alt_gps, y=alt_DS))

dat_comp<-rbind(data.frame(method='Dynamic soaring', Altitude=dat_flying$alt_DS, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID) ,
                data.frame(method='Satellite ocean', Altitude=dat_flying$alt_SO, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID),
                data.frame(method='GPS', Altitude=dat_flying$alt_gps, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID))

#summarise
dat_comp%>%group_by(method)%>%summarise(mn_alt=mean(Altitude), sd_alt=sd(Altitude), median=median(Altitude),
                                    min=min(Altitude), max=max(Altitude),
                                    q25=quantile(Altitude, 0.25), q75=quantile(Altitude, 0.75), skew=skewness(Altitude))
# make plot
cols <- c("#000000",'#00A9FF','#E68613')

cols.alpha<-c(grDevices::adjustcolor(cols[1], alpha.f = 0.75),
        grDevices::adjustcolor(cols[2], alpha.f = 0.75),
        grDevices::adjustcolor(cols[3], alpha.f = 0.75))

ggplot(data=dat_comp)+geom_density(aes(x=Altitude, colour=method), fill=NA, size=2)+
  theme_bw()+geom_vline(xintercept = 0, linetype='dotted')+scale_x_continuous(breaks=seq(-60,60,2))+
  scale_colour_manual(values = cols.alpha)+coord_cartesian(xlim=c(-20, 40))+
  theme(legend.position= c(0.8,0.8), axis.text=element_text(size=10),axis.title=element_text(size=12),
        legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  scale_colour_manual("Altitude estimation method", values=cols.alpha, labels=c("Dynamic soaring calibrated altimeter", 
  "GPS Altitude", "Sitting satellite calibrated altimeter"))+labs(x="Altitude (m)", y="Density")

# Now run stats on difference data to keep things normal

dat_diff<-rbind(data.frame(method='Dynamic soaring - GPS', Altitude=dat_flying$alt_DS-dat_flying$alt_gps, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID) ,
                data.frame(method='Dynamic soaring - Satellite ocean', Altitude=dat_flying$alt_DS-dat_flying$alt_SO, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID),
                data.frame(method='GPS - Satellite ocean', Altitude=dat_flying$alt_gps-dat_flying$alt_SO, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID))


m1<-lme(Altitude~method, random=~1|burstID, weights=varIdent(form=~1|method), data=dat_diff)
boxplot(residuals(m1, type='pearson')~dat_diff$method)
m2<-lme(Altitude~method, random=~1|Logger, weights=varIdent(form=~1|method), data=dat_diff)

AIC(m1, m2) # logger ID as RE not as good
resid_panel(m1)
summary(m1)

anova(m1)
em1<-emmeans(m1, specs='method')
em1
test(em1, adjust="bonferroni")

#### Make prop time in 1m band fig and table ####

min(dat_flying$alt_DS) # get offset to make all vals above zero for lognormal and Gamma
#4.166551

lognormal.SH = fitdistr(dat_flying$alt_DS+4.166551, "lognormal")
lognormal.SH.Fit = dlnorm(seq(0,30,0.1), lognormal.SH$estimate[1], lognormal.SH$estimate[2])
plot(seq(0,30,0.1), lognormal.SH.Fit, type = "l", ylab = "Prop. at height", xlab = "Height above Sea-level")

# dont do this, do own bootstrap with sample and refit and predict the dist. 200 times (apparently) can alos use these to reoport 
#ucl and lcl in the table on the paper


set.seed(123)
lnorm.boot = data.frame(matrix(data = 0, nrow = 3001, ncol = 200))
for(i in 1:200){
  balanced_ID_alts<-c(sample(dat_flying[dat_flying$ID==8611649,]$alt_DS, round(nrow(dat_flying)/3),replace=TRUE),
                 sample(dat_flying[dat_flying$ID==8611854,]$alt_DS, round(nrow(dat_flying)/3),replace=TRUE),
                 sample(dat_flying[dat_flying$ID==41490936,]$alt_DS, round(nrow(dat_flying)/3),replace=TRUE))
  lognormal.SH.boot = fitdistr(balanced_ID_alts+4.166551, "lognormal")
  lnorm.boot[,i] = dlnorm(seq(0,300,0.1), lognormal.SH.boot$estimate[1], lognormal.SH.boot$estimate[2])
}

lnorm.boot[,201]<-rowMeans(lnorm.boot)

dimnames(lnorm.boot)[[1]]<-(seq(0,300,0.1)-4.166551)+1.15
dimnames(lnorm.boot)[[2]]<-c(paste0('bootId_', 1:200), "mean")

# sd col now then get 95 CI or 97.5 actually

boot_for_plot<-data.frame(ht=as.numeric(row.names(lnorm.boot[as.numeric(row.names(lnorm.boot))<max(dat_flying$alt_DS+1.15),])),
                             prop = lnorm.boot[as.numeric(row.names(lnorm.boot))<max(dat_flying$alt_DS+1.15),"mean"])

p1<-ggplot()+
  geom_vline(xintercept=0, colour='blue', size=1)+
  geom_vline(xintercept=10, colour='blue', size=0.5)+
  geom_vline(xintercept=20, colour='blue', size=0.5)+
  geom_rect(aes(xmin=30, xmax=31, ymin=0, ymax=0.17), fill='red', size=0.5, alpha=0.3)+
  geom_histogram(data=dat_flying, aes(x=alt_DS+1.15, after_stat(density)), fill='grey', colour='darkgrey', binwidth = 1)+
  geom_line(data=boot_for_plot, aes(x=ht, y=prop), size=1)+
  scale_x_continuous(breaks=seq(-5, 50, 5),minor_breaks=seq(-5, 50, 1), limits=c(-5, 31), expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(x="Height above sea level (m)", y="Proportion at height")+theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))


# make proportion table, per 1m
propt<-lnorm.boot%>%mutate(ht_1m=cut(as.numeric(row.names(lnorm.boot)), breaks = c(-5, 1:297)))%>%group_by(ht_1m)%>%summarise_all(sum)

propt[,2:202]<-propt[,2:202]/10

propt[,203]<-apply(propt[,2:201], 1,  function(x){confint(lm(x~1), level=0.95)[1]})
propt[,204]<-apply(propt[,2:201], 1,  function(x){confint(lm(x~1), level=0.95)[2]})
colnames(propt)[203:204]<-c("lci", "uci")

propt$ht_1m<-gsub( ",", "-",propt$ht_1m)

#write out
#write.csv(propt, "C:/Users/mmil0049/OneDrive - Monash University/projects/02 flight heights/outputs/flight_height_boots.csv", quote=F, row.names=F)

# make table for plot up to 25m

table_dat<-propt[1:25, c("ht_1m","mean", "lci", "uci")]


specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

table_dat$mean<-specify_decimal(table_dat$mean, 5)
table_dat$lci<-specify_decimal(table_dat$lci, 5)
table_dat$uci<-specify_decimal(table_dat$uci, 5)
table_dat$Altitude=c("<1m", "1-2m",   "2-3m",   "3-4m",   "4-5m",   "5-6m",   "6-7m",   "7-8m",   "8-9m",   "9-10m",
"10-11m", "11-12m", "12-13m", "13-14m", "14-15m", "15-16m", "16-17m", "17-18m", "18-19m",
"19-20m", "20-21m", "21-22m", "22-23m", "23-24m", "24-25m")
names(table_dat)[2:4]<-c("Proportion", "LCI", "UCL")
sum(as.numeric(table_dat$Proportion))

p1 + gridExtra::tableGrob(table_dat[c('Altitude', 'Proportion', "LCI", "UCL")], row=NULL)

#### ^^ ####

ggplot(data=dat)+geom_point(aes(x=DateTime_AEDT, y=alt, colour=class))+
  facet_wrap(~ID, scales="free", nrow=3)
# do overall distribution for flying birds

mean(dat%>%filter(class %in% c("T", "L") & sit_fly=="fly")%>%pull(pres_alt), na.rm=T)

ggplot()+
  geom_vline(xintercept=0, colour='blue', size=1)+
  geom_vline(xintercept=10, colour='blue', size=0.5)+
  geom_vline(xintercept=20, colour='blue', size=0.5)+
  geom_vline(xintercept=4.209, colour='red',linetype='dotted', size=1)+
  geom_label(aes(x=6.9, y=4200, label="Mean\n4.21 m"),size=5, col='red')+
  geom_rect(aes(xmin=30, xmax=50, ymin=0, ymax=4300), fill='red', size=0.5, alpha=0.3)+
  geom_histogram(data=dat%>%filter(class %in% c("T", "L") & sit_fly=="fly"), aes(x=alt_DS),colour=1, binwidth = 1)+
  scale_x_continuous(breaks=seq(-5, 50, 5),minor_breaks=seq(-5, 50, 1), limits=c(-5, 50))+
  labs(x="Altitude (m)", y="Number 3D datapoints (Lat,Lon,Pressure)")+theme_bw()+
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold"))

tl_dat<-dat%>%filter(class %in% c("T", "L") & sit_fly=="fly" & burstID!="-1_93")

mean(tl_dat[tl_dat$pres_alt>0,]$pres_alt) # exclude negative vals
  
nrow(tl_dat[tl_dat$pres_alt<20,])/nrow(tl_dat)  
nrow(tl_dat[tl_dat$pres_alt<30,])/nrow(tl_dat)
nrow(tl_dat[tl_dat$pres_alt<10,])/nrow(tl_dat)
nrow(tl_dat[tl_dat$pres_alt<5,])/nrow(tl_dat)

# write out kmz of bird that crosses the island - removed prior to analyes now
#sf3d<-dat%>%filter(burstID=="-1_93")%>%st_as_sf(coords = c("Longitude", "Latitude", "pres_alt"), crs = 4326, dim = "XYZ")
#st_write(sf3d, "analyses/GIS/over_the_island_burst.kml")

# test GPS difference

sensor_comp<-rbind(data.frame(Sensor="gps", Altitude=tl_dat[tl_dat$alt<2000,]$alt),
      data.frame(Sensor="pressure", Altitude=tl_dat[tl_dat$alt<2000,]$pres_alt))

sensor_comp%>%group_by(Sensor)%>%summarise(mean_alt=mean(Altitude))

ggplot(data=sensor_comp%>%filter(Altitude<90))+
  geom_hline(yintercept=0, colour='blue', size=1)+
  geom_jitter(aes(x=Sensor, y=Altitude), height=0, width=0.4, alpha=0.1, shape=16, size=1.5)+
  geom_violin(aes(x=Sensor, y=Altitude, colour=Sensor), fill=NA, size=1)+
  scale_y_continuous(breaks=seq(-50, 90, 10),minor_breaks=seq(-50, 90, 5))+
  labs(y="Altitude (m)")+theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),legend.position="none")
  
# test windspeed

ggplot(data=tl_dat, aes(x=wind_speed, y=pres_alt))+
  geom_point(alpha=0.1)+geom_smooth(method='lm')+facet_wrap(~b_w_class)

av_dat<-tl_dat%>%group_by(burstID)%>%summarise(mn_alt=mean(pres_alt), mn_w_speed=mean(wind_speed, na.rm=T))

ggplot()+
  geom_hline(yintercept=0, colour='blue', size=1)+
  geom_point(data=tl_dat, aes(x=wind_speed, y=pres_alt), alpha=0.05)+
  geom_smooth(data=av_dat, aes(x=mn_w_speed, y=mn_alt), method='lm', colour='red', fill='pink')+
  scale_y_continuous(breaks=seq(-4,24,2), limits=c(-4, 24))+
  scale_x_continuous(breaks=seq(0,15,1))+
  labs(y="Altitude (m)", x="Windspeed (m/s)")+theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),legend.position="none")

 mjyuy# compare pressure and GPS altitude

ggplot(data=dat%>%filter(class %in% c("T", "L")& sit_fly=="fly"))+geom_density(aes(y=alt, x=pres_alt))



### example plots for talk

# pressure + dyn soaring
ggplot(data=dat[dat$burstID=="08611854_02_71",])+geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa), size=1)+scale_y_continuous(breaks=c(100800, 100850, 100900, 100950, 101000, 101050))+
  labs(y="Pressure (mb)", x="Time")+theme(axis.text=element_text(size=12),
                                                                                                                                                                           axis.title=element_text(size=14,face="bold"))

# 3d plot with rayshader - example track ALSO check 08611854_06_146 for loopy example and "08611854_04_120" for Alight
p1<-ggplot(data=dat[dat$burstID=="41490936_01_24",])+
  geom_point(aes(x=Longitude, y=Latitude, colour=pres_alt))+scale_color_viridis()+
  labs(x="Longitude", y="Latitude", colour="   Altitude (m)")

p4<-render_snapshot(plot_gg(p1, height=3, width=8, pointcontract = 0.7))

#comparison with GPS
p1<-ggplot(data=dat[dat$burstID=="41490936_01_24",])+
  geom_point(aes(x=Longitude, y=Latitude, colour=alt))+scale_color_viridis()+
  labs(x="Longitude", y="Latitude", colour="   Altitude (m)")

plot_gg(p1, height=3, width=8, pointcontract = 0.7)

# showing altitude against GPS speed - not manually calculated speed!

p2<-ggplot(data=dat[dat$burstID=="41490936_01_24",])+
  geom_point(aes(x=Longitude, y=Latitude, colour=speed))+scale_color_viridis()
p1<-ggplot(data=dat[dat$burstID=="41490936_01_24",])+
  geom_point(aes(x=Longitude, y=Latitude, colour=pres_alt))+scale_color_viridis()

p1+p2

p2<-ggplot(data=dat[dat$burstID=="08611854_06_146",])+
  geom_point(aes(x=Longitude, y=Latitude, colour=speed))+scale_color_viridis()
p1<-ggplot(data=dat[dat$burstID=="08611854_06_146",])+
  geom_point(aes(x=Longitude, y=Latitude, colour=pres_alt))+scale_color_viridis()

p1+p2
# 3d plot to google earth

temp1<-dat
temp1$DESCRIPTION=temp1$class
temp1$TIME<-ymd_hms(temp1$DateTime_AEDT)

sf4d<-temp1[c("DESCRIPTION","Longitude", "Latitude", "pres_alt", "TIME")]%>%
  st_as_sf(coords = c("Longitude", "Latitude", "pres_alt", "TIME"), crs = 4326, dim = "XYZM")

sf3d<-temp1[c("DESCRIPTION","Longitude", "Latitude", "pres_alt", "TIME")]%>%
  st_as_sf(coords = c("Longitude", "Latitude", "pres_alt"), crs = 4326, dim = "XYZ")

#st_write(sf3d, "analyses/GIS/google_earth_3d_vis.kml", driver='kml')

#st_write(sf3d%>%filter(DESCRIPTION=="T"), "analyses/GIS/google_earth_3d_vis_transit.kml", driver='kml')
#st_write(sf3d%>%filter(DESCRIPTION=="L"), "analyses/GIS/google_earth_3d_vis.kml_loop.kml", driver='kml')
#st_write(sf3d%>%filter(DESCRIPTION=="A"), "analyses/GIS/google_earth_3d_vis.kml_alight.kml", driver='kml')

# find good bursts

T_b<-dat%>%filter(class=='T')%>%pull(burstID)%>%unique()

for(j in T_b){
p2<-ggplot(data=dat[dat$burstID==j,])+geom_line(aes(x=tim_UTC, y=pres_pa, group=1))+
  geom_point(aes(x=tim_UTC, y=pres_pa), size=1)+scale_y_reverse()
p3<-ggplot(data=dat[dat$burstID==j,])+geom_line(aes(x=X, y=Y, group=1))+
  geom_point(aes(x=X, y=Y, colour=pres_alt))
p2/p3             
print(p2/p3)
print(j)
readline("")
}

#"08611854_01_64"
"08611854_02_68"
"08611854_02_69"
"08611854_02_70"
"08611854_02_71"

p2<-ggplot(data=dat[dat$burstID==j,])+geom_line(aes(x=tim_UTC, y=pres_pa, group=1))+
  geom_point(aes(x=tim_UTC, y=pres_pa), size=1)+geom_hline(yintercept=upper_95, col='red')+scale_y_reverse()

# test export as kmz

temp1<-dat
temp1$tim_UTC<-ymd_hms(temp1$tim_UTC)

#sf4d<-temp1%>%st_as_sf(coords = c("X", "Y", "pres_alt", "tim_UTC"), crs = 4326, dim = "XYZM")

sf3d<-temp1%>%st_as_sf(coords = c("X", "Y", "pres_alt"), crs = 4326, dim = "XYZ")

st_write(sf3d, "C:/Users/mmil0049/Downloads/temp.kml")

#### Investigating 'A' class alighting/landing bursts - OLD CODE ####

p1<-ggplot(data=dat[dat$class=="A",])+
  geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa, colour=sit_fly), size=1)+geom_line(aes(x=DateTime_AEDT, y=p0), col='red')+
  scale_y_reverse()+facet_wrap(~burstID, scales="free")

p2<-ggplot(data=dat[dat$class=="A",])+
  geom_line(aes(x=DateTime_AEDT, y=alt, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=alt), size=1)+facet_wrap(~burstID, scales="free")

p1/p2

ggplot(data=dat[dat$class=="A",])+
  geom_line(aes(x=DateTime_AEDT, y=temp, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=temp), size=1)+
  facet_wrap(~burstID, scales="free")

ggplot(data=dat[dat$class=="A",])+
  geom_line(aes(x=DateTime_AEDT, y=alt_DS, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=alt_DS, colour=wave_height), size=1)+
  facet_wrap(~burstID, scales="free")+scale_colour_viridis() # high waves can look like DS - remember bird speed!

#### ^^ ####


## have a look at cleasby gannet data

gps1<-read.csv("sourced_data/Cleasby (2015) data doi_10.5061_dryad.1ds1q__v1/GPS/1446280_260611.csv", h=T)

gps1$TIME<-dmy_hms(gps1$TIME)
gps1$day<-day(gps1$TIME)
gps1$hour<-hour(gps1$TIME)
gps1$colD_class<-round((gps1$COL_DIST/10),0)

ggplot(data=gps1)+geom_point(aes(x=TIME, y=COL_DIST, colour=STATUS1))

ggplot(data=gps1%>%filter(day==26&STATUS1!='DIVE'&P<1025))+geom_line(aes(x=TIME, y=P, group=1))+
  geom_point(aes(x=TIME, y=P, colour=STATUS1), size=1)

ggplot(data=gps1%>%filter(day==26&STATUS1!='DIVE'&P<1025))+geom_line(aes(x=TIME, y=P, group=1))+
  geom_point(aes(x=TIME, y=P, colour=STATUS1), size=1)+facet_wrap(~hour+colD_class, scales="free")

# seems like pressure already in 1 sec GPS data
gps1$key<-paste(gps1$day, gps1$hour,gps1$colD_class)

for(i in unique(gps1$key))
{
  print(ggplot(data=gps1%>%filter(STATUS1!='DIVE'&P<1025&key==i))+geom_line(aes(x=TIME, y=P, group=1))+
    geom_point(aes(x=TIME, y=P, colour=STATUS1), size=1.5)+scale_y_reverse()+labs(main=i))
  readline("")
}

# ok looks interesting, there is some 'dynamic soaring' at low latitude but clear flapping to gain altitude. Few obvious waves while sitting
