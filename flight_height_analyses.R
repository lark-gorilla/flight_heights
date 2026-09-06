# Code to conduct paper analyses.

library(ggplot2)
library(dplyr)
library(tidyr)
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
library(car)
library(zoo)
library(pracma)
library(mgcv)

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

###~ GPS correction and tweak ~###

#Add -9m correction to GPS elevation 
dat$alt_gps<-dat$alt-9

#remove first GPS fix of each burst as higher error
dat<-dat %>% group_by(burstID) %>%
  filter(row_number()!=1)%>%ungroup()%>%as.data.frame()
### ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

#make flying subset for flight height analyses
dat_flying<-dat%>%filter(class %in% c("T", "L"))%>%ungroup()

#### ^^ ####

#### dynamic soaring cycle segmentation ####
dat_flying$ds_seg_pressure<-NA
dat_flying$ds_seg_gps<-NA
for(i in unique(dat_flying$burstID))
{
i_burst<-dat_flying[dat_flying$burstID==i,]

# don't detrend before
#loess_fit <- loess(alt_gps ~ as.numeric(DateTime_AEDT), data = i_burst, span = 0.75)
#detrended <- i_burst$alt_gps - predict(loess_fit)

pres_smth <- rollmean(i_burst$pres_pa, k = 3, fill = NA) # apply moving window 3-pt
gps_smth <- rollmean(i_burst$alt_gps, k = 3, fill = NA)

# find valleys - could possibly be optimised including zero-crossing too
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

#### Sensitivity tests for p0 method ####

# use absolute maximum
dat_flying<-dat_flying%>%group_by(burstID)%>%mutate(p0_mx=max(pres_pa)) 

# per dynamic soaring segment (above loop)
temp_d<-dat_flying %>% group_by(burstID) %>% slice(c(1,1:(n()-1)))%>%ungroup()%>%arrange(ID, burstID,ds_seg_pressure)
dat_flying<-dat_flying%>%ungroup()%>%arrange(ID, burstID,ds_seg_pressure) %>%mutate(pres_pa_1=temp_d$pres_pa)
dat_flying<-dat_flying%>%group_by(burstID, ds_seg_pressure)%>%mutate(p0_ds_segSRT=max(pres_pa_1), p0_ds_segEND=max(pres_pa),
                                                                     p0_ds_seg=max(p0_ds_segSRT, p0_ds_segEND) ) 
#calculate most likely p0
dat_flying$p0_diff<-dat_flying$p0_mx-dat_flying$p0_ds_seg
dat_flying<-dat_flying%>%group_by(burstID)%>%mutate(p0_gam=
                               fitted(gam((p0_mx-(p0_diff/2))~s(DateTime_AEDT%>%as.numeric(), k=7))))%>%ungroup
#make sure no predictions higher than p0_max
dat_flying$p0_gam<-ifelse(dat_flying$p0_gam>dat_flying$p0_mx, dat_flying$p0_mx, dat_flying$p0_gam)

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
dat_flying$alt_p0_mx<-(-1*((k*(dat_flying$temp+273.15))/(m*g))*log(dat_flying$pres_pa/dat_flying$p0_mx))
dat_flying$alt_p0_ds_seg<-(-1*((k*(dat_flying$temp+273.15))/(m*g))*log(dat_flying$pres_pa/dat_flying$p0_ds_seg))
dat_flying$alt_p0_gam<-(-1*((k*(dat_flying$temp+273.15))/(m*g))*log(dat_flying$pres_pa/dat_flying$p0_gam))
dat_flying$alt_p0_error<-(-1*((k*(dat_flying$temp+273.15))/(m*g))*log(dat_flying$p0_ds_seg/dat_flying$p0_mx))

# view error distribution

ggplot(data=dat_flying)+geom_histogram(aes(x=alt_p0_error), binwidth=0.5)+
  geom_vline(data=dat_flying%>%summarise(med_er=median(alt_p0_error)),aes(xintercept=med_er), col=3)+
  scale_x_continuous(breaks=0:11)+theme_bw()+xlab("Error in altimeter flight height estimates (m)")+ylab("Count of datapoints")

ggplot(data=dat_flying)+geom_histogram(aes(x=alt_p0_error), binwidth=0.5)+
  geom_vline(data=dat_flying%>%group_by(ID)%>%summarise(med_er=median(alt_p0_error)),aes(xintercept=med_er), col=3)+
  scale_x_continuous(breaks=0:11)+ theme_bw()+
  facet_wrap(~ID, scales='free_y')+xlab("Error in altimeter flight height estimates (m)")+ylab("Count of datapoints")

summary(dat_flying$alt_p0_error)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   1.085   2.087   2.376   3.270  10.842
dat_flying%>%group_by(ID)%>%summarise(med_er=median(alt_p0_error))

#### ^^ ####

#### Summarise altitude from the three DS methods and GPS, and compare  ####

dat_comp<-rbind(data.frame(method='Dynamic soaring - lwr bound', Altitude=dat_flying$alt_p0_ds_seg, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID),
                data.frame(method='Dynamic soaring - upr bound', Altitude=dat_flying$alt_p0_mx, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID),
                data.frame(method='Dynamic soaring - mean', Altitude=dat_flying$alt_p0_gam, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID),
                data.frame(method='GPS', Altitude=dat_flying$alt_gps, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID))

#summarise
dat_comp%>%group_by(method)%>%summarise(mn_alt=mean(Altitude), sd_alt=sd(Altitude), median=median(Altitude),
                                        min=min(Altitude), max=max(Altitude),
                                        q25=quantile(Altitude, 0.25), q75=quantile(Altitude, 0.75), skew=skewness(Altitude))

#method                      mn_alt sd_alt median    min   max   q25   q75   skew
#<chr>                        <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>  <dbl>
#1 Dynamic soaring - lwr bound   3.23   2.70   2.67   0     26.1  1.26  4.55 1.46  
#2 Dynamic soaring - mean        4.42   2.89   3.86  -1.61  25.9  2.36  5.91 1.20  
#3 Dynamic soaring - upr bound   5.61   3.24   4.97   0     26.1  3.28  7.32 1.09  
#4 GPS                           5.10   9.72   4    -63     89    0     9    0.0607

# skewness stats
test.skew(dat_comp%>%filter(method=='GPS')%>%pull(Altitude))


# make plot
cols <- c('#dc267f','#648fff','#fe6100','#ffb000')

cols.alpha<-c(grDevices::adjustcolor(cols[1], alpha.f = 1),
              grDevices::adjustcolor(cols[2], alpha.f = 1),
              grDevices::adjustcolor(cols[3], alpha.f = 1),
              grDevices::adjustcolor(cols[4], alpha.f = 1))

ggplot(data=dat_comp)+geom_density(aes(x=Altitude, colour=method), fill=NA, size=1.5)+
  theme_bw()+geom_vline(xintercept = 0, linetype='dotted')+geom_hline(yintercept = 0,size=1.5)+
  scale_x_continuous(breaks=seq(-60,60,2))+
  scale_colour_manual(values = cols.alpha)+coord_cartesian(xlim=c(-20, 40))+
  theme(legend.position= c(0.8,0.8), axis.text=element_text(size=10),axis.title=element_text(size=12),
        legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  scale_colour_manual("Flight height estimation method", values=cols.alpha, labels=c("Altimeters (lower scenario)",
                                                                                     "Altimeters (most likely scenario)",
                                                                                     "Altimeters (upper scenario)",
                                                                                     "GPS Altitude"))+labs(x="Flight height (m)", y="Density")


wilcox.test(x= dat_flying$alt_gps,
            y=dat_flying$alt_p0_gam, paired=T)

#Wilcoxon signed rank test with continuity correction

#data:  dat_flying$alt_gps and dat_flying$alt_p0_gam
#V = 169251078, p-value = 1.703e-08

wilcox.test(x= dat_flying$alt_gps,
            y=dat_flying$alt_p0_mx, paired=T)

#data:  dat_flying$alt_gps and dat_flying$alt_p0_mx
#V = 138154406, p-value < 2.2e-16

wilcox.test(x= dat_flying$alt_gps,
            y=dat_flying$alt_p0_ds_seg, paired=T)

#data:  dat_flying$alt_gps and dat_flying$alt_p0_ds_seg
#V = 198704999, p-value < 2.2e-16

#### ^^ ####


#### Check suitability floating-satellite method (Johnston et al 2023) ####

dat<-dat%>%group_by(burstID)%>%mutate(mean_sea_level_pressure =mean(mean_sea_level_pressure ,na.rm = T))%>%
  ungroup()%>%as.data.frame() # average satellite pressure dat per burst

#get sitting pressure diff 
dat$sat_sit_pdiff<-NA
dat[dat$class %in% c('A', 'S') & dat$sit_fly=='sit',]$sat_sit_pdiff<-
 (dat[dat$class %in% c('A', 'S') & dat$sit_fly=='sit',]$pres_pa-
      dat[dat$class %in% c('A', 'S') & dat$sit_fly=='sit',]$mean_sea_level_pressure) 

# summarise per burst
dat<-dat%>%group_by(burstID)%>%
  mutate(burstID_sat_sit_pdiff=mean(sat_sit_pdiff,na.rm = T))%>%ungroup()%>%as.data.frame()

sumr<-NULL
for(i in unique(dat$burstID))
{
  dtemp<-dat%>%filter(burstID==i)
  IDtemp<-dat%>%filter(ID==unique(dtemp$ID)) # get nearest from same logger
  sit_burst<-IDtemp[! is.na(IDtemp$burstID_sat_sit_pdiff),] # only bursts with sitting diffs
  
  out1<-data.frame(burstID=i, ID=unique(dtemp$ID), class=unique(dtemp$class), 
                   min_t=min(abs((sit_burst$DateTime_AEDT- 
              median(dtemp$DateTime_AEDT)))))
  sumr<-rbind(sumr, out1)
}

sumr$min_t<-as.numeric(sumr$min_t)
ggplot(data=sumr%>%filter(class%in%c('T', 'L')))+geom_histogram(aes(x=min_t/3600))+facet_wrap(~ID)
summary(sumr[sumr$class%in%c('T', 'L'),]$min_t/3600)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9339  0.9644  1.9719  6.7607 14.8675 23.9706 - not suitable to large time gap
sd(sumr[sumr$class%in%c('T', 'L'),]$min_t/3600) 
# 7.652972
#### ^^^ ####


#### Zero-crossing analysis  ####

# Using Zero-crossing method in oceanwaves package

# Use floating data only
dat_floating<-dat%>%filter(class=="S")

# Use max p0 to calc wave height from pressure data - zero crossing doesn't care which p0 is used, 
# using barometric equation just converts units from pascals to metres to allow comparison with GPS
dat_floating<-dat_floating%>%group_by(burstID)%>%mutate(p0_mx=max(pres_pa)) 
k=8.31432
m=0.0289644
g=9.80665
dat_floating$alt_p0_mx<-(-1*((k*(dat_floating$temp+273.15))/(m*g))*log(dat_floating$pres_pa/dat_floating$p0_mx))

zc_summary<-NULL
for ( i in unique(dat_floating$burstID))
{
 tout<-data.frame(class=unique(dat_floating[dat_floating$burstID==i,]$class), burstID=i, 
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
geom_point(aes(x=Longitude, y=Latitude, colour=pres_pa))+scale_color_viridis(trans="reverse")+
labs(x="Longitude", y="Latitude", colour="Pressure\n(Pa)\n\n")+theme_bw()
plot_gg(p1.5, height=4, width=8, pointcontract = 0.5, sunangle = 40)

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
            aes(xmin=xmin, xmax=xmax, ymin=-10, ymax=25),fill='khaki2', alpha=0.5)+
  geom_line(aes(x=index2, y=alt_DS), colour='black', linewidth=0.1)+
  geom_line(aes(x=index2, y=alt_gps), col='#00A9FF', alpha=0.75, linewidth=0.1)+
  geom_line(aes(x=index2, y=alt_SO), colour='#E68613', alpha=0.75, linewidth=0.1)+
  geom_text(aes(x=130, y=21), label="a)", size=8)+
  scale_y_continuous(limits=c(-10, 25), breaks=c(-10,-5,0,5,10,15,20,25),  expand = c(0,0))+labs(y='Altitude (m)', x='5 minute burst index')+
  scale_x_continuous(minor_breaks=NULL,breaks = fig_dat%>%group_by(burstID)%>%summarise(min_i=min(index2))%>%arrange(min_i)%>%pull(min_i), 
                     labels=c("               Sit1", "               Sit2", "               Fly1", "               Fly2",
                              "               Fly3", "                 Land1", "               Fly4", "               Fly5",
                              "               Sit3", "               Sit4", "               Sit5", "                 Land2",
                              "                 Land3","               Sit6", "               Land4", "               Fly7",
                              "          Fly8"),expand = c(0,0))+
  theme_bw()+ theme( axis.text=element_text(size=10),axis.title=element_text(size=12))+
  geom_rect(xmin=4310, xmax=4620 , ymin=-9.8, ymax=24.8,colour='purple',fill=NA, linewidth=1)
# may need to tweak labels but OK for now

p2<-ggplot(data=fig_dat[fig_dat$burstID=="41490936_01_17",])+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa), colour='grey', size=0.5, alpha=0.3)+
  geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1), colour="darkgrey")+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa, colour=sit_fly), size=1)+
  scale_y_reverse(name='Pressure (Pa) - reversed')+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)',  expand = c(0,0))+
  theme_bw() +
  geom_text(aes(x=ymd_hms("2023-04-03 09:48:00", tz="Australia/Sydney"), y=102150), label="b)", size=8)+
  theme(legend.position= c(0.8,0.7), legend.text=element_text(size=10), axis.text=element_text(size=10),axis.title=element_text(size=12),
        plot.background = element_rect(color = "purple", size = 1), legend.box.background = element_blank(),
        legend.background = element_blank())+
  scale_colour_manual("Behaviour", values=c("red", "blue"),labels=c("flying", "sitting"))

colors <- c("alt_DS" = "black", "alt_gps" = "#00A9FF", "alt_SO" = "#E68613")

p3<-ggplot(data=fig_dat[fig_dat$burstID=="41490936_01_17",])+
    geom_line(aes(x=DateTime_AEDT, y=alt_DS, color="alt_DS"))+
  geom_line(aes(x=DateTime_AEDT, y=alt_gps, color="alt_gps"))+
  geom_line(aes(x=DateTime_AEDT, y=alt_SO, color="alt_SO"))+
  scale_y_continuous(limits=c(-3, 18), breaks=c(seq(-3, 18, 3)),  expand = c(0,0))+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)',  expand = c(0,0))+
  theme_bw() +
  labs(y="Flight height (m)")+
  geom_text(aes(x=ymd_hms("2023-04-03 09:48:00", tz="Australia/Sydney"), y=15), label="c)", size=8)+
  theme(legend.position= c(0.8,0.7), legend.text=element_text(size=10), axis.text=element_text(size=10),axis.title=element_text(size=12),
        plot.background = element_rect(color = "purple", size = 1), legend.box.background = element_blank(),
        legend.background = element_blank())+
  scale_colour_manual("Flight height estimator", values=colors,labels=c("alt_DS" ="Altimeters zeroed with dynamic soaring", 
                                                                                 "alt_gps" = "GPS Altitude",
                                                                                 "alt_SO"="Altimeters zeroed from sitting at sea or satellite"))

#make fig 4
p1/p2/p3  

#### ^^ ####
  
#### Measuring wave height w/ altimeters and making figure 5  ####

#summarise first
dat_floating%>%filter(wave_height!="NA")%>%
  summarise(mn_gps=mean(alt_gps), sd_gps=sd(alt_gps), mn_ds=mean(alt_p0_mx), sd_ds=sd(alt_p0_mx), 
                mn_wh=mean(wave_height), sd_wh=sd(wave_height)) 

ggplot(data=dat_floating%>%filter(wave_height!="NA"))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa, colour=wave_height))+geom_line(aes(x=DateTime_AEDT, y=pres_pa))+facet_wrap(~burstID, scales="free")+scale_colour_viridis()

ggplot(data=dat_floating%>%filter(wave_height!="NA"))+
  geom_point(aes(x=DateTime_AEDT, y=alt_p0_mx), colour='red')+geom_line(aes(x=DateTime_AEDT, y=alt_p0_mx), colour='red')+
  geom_point(aes(x=DateTime_AEDT, y=alt_gps), colour='green')+geom_line(aes(x=DateTime_AEDT, y=alt_gps), colour='green')+  
  facet_wrap(~burstID, scales="free")

# FYI Significant wave height = the average wave height of the top one-third highest waves

wave_temp<-dat_floating%>%group_by(burstID)%>%summarise(w_height=mean(wave_height),
                                                        w_period=mean(wave_period)) 

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


# for main fig - not used anymore. Looks better without pred line
w1<-lm(w_height~ds_hsig, data=wave_sum[-c(13, 17,18,32),])
w3<-lm(w_period~ds_tmean, data=wave_sum[-c(13, 17,18,32),])

w2<-lm(w_height~ds_hsig,data=wave_sum[wave_sum$gps_hsig<6,])

# altimeter height

p4<-ggplot()+geom_point(data=wave_sum, aes(y=w_height, x=ds_hsig))+
  scale_x_continuous(limits=c(1, 9), breaks=1:9, expand = c(0,0))+
  scale_y_continuous(limits=c(1, 4), expand = c(0,0))+
  labs(x='Wave height from altimeters (m)',y='Wave height from satellite (m)', size=4)+
  geom_text(aes(x=6.8, y=1.25), label=expression(italic(r)*" = "*"0.58, "* italic(p) < 0.001), size=4)+
  theme_bw()+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=11))

# altimeter periodicity
p5<-ggplot()+geom_point(data=wave_sum[-c(13, 17,18,32),], aes(y=w_period, x=ds_tmean))+
  scale_x_continuous(limits=c(4, 13), breaks=4:13, expand = c(0,0))+scale_y_continuous(limits=c(6, 12), breaks=6:12, expand = c(0,0))+
  labs(x='Wave period from altimeters (s)',y='Wave period from satellite (s)', size=4)+
  geom_text(aes(x=10.5, y=6.5), label=expression(italic(r)*" = "*"0.86, "* italic(p) < 0.001), size=4)+
  theme_bw()+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=11))

# GPS height
p6<-ggplot()+geom_point(data=wave_sum, aes(y=w_height, x=gps_hsig))+
   scale_x_continuous(limits=c(1, 9), breaks=1:9, expand = c(0,0))+
  scale_y_continuous(limits=c(1, 4), expand = c(0,0))+
  labs(x='Wave height from GPS (m)',y='Wave height from satellite (m)', size=4)+
  geom_text(aes(x=7, y=1.25), label=expression(italic(r)*" = "*"0.37, "* italic(p)*" = "*0.03), size=4)+
  theme_bw()+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=11))

# GPS periodicity
#ns
p7<-ggplot()+geom_point(data=wave_sum[wave_sum$gps_hsig<6,], aes(y=w_period, x=gps_tmean))+
   scale_y_continuous(limits=c(6, 12), breaks=6:12, expand = c(0,0))+
  labs(x='Wave period from GPS (s)',y='Wave period from satellite (s)', size=4)+
  geom_text(aes(x=42, y=6.5), label=expression(italic(p)*" = NS"), size=4)+
  theme_bw()+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=11))

(p4+p6)/(p5+p7)

# Make mega plot!! - Not used not, split into two figures
#areas <- c(patchwork::area(1, 1, 1, 3),patchwork::area(2, 1, 2, 2), patchwork::area(2, 3, 2,3))
#p1 + (p2/p3) + (p4/p5) + plot_layout(design = areas)

#### ^^ ####


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
