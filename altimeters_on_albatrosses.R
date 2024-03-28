# Code to conduct analyses for manuscript: Altimeters on albatrosses quantifying flight heights for dynamic soaring seabirds

#### Housekeeping and libraries ####

rm(list=ls())

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

setwd("C:/myloc")

#### ^^ ####

#### Read in data ####

dat<-read.csv("altimeters_on_albatrosses_public_data.csv", h=T)

dat$DateTime_AEDT<-ymd_hms(dat$DateTime_AEDT, tz="Australia/Sydney")

#Add -9m correction to GPS elevation due to geoid model mismatch
dat$alt_gps<-dat$alt-9

dat%>%group_by(class)%>%summarise(n(unique(burstID)))
# S=sitting, A=takeoff/landing, L&T=flying, C=colony

#### ^^ ####

#### Summary table of env variables (Table 2) ####
tabl1<-dat%>%filter(class %in% c("T", "L", "A"))%>%dplyr::select(ColDist , dist2coast, wind_speed, wind_dir,
                                                          chla, sst, wave_height, wave_period)%>%
                            summarise_all(list(mean=mean, sd=sd, qz=quantile), na.rm=T)

#write_xlsx(tabl1, 'analyses/env_summary_table.xlsx')

#### ^^ ####

#### Calculation of flight height using dynamic soaring altimeter calibration method ####

# barometric formula (Berberan Santos et al. 1997)
#h=((k*T)/(m*g))*ln(p/p0)
k=8.31432
m=0.0289644
g=9.80665

dat$index<-1:nrow(dat)
dat$p0<-0

dat$alt_DS<-NA
for (i in unique(dat$burstID))
{
  # 95% upper quantile of pressure to set p0 for entire burst
  dat[dat$burstID==i,]$p0<-quantile(dat[dat$burstID==i,]$pres_pa, probs=0.95)
}
  
dat$alt_DS<-(-1*  # *-1 flips negative/positive values
             ((k*(dat$temp+273.15))/(m*g))*log(dat$pres_pa/dat$p0))

#### ^^^ ####

#### Calculation of flight height using 'sitting satellite' altimeter calibration (sensu Johnston et al 2023) ####

# Work out difference between pressure of sitting bursts and ECMWF.ERA5.SL.Mean.Sea.Level.Pressure
# then use value to calibrate satellite data to 'true' surface pressure (p0). Apply 'true'
# p0 value sitting bursts to flying bursts within 1 day. Nearest (in time) sitting burst has priority
# analysis run per logger as each will have unique sensor calibration

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

# check outputs - reproduce Johnston et al (2023) Fig 2
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

#calculate p0 for 'sitting satellite' method, than apply barometric formula
dat$p0_SO<-ifelse(is.na(dat$sat_sit_pdiff),dat$mean_sea_level_pressure+dat$nearest_sat_sit_pdiff,dat$mean_sea_level_pressure+dat$sat_sit_pdiff)
dat$alt_SO<-(-1*  # *-1 flips negative/positive values
               ((k*(dat$temp+273.15))/(m*g))*log(dat$pres_pa/dat$p0_SO))

#### ^^^ ####


#### Zero-crossing analysis to get periodicity of dynamic soaring and wave metrics when sitting water ####

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
    error=function(e) e) # catches one or two NA bursts
  
  if(!inherits(possibleError, "error")){
  g1<-waveStatsZC(dat[dat$burstID==i,]%>%filter(!alt_gps %in% boxplot(alt_gps)$out)%>%pull(alt_gps), 1,) 
  tout$gps_hsig=g1$Hsig
  tout$gps_hmean=g1$Hmean
  tout$gps_tmean=g1$Tmean
  tout$gps_tsig=g1$Tsig
  }else{}
  
  zc_summary<-rbind(zc_summary, tout)
}
# two altimeter methods identical; sitting satellite method not run

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


#### Make  Fig 2 ####

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

render_snapshot("projects/02 flight heights/writeup/3dplot.png", clear = T)

p2<-ggplot(data=dat[dat$burstID=="08611854_04_122",])+geom_line(aes(x=DateTime_AEDT, y=alt_DS, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=alt_DS), size=1)+
geom_line(aes(x=DateTime_AEDT, y=alt_gps, group=1), col='#00A9FF')+
  geom_point(aes(x=DateTime_AEDT, y=alt_gps), size=1, col='#00A9FF')+
  labs(y="Altitude (m)")+
  geom_hline(yintercept=0, linetype='dotted')+theme_bw()+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  scale_y_continuous(limits=c(-2, 26), breaks=seq(-2,26,2), minor_breaks = NULL)+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)')

p1.6 <- fig("02 flight heights/writeup/3dplot_zoom.png")
  
wrap_plots(p1, p1.6, p2, nrow=3)
p1+p1.6+p2 + plot_layout(nrow=3, heights = c(1,2))+ plot_annotation(tag_levels = 'a',tag_suffix = ')')+theme(plot.tag = element_text(size = 14))

p1/p1.6/p2

#### ^^^ ####
  
#### Compare wave height and periodicity measured from albatrosses and satellite  ####

#summarise first
dat%>%filter(class=="S"& wave_height!="NA")%>%
  summarise(mn_gps=mean(alt_gps), sd_gps=sd(alt_gps), mn_ds=mean(alt_DS), sd_ds=sd(alt_DS), 
                mn_wh=mean(wave_height), sd_wh=sd(wave_height)) 

ggplot(data=dat%>%filter(class=="S"& wave_height!="NA"))+
  geom_point(aes(x=DateTime_AEDT, y=alt_DS), colour='red')+geom_line(aes(x=DateTime_AEDT, y=alt_DS), colour='red')+
  geom_point(aes(x=DateTime_AEDT, y=alt_gps), colour='green')+geom_line(aes(x=DateTime_AEDT, y=alt_gps), colour='green')+  
  facet_wrap(~burstID, scales="free") # gps error makes it hard to see

# Significant wave height can be shown to correspond to the average wave height of the top one-third highest waves
wave_temp<-dat%>%filter(class=="S")%>%group_by(burstID)%>%summarise(w_height=mean(wave_height),
                                                        w_period=mean(wave_period)) 

wave_sum<-left_join(zc_summary%>%filter(class=="S"), wave_temp, by="burstID")

# lms tell us which outliers to remove before pearsons correlation
w1<-lm(w_height~ds_hsig, data=wave_sum)
w2<-lm(w_height~gps_hsig, data=wave_sum)
w3<-lm(w_period~ds_tmean, data=wave_sum)
w4<-lm(w_period~gps_tmean, data=wave_sum)

check_model(w1)
check_model(w2)
check_model(w3)
check_model(w4)
check_model(w5)

# run correlations with few outlier bursts removed (these are bursts with large error due to GPS error or pressure sensor 'drift' while sitting)
cor.test(x=wave_sum[-c(13, 17,18,32),]$w_height, y=wave_sum[-c(13, 17,18,32),]$ds_hsig, method='pearson', na.action=na.omit) 
cor.test(x=wave_sum[wave_sum$gps_hsig<6,]$w_height, y=wave_sum[wave_sum$gps_hsig<6,]$gps_hsig, method='pearson', na.acmtion=na.omit)

cor.test(x=wave_sum[-c(13, 17,18,32),]$w_period, y=wave_sum[-c(13, 17,18,32),]$ds_tmean, method='pearson', na.action=na.omit)
cor.test(x=wave_sum[wave_sum$gps_hsig<6,]$w_period, y=wave_sum[wave_sum$gps_hsig<6,]$gps_tmean, method='pearson', na.acmtion=na.omit)

#Make plots for supplementary material
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

(wp1+wp2)/(wp3+wp4)

#### ^^ ####

#### Make 4-plot Fig 4 ####

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


# refit wave height and periodicity models and predict relationship
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

areas <- c(patchwork::area(1, 1, 1, 3),patchwork::area(2, 1, 2, 2), patchwork::area(2, 3, 2,3))
p1 + p2 + (p3/p4) + plot_layout(design = areas)

#### ^^ ####

#### Summarise and compare altitude from dynamic soaring calibrated altimeters,
#sitting satellite calibrated altimeters and GPS altitude  ####

#remove first GPS fix of each burst as v high error
dat<-dat %>% group_by(burstID) %>%
  filter(row_number()!=1)%>%ungroup()%>%as.data.frame()

dat_flying<-dat%>%filter(class %in% c('T', 'L') & sit_fly=='fly') # select flying only data


dat_comp<-rbind(data.frame(method='Dynamic soaring', Altitude=dat_flying$alt_DS, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID) ,
                data.frame(method='Satellite ocean', Altitude=dat_flying$alt_SO, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID),
                data.frame(method='GPS', Altitude=dat_flying$alt_gps, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID))

#summarise
dat_comp%>%group_by(method)%>%summarise(mn_alt=mean(Altitude), sd_alt=sd(Altitude), median=median(Altitude),
                                    min=min(Altitude), max=max(Altitude),
                                    q25=quantile(Altitude, 0.25), q75=quantile(Altitude, 0.75), skew=skewness(Altitude))
# make Fig 3
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

# Now run stats on difference in altitude between methods rather than actual altitude values to better respect normality
dat_diff<-rbind(data.frame(method='Dynamic soaring - GPS', Altitude=dat_flying$alt_DS-dat_flying$alt_gps, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID) ,
                data.frame(method='Dynamic soaring - Satellite ocean', Altitude=dat_flying$alt_DS-dat_flying$alt_SO, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID),
                data.frame(method='GPS - Satellite ocean', Altitude=dat_flying$alt_gps-dat_flying$alt_SO, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID))

m1<-lme(Altitude~method, random=~1|burstID, weights=varIdent(form=~1|method), data=dat_diff)
boxplot(residuals(m1, type='pearson')~dat_diff$method)
m2<-lme(Altitude~method, random=~1|Logger, weights=varIdent(form=~1|method), data=dat_diff)

AIC(m1, m2) # use BurstID as RE as logger not as good
resid_panel(m1)
summary(m1)

anova(m1)
em1<-emmeans(m1, specs='method')
em1
test(em1, adjust="bonferroni")

#### Proportion of time in 1m band calculation and Fig 6 ####

min(dat_flying$alt_DS) # get offset to make all vals above zero for lognormal distribution
#4.166551

lognormal.SH = fitdistr(dat_flying$alt_DS+4.166551, "lognormal") # add offset and fit model
lognormal.SH.Fit = dlnorm(seq(0,30,0.1), lognormal.SH$estimate[1], lognormal.SH$estimate[2]) # predict over 30m
plot(seq(0,30,0.1), lognormal.SH.Fit, type = "l", ylab = "Prop. at height", xlab = "Height above Sea-level")

#bootstrap as per Johnston et al (2014), but replacing site iwth loggerID: "using the site as the bootstrap unit, 200 bootstrap samples were produced,
#with a balanced design, such that each site appeared 200 times across all bootstraps""
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

dimnames(lnorm.boot)[[1]]<-(seq(0,300,0.1)-4.166551)+1.15 # add +1.15 wing length offset
dimnames(lnorm.boot)[[2]]<-c(paste0('bootId_', 1:200), "mean")

# make proportion table, up to 300m, and with all 200 boots at 1m intervals
propt<-lnorm.boot%>%mutate(ht_1m=cut(as.numeric(row.names(lnorm.boot)), breaks = c(-5, 1:297)))%>%group_by(ht_1m)%>%summarise_all(sum)

propt[,2:202]<-propt[,2:202]/10

propt[,203]<-apply(propt[,2:201], 1,  function(x){confint(lm(x~1), level=0.95)[1]})
propt[,204]<-apply(propt[,2:201], 1,  function(x){confint(lm(x~1), level=0.95)[2]})
colnames(propt)[203:204]<-c("lci", "uci")
propt$ht_1m<-gsub( ",", "-",propt$ht_1m)

#Make Fig 6 plot
boot_for_plot<-data.frame(ht=as.numeric(row.names(lnorm.boot[as.numeric(row.names(lnorm.boot))<max(dat_flying$alt_DS+1.15),])),
                          prop = lnorm.boot[as.numeric(row.names(lnorm.boot))<max(dat_flying$alt_DS+1.15),"mean"]) # add +1.15 wing length offset

p1<-ggplot()+
  geom_vline(xintercept=0, colour='blue', size=1)+
  geom_vline(xintercept=10, colour='blue', size=0.5)+
  geom_vline(xintercept=20, colour='blue', size=0.5)+
  geom_rect(aes(xmin=30, xmax=31, ymin=0, ymax=0.17), fill='red', size=0.5, alpha=0.3)+
  geom_histogram(data=dat_flying, aes(x=alt_DS+1.15, after_stat(density)), fill='grey', colour='darkgrey', binwidth = 1)+ # add +1.15 wing length offset
  geom_line(data=boot_for_plot, aes(x=ht, y=prop), size=1)+
  scale_x_continuous(breaks=seq(-5, 50, 5),minor_breaks=seq(-5, 50, 1), limits=c(-5, 31), expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(x="Height above sea level (m)", y="Proportion at height")+theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

# Fig 6 table (up to 25m)

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