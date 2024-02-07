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

#### Detect loops #### https://stackoverflow.com/questions/73404664/detecting-looping-behavior-in-track-data

sf_use_s2(FALSE)# needs to be set otherwise get errors from overlapping vertices in polys

dat$loop=F

for(i in unique(dat$burstID))
{
  track_dat <- dat[dat$burstID==i,]%>%
  st_as_sf(coords = c("Longitude","Latitude"), crs = 4326)

track_ls<-track_dat%>%st_combine() |> st_cast("LINESTRING")

# could use this code to remove either v small or large loops track_poly_df <- 
#track_polys |>  st_as_sf(crs = 4326) |> st_set_geometry("geometry") |> 
#  mutate( size = st_area(geometry) ) |> arrange(desc(size)) 

track_polys <- st_intersection(track_ls,track_ls) |>st_polygonize() |> st_cast() %>%st_union()

int1<-st_intersects(track_dat, track_polys, sparse=F)

if(!TRUE %in% int1){next} # if not loops skip to next

dat[dat$burstID==i,][unlist(int1),]$loop<-"T"
}

# finally set all sitting points that might be looped over as not in the loop

#ggplot()+geom_sf(data=track_polys, fill='red')+geom_sf(data=track_ls)+geom_sf(data=track_dat[unlist(int1),], col='green', aes(shape=sit_fly))

#### ^^ ####

#### Calculation of flight height using dynamic soaring method ####

# barometric formula (Berberan Santos et al. 1997)
#h=((k*T)/(m*g))*ln(p/p0)
k=8.31432
m=0.0289644
g=9.80665

# for p0 per burst we want to identify pressure level at sea level, we can assume waves
# so finding upper 95th? quantile could represent burst mean sea level, negative values 
# in resultant flight height would be wave troughs. A more precise method would find the 
# max pressure point of each osscillation (after a short window smooth to remove error),
# and then take the mean of these points, you could even run a gam through these if we
# expect local pressure changes to change over the burst - possible on transiting bursts.
# Question we need answered if birds always meet the surface during their soaring oscillation
# Anchoring each oscillaition to the sea surface removes negative values but probably 
# overestimates.. hmm also look at a small moving window to clean up 

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

#### Calculation of flight height using satellite ocean data (Johnston et al 2023) + 9m GPS offset (final line of code) ####

# Work out difference between pressure of sitting bursts and ECMWF.ERA5.SL.Mean.Sea.Level.Pressure
# then use value to calibrate satellite data to 'true' surface pressure (p0). Apply 'true'
# p0 value sitting bursts to flying bursts within 1 day. Nearest (in time) sitting burst has priority

# analysis run per logger as each will have unique sensor calibration


ggplot(data=dat)+geom_point(aes(x=pres_pa, y=mean_sea_level_pressure))+
  geom_abline()+facet_wrap(~class, scales="free")

ggplot(data=dat%>%filter(class=="S")%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_point(aes(x=index, y=pres_pa), col='red')+
  facet_wrap(~ID, scales="free")

# check each individually
p1<-ggplot(data=dat%>%filter(ID==08611649)%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_point(aes(x=index, y=pres_pa, colour=embc, shape=class))
p2<-ggplot(data=dat%>%filter(ID==08611649)%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_point(aes(x=index, y=pres_pa, colour=sit_fly, shape=class))

p1/p2

p1<-ggplot(data=dat%>%filter(ID==8611854 )%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_point(aes(x=index, y=pres_pa, colour=embc, shape=class))
p2<-ggplot(data=dat%>%filter(ID==8611854 )%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_point(aes(x=index, y=pres_pa, colour=sit_fly, shape=class))

p1/p2

p1<-ggplot(data=dat%>%filter(ID==41490936   )%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_point(aes(x=index, y=pres_pa, colour=embc, shape=class))
p2<-ggplot(data=dat%>%filter(ID==41490936   )%>%mutate(index=1:nrow(.)))+geom_point(aes(x=index, y=mean_sea_level_pressure), col="green")+
  geom_point(aes(x=index, y=pres_pa, colour=sit_fly, shape=class))

p1/p2

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

#Add -9m correction to GPS elevation
dat$alt_gps<-dat$alt-9

#### ^^^ ####


#### Make  Fig 1 ####

p1<-ggplot(data=dat[dat$burstID=="08611854_02_71",])+geom_line(aes(x=DateTime_AEDT, y=pres_pa, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa), size=1)+scale_y_reverse(breaks=c(100800, 100850, 100900, 100950, 101000, 101050))+
  labs(y="Pressure (mb) - reversed", x="Time")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  geom_hline(yintercept=unique(dat[dat$burstID=="08611854_02_71",]$p0), colour='red')+theme_bw()+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)')


p2<-ggplot(data=dat[dat$burstID=="08611854_02_71",])+geom_line(aes(x=DateTime_AEDT, y=alt_DS, group=1))+
  geom_point(aes(x=DateTime_AEDT, y=alt_DS), size=1)+
geom_line(aes(x=DateTime_AEDT, y=alt_gps, group=1), col='green')+
  geom_point(aes(x=DateTime_AEDT, y=alt_gps), size=1, col='green')+
  labs(y="Altitude", x="Time")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)')


#### ^^^ ####

#### Make  plots 1 and 2 of 3 panel figure ####

ggplot(data=dat%>%group_by(ID)%>%mutate(index=1:n())%>%ungroup())+
  geom_line(aes(x=index, y=alt_DS), colour='black')+
  geom_line(aes(x=index, y=alt_gps), colour="green", alpha=0.5)+
  geom_line(aes(x=index, y=alt_SO), colour='orange', alpha=0.5)+
  facet_wrap(~ID, nrow=3, scales='free')+labs(y='Altitude (m)', x='Index') # ok looks good

# make fig
fig_dat<-dat%>%filter(ID==41490936)%>%mutate(index2=1:n())%>%filter(index2<5109)
p1<-ggplot(data=fig_dat)+
  geom_rect(data=fig_dat%>%filter(class=='S')%>%group_by(burstID)%>%
              summarise(xmin=min(index2), xmax=max(index2)),
            aes(xmin=xmin, xmax=xmax, ymin=-10, ymax=25),fill='grey', alpha=0.5)+
  geom_rect(xmin=4304, xmax=4614 , ymin=-10, ymax=25,colour='purple',fill=NA)+
  geom_line(aes(x=index2, y=alt_DS), colour='black', linewidth=0.1)+
  geom_line(aes(x=index2, y=alt_gps), col="green", alpha=0.5, linewidth=0.1)+
  geom_line(aes(x=index2, y=alt_SO), colour='orange', alpha=0.5, linewidth=0.1)+
  geom_point(aes(x=index2, y=alt_DS), colour='black', size=0.5)+
  geom_point(aes(x=index2, y=alt_gps), col="green", alpha=0.5, size=0.5)+
  geom_point(aes(x=index2, y=alt_SO), colour='orange', alpha=0.5, size=0.5)+
  scale_y_continuous(limits=c(-10, 25), breaks=c(-10,-5,0,5,10,15,20,25))+labs(y='Altitude (m)', x='5 minute burst index')+
  scale_x_continuous(minor_breaks=NULL,breaks = fig_dat%>%group_by(burstID)%>%summarise(min_i=min(index2))%>%arrange(min_i)%>%pull(min_i))+
  theme_bw()+ theme( axis.text=element_text(size=12),axis.title=element_text(size=14))

p2<-ggplot(data=fig_dat[fig_dat$burstID=="41490936_01_17",])+
  geom_hline(aes(yintercept=8), linetype='dotted')+
  geom_hline(aes(yintercept=0), linetype='dotted')+
geom_line(aes(x=DateTime_AEDT, y=alt_DS, group=1))+
geom_point(aes(x=DateTime_AEDT, y=alt_DS, colour=sit_fly))+
  scale_y_continuous(name='Ocean satellite date calibrated altitude (m)', breaks=seq(0, 20, 2),
  sec.axis = sec_axis(~.-8, name="Dynamic soaring calibrated altitude (m)",
breaks=seq(0, 10, 2),labels = function(x) {ifelse(x>-1, x, "")}))+
  scale_x_datetime(date_breaks = "1 min", date_labels= '%H:%M:%S', name='Burst time (AEDT)')+
theme_bw() +
  theme(legend.position= c(0.8,0.8), axis.text=element_text(size=12),axis.title=element_text(size=14),
        plot.background = element_rect(color = "purple", size = 1))+
  scale_colour_manual("Behaviour", values=c("red", "blue"),labels=c("flying", "sitting"))
#### ^^ ####
  
#### Measuring wave height w/ altimeters and making plot 3 of 3 panel figure  ####

ggplot(data=dat%>%filter(class=="S"& wave_height!="NA"))+
  geom_point(aes(x=DateTime_AEDT, y=pres_pa, colour=wave_height))+geom_line(aes(x=DateTime_AEDT, y=pres_pa))+facet_wrap(~burstID, scales="free")+scale_colour_viridis()

# ok nice, so does variance increase with wave_height
# very rough calc, need to tweak variance/mean etc for some dodgy bursts
# Significant wave height can be shown to correspond to the average wave height of the top one-third highest waves
wave_sum<-dat%>%filter(class=="S")%>%group_by(burstID)%>%summarise(p_var=var(pres_pa), w_height=mean(wave_height), 
                                                                   pres_h_0.66=quantile(alt_DS, probs=0.66)) # could do 0.6 as we already take 95% percentile

ggplot(data=wave_sum)+
  geom_point(aes(x=w_height, y=p_var)) # looks like it

ggplot(data=wave_sum)+
  geom_point(aes(x=w_height, y=pres_h_0.66))
cor.test(x=wave_sum$w_height, y=wave_sum$pres_h_0.66, method='pearson', na.action=na.omit) # looks like it

# do some stats
w1<-lm(w_height~pres_h_0.66, data=wave_sum)
w2<-lm(log(w_height)~pres_h_0.66, data=wave_sum)
w3<-lm(w_height~sqrt(pres_h_0.66), data=wave_sum)
w4<-lm(w_height~pres_h_0.66, data=wave_sum[wave_sum$w_height<3 & wave_sum$pres_h_0.66<8,])
w5<-lm(w_height~pres_h_0.66, data=wave_sum[wave_sum$pres_h_0.66<5,])

check_model(w1)
check_model(w2)
check_model(w3)
check_model(w4)
check_model(w5)

AIC(w1, w2, w3, w4, w5)

anova(w1)
anova(w2)
anova(w3) 
anova(w4) # thats the one
anova(w5) # thats the one

summary(w5)

new_d<-data.frame(pres_h_0.66=seq(0, 10, 0.2))
new_d<-cbind(new_d, predict(w5, new_d, se.fit = T)[c('fit', 'se.fit')])
new_d$lci<-new_d$fit-(new_d$se.fit*1.96)
new_d$uci<-new_d$fit+(new_d$se.fit*1.96)

p3<-ggplot()+
  geom_point(data=wave_sum, aes(y=w_height, x=pres_h_0.66))+
  geom_line(data=new_d, aes(x=pres_h_0.66, y=fit),colour='red', linewidth=1)+
  geom_ribbon(data=new_d, aes(x=pres_h_0.66, ymin=lci, ymax=uci), alpha=0.5, fill='grey')+
  theme_bw()+
  scale_x_continuous(limits=c(0, 4), breaks=0:9, expand = c(0,0))+scale_y_continuous(limits=c(0, 4), expand = c(0,0))+
  labs(x='Wave height from albatross altimeters (m)',y='Wave height from satellite (m)', size=5)+
  geom_text(aes(x=2, y=0.5), label=expression("Y = 0.149 * X + 1.55 (F=9.79"["1,34"]*", "* italic(p) < 0.001* ")"), size=4)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

# Make mega plot!!

areas <- c(area(1, 1, 1, 3),area(2, 1, 2, 2), area(2, 3, 2,3))
p1 + p2 + p3 + plot_layout(design = areas)

#### ^^ ####

#### Summarise altitude from the three methods and compare  ####

#remove first GPS fix of each burst as higher error
dat<-dat %>% group_by(burstID) %>%
  filter(row_number()!=1)%>%ungroup()%>%as.data.frame()

# compare differences between three methods

# format dataset for comparison # not added 1000 to all alts to make positive for Gamma

dat_flying<-dat%>%filter(class %in% c('T', 'L') & sit_fly=='fly')

dat_comp<-rbind(data.frame(method='Dynamic soaring', Altitude=dat_flying$alt_DS, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID) ,
                data.frame(method='Satellite ocean', Altitude=dat_flying$alt_SO, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID),
                data.frame(method='GPS', Altitude=dat_flying$alt_gps, Logger=as.character(dat_flying$ID), burstID=dat_flying$burstID))

#summarise
dat_comp%>%group_by(method)%>%summarise(mn_alt=mean(Altitude), sd_alt=sd(Altitude), median=median(Altitude),
                                    min=min(Altitude), max=max(Altitude),
                                    q25=quantile(Altitude, 0.25), q75=quantile(Altitude, 0.75))
# make plot
cols <- scales::hue_pal()(3)
cols.alpha<-c(grDevices::adjustcolor(cols[1], alpha.f = 0.5),
        grDevices::adjustcolor(cols[2], alpha.f = 0.5),
        grDevices::adjustcolor(cols[3], alpha.f = 0.5))

ggplot(data=dat_comp)+geom_density(aes(x=Altitude, colour=method), fill=NA, size=2)+
  theme_bw()+geom_vline(xintercept = 0, linetype='dotted')+scale_x_continuous(breaks=seq(-60,60,2))+
  scale_colour_manual(values = cols.alpha)+coord_cartesian(xlim=c(-20, 40))+
  theme(legend.position= c(0.8,0.8), axis.text=element_text(size=10),axis.title=element_text(size=12),
        legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  scale_colour_manual("Altitude estimation method", values=cols.alpha, labels=c("Dynamic soaring calibrated altimeter", 
  "GPS Altitude", "Satellite ocean data calibrated altimeter"))+labs(x="Altitude (m)", y="Density")

#m1_altD<-glmer(Altitude~method+(1|burstID:Logger), data=dat_comp, family=Gamma(link='log')) 
#initially tried gamma with +1000 added to alt, but distirbution more similar to normal so went with LMM 

m1_altD<-lmer(Altitude~method+(1|burstID:Logger), data=dat_comp) # need to formulate with nlme

m1<-lme(Altitude~method, random=~1|burstID, weights=varIdent(form=~1|method), data=dat_comp)
#resid_panel(m1)
summary(m1)

em1<-emmeans(m1, specs='method')
pairs(em1)
plot(em1, comparisons = TRUE)


anova(m1_altD)

rg1<-ref_grid(m1_altD, specs='method')
emmeans(rg1, specs='method')
contrast(rg1, method='pairwise')

# model diags
check_model(m1)



#### Make DS fig and prop time in 1m band fig and table ####

p1<-ggplot()+
  geom_vline(xintercept=0, colour='blue', size=1)+
  geom_vline(xintercept=10, colour='blue', size=0.5)+
  geom_vline(xintercept=20, colour='blue', size=0.5)+
  geom_rect(aes(xmin=30, xmax=31, ymin=0, ymax=4300), fill='red', size=0.5, alpha=0.3)+
  geom_histogram(data=dat_flying, aes(x=alt_DS),colour=1, binwidth = 1)+
  scale_x_continuous(breaks=seq(-5, 50, 5),minor_breaks=seq(-5, 50, 1), limits=c(-5, 31))+
  labs(x="Altitude (m)", y="Count of altitude datapoints")+theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

prop_t<-data.frame(cut(dat_flying$alt_DS, breaks = c(-5, 1:24)) %>% table)%>%
  mutate(Proportion=Freq/length(dat_flying$alt_DS))

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

prop_t$Proportion<-specify_decimal(prop_t$Proportion, 5)
prop_t$Altitude=c("<1m", "1-2m",   "2-3m",   "3-4m",   "4-5m",   "5-6m",   "6-7m",   "7-8m",   "8-9m",   "9-10m",
"10-11m", "11-12m", "12-13m", "13-14m", "14-15m", "15-16m", "16-17m", "17-18m", "18-19m", "19-20m", "20-21m", "21-22m", "22-23m", "23-24m")
sum(prop_t$Proportion)

library(gridExtra)
p1 + gridExtra::tableGrob(prop_t[c('Altitude', 'Proportion')], row=NULL)

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
