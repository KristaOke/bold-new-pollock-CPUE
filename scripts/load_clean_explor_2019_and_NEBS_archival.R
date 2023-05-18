#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Load, clean, explore 2021 data and 1982-2019 online NEBS csv  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Notes: 
#found out that the publicly available 1982-2019 csv for the NEBS has more years
#and has 2018 stations (and 80s-90s stations?) automatically renamed to match standard stations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


library("rnaturalearth")
library("rnaturalearthdata")
library( "ggspatial" )
library("sf")
library(ggplot2)
library(tidyverse)


#load data-------

wd <- getwd()
nebs2021dat <- read.csv(paste(wd,"/data/survey data/GAP_BottomTrawl_SurveyData_2021/NBS2021.csv",sep=""))
sebs2021dat <- read.csv(paste(wd,"/data/survey data/GAP_BottomTrawl_SurveyData_2021/EBS2021.csv",sep=""))

nebs82_19 <- read.csv(paste(wd,"/data/survey data/nbs1982_2019.csv",sep=""))

#all same column names which is nice

newdat <- rbind(nebs2021dat, sebs2021dat, nebs82_19)

#will limit to SID 21740 but not yet

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  # scale_colour_gradient2(low="blue", high="red", guide="colorbar") + 
  geom_point(aes(LONGITUDE, LATITUDE, 
                 col=as.factor(STRATUM)), data=newdat) + theme_bw() 

#look at each year
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, 
                 col=as.factor(STRATUM)), data=newdat) + theme_bw() + facet_wrap(~YEAR)

#look at only early year
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, 
                 col=as.factor(STRATUM)), data=newdat[which(newdat$YEAR<1992),]) + theme_bw() + facet_wrap(~YEAR)

#1985 and 1991 stations need stratum assigned, 1988 is missing 1 or 2 stratums
#Norton Sound stations are 0 in 1982, no stratum otherwise, assign?

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(60, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_text(aes(LONGITUDE, LATITUDE, 
                 col=as.factor(STRATUM),
                label=STATION), data=newdat[which(newdat$YEAR<1992),]) + theme_bw() + facet_wrap(~YEAR)


ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(60, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_text(aes(LONGITUDE, LATITUDE, 
                col=as.factor(STRATUM),
                label=STATION), data=newdat[which(newdat$YEAR==2018),]) + theme_bw() + facet_wrap(~YEAR)

newdat$STRATUM[which(newdat$STATION=="ZZ-23" & newdat$YEAR<1992|
                       newdat$STATION=="Y-24" & newdat$YEAR<1992| 
                       newdat$STATION=="X-25" & newdat$YEAR<1992| 
                       newdat$STATION=="X-23" & newdat$YEAR<1992| 
                       newdat$STATION=="T-23" & newdat$YEAR<1992| 
                       newdat$STATION=="T-21" & newdat$YEAR<1992| 
                       newdat$STATION=="R-21" & newdat$YEAR<1992| 
                       newdat$STATION=="S-21" & newdat$YEAR<1992| 
                       newdat$STATION=="W-26" & newdat$YEAR<1992| 
                       newdat$STATION=="W-24" & newdat$YEAR<1992| 
                       newdat$STATION=="V-23" & newdat$YEAR<1992| 
                       newdat$STATION=="U-22" & newdat$YEAR<1992| 
                       newdat$STATION=="AA-23" & newdat$YEAR<1992| 
                       newdat$STATION=="U-24" & newdat$YEAR<1992 )] <- "81"

newdat$STRATUM[which(newdat$STATION=="Y-22" & newdat$YEAR<1992|
                       newdat$STATION=="Y-20" & newdat$YEAR<1992| 
                       newdat$STATION=="X-21" & newdat$YEAR<1992| 
                       newdat$STATION=="X-19" & newdat$YEAR<1992| 
                       newdat$STATION=="W-22" & newdat$YEAR<1992| 
                       newdat$STATION=="W-20" & newdat$YEAR<1992| 
                       newdat$STATION=="V-21" & newdat$YEAR<1992| 
                       newdat$STATION=="V-19" & newdat$YEAR<1992| 
                       newdat$STATION=="U-20" & newdat$YEAR<1992| 
                       newdat$STATION=="U-18" & newdat$YEAR<1992| 
                       newdat$STATION=="T-19" & newdat$YEAR<1992| 
                       newdat$STATION=="U-02" & newdat$YEAR<1992|
                       newdat$STATION=="T-01" & newdat$YEAR<1992| 
                       newdat$STATION=="T-03" & newdat$YEAR<1992| 
                       newdat$STATION=="S-20" & newdat$YEAR<1992| 
                       newdat$STATION=="S-18" & newdat$YEAR<1992| 
                       newdat$STATION=="R-19" & newdat$YEAR<1992| 
                       newdat$STATION=="R-01" & newdat$YEAR<1992| 
                       newdat$STATION=="R-03" & newdat$YEAR<1992| 
                       newdat$STATION=="V-03" & newdat$YEAR<1992| 
                       newdat$STATION=="W-02" & newdat$YEAR<1992| 
                       newdat$STATION=="X-03" & newdat$YEAR<1992| 
                       newdat$STATION=="Y-18" & newdat$YEAR<1992| 
                       newdat$STATION=="X-01" & newdat$YEAR<1992| 
                       newdat$STATION=="W-18" & newdat$YEAR<1992| 
                       newdat$STATION=="V-01" & newdat$YEAR<1992| 
                       newdat$STATION=="S-02" & newdat$YEAR<1992| 
                       newdat$STATION=="Y-02" & newdat$YEAR<1992 )] <- "70"

newdat$STRATUM[which(newdat$STATION=="DD-19" & newdat$YEAR==2018|
                       newdat$STATION=="CC-19" & newdat$YEAR==2018| 
                       newdat$STATION=="AA-19" & newdat$YEAR==2018)] <- "71"


newdat$BOT_TEMP[which(newdat$BOT_TEMP=="-9999")]<-NA
newdat$SURF_TEMP[which(newdat$SURF_TEMP=="-9999")]<-NA
newdat$WTCPUE[which(newdat$WTCPUE=="-9999")]<-NA
newdat$NUMCPUE[which(newdat$NUMCPUE=="-9999")]<-NA

#look at 2021 wider area

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(55, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_text(aes(LONGITUDE, LATITUDE, 
                col=as.factor(STRATUM),
                label=STATION), data=newdat[which(newdat$YEAR==2021),]) + theme_bw() + facet_wrap(~YEAR)

table(newdat$YEAR)

#look at other years one at a time
#82,85,88,91
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(60, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_text(aes(LONGITUDE, LATITUDE, 
                col=as.factor(STRATUM),
                label=STATION), data=newdat[which(newdat$YEAR==2019),]) + theme_bw() + facet_wrap(~YEAR)


#lets look at CPUE
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(60, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, 
                col=NUMCPUE,
                label=STATION), data=newdat[which(newdat$SID==21740),]) + theme_bw() + facet_wrap(~YEAR)

newdat$logCPUE <- log(newdat$WTCPUE + 1)

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(60, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, 
                 col=logCPUE,
                 label=STATION), data=newdat[which(newdat$SID==21740),]) + 
  scale_colour_distiller(palette = "Spectral") + theme_bw() + facet_wrap(~YEAR)

#and temp
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(60, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, 
                 col=BOT_TEMP,
                 label=STATION), data=newdat[which(newdat$SID==21740),]) + 
  scale_colour_distiller(palette = "Spectral") + theme_bw() + facet_wrap(~YEAR)

range(newdat$BOT_TEMP)

#need julian day
library(lubridate)
#split date and time
newdat <- separate(newdat, DATETIME, into = c("Date", "Time"), sep = " ")

newdat$julian <- yday(parse_date_time(newdat$Date, orders = "mdy"))

#check for duplicates
unique(duplicated(newdat)) #doesn't find any duplicates

#quick look
ggplot(newdat, aes(julian, logCPUE)) + geom_point() +
  facet_wrap(~YEAR)

ggplot(newdat, aes(BOT_DEPTH, logCPUE)) + geom_point() +
  facet_wrap(~YEAR)

ggplot(newdat, aes(BOT_TEMP, logCPUE)) + geom_point() +
  facet_wrap(~YEAR)

ggplot(newdat, aes(SURF_TEMP, logCPUE)) + geom_point() +
  facet_wrap(~YEAR)

#write csv for cleaned new data 
wd <- getwd()
write_csv(newdat, file=paste(wd,"/data/survey data/cleaned_nebs_and_2021_bot_trawl_data.csv", sep=""))




