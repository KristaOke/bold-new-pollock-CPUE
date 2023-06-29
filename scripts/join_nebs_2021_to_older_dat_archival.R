#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Join new NEBS and 2021 data to older sebs data

#Krista, May 2022
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Notes:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(lubridate)
library(tidyverse)
library(datatable)
library(mapproj)

#load data-----
#load cleaned new data 
wd <- getwd()
newdat <- read.csv(paste(wd,"/data/survey data/cleaned_nebs_and_2021_bot_trawl_data.csv", sep=""))

#load other data also cleaned
#periods_analysis_dat <- read.csv(paste(wd,"/data/processed_periods_analysis_data.csv", sep=""), row.names = 1)

#long format older data
sel.trawl.dat <- read.csv(paste(wd,"/data/select_trawl_dat.csv", sep=""), row.names = 1)

#cleaning  -----
sel.trawl.dat$YEAR_factor <- as.factor(sel.trawl.dat$YEAR)
sel.trawl.dat$BOT_TEMP[which(sel.trawl.dat$BOT_TEMP=="-9999")]<-NA
sel.trawl.dat$SURF_TEMP[which(sel.trawl.dat$SURF_TEMP=="-9999")]<-NA
sel.trawl.dat$WTCPUE[which(sel.trawl.dat$WTCPUE=="-9999")]<-NA
sel.trawl.dat$NUMCPUE[which(sel.trawl.dat$NUMCPUE=="-9999")]<-NA

sel.trawl.dat$LATITUDE  <- as.numeric(sel.trawl.dat$LATITUDE)
sel.trawl.dat$LONGITUDE <- as.numeric(sel.trawl.dat$LONGITUDE)
sel.trawl.dat$YEAR <- as.numeric(sel.trawl.dat$YEAR )
sel.trawl.dat$WTCPUE <- as.numeric(sel.trawl.dat$WTCPUE)
sel.trawl.dat$NUMCPUE <- as.numeric(sel.trawl.dat$NUMCPUE)
sel.trawl.dat$BOT_DEPTH <- as.numeric(sel.trawl.dat$BOT_DEPTH )
sel.trawl.dat$BOT_TEMP <- as.numeric(sel.trawl.dat$BOT_TEMP)
sel.trawl.dat$SURF_TEMP <- as.numeric(sel.trawl.dat$SURF_TEMP)


sel.trawl.dat$logCPUE <- log(sel.trawl.dat$WTCPUE + 1)

#fill in columns missing from one or the other - - - - - 
#add date, time, and julian to sel trawl dat

#split date and time
#2019 date time doesn't match!
sel.trawl.dat <- separate(sel.trawl.dat, DATETIME, into = c("Date", "Time"), sep = " ")
#this throws an error because many 2019 times are missing, but date is ok (if different format that other years)

#doing this seperately b/c there are two different date formats (some in 2019 and everything else)
sel.trawl.dat$julian[is.na(sel.trawl.dat$Time)==FALSE] <- yday(parse_date_time(sel.trawl.dat$Date[is.na(sel.trawl.dat$Time)==FALSE], orders = "mdy"))
sel.trawl.dat$julian[is.na(sel.trawl.dat$Time)==TRUE] <-  yday(parse_date_time(sel.trawl.dat$Date[is.na(sel.trawl.dat$Time)==TRUE], orders = "dmy")) 
#yep looks good

#add Year_factor to new dat
newdat$YEAR_factor <- as.factor(newdat$YEAR)

#check station assignments in sel.trawl.dat

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(55, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_text(aes(LONGITUDE, LATITUDE, 
                col=as.factor(STRATUM),
                label=STATION), data=sel.trawl.dat[which(sel.trawl.dat$YEAR>2010 & sel.trawl.dat$YEAR<2020),]) + theme_bw() + facet_wrap(~YEAR)
#station assignments look ok not seeing any issues

#remove nebs from sel.trawl.dat for join
#these are the data I am replacing
trawl_sebs_only <- sel.trawl.dat[which(sel.trawl.dat$STRATUM!="70"&
                                         sel.trawl.dat$STRATUM!="71"&
                                         sel.trawl.dat$STRATUM!="81"),]

#join----
northsouthdata_all <- rbind(trawl_sebs_only, newdat) 

#check for & remove duplicates
unique(duplicated(northsouthdata_all)) #NO duplicates


#DO THIS AFTER JOIN
#convert latitude & longitude to x/y coordinates that reflect 
# actual distances using the equal-distance Albers projection:
# (I chose the parameters based on a 'rule of thumb' in the help file)
x <- mapproject(northsouthdata_all$LONGITUDE, northsouthdata_all$LATITUDE, "albers", param=c(55.9, 60.8))
northsouthdata_all$long_albers <- x$x
northsouthdata_all$lat_albers <- x$y

#save data----
wd <- getwd()
write_csv(northsouthdata_all, file=paste(wd,"/data/survey data/combined_cleaned_north-south_1982-2021_bot_trawl_data.csv", sep=""))


                   
                   





















