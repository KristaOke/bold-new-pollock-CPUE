#try plotting depth * mean_station_bottemp for NEBS

#using joinboth from trawl_biomass_GAM_explor.R
#hmm still seems to remove NEBS
#use pol.trawl.dat instead from data_explor_survey.R

wd <- getwd()
# sel.trawl.dat <- read.csv(paste(wd,"/data/select_trawl_dat.csv",sep=""), row.names = 1)

#updating to use new nebs data
newdat <- read.csv(file=paste(wd,"/data/survey data/cleaned_nebs_and_2021_bot_trawl_data.csv", sep=""))



#to do this I am going to need to:
# 1 - get mean station temp for NEBS
# 2 - figure out which stations are NEBS!

#ggplot(joinboth, aes(BOT_DEPTH))

library("rnaturalearth")
library("rnaturalearthdata")
library( "ggspatial" )
library("sf")
library(mapproj)

# 
# 
# pol.trawl.dat$YEAR_factor <- as.factor(pol.trawl.dat$YEAR)
# pol.trawl.dat$BOT_TEMP[which(pol.trawl.dat$BOT_TEMP=="-9999")]<-NA
# pol.trawl.dat$SURF_TEMP[which(pol.trawl.dat$SURF_TEMP=="-9999")]<-NA
# pol.trawl.dat$WTCPUE[which(pol.trawl.dat$WTCPUE=="-9999")]<-NA
# pol.trawl.dat$NUMCPUE[which(pol.trawl.dat$NUMCPUE=="-9999")]<-NA
# 
# pol.trawl.dat$LATITUDE  <- as.numeric(pol.trawl.dat$LATITUDE)
# pol.trawl.dat$LONGITUDE <- as.numeric(pol.trawl.dat$LONGITUDE)
# pol.trawl.dat$YEAR <- as.numeric(pol.trawl.dat$YEAR )
# pol.trawl.dat$WTCPUE <- as.numeric(pol.trawl.dat$WTCPUE)
# pol.trawl.dat$NUMCPUE <- as.numeric(pol.trawl.dat$NUMCPUE)
# pol.trawl.dat$BOT_DEPTH <- as.numeric(pol.trawl.dat$BOT_DEPTH )
# pol.trawl.dat$BOT_TEMP <- as.numeric(pol.trawl.dat$BOT_TEMP)
# pol.trawl.dat$SURF_TEMP <- as.numeric(pol.trawl.dat$SURF_TEMP)
# 
# 
# pol.trawl.dat$logCPUE <- log(pol.trawl.dat$WTCPUE + 1)


world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE), data=newdat) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar")  

x <- mapproject(newdat$LONGITUDE, newdat$LATITUDE, "albers", param=c(55.9, 60.8))
newdat$long_albers <- x$x
newdat$lat_albers <- x$y



#do I need the section 'deal w missing data'?
#going to guess no so that it doesn't fill in NEBS stations for early period

#calculate anomalies and means=====================================================================

#annual 'global' (w/in dataset) mean
# poltempmeans <- pol.trawl.dat %>% group_by(STATION) %>% summarize(mean_station_bottemp=mean(BOT_TEMP, na.rm=TRUE))
# poltempmeans <- left_join(pol.trawl.dat, poltempmeans )
# 
# poltempmeans$bottemp_anom <- poltempmeans$BOT_TEMP - poltempmeans$mean_station_bottemp
# 
# all_analysis_dat <- poltempmeans



ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") 

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  # scale_colour_gradient2(low="blue", high="red", guide="colorbar") + 
  geom_text(aes(LONGITUDE, LATITUDE, 
                         label=STATION), data=all_analysis_dat)

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  # scale_colour_gradient2(low="blue", high="red", guide="colorbar") + 
  geom_text(aes(LONGITUDE, LATITUDE, 
                label=STRATUM), data=all_analysis_dat)

post2010 <- all_analysis_dat[which(all_analysis_dat$YEAR>2009),]

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  # scale_colour_gradient2(low="blue", high="red", guide="colorbar") + 
  geom_text(aes(LONGITUDE, LATITUDE, 
                label=STATION), data=post2010[which(post2010$YEAR=="2010"|
                                                      post2010$YEAR=="2017"|
                                                      post2010$YEAR=="2018"|
                                                      post2010$YEAR=="2019"),]) + facet_wrap(~YEAR)


#sort into NEBS and SEBS

newdat$region <- "SEBS"
newdat$region[which(newdat$STRATUM==81 | 
                                newdat$STRATUM==70 |
                                newdat$STRATUM==71)] <- "NEBS"

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  # scale_colour_gradient2(low="blue", high="red", guide="colorbar") + 
  geom_point(aes(LONGITUDE, LATITUDE, 
                col=region), data=newdat) 

#add shelves 
#based on table 1 in Laurth et al 2019 NOAA Technical Memorandum NMFS-AFSC-396

newdat$shelf <- NA
newdat$shelf[which(newdat$STRATUM==81 | 
                     newdat$STRATUM==70 |
                     newdat$STRATUM==71)] <- "NEBS"
newdat$shelf[which(newdat$STRATUM==10 | 
                     newdat$STRATUM==20)] <- "EBS_inner"
newdat$shelf[which(newdat$STRATUM==31 | 
                     newdat$STRATUM==32 | 
                     newdat$STRATUM==41 | 
                     newdat$STRATUM==42 | 
                     newdat$STRATUM==43 | 
                     newdat$STRATUM==82)] <- "EBS_middle"
newdat$shelf[which(newdat$STRATUM==50 | 
                     newdat$STRATUM==61 | 
                     newdat$STRATUM==62 | 
                     newdat$STRATUM==90)] <- "EBS_outer"

load(paste(wd,"/data/Bering_Chukchi_contours.RData",sep=""))

contour.lines.adj <- contour.lines
contour.lines.adj$z20$x_adj <- contour.lines.adj$z20$x - 360
contour.lines.adj$z50$x_adj <- contour.lines.adj$z50$x - 360
contour.lines.adj$z100$x_adj <- contour.lines.adj$z100$x - 360
contour.lines.adj$z200$x_adj <- contour.lines.adj$z200$x - 360
contour.lines.adj$z1000$x_adj <- contour.lines.adj$z1000$x - 360

#Manuscript
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  # scale_colour_gradient2(low="blue", high="red", guide="colorbar") + 
  geom_point(aes(LONGITUDE, LATITUDE, 
                 col=region), data=newdat) + theme_bw() + 
  scale_color_manual(name="Region", 
   labels = c("NEBS", "SEBS"),
   values = c("NEBS"="#d8b365", "SEBS"="#5ab4ac")) +
  theme( legend.position = c(0.87, 0.85), legend.key = element_blank(),
         legend.background=element_blank()) + #geom_path(aes(x_adj,y), data=contour.lines.adj$z20, col="dark blue") +
  geom_path(aes(x_adj,y), data=contour.lines.adj$z50, col="#9ecae1") + geom_path(aes(x_adj,y), data=contour.lines.adj$z100, col="#3182bd") +
  geom_path(aes(x_adj,y), data=contour.lines.adj$z200, col="navy blue") #+ geom_path(aes(x_adj,y), data=contour.lines.adj$z1000, col="dark green") 


#poster
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  # scale_colour_gradient2(low="blue", high="red", guide="colorbar") + 
  geom_point(aes(LONGITUDE, LATITUDE, 
                 col=region), data=newdat[which(newdat$STRATUM!="0"),], size=0.5) + theme_bw() + 
  scale_color_manual(name="Region", 
                     labels = c("NEBS", "SEBS"),
                     values = c("NEBS"="#d8b365", "SEBS"="#5ab4ac")) +
  theme( legend.position = c(0.85, 0.85), legend.key = element_blank(),
         legend.background=element_blank()) + #geom_path(aes(x_adj,y), data=contour.lines.adj$z20, col="dark blue") +
  geom_path(aes(x_adj,y), data=contour.lines.adj$z50, col="#9ecae1") + geom_path(aes(x_adj,y), data=contour.lines.adj$z100, col="#3182bd") +
  geom_path(aes(x_adj,y), data=contour.lines.adj$z200, col="navy blue") #+ geom_path(aes(x_adj,y), data=contour.lines.adj$z1000, col="dark green") 

#area of detail map
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-177.5, -120), ylim = c(45, 69), expand = TRUE) + theme_bw() +
  annotate("rect", xmin = -178, xmax = -155, ymin = 53, ymax = 65,
           alpha = .1,fill = "blue") + annotate("text", x=-151, y=63, label= "Alaska", size=9) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


newdat$period <- NA

newdat$period[which(newdat$YEAR<2014)] <- "early"
newdat$period[which(newdat$YEAR>2013)] <- "late"

x <- mapproject(newdat$LONGITUDE, newdat$LATITUDE, "albers", param=c(55.9, 60.8))
newdat$long_albers <- x$x
newdat$lat_albers <- x$y

#plot temps * depth=================

r1 <- ggplot(newdat, aes(mean_station_bottemp, BOT_DEPTH))
r1 + geom_point() + facet_wrap(~region)


r2 <- ggplot(newdat, aes(mean_station_bottemp, BOT_DEPTH, col=region))
r2 + geom_point()

r3 <- ggplot(newdat, aes(mean_station_bottemp, BOT_DEPTH))
r3 + geom_point() + facet_wrap(~shelf)

r4 <- ggplot(newdat, aes(mean_station_bottemp, BOT_DEPTH, col=shelf))
r4 + geom_point() 

