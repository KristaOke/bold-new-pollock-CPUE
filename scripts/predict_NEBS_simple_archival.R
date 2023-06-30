#====================================================================================================================================
# Predict into NEBS with new simplified model
#
#Krista, Jun 2022
#====================================================================================================================================
#Notes:
#====================================================================================================================================
#
library("rnaturalearth")
library("rnaturalearthdata")
library( "ggspatial" )
library("sf")
library(tidyverse)
library(mgcv)
library(cowplot)


#import the model fit for the model with a linear interaction
cmod1_noint <- read_rds(paste(wd,"/scripts/new_no_int.RDS", sep=""))
summary(cmod1_noint)

#import the data

wd <- getwd()
northsouthdata_all <- read.csv(paste(wd,"/data/survey data/combined_cleaned_north-south_1982-2021_bot_trawl_data.csv",sep=""))

#limit to pollock data
NS_pollock <- northsouthdata_all[which(northsouthdata_all$SID=="21740"),]

#add period
NS_pollock$period <- NA
NS_pollock$period[which(NS_pollock$YEAR<2014)] <- "early"
NS_pollock$period[which(NS_pollock$YEAR>2013)] <- "late"
NS_pollock$period <- as.factor(NS_pollock$period)


pdat <- NS_pollock

#----
#sort into NEBS and SEBS

pdat$region <- "SEBS"
pdat$region[which(pdat$STRATUM==81 | 
                    pdat$STRATUM==70 |
                    pdat$STRATUM==71)] <- "NEBS"


#add shelves 
#based on table 1 in Laurth et al 2019 NOAA Technical Memorandum NMFS-AFSC-396

pdat$shelf <- NA
pdat$shelf[which(pdat$STRATUM==81 | 
                   pdat$STRATUM==70 |
                   pdat$STRATUM==71)] <- "NEBS"
pdat$shelf[which(pdat$STRATUM==10 | 
                   pdat$STRATUM==20)] <- "EBS_inner"
pdat$shelf[which(pdat$STRATUM==31 | 
                   pdat$STRATUM==32 | 
                   pdat$STRATUM==41 | 
                   pdat$STRATUM==42 | 
                   pdat$STRATUM==43 | 
                   pdat$STRATUM==82)] <- "EBS_middle"
pdat$shelf[which(pdat$STRATUM==50 | 
                   pdat$STRATUM==61 | 
                   pdat$STRATUM==62 | 
                   pdat$STRATUM==90)] <- "EBS_outer"
#----

#predict-------


pdat_NEBS <- pdat[which(pdat$region=="NEBS"),]

nebs_sel<- pdat_NEBS[,c(5, 13, 14, 19, 20:26)]
names(nebs_sel)

NEBSpred3 <- predict.gam(cmod1_noint$gam, newdata = nebs_sel)
length(NEBSpred3)
length(nebs_sel$BOT_DEPTH) #same length

pdat_NEBS$predicted <- NEBSpred3


NEBSpred4 <- predict.gam(cmod1_noint$gam, newdata = nebs_sel, type="response")
#seems same, still neg


#get difference
pdat_NEBS$difference <- pdat_NEBS$predicted - pdat_NEBS$logCPUE

#get RMSE----
vanilla_mod_rsme <- sqrt(mean((pdat_NEBS$logCPUE - pdat_NEBS$predicted)^2, na.rm=TRUE))

pdat1982 <- pdat_NEBS[which(pdat_NEBS$YEAR=="1982"),]
pdat1985 <- pdat_NEBS[which(pdat_NEBS$YEAR=="1985"),]
pdat1988 <- pdat_NEBS[which(pdat_NEBS$YEAR=="1988"),]
pdat1991 <- pdat_NEBS[which(pdat_NEBS$YEAR=="1991"),]

pdat2010 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2010"),]
pdat2017 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2017"),]
pdat2018 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2018"),]
pdat2019 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2019"),]
pdat2021 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2021"),]

vanilla_1982_rsme <- sqrt(mean((pdat1982$logCPUE - pdat1982$predicted)^2, na.rm=TRUE))
vanilla_1985_rsme <- sqrt(mean((pdat1985$logCPUE - pdat1985$predicted)^2, na.rm=TRUE))
vanilla_1988_rsme <- sqrt(mean((pdat1988$logCPUE - pdat1988$predicted)^2, na.rm=TRUE))
vanilla_1991_rsme <- sqrt(mean((pdat1991$logCPUE - pdat1991$predicted)^2, na.rm=TRUE))
vanilla_2010_rsme <- sqrt(mean((pdat2010$logCPUE - pdat2010$predicted)^2, na.rm=TRUE))
vanilla_2017_rsme <- sqrt(mean((pdat2017$logCPUE - pdat2017$predicted)^2, na.rm=TRUE))
vanilla_2018_rsme <- sqrt(mean((pdat2018$logCPUE - pdat2018$predicted)^2, na.rm=TRUE))
vanilla_2019_rsme <- sqrt(mean((pdat2019$logCPUE - pdat2019$predicted)^2, na.rm=TRUE))
vanilla_2021_rsme <- sqrt(mean((pdat2021$logCPUE - pdat2021$predicted)^2, na.rm=TRUE))

vanilla_1982_rsme
vanilla_1985_rsme
vanilla_1988_rsme
vanilla_1991_rsme
vanilla_2010_rsme
vanilla_2017_rsme
vanilla_2018_rsme 
vanilla_2019_rsme
vanilla_2021_rsme

#plot----
#pivot longer so that can plot on same scale!

plot_pred_dat <- pdat_NEBS[,c(1:3, 5,  19, 24:28)] %>% pivot_longer(!c(LATITUDE, LONGITUDE, STATION, YEAR,  region, period, shelf), 
                                                                                    names_to="response_type", values_to="value")



names(plot_pred_dat)
table(plot_pred_dat$response_type) #three coloumns, difference, logCPUE... and predicted


#plot----
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  stat_summary_2d(aes(LONGITUDE,LATITUDE,  z=value), bins = 20, fun = mean, data=plot_pred_dat) + 
  facet_wrap(~interaction( YEAR, response_type), nrow=3)  +
  scale_fill_distiller(palette = "Spectral")

#difference
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  stat_summary_2d(aes(LONGITUDE,LATITUDE,  z=value), bins = 20, fun = mean, data=plot_pred_dat[which(plot_pred_dat$response_type=="difference"),]) + 
  facet_wrap(~interaction( YEAR, response_type), nrow=3)  +
  scale_fill_distiller(palette = "Spectral")



#what about as points?

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=plot_pred_dat) + 
  facet_wrap(response_type~YEAR, nrow=3)  +
  scale_colour_distiller(palette = "Spectral") + theme_bw() +
  theme( legend.position = c(0.97, 0.25), legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank()) 

#just difference
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=plot_pred_dat[which(plot_pred_dat$response_type=="difference"),]) + 
  facet_wrap(response_type~YEAR, nrow=3)  +
  scale_colour_distiller(palette = "Spectral") + theme_bw() +
  theme( legend.position = c(0.97, 0.25), legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank()) 

#look at bottom temp
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=BOT_TEMP), 
             data=pdat_NEBS) + 
  facet_wrap(~YEAR, nrow=3)  +
  scale_colour_distiller(palette = "Spectral") + theme_bw() +
  theme( legend.position = c(0.97, 0.25), legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank()) 


#look at bottom temp across both nebs and sebs
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=BOT_TEMP), 
             data=NS_pollock) + 
  facet_wrap(~YEAR)  +
  scale_colour_distiller(palette = "Spectral") + theme_bw() +
  theme( legend.position = c(0.97, 0.25), legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank()) 


#one-to-one plots========

ggplot(pdat_NEBS, aes(predicted, logCPUE)) + geom_point() + geom_smooth(method="lm") + geom_abline(intercept=0)

ggplot(pdat_NEBS, aes(predicted, logCPUE, col=as.factor(YEAR))) + geom_point() + geom_smooth(method="lm") + geom_abline(intercept=0)
#OH this is interesting!!!!!

ggplot(pdat_NEBS, aes(predicted, logCPUE, col=BOT_TEMP)) + geom_point() + geom_smooth(method="lm") + geom_abline(intercept=0) +
  scale_colour_distiller(palette = "Spectral")

ggplot(pdat_NEBS, aes(predicted, logCPUE, col=BOT_DEPTH)) + geom_point() + geom_smooth(method="lm") + geom_abline(intercept=0)+
  scale_colour_distiller(palette = "Spectral")

ggplot(pdat_NEBS, aes(predicted, logCPUE, col=long_albers)) + geom_point() + geom_smooth(method="lm") + geom_abline(intercept=0)+
  scale_colour_distiller(palette = "Spectral") #

ggplot(pdat_NEBS, aes(predicted, logCPUE, col=lat_albers)) + geom_point() + geom_smooth(method="lm") + geom_abline(intercept=0)+
  scale_colour_distiller(palette = "Spectral") #

ggplot(pdat_NEBS, aes(predicted, logCPUE, col=as.factor(YEAR))) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
 facet_wrap(~YEAR)
  # scale_color_manual(values=c("#b2df8a", "#66c2a5", "#fc8d62", "#8da0cb"))

#adding in random effects======================================

#following examples from gamm() documentation

refa <- ranef(cmod1_noint$lme,level=3) #extract random effect for year
rownames(refa) <- substr(rownames(refa),start=5,stop=8) #rename names

## make a prediction, with random effects zero...
#p0 <- predict(b$gam,data.frame(x0=.3,x1=.6,x2=.98,x3=.77))
p0_21 <- predict(cmod1_noint$gam, newdata = nebs_sel[which(nebs_sel$YEAR=="2021"),])

## add in effect for fa = "2" and fb="2/4"...
p <- p0_21 + refa["2021",1] 

refa$year <- rownames(refa)

#and loop through the years that we have nebs data for
yrs <- unique(nebs_sel$YEAR)
output_df <- data.frame(matrix(ncol = 11, nrow = 0))
nms <- c("ptemp_wrand", "year", "predicted_not_conditional",
"long_albers", "lat_albers", "logCPUE", "BOT_DEPTH", "BOT_TEMP",
"julian", "region", "shelf")
colnames(output_df) <- nms
i <- 1
for(i in 1:length(yrs)){
  temp_yr <- yrs[i]
  yr_dat <- nebs_sel[which(nebs_sel$YEAR==temp_yr),]
  
  ## make a prediction, with random effects zero...
  ptemp <- predict(cmod1_noint$gam, newdata = yr_dat)
  
  ## add in effect for fa = "2" and fb="2/4"...
  ptemp_wrand <- ptemp + refa[which(refa$year==temp_yr),1] 
  
  dftemp <- as.data.frame(ptemp_wrand)
  dftemp$year <- temp_yr
  dftemp$predicted_not_conditional <- ptemp
  dftemp$long_albers <- yr_dat$long_albers
  dftemp$lat_albers <- yr_dat$lat_albers
  dftemp$logCPUE <- yr_dat$logCPUE
  dftemp$BOT_DEPTH <- yr_dat$BOT_DEPTH
  dftemp$BOT_TEMP <- yr_dat$BOT_TEMP
  dftemp$julian <- yr_dat$julian
  dftemp$region <- yr_dat$region
  dftemp$shelf <- yr_dat$shelf
  
  output_df <- rbind(output_df, dftemp)
}

output_df <- rename(output_df, predicted_cond_on_random = ptemp_wrand)

#ok nice the loop output seems to be working
#need to double check
#then need to add in metadata to loop df so that it can be used to plot etc

#look at 1:1 plot----
#not conditional
ggplot(output_df, aes(predicted_not_conditional, logCPUE, col=as.factor(year))) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~year)

#conditional on random year
ggplot(output_df, aes(predicted_cond_on_random, logCPUE, col=as.factor(year))) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~year)

#manuscript
ggplot(output_df, aes(predicted_cond_on_random, logCPUE)) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~year)

#manuscript WITH correlation (calculated in section 'try to calc cor')
ggplot(output_df, aes(predicted_cond_on_random, logCPUE)) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~year) + #geom_text(aes(res_cor_vec), data=results_cor_output)
  geom_label(data = results_cor_output, aes(x = 0.25, y = 5.5, 
                                            label = round(res_cor_vec, digits=3),
                                            label.size = NA))

#above isn't working so maybe join to other database
plotdf <- left_join(output_df, results_cor_output, by=c("year"="res_yr_vec"))

ggplot(plotdf, aes(predicted_cond_on_random, logCPUE)) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~year) + #geom_text(aes(res_cor_vec), data=results_cor_output)
  geom_label(aes(x = 0.25, y = 5.5,   label = round(res_cor_vec, digits=3)))

#plot----
#pivot longer so that can plot on same scale!

plot_out_dat <- output_df[,-c(7:9)] %>% pivot_longer(!c(lat_albers, long_albers, year,  region,  shelf), 
                                              names_to="response_type", values_to="value")


names(plot_out_dat)
table(plot_out_dat$response_type) #three coloumns, difference, logCPUE... and predicted

#match lat/long albers to lat long to plot
lats_longs <- pdat_NEBS[,c(1,2,22,23)]

pred_map_plot_dat <- left_join(plot_out_dat, lats_longs)
length(plot_out_dat$long_albers)
length(pred_map_plot_dat$long_albers) #same length

#plot----
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dat) + 
  facet_wrap(response_type~year, nrow=3)  +
  scale_colour_distiller(palette = "Spectral") + theme_bw() +
  theme( legend.position = c(0.97, 0.25), legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank()) 

#for ICES poster
pred_map_plot_dat$response_type2 <- pred_map_plot_dat$response_type
# 
pred_map_plot_dat$response_type2[which(pred_map_plot_dat$response_type=="logCPUE")] <- "Actual"
pred_map_plot_dat$response_type2[which(pred_map_plot_dat$response_type=="predicted_cond_on_random")] <- "Predicted"


ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dat[which(pred_map_plot_dat$response_type2!="predicted_not_conditional"),]) + 
  facet_grid(response_type2~year)  +
  scale_colour_distiller(palette = "Spectral") + theme_bw() +
  theme( legend.position = c(0.97, 0.25), legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank())

#try two sets of columns
library(cowplots)

#NEED TO SET TO SAME SCALE IF PLOTTING LIKE THIS

p1 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dat[which(pred_map_plot_dat$response_type2!="predicted_not_conditional"&
                                         pred_map_plot_dat$year<2000 ),]) + 
  facet_grid(year~response_type2)  +
  scale_colour_distiller(palette = "Spectral", 
  limits = c(min(pred_map_plot_dat$value, na.rm=TRUE), max(pred_map_plot_dat$value, na.rm=TRUE))) + theme_bw() +
  theme( legend.position = "right", legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank())

p2 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dat[which(pred_map_plot_dat$response_type2!="predicted_not_conditional"&
                                            pred_map_plot_dat$year>2000 ),]) + 
  facet_grid(year~response_type2)  +
  scale_colour_distiller(palette = "Spectral", 
limits = c(min(pred_map_plot_dat$value, na.rm=TRUE), max(pred_map_plot_dat$value, na.rm=TRUE))) + theme_bw() +
  theme( legend.position = "right", legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank())

plot_grid(p1, p2)


#let's do only recent years for poster
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dat[which(pred_map_plot_dat$response_type2!="predicted_not_conditional"&
                                            pred_map_plot_dat$year>2000 ),]) + 
  facet_grid(response_type2~year)  +
  scale_colour_distiller(palette = "Spectral") + theme_bw() +
  theme( legend.position = c(0.9, 0.07), legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank(),
         legend.direction = "horizontal")




#get RMSE----
cond_mod_rsme <- sqrt(mean((output_df$logCPUE - output_df$predicted_cond_on_random)^2, na.rm=TRUE))

out1982 <- output_df[which(output_df$year=="1982"),]
out1985 <- output_df[which(output_df$year=="1985"),]
out1988 <- output_df[which(output_df$year=="1988"),]
out1991 <- output_df[which(output_df$year=="1991"),]

out2010 <- output_df[which(output_df$year=="2010"),]
out2017 <- output_df[which(output_df$year=="2017"),]
out2018 <- output_df[which(output_df$year=="2018"),]
out2019 <- output_df[which(output_df$year=="2019"),]
out2021 <- output_df[which(output_df$year=="2021"),]

cond_1982_rsme <- sqrt(mean((out1982$logCPUE - out1982$predicted_cond_on_random)^2, na.rm=TRUE))
cond_1985_rsme <- sqrt(mean((out1985$logCPUE - out1985$predicted_cond_on_random)^2, na.rm=TRUE))
cond_1988_rsme <- sqrt(mean((out1988$logCPUE - out1988$predicted_cond_on_random)^2, na.rm=TRUE))
cond_1991_rsme <- sqrt(mean((out1991$logCPUE - out1991$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2010_rsme <- sqrt(mean((out2010$logCPUE - out2010$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2017_rsme <- sqrt(mean((out2017$logCPUE - out2017$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2018_rsme <- sqrt(mean((out2018$logCPUE - out2018$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2019_rsme <- sqrt(mean((out2019$logCPUE - out2019$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2021_rsme <- sqrt(mean((out2021$logCPUE - out2021$predicted_cond_on_random)^2, na.rm=TRUE))

cond_1982_rsme
cond_1985_rsme
cond_1988_rsme
cond_1991_rsme
cond_2010_rsme
cond_2017_rsme
cond_2018_rsme 
cond_2019_rsme
cond_2021_rsme


#try to calc cor-------
#from other script
#tempmat <- output_df[,c(1,6)]

#tempcor <- cor(tempmat, method = "pearson", use = "complete.obs")

#copied below
yrs <- unique(output_df$year)
#loop through species, loop through each year, get correlation and store

res_yr_vec <- vector(mode="numeric", length=length(yrs))
res_sps_vec <- vector(mode="numeric", length=length(yrs))
res_cor_vec <- vector(mode="numeric", length=length(yrs))

#yrs <- unique(full_wide$YEAR)
i <- 1
counter <- 1
for(i in 1:length(yrs)){
  temp_yr <- yrs[i]
  output_slice <- output_df[which(output_df$year==temp_yr),]
  tempmat <- output_slice[,c(1,6)]
  
  tempcor <- cor(tempmat, method = "pearson", use = "complete.obs")

  
  #save output
  #for(l in 1:10){
    res_yr_vec[i] <- temp_yr      
    #res_sps_vec[counter] <- rownames(tempcor)[l]
    res_cor_vec[i] <- tempcor[1,2]
    
   # counter <- counter + 1
  #}
}

results_cor_output <- data.frame(res_yr_vec, #res_sps_vec, 
                                 res_cor_vec)


#re-run with time-varying model==============================================================

#tvmod2_REML is better than previous best model!

tvmod2_REML <- readRDS(paste(wd,"/scripts/time_var2_REML_lessthan180.RDS", sep=""))

#predict-------

names(nebs_sel)

NEBSpred_tv <- predict.gam(tvmod2_REML$gam, newdata = nebs_sel)
length(NEBSpred_tv)
length(nebs_sel$BOT_DEPTH) #same length

pdat_NEBS$predicted_time_varying <- NEBSpred_tv


NEBSpred4 <- predict.gam(cmod1_noint$gam, newdata = nebs_sel, type="response")
#


#get difference
#pdat_NEBS$difference <- pdat_NEBS$predicted - pdat_NEBS$logCPUE

#get RMSE----
tv_mod_rsme <- sqrt(mean((pdat_NEBS$logCPUE - pdat_NEBS$predicted_time_varying)^2, na.rm=TRUE))

pdat1982 <- pdat_NEBS[which(pdat_NEBS$YEAR=="1982"),]
pdat1985 <- pdat_NEBS[which(pdat_NEBS$YEAR=="1985"),]
pdat1988 <- pdat_NEBS[which(pdat_NEBS$YEAR=="1988"),]
pdat1991 <- pdat_NEBS[which(pdat_NEBS$YEAR=="1991"),]

pdat2010 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2010"),]
pdat2017 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2017"),]
pdat2018 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2018"),]
pdat2019 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2019"),]
pdat2021 <- pdat_NEBS[which(pdat_NEBS$YEAR=="2021"),]

tv_1982_rsme <- sqrt(mean((pdat1982$logCPUE - pdat1982$predicted_time_varying)^2, na.rm=TRUE))
tv_1985_rsme <- sqrt(mean((pdat1985$logCPUE - pdat1985$predicted_time_varying)^2, na.rm=TRUE))
tv_1988_rsme <- sqrt(mean((pdat1988$logCPUE - pdat1988$predicted_time_varying)^2, na.rm=TRUE))
tv_1991_rsme <- sqrt(mean((pdat1991$logCPUE - pdat1991$predicted_time_varying)^2, na.rm=TRUE))
tv_2010_rsme <- sqrt(mean((pdat2010$logCPUE - pdat2010$predicted_time_varying)^2, na.rm=TRUE))
tv_2017_rsme <- sqrt(mean((pdat2017$logCPUE - pdat2017$predicted_time_varying)^2, na.rm=TRUE))
tv_2018_rsme <- sqrt(mean((pdat2018$logCPUE - pdat2018$predicted_time_varying)^2, na.rm=TRUE))
tv_2019_rsme <- sqrt(mean((pdat2019$logCPUE - pdat2019$predicted_time_varying)^2, na.rm=TRUE))
tv_2021_rsme <- sqrt(mean((pdat2021$logCPUE - pdat2021$predicted_time_varying)^2, na.rm=TRUE))

tv_1982_rsme
tv_1985_rsme
tv_1988_rsme
tv_1991_rsme
tv_2010_rsme
tv_2017_rsme
tv_2018_rsme 
tv_2019_rsme
tv_2021_rsme

#plot----
#pivot longer so that can plot on same scale!

plot_pred_dat <- pdat_NEBS[,c(1:3, 5,  19, 24:29)] %>% pivot_longer(!c(LATITUDE, LONGITUDE, STATION, YEAR,  region, period, shelf), 
                                                                    names_to="response_type", values_to="value")



names(plot_pred_dat)
table(plot_pred_dat$response_type) #now four coloumns, difference, logCPUE... and predicted, and predicted_time_varying


#one-to-one plots========

ggplot(pdat_NEBS, aes(predicted_time_varying, logCPUE)) + geom_point() + geom_smooth(method="lm") + geom_abline(intercept=0)

#NEW
ggplot(pdat_NEBS, aes(predicted_time_varying, logCPUE, col=as.factor(YEAR))) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~YEAR)

#VS OLD
ggplot(pdat_NEBS, aes(predicted, logCPUE, col=as.factor(YEAR))) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~YEAR)
# scale_color_manual(values=c("#b2df8a", "#66c2a5", "#fc8d62", "#8da0cb"))

#adding in random effects======================================

#following examples from gamm() documentation

refatv <- ranef(tvmod2_REML$lme,level=5) #extract random effect for year
#updated here, seems to be different level?
rownames(refatv) <- substr(rownames(refatv),start=9,stop=12) #rename names

## make a prediction, with random effects zero...
#p0 <- predict(b$gam,data.frame(x0=.3,x1=.6,x2=.98,x3=.77))
p0_21tv <- predict(tvmod2_REML$gam, newdata = nebs_sel[which(nebs_sel$YEAR=="2021"),])

## add in effect for fa = "2" and fb="2/4"...
ptv <- p0_21tv + refatv["2021",1] 

refatv$year <- rownames(refatv)

#and loop through the years that we have nebs data for
yrs <- unique(nebs_sel$YEAR)
output_dftv <- data.frame(matrix(ncol = 11, nrow = 0))
nms <- c("ptemp_wrand", "year", "predicted_not_conditional",
         "long_albers", "lat_albers", "logCPUE", "BOT_DEPTH", "BOT_TEMP",
         "julian", "region", "shelf")
colnames(output_dftv) <- nms
i <- 1
for(i in 1:length(yrs)){
  temp_yr <- yrs[i]
  yr_dat <- nebs_sel[which(nebs_sel$YEAR==temp_yr),]
  
  ## make a prediction, with random effects zero...
  ptemp <- predict(tvmod2_REML$gam, newdata = yr_dat)
  
  ## add in effect for fa = "2" and fb="2/4"...
  ptemp_wrand <- ptemp + refatv[which(refatv$year==temp_yr),1] 
  
  dftemp <- as.data.frame(ptemp_wrand)
  dftemp$year <- temp_yr
  dftemp$predicted_not_conditional <- ptemp
  dftemp$long_albers <- yr_dat$long_albers
  dftemp$lat_albers <- yr_dat$lat_albers
  dftemp$logCPUE <- yr_dat$logCPUE
  dftemp$BOT_DEPTH <- yr_dat$BOT_DEPTH
  dftemp$BOT_TEMP <- yr_dat$BOT_TEMP
  dftemp$julian <- yr_dat$julian
  dftemp$region <- yr_dat$region
  dftemp$shelf <- yr_dat$shelf
  
  output_dftv <- rbind(output_dftv, dftemp)
}

output_dftv <- rename(output_dftv, predicted_cond_on_random = ptemp_wrand)

#ok nice the loop output seems to be working
#need to double check
#then need to add in metadata to loop df so that it can be used to plot etc

#look at 1:1 plot----
#not conditional
ggplot(output_dftv, aes(predicted_not_conditional, logCPUE, col=as.factor(year))) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~year)

#conditional on random year
ggplot(output_dftv, aes(predicted_cond_on_random, logCPUE, col=as.factor(year))) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~year)

#manuscript
ggplot(output_dftv, aes(predicted_cond_on_random, logCPUE)) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~year)

#manuscript WITH correlation (calculated in section 'try to calc cor')
# ggplot(output_dftv, aes(predicted_cond_on_random, logCPUE)) + geom_point() + 
#   geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
#   facet_wrap(~year) + #geom_text(aes(res_cor_vec), data=results_cor_output)
#   geom_label(data = results_cor_outputtv, aes(x = 0.25, y = 5.5, 
#                                             label = round(res_cor_vec, digits=3),
#                                             label.size = NA))

#above isn't working so maybe join to other database
plotdftv <- left_join(output_dftv, results_cor_outputtv, by=c("year"="res_yr_vec"))

ggplot(plotdftv, aes(predicted_cond_on_random, logCPUE)) + geom_point() + 
  geom_smooth(method="lm") + geom_abline(intercept=0) + theme_bw() + ylab("log(CPUE+1)") + xlab("Predicted log(CPUE+1)") +
  facet_wrap(~year) + #geom_text(aes(res_cor_vec), data=results_cor_outputtv)
  geom_label(aes(x = 0.25, y = 5.5,   label = round(res_cor_vec, digits=3)))

#plot----
#pivot longer so that can plot on same scale!

plot_out_dattv <- output_dftv[,-c(7:9)] %>% pivot_longer(!c(lat_albers, long_albers, year,  region,  shelf), 
                                                     names_to="response_type", values_to="value")


names(plot_out_dattv)
table(plot_out_dattv$response_type) #three coloumns, difference, logCPUE... and predicted

#match lat/long albers to lat long to plot
lats_longs <- pdat_NEBS[,c(1,2,22,23)]

pred_map_plot_dattv <- left_join(plot_out_dattv, lats_longs)
length(plot_out_dattv$long_albers)
length(pred_map_plot_dattv$long_albers) #same length

#plot----
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dattv) + 
  facet_wrap(response_type~year, nrow=3)  +
  scale_colour_distiller(palette = "Spectral") + theme_bw() +
  theme( legend.position = c(0.97, 0.25), legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank()) 


#try two sets of columns
library(cowplots)

#NEED TO SET TO SAME SCALE IF PLOTTING LIKE THIS

pred_map_plot_dattv$response_type2 <- pred_map_plot_dattv$response_type
# 
pred_map_plot_dattv$response_type2[which(pred_map_plot_dattv$response_type=="logCPUE")] <- "Actual"
pred_map_plot_dattv$response_type2[which(pred_map_plot_dattv$response_type=="predicted_cond_on_random")] <- "Predicted"


p1 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dattv[which(pred_map_plot_dattv$response_type2!="predicted_not_conditional"&
                                            pred_map_plot_dattv$year<2000 ),]) + 
  facet_grid(year~response_type2)  +
  scale_colour_distiller(palette = "Spectral", 
                         limits = c(min(pred_map_plot_dattv$value, na.rm=TRUE), max(pred_map_plot_dattv$value, na.rm=TRUE))) + theme_bw() +
  theme( legend.position = "right", legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank())
#switched this a bit for presentation
p2 <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dattv[which(pred_map_plot_dattv$response_type2!="predicted_not_conditional"&
                                            pred_map_plot_dattv$year>2000 ),]) + 
  facet_grid(response_type2~year)  +
  scale_colour_distiller(palette = "Spectral", 
                         limits = c(min(pred_map_plot_dattv$value, na.rm=TRUE), max(pred_map_plot_dattv$value, na.rm=TRUE))) + theme_bw() +
  theme( legend.position = "right", legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank())

plot_grid(p1, p2)

#can I reorder these?
ptest <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dattv[which(pred_map_plot_dattv$response_type2!="predicted_not_conditional"),]) + 
  facet_wrap(year~response_type2, ncol = 6)  +
  scale_colour_distiller(palette = "Spectral", 
                         limits = c(min(pred_map_plot_dattv$value, na.rm=TRUE), max(pred_map_plot_dattv$value, na.rm=TRUE))) + theme_bw() +
  theme( legend.position = "right", legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank())
ptest

ptest <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(58, 66), expand = TRUE) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE,LATITUDE,  col=value), 
             data=pred_map_plot_dattv[which(pred_map_plot_dattv$response_type2!="predicted_not_conditional"),]) + 
  facet_grid(vars(year), vars(response_type2))  +
  scale_colour_distiller(palette = "Spectral", 
                         limits = c(min(pred_map_plot_dattv$value, na.rm=TRUE), max(pred_map_plot_dattv$value, na.rm=TRUE))) + theme_bw() +
  theme( legend.position = "right", legend.key = element_blank(),
         legend.background=element_blank(), legend.title = element_blank(),
         axis.text = element_text(size=8)) 
ptest

ggsave(filename = "figs/2022_03_31_predict_small_text.png", ptest,
               width=10, height=13, dpi=600, units="in")

#get RMSE----
cond_mod_rsmetv <- sqrt(mean((output_dftv$logCPUE - output_dftv$predicted_cond_on_random)^2, na.rm=TRUE))

out1982tv <- output_dftv[which(output_dftv$year=="1982"),]
out1985tv <- output_dftv[which(output_dftv$year=="1985"),]
out1988tv <- output_dftv[which(output_dftv$year=="1988"),]
out1991tv <- output_dftv[which(output_dftv$year=="1991"),]

out2010tv <- output_dftv[which(output_dftv$year=="2010"),]
out2017tv <- output_dftv[which(output_dftv$year=="2017"),]
out2018tv <- output_dftv[which(output_dftv$year=="2018"),]
out2019tv <- output_dftv[which(output_dftv$year=="2019"),]
out2021tv <- output_dftv[which(output_dftv$year=="2021"),]

cond_1982_rsmetv <- sqrt(mean((out1982tv$logCPUE - out1982tv$predicted_cond_on_random)^2, na.rm=TRUE))
cond_1985_rsmetv <- sqrt(mean((out1985tv$logCPUE - out1985tv$predicted_cond_on_random)^2, na.rm=TRUE))
cond_1988_rsmetv <- sqrt(mean((out1988tv$logCPUE - out1988tv$predicted_cond_on_random)^2, na.rm=TRUE))
cond_1991_rsmetv <- sqrt(mean((out1991tv$logCPUE - out1991tv$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2010_rsmetv <- sqrt(mean((out2010tv$logCPUE - out2010tv$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2017_rsmetv <- sqrt(mean((out2017tv$logCPUE - out2017tv$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2018_rsmetv <- sqrt(mean((out2018tv$logCPUE - out2018tv$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2019_rsmetv <- sqrt(mean((out2019tv$logCPUE - out2019tv$predicted_cond_on_random)^2, na.rm=TRUE))
cond_2021_rsmetv <- sqrt(mean((out2021tv$logCPUE - out2021tv$predicted_cond_on_random)^2, na.rm=TRUE))

cond_1982_rsmetv
cond_1985_rsmetv
cond_1988_rsmetv
cond_1991_rsmetv
cond_2010_rsmetv
cond_2017_rsmetv
cond_2018_rsmetv
cond_2019_rsmetv
cond_2021_rsmetv


#try to calc cor-------
#from other script
#tempmat <- output_df[,c(1,6)]

#tempcor <- cor(tempmat, method = "pearson", use = "complete.obs")

#copied below
yrs <- unique(output_dftv$year)
#loop through species, loop through each year, get correlation and store

res_yr_vec <- vector(mode="numeric", length=length(yrs))
res_sps_vec <- vector(mode="numeric", length=length(yrs))
res_cor_vec <- vector(mode="numeric", length=length(yrs))

#yrs <- unique(full_wide$YEAR)
i <- 1
counter <- 1
for(i in 1:length(yrs)){
  temp_yr <- yrs[i]
  output_slice <- output_dftv[which(output_dftv$year==temp_yr),]
  tempmat <- output_slice[,c(1,6)]
  
  tempcor <- cor(tempmat, method = "pearson", use = "complete.obs")
  
  
  #save output
  #for(l in 1:10){
  res_yr_vec[i] <- temp_yr      
  #res_sps_vec[counter] <- rownames(tempcor)[l]
  res_cor_vec[i] <- tempcor[1,2]
  
  # counter <- counter + 1
  #}
}

results_cor_outputtv <- data.frame(res_yr_vec, #res_sps_vec, 
                                 res_cor_vec)




