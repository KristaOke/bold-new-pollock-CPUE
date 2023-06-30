#====================================================================================================================
#Looping through to test different windows and breakpoints
#Again with simple models and SEBS only

#by Krista, July-2022
#====================================================================================================================
#Notes: editted Nov 2020 to update models to te
#====================================================================================================================

#read in periods_analysis_data from trawl_biomass_GAM_explor but rename
wd <- getwd()
#loops_dat <- read.csv(file=paste(wd,"/data/processed_periods_analysis_data.csv", sep=""), row.names = 1)
names(loops_dat)
loops_dat <- read.csv( file=paste(wd,"/data/survey data/south_1982-2021_bot_trawl_data_for_sebs_models.csv", sep=""))


yrs <- unique(loops_dat$YEAR)
n_yrs <- length(yrs)
possible_breaks <- yrs[which(yrs>(1982+5) & yrs<(2021-5))]


#starting with breakpoint
#will try all possible breakpoints with at least 6 yrs data on either side 
#(b/c don't want a smaller # yrs in either levels than the 2014 breakpoint)

breakpoint_vec <- vector(mode="numeric", length=length(possible_breaks))
AIC_split_vec <- vector(mode="numeric", length=length(possible_breaks))
AIC_drop_vec <- vector(mode="numeric", length=length(possible_breaks))

loops_dat <- loops_dat[which(loops_dat$BOT_DEPTH<180),]

i<-1
for(i in 1:length(possible_breaks)){
  #set up data with new breakpoint
  print(possible_breaks[i])
  temp_break <- possible_breaks[i]
  temp_dat <- loops_dat
  temp_dat$split <- NA
  temp_dat$split[which(temp_dat$YEAR>(temp_break-1))] <- "post"
  temp_dat$split[which(temp_dat$YEAR<(temp_break))] <- "pre"
  
  temp_dat$split <- as.factor(temp_dat$split)
  
  #run model
  # temp_mod <- gamm(logCPUE_Gadus_chalcogrammus ~ bottemp_anom:split + te(mean_station_bottemp, BOT_DEPTH, k=29),
  #                  random=list(YEAR_factor=~1), 
  #                  correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE),
  #                  data=temp_dat)
  
  temp_mod <- gamm(logCPUE ~  s(BOT_DEPTH) +
                         s(BOT_TEMP, by=split, bs="fs"),  random=list(YEAR_factor=~1), 
                       correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), REML=TRUE,
                       data=temp_dat)
  
  #save output
  breakpoint_vec[i] <- temp_break      
  AIC_split_vec[i] <- AIC(temp_mod[[1]])
  saveRDS(temp_mod, file=paste(wd,"/data/breakpoint_model_simple", temp_break, ".csv", sep=""))
  
  
  #run drop model
  # temp_drop <- gamm(logCPUE_Gadus_chalcogrammus ~ s(bottemp_anom, k=4) + te(mean_station_bottemp, BOT_DEPTH, k=29) ,
  #                   random=list(YEAR_factor=~1), 
  #                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE),
  #                   data=temp_dat)
  temp_drop <-   gamm(logCPUE ~  s(BOT_DEPTH) +
         s(BOT_TEMP),  random=list(YEAR_factor=~1), 
       correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), REML=TRUE,
       data=temp_dat)
  
  #save output
  AIC_drop_vec[i] <- AIC(temp_drop[[1]])
  saveRDS(temp_drop, file=paste(wd,"/data/breakpoint_drop_model_simple", temp_break, ".csv", sep=""))
  
}

outputdf <- as.data.frame(cbind(breakpoint_vec, AIC_split_vec, AIC_drop_vec))

outputdf$difference_AIC <- outputdf$AIC_split_vec - outputdf$AIC_drop_vec

write_csv(outputdf, file=paste(wd,"/data/breakpoint_loop1_output_timevarying",  ".csv", sep=""))
#outputdf <- read_csv(file=paste(wd,"/data/breakpoint_loop1_output",  ".csv", sep=""))


#next cut off data six yrs after breakpoint and try again



breakpoint_six_vec <- vector(mode="numeric", length=length(possible_breaks))
AIC_split_six_vec <- vector(mode="numeric", length=length(possible_breaks))
AIC_drop_six_vec <- vector(mode="numeric", length=length(possible_breaks))


i<-1
for(i in 1:length(possible_breaks)){
  #set up data with new breakpoint
  print(possible_breaks[i])
  temp_break <- possible_breaks[i]
  temp_dat <- loops_dat[which(loops_dat$YEAR<(temp_break+5)),]
  temp_dat$split <- NA
  temp_dat$split[which(temp_dat$YEAR>(temp_break-1))] <- "post"
  temp_dat$split[which(temp_dat$YEAR<(temp_break))] <- "pre"
  
  temp_dat$split <- as.factor(temp_dat$split)
  
  #run model
  temp_mod <- gamm(logCPUE ~  s(BOT_DEPTH) +
                     s(BOT_TEMP, by=split, bs="fs"),  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), REML=TRUE,
                   data=temp_dat)
  
  #save output
  breakpoint_six_vec[i] <- temp_break      
  AIC_split_six_vec[i] <- AIC(temp_mod[[1]])
  saveRDS(temp_mod, file=paste(wd,"/data/breakpoint_six_model", temp_break, ".csv", sep=""))
  
  
  #run drop model
  temp_drop <-  gamm(logCPUE ~  s(BOT_DEPTH) +
         s(BOT_TEMP),  random=list(YEAR_factor=~1), 
       correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), REML=TRUE,
       data=temp_dat)
  
  #save output
  AIC_drop_six_vec[i] <- AIC(temp_drop[[1]])
  saveRDS(temp_drop, file=paste(wd,"/data/breakpoint_six_drop_model", temp_break, ".csv", sep=""))
  
}

outputdf_six <- as.data.frame(cbind(breakpoint_six_vec, AIC_split_six_vec, AIC_drop_six_vec))

outputdf_six$difference_AIC <- outputdf_six$AIC_split_six_vec - outputdf_six$AIC_drop_six_vec

write_csv(outputdf_six, file=paste(wd,"/data/breakpoint_loop2_output_timevarying",  ".csv", sep=""))
#outputdf_six <- read_csv(file=paste(wd,"/data/breakpoint_loop2_output",  ".csv", sep=""))

AIC(tmod1E.1L[[1]]) - AIC(tmod1E.1Ldrop[[1]]) #

#let's look at models with significant period
mod05 <- readRDS(paste(wd,"/data/breakpoint_six_model", "2005", ".csv", sep=""))
plot_model(mod05[[2]], type="int") #


mod06 <- readRDS(paste(wd,"/data/breakpoint_six_model", "2006", ".csv", sep=""))
plot_model(mod06[[2]], type="int") #


mod07 <- readRDS(paste(wd,"/data/breakpoint_six_model", "2007", ".csv", sep=""))
plot_model(mod07[[2]], type="int") #


mod08 <- readRDS(paste(wd,"/data/breakpoint_six_model", "2008", ".csv", sep=""))
plot_model(mod08[[2]], type="int") #


mod09 <- readRDS(paste(wd,"/data/breakpoint_six_model", "2009", ".csv", sep=""))
plot_model(mod09[[2]], type="int") #


mod13 <- readRDS(paste(wd,"/data/breakpoint_six_model", "2013", ".csv", sep=""))
plot_model(mod13[[2]], type="int") #

#smaller AIC differences

mod94 <- readRDS(paste(wd,"/data/breakpoint_six_model", "1994", ".csv", sep=""))
plot_model(mod94[[2]], type="int") #

mod95 <- readRDS(paste(wd,"/data/breakpoint_six_model", "1995", ".csv", sep=""))
plot_model(mod95[[2]], type="int") #

mod97 <- readRDS(paste(wd,"/data/breakpoint_six_model", "1997", ".csv", sep=""))
plot_model(mod97[[2]], type="int") #

mod04 <- readRDS(paste(wd,"/data/breakpoint_six_model", "2004", ".csv", sep=""))
plot_model(mod04[[2]], type="int") #


ggplot(loops_dat, aes(YEAR, bottemp_anom)) + geom_line(aes(col=STATION)) + theme(legend.position = "none") +
  geom_smooth()




#now with a 6 yr window

#next cut off data six yrs after breakpoint and try again



breakpoint_window_vec <- vector(mode="numeric", length=length(possible_breaks))
AIC_split_window_vec <- vector(mode="numeric", length=length(possible_breaks))
AIC_drop_window_vec <- vector(mode="numeric", length=length(possible_breaks))


i<-1
for(i in 1:length(possible_breaks)){
  #set up data with new breakpoint
  print(possible_breaks[i])
  temp_break <- possible_breaks[i]
  temp_dat <- loops_dat[which(loops_dat$YEAR<(temp_break+5)),]
  temp_dat$split <- NA
  temp_dat$split[which(temp_dat$YEAR>((temp_break-1)+6))] <- "post"
  temp_dat$split[which(temp_dat$YEAR<(temp_break+6) & temp_dat$YEAR>(temp_break-1))] <- "window"
  temp_dat$split[which(temp_dat$YEAR<(temp_break))] <- "pre"
  
  temp_dat$split <- as.factor(temp_dat$split)
  
  #run model
  temp_mod <- gamm(logCPUE_Gadus_chalcogrammus ~ bottemp_anom:split + te(mean_station_bottemp, BOT_DEPTH, k=29),
                   random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE),
                   data=temp_dat)
  
  #save output
  breakpoint_window_vec[i] <- temp_break      
  AIC_split_window_vec[i] <- AIC(temp_mod[[1]])
  saveRDS(temp_mod, file=paste(wd,"/data/breakpoint_window_model", temp_break, ".csv", sep=""))
  
  
  #run drop model
  temp_drop <- gamm(logCPUE_Gadus_chalcogrammus ~ s(bottemp_anom, k=4) + te(mean_station_bottemp, BOT_DEPTH, k=29),
                    random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE),
                    data=temp_dat)
  
  #save output
  AIC_drop_window_vec[i] <- AIC(temp_drop[[1]])
  saveRDS(temp_drop, file=paste(wd,"/data/breakpoint_window_drop_model", temp_break, ".csv", sep=""))
  
}
#stopped at 2009
outputdf_window <- as.data.frame(cbind(breakpoint_window_vec, AIC_split_window_vec, AIC_drop_window_vec))

outputdf_window$difference_AIC <- outputdf_window$AIC_split_window_vec - outputdf_window$AIC_drop_window_vec

