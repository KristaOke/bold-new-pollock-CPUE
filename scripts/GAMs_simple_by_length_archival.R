#stats on data binned by length
#REDONE by Krista June 2022

#using binmeta from join_CPUEdat_w_sizeCPUEdat.R
library(ggplot2)
library(sjPlot)
library(tidyverse)
library(mgcv)
library(visreg)
library(mgcViz)
library("rnaturalearth")
library("rnaturalearthdata")
library( "ggspatial" )
library("sf")

wd <- getwd()
binmeta_clean_anom<- read.csv(file=paste(wd,"/data/clean_binned_anom_data.csv", sep=""), row.names = 1)



binmeta_clean_anom$log_sum_WGTCPUE_LEN <- log(binmeta_clean_anom$sum_wgtCPUE_len + 1)


binmeta_clean_anom$STATION <- as.factor(binmeta_clean_anom$STATION)
binmeta_clean_anom$VESSEL <- as.factor(binmeta_clean_anom$VESSEL)
binmeta_clean_anom$CRUISE <- as.factor(binmeta_clean_anom$CRUISE)
binmeta_clean_anom$HAUL <- as.factor(binmeta_clean_anom$HAUL)
binmeta_clean_anom$bin <- as.factor(binmeta_clean_anom$bin)

binmeta_clean_anom <- binmeta_clean_anom[which( binmeta_clean_anom$STRATUM!="70" &
                                                  binmeta_clean_anom$STRATUM!="71" &
                                                  binmeta_clean_anom$STRATUM!="81" ),]

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  # scale_colour_gradient2(low="blue", high="red", guide="colorbar") + 
  geom_point(aes(LONGITUDE, LATITUDE, 
                 col=as.factor(STRATUM)), data=binmeta_clean_anom) + theme_bw() 

binmeta_clean_anom$period <- as.factor(binmeta_clean_anom$period)
binmeta_clean_anom <- binmeta_clean_anom[which(binmeta_clean_anom$BOT_DEPTH<180),]

##exclusion criteria===============

#how to exclude stations w few obs? Maybe n in a single bin? yeah do once split out


bin1dat <- binmeta_clean_anom[which(binmeta_clean_anom$bin=="0-200"),]
bin2dat <- binmeta_clean_anom[which(binmeta_clean_anom$bin=="200-300"),]
bin3dat <- binmeta_clean_anom[which(binmeta_clean_anom$bin=="300-400"),]
bin4dat <- binmeta_clean_anom[which(binmeta_clean_anom$bin=="400-500"),]
bin5dat <- binmeta_clean_anom[which(binmeta_clean_anom$bin=="500+"),]

#exclude by bin

#bin 1
station_bin1 <- bin1dat %>% group_by(STATION) %>%
  summarize(n_yrs=n()) #

join1stat <- left_join(bin1dat, station_bin1)

bin1_analysis_dat <- join1stat[which(join1stat$n_yrs>5),] #none even close to 5


#bin 2
station_bin2 <- bin2dat %>% group_by(STATION) %>%
  summarize(n_yrs=n()) #

join2stat <- left_join(bin2dat, station_bin2)

bin2_analysis_dat <- join2stat[which(join2stat$n_yrs>5),] #none even close to 5


#bin 3
station_bin3 <- bin3dat %>% group_by(STATION) %>%
  summarize(n_yrs=n()) #

join3stat <- left_join(bin3dat, station_bin3)

bin3_analysis_dat <- join3stat[which(join3stat$n_yrs>5),] #none even close to 5


#bin 4
station_bin4 <- bin4dat %>% group_by(STATION) %>%
  summarize(n_yrs=n()) #

join4stat <- left_join(bin4dat, station_bin4)

bin4_analysis_dat <- join4stat[which(join4stat$n_yrs>5),] #none even close to 5


#bin 5
station_bin5 <- bin5dat %>% group_by(STATION) %>%
  summarize(n_yrs=n()) #

join5stat <- left_join(bin5dat, station_bin5)

bin5_analysis_dat <- join5stat[which(join5stat$n_yrs>5),] #none even close to 5





#models by bin======

#bin 1======
smint_bin1_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                  s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                data=bin1_analysis_dat, REML=FALSE)
saveRDS(smint_bin1_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin1_sebs_model180.RDS")
smint_bin1_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin1_sebs_model180.RDS")
gam.check(smint_bin1_ML[[2]]) #oh no
summary(smint_bin1_ML[[1]]) #  
summary(smint_bin1_ML[[2]]) #
smint1_aic <- AIC(smint_bin1_ML$lme)

gammodsmintbin1 <- smint_bin1_ML$gam

draw(gammodsmintbin1, select = 1)
draw(gammodsmintbin1, select = 1, dist=0.05)
draw(gammodsmintbin1, select = 1, dist=0.01)

draw(gammodsmintbin1, select = 2)
draw(gammodsmintbin1, select = 3)

appraise(gammodsmintbin1)

anova(smint_bin1_ML[[2]])
plot(smint_bin1_ML[[2]])


binviz1 <- getViz(smint_bin1_ML$gam)
plot(sm(binviz1 , 1))
plot(sm(binviz1 , 2))
plot(sm(binviz1 , 3))

library(sjPlot)
library(sjmisc)

# plot_model(modbin1[[2]], type="int") #conditioned on fixed effects
# plot_model(modbin1[[2]], type="int", pred.type = "re") #conditioned on random effects
# plot_model(modbin1[[2]], type="int", pred.type = "re",
#            show.data = TRUE) #conditioned on random effects
# plot_model(modbin1[[2]], type="int", pred.type = "re",
#            show.values = TRUE) 
# plot_model(modbin1[[2]], type="resid")
# 
# visreg(modbin1$gam, "bottemp_anom", by="period", data=bin1_analysis_dat,
#        overlay=TRUE, partial=FALSE, rug=FALSE, ylim=c(-0.5,5))

#bin 1 no interaction
dropbin1 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin1_analysis_dat, REML=FALSE)
saveRDS(dropbin1, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin1_model180.RDS")
dropbin1 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin1_model180.RDS")

dropbin1REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin1_analysis_dat, REML=TRUE)
saveRDS(dropbin1REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin1_model180REML.RDS")
dropbin1REML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin1_model180REML.RDS")

summary(dropbin1REML$gam)

gam.check(dropbin1[[2]]) 
summary(dropbin1[[1]]) #  
summary(dropbin1[[2]]) #rsq 0.
smooth1_aic <- AIC(dropbin1$lme)

gamdropbin1 <- dropbin1$gam

draw(gamdropbin1, select = 1)
draw(gamdropbin1, select = 1, dist=0.05)
draw(gamdropbin1, select = 1, dist=0.01)
draw(gamdropbin1, select = 2)

appraise(dropdropbin1$gam)

anova(dropdropbin1[[2]])
plot(dropdropbin1[[2]])



#bin 1 linear no interaction
droplinbin1 <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin1_analysis_dat, REML=FALSE)
saveRDS(droplinbin1, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin1_model180.RDS")
droplinbin1 <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin1_model180.RDS")
gam.check(droplinbin1[[2]]) 
summary(droplinbin1[[1]]) #   
summary(droplinbin1[[2]]) #rsq 0.
null1_aic <- AIC(droplinbin1$lme)

droplinbin1REML <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin1_analysis_dat, REML=TRUE)
saveRDS(droplinbin1REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin1_model180REML.RDS")
summary(droplinbin1REML$gam)

gamdroplinbin1 <- droplinbin1$gam

draw(gamdroplinbin1, select = 1)
draw(gamdroplinbin1, select = 1, dist=0.05)
draw(gamdroplinbin1, select = 1, dist=0.01)
draw(gamdroplinbin1, select = 2)

appraise(droplinbin1$gam)

anova(droplinbin1[[2]])
plot(droplinbin1[[2]])

linintbin1 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                  data=bin1_analysis_dat, REML=FALSE)
saveRDS(linintbin1, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin1_model180.RDS")
linintbin1 <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin1_model180.RDS")
gam.check(linintbin1[[2]]) #quite bad!
summary(linintbin1[[1]]) #   
summary(linintbin1[[2]]) #rsq 0.0254
linint1_aic <- AIC(linintbin1$lme)

gamlinintbin1 <- linintbin1$gam

linintbin1REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin1_analysis_dat, REML=TRUE)
saveRDS(linintbin1REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin1_model180REML.RDS")
summary(linintbin1REML$gam)
gam.check(linintbin1REML$gam) #also bad

draw(gamlinintbin1, select = 1)
draw(gamlinintbin1, select = 1, dist=0.05)
draw(gamlinintbin1, select = 1, dist=0.01)
draw(gamlinintbin1, select = 2)

smint_bin1_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin1_analysis_dat, REML=FALSE)
saveRDS(smint_bin1_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin1_sebs_model180.RDS")
smint_bin1_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin1_sebs_model180.RDS")
gam.check(smint_bin1_ML[[2]]) #oh no
summary(smint_bin1_ML[[1]]) #  
summary(smint_bin1_ML[[2]]) #R2 = 0.05
smint1_aic <- AIC(smint_bin1_ML$lme)
plot(smint_bin1_ML$gam)

smint_bin1_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin1_analysis_dat, REML=TRUE)
saveRDS(smint_bin1_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin1_sebs_model180REML.RDS")
summary(smint_bin1_REML$gam) #R2 = 0.05

tv_bin1_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                     t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                        data=bin1_analysis_dat, REML=FALSE)
saveRDS(tv_bin1_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin1_sebs_model180ML.RDS")
summary(tv_bin1_ML$gam) #R2 = 0
gam.check(tv_bin1_ML[[2]]) #not great
summary(tv_bin1_ML[[1]]) #  
summary(tv_bin1_ML[[2]]) #R2 = 0.06
smtv1_aic <- AIC(tv_bin1_ML$lme)
plot(tv_bin1_ML$gam)

tv_bin1_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                       t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin1_analysis_dat, REML=TRUE)
saveRDS(tv_bin1_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin1_sebs_model180REML.RDS")
tv_bin1_REML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin1_sebs_model180REML.RDS")
summary(tv_bin1_REML$gam) #R2 = 0
plot(tv_bin1_REML$gam)
draw(tv_bin1_REML$gam)


#delta aic bin 1

#REDO HERE
print(c(null1_aic, smooth1_aic, linint1_aic, smint_1_aic, smtv1_aic ))
AIC(smint_bin1_ML$lme, droplinbin1$lme,  dropbin1$lme, linintbin1$lme, tv_bin1_ML$lme) #dropbin1 is still best, that's single smooth


#null is lowest
 delta_dropbin1 <- 0
# delta_smooth1 <- smooth1_aic - smooth1_aic
 delta_linint1 <- linint1_aic - smooth1_aic
 delta_smoothint1 <- smint1_aic - smooth1_aic
 delta_tv1 <- smtv1_aic - smooth1_aic


#bin 1 akaike weights

 sumbin1aic <- sum(exp(-0.5*delta_smoothint1), exp(-0.5*delta_dropbin1), exp(-0.5*delta_linint1), exp(-0.5*delta_tv1)) 
# 

 aw_smooth1 <- exp(-0.5*delta_dropbin1)/sumbin1aic
 aw_smint1 <- exp(-0.5*delta_smoothint1)/sumbin1aic
 aw_linint1 <- exp(-0.5*delta_linint1)/sumbin1aic
 aw_tv1 <- exp(-0.5*delta_tv1)/sumbin1aic


#bin 2----

smint_bin2_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin2_analysis_dat, REML=FALSE)
saveRDS(smint_bin2_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin2_sebs_model180.RDS")
smint_bin2_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin2_sebs_model180.RDS")
gam.check(smint_bin2_ML[[2]]) #oh no
summary(smint_bin2_ML[[1]]) #  
summary(smint_bin2_ML[[2]]) #
smint2_aic <- AIC(smint_bin2_ML$lme)
plot(smint_bin2_ML$gam)

smint_bin2_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin2_analysis_dat, REML=TRUE)
saveRDS(smint_bin2_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin2_sebs_model180REML.RDS")
summary(smint_bin2_REML$gam)


#bin 2 no interaction
dropbin2 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin2_analysis_dat, REML=FALSE)
saveRDS(dropbin2, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin2_model180.RDS")
dropbin2 <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin2_model.RDS")
gam.check(dropbin2[[2]]) 
summary(dropbin2[[1]]) #  
summary(dropbin2[[2]]) #rsq 0.
smooth2_aic <- AIC(dropbin2$lme)
plot(dropbin2$gam)

dropbin2REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin2_analysis_dat, REML=TRUE)
saveRDS(dropbin2REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin2_model180REML.RDS")
summary(dropbin2REML$gam)


#bin 2 linear no interaction
droplinbin2 <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin2_analysis_dat, REML=FALSE)
saveRDS(droplinbin2, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin2_model180.RDS")
droplinbin2 <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin2_model.RDS")
gam.check(droplinbin2[[2]]) 
summary(droplinbin2[[1]]) #   
summary(droplinbin2[[2]]) #rsq 0.
null2_aic <- AIC(droplinbin2$lme)

droplinbin2REML <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin2_analysis_dat, REML=TRUE)
saveRDS(droplinbin2REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin2_model180REML.RDS")




linintbin2 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin2_analysis_dat, REML=FALSE)
saveRDS(linintbin2, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin2_model180.RDS")
linintbin2 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin2_model.RDS")
gam.check(linintbin2[[2]]) 
summary(linintbin2[[1]]) #   did not converge
summary(linintbin2[[2]]) #rsq 0.
linint2_aic <- AIC(linintbin2$lme)

linintbin2REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin2_analysis_dat, REML=TRUE)
saveRDS(linintbin2REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin2_model180REML.RDS")
summary(linintbin2REML$gam)


tv_bin2_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                     t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin2_analysis_dat, REML=FALSE)
saveRDS(tv_bin2_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin2_sebs_model180ML.RDS")
summary(tv_bin2_ML$gam) #R2 = 0.149
gam.check(tv_bin2_ML[[2]]) #not great
summary(tv_bin2_ML[[1]]) #  
summary(tv_bin2_ML[[2]]) #
smtv2_aic <- AIC(tv_bin2_ML$lme)
plot(tv_bin2_ML$gam)

tv_bin2_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                       t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                     correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                     data=bin2_analysis_dat, REML=TRUE)
saveRDS(tv_bin2_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin2_sebs_model180REML.RDS")
summary(tv_bin2_REML$gam) #R2 = 0.149
draw(tv_bin2_REML$gam)

#delta aic bin 2


print(c(null2_aic, smooth2_aic, linint2_aic, smint2_aic, smtv2_aic))
AIC(smint_bin2_ML$lme, droplinbin2$lme, linintbin2$lme, dropbin2$lme, tv_bin2_ML$lme) #tv and smooth int best
#they are close, maybe no int better after accounting for parameters?

#tv is lowest
 # delta_null2 <- 0
delta_tv2 <-  0
 #delta_null2 <- null2_aic - smint2_aic
 delta_linint2 <- linint2_aic - smtv2_aic
 delta_smooth2 <- smooth2_aic - smtv2_aic
 delta_smint2 <- smint2_aic - smtv2_aic



#bin 2 akaike weights

 sumbin2aic <- sum(exp(-0.5*delta_smoothint2), exp(-0.5*delta_smooth2), exp(-0.5*delta_linint2), exp(-0.5*delta_tv2)) 
# 
 aw_smoothint2 <- exp(-0.5*delta_smoothint2)/sumbin2aic
 aw_smooth2 <- exp(-0.5*delta_smooth2)/sumbin2aic
 aw_linint2 <- exp(-0.5*delta_linint2)/sumbin2aic
 aw_tv2 <- exp(-0.5*delta_tv2)/sumbin2aic






#bin 3----

smint_bin3_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin3_analysis_dat, REML=FALSE)
saveRDS(smint_bin3_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin3_sebs_model180.RDS")
smint_bin3_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin3_sebs_model180.RDS")
gam.check(smint_bin3_ML[[2]]) #oh no
summary(smint_bin3_ML[[1]]) #  
summary(smint_bin3_ML[[2]]) #
smint3_aic <- AIC(smint_bin3_ML$lme)

smint_bin3_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin3_analysis_dat, REML=TRUE)
saveRDS(smint_bin3_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin3_sebs_model180REML.RDS")
summary(smint_bin3_REML$gam)


#bin 3 no interaction
dropbin3 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin3_analysis_dat, REML=FALSE)
saveRDS(dropbin3, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin3_model180.RDS")
dropbin3 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin3_model.RDS")
gam.check(dropbin3[[2]]) 
summary(dropbin3[[1]]) #  
summary(dropbin3[[2]]) #rsq 0.
smooth3_aic <- AIC(dropbin3$lme)
plot(dropbin3$gam)

dropbin3REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin3_analysis_dat, REML=TRUE)
saveRDS(dropbin3REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin3_model180REML.RDS")
summary(dropbin3REML$gam)


#bin 3 linear no interaction
droplinbin3 <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin3_analysis_dat, REML=FALSE)
saveRDS(droplinbin3, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin3_model180.RDS")
droplinbin3 <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin3_model.RDS")
gam.check(droplinbin3[[2]]) 
summary(droplinbin3[[1]]) #   
summary(droplinbin3[[2]]) #rsq 0.
null3_aic <- AIC(droplinbin3$lme)

droplinbin3REML <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin3_analysis_dat, REML=TRUE)
saveRDS(droplinbin3REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin3_model180REML.RDS")



linintbin3 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin3_analysis_dat, REML=FALSE)
saveRDS(linintbin3, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin3_model180.RDS")
linintbin3 <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin3_model.RDS")
gam.check(linintbin3[[2]]) 
summary(linintbin3[[1]]) #   did not converge
summary(linintbin3[[2]]) #rsq 0.
linint3_aic <- AIC(linintbin3$lme)

linintbin3REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin3_analysis_dat, REML=TRUE)
saveRDS(linintbin3REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin3_model180REML.RDS")
summary(linintbin3REML$gam)


tv_bin3_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                     t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin3_analysis_dat, REML=FALSE)
saveRDS(tv_bin3_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin3_sebs_model180ML.RDS")
summary(tv_bin3_ML$gam) #R2 = 0.242
gam.check(tv_bin3_ML[[2]]) #not great
summary(tv_bin3_ML[[1]]) #  
summary(tv_bin3_ML[[2]]) #
smtv3_aic <- AIC(tv_bin3_ML$lme)
plot(tv_bin3_ML$gam)


tv_bin3_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                       t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                     correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                     data=bin3_analysis_dat, REML=TRUE)
saveRDS(tv_bin3_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin3_sebs_model180REML.RDS")
summary(tv_bin3_REML$gam) #R2 = 0.242
draw(tv_bin3_ML$gam)


#delta aic bin 3

#REDO HERE
print(c(null3_aic, smooth3_aic, linint3_aic, smint3_aic, smtv3_aic))
AIC(smint_bin3_ML$lme, droplinbin3$lme, linintbin3$lme, dropbin3$lme, tv_bin3_ML$lme) #tv best
#

#tv is lowest
# delta_null3 <- 0
delta_tv3 <- 0
 delta_smooth3 <- smooth3_aic - smtv3_aic
 delta_smoothint3  <- smint3_aic - smtv3_aic
 delta_linint3 <- linint3_aic - smtv3_aic
 
# 


#bin 3 akaike weights

 sumbin3aic <- sum(exp(-0.5*delta_smoothint3), exp(-0.5*delta_smooth3), exp(-0.5*delta_linint3), exp(-0.5*delta_tv3)) 
# 
 aw_smint3 <- exp(-0.5*delta_smoothint3)/sumbin3aic
 aw_smooth3 <- exp(-0.5*delta_smooth3)/sumbin3aic
 aw_linint3 <- exp(-0.5*delta_linint3)/sumbin3aic
 aw_tv3 <- exp(-0.5*delta_tv3)/sumbin3aic





#bin 4----

smint_bin4_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin4_analysis_dat, REML=FALSE)
saveRDS(smint_bin4_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin4_sebs_model180.RDS")
smint_bin4_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin4_sebs_model180.RDS")
gam.check(smint_bin4_ML[[2]]) #oh no
summary(smint_bin4_ML[[1]]) #  
summary(smint_bin4_ML[[2]]) #
smint4_aic <- AIC(smint_bin4_ML$lme)
plot(smint_bin4_ML$gam)

smint_bin4_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin4_analysis_dat, REML=TRUE)
saveRDS(smint_bin4_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin4_sebs_model180REML.RDS")
summary(smint_bin4_REML$gam)


#bin 4 no interaction
dropbin4 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin4_analysis_dat, REML=FALSE)
saveRDS(dropbin4, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin4_model180.RDS")
dropbin4 <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin4_model.RDS")
gam.check(dropbin4[[2]]) 
summary(dropbin4[[1]]) #  
summary(dropbin4[[2]]) #rsq 0.
smooth4_aic <- AIC(dropbin4$lme)
plot(dropbin4$gam)

dropbin4REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin4_analysis_dat, REML=TRUE)
saveRDS(dropbin4REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin4_model180REML.RDS")
summary(dropbin4REML$gam)

#bin 4 linear no interaction
droplinbin4 <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin4_analysis_dat, REML=FALSE)
saveRDS(droplinbin4, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin4_model180.RDS")
droplinbin4 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin4_model.RDS")
gam.check(droplinbin4[[2]]) 
summary(droplinbin4[[1]]) #   
summary(droplinbin4[[2]]) #rsq 0.
null4_aic <- AIC(droplinbin4$lme)

droplinbin4REML <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin4_analysis_dat, REML=TRUE)
saveRDS(droplinbin4REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin4_model180REML.RDS")
summary(droplinbin4REML$gam)


linintbin4 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin4_analysis_dat, REML=FALSE)
saveRDS(linintbin4, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin4_model180.RDS")
linintbin4 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin4_model.RDS")
gam.check(linintbin4[[2]]) 
summary(linintbin4[[1]]) #   did not converge
summary(linintbin4[[2]]) #rsq 0.
linint4_aic <- AIC(linintbin4$lme)

linintbin4REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin4_analysis_dat, REML=TRUE)
saveRDS(linintbin4REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin4_model180REML.RDS")
summary(linintbin4REML$gam)


tv_bin4_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                     t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin4_analysis_dat, REML=FALSE)
saveRDS(tv_bin4_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin4_sebs_model180ML.RDS")
summary(tv_bin4_ML$gam) #R2 = 0.399
gam.check(tv_bin4_ML[[2]]) #not as bad
summary(tv_bin4_ML[[1]]) #  
summary(tv_bin4_ML[[2]]) #
smtv4_aic <- AIC(tv_bin4_ML$lme)
plot(tv_bin4_ML$gam)

tv_bin4_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                       t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                     correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                     data=bin4_analysis_dat, REML=TRUE)
saveRDS(tv_bin4_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin4_sebs_model180REML.RDS")
summary(tv_bin4_REML$gam) #R2 = 0.399
draw(tv_bin4_REML$gam)


#delta aic bin 4

#REDO HERE
print(c(null4_aic, smooth4_aic, linint4_aic, smint4_aic, smtv4_aic))
AIC(smint_bin4_ML$lme, droplinbin4$lme, linintbin4$lme, dropbin4$lme, tv_bin4_ML$lme) #
#tv best

#tv is lowest
# delta_null4 <- 0
#delta_smoothint4 <- 0
delta_tv4 <- 0
 delta_smooth4 <- smooth4_aic - smtv4_aic
 delta_linint4 <- linint4_aic - smtv4_aic
 delta_smint4 <- smint4_aic - smtv4_aic



#bin 4 akaike weights

 sumbin4aic <- sum(exp(-0.5*delta_smoothint4), exp(-0.5*delta_smooth4), exp(-0.5*delta_linint4), exp(-0.5*delta_tv4)) 
# 
 aw_smint4 <- exp(-0.5*delta_smoothint4)/sumbin4aic
 aw_smooth4 <- exp(-0.5*delta_smooth4)/sumbin4aic
 aw_linint4 <- exp(-0.5*delta_linint4)/sumbin4aic
 aw_tv4 <- exp(-0.5*delta_tv4)/sumbin4aic





#bin 5----

smint_bin5_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin5_analysis_dat, REML=FALSE)
saveRDS(smint_bin5_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin5_sebs_model180.RDS")
smint_bin5_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin5_sebs_model180.RDS")
gam.check(smint_bin5_ML[[2]]) #not nearly as bad!
summary(smint_bin5_ML[[1]]) #  
summary(smint_bin5_ML[[2]]) #
smint5_aic <- AIC(smint_bin5_ML$lme)
plot(smint_bin5_ML$gam)

smint_bin5_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                      correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                      data=bin5_analysis_dat, REML=TRUE)
saveRDS(smint_bin5_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_bin5_sebs_model180REML.RDS")
summary(smint_bin5_REML$gam)


#bin 5 no interaction
dropbin5 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin5_analysis_dat, REML=FALSE)
saveRDS(dropbin5, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin5_model180.RDS")
dropbin5 <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin5_model.RDS")
gam.check(dropbin5[[2]]) 
summary(dropbin5[[1]]) #  
summary(dropbin5[[2]]) #rsq 0.
smooth5_aic <- AIC(dropbin5$lme)
plot(dropbin5$gam)

dropbin5REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                   s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                 data=bin5_analysis_dat, REML=TRUE)
saveRDS(dropbin5REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothtemp_bin5_model180REML.RDS")
summary(dropbin5REML$gam)


#bin 5 linear no interaction
droplinbin5 <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin5_analysis_dat, REML=FALSE)
saveRDS(droplinbin5, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin5_model180.RDS")
droplinbin5 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin5_model.RDS")
gam.check(droplinbin5[[2]]) 
summary(droplinbin5[[1]]) #   
summary(droplinbin5[[2]]) #rsq 0.
null5_aic <- AIC(droplinbin5$lme)

droplinbin5REML <- gamm(log_sum_WGTCPUE_LEN ~  s(BOT_DEPTH) +
                      BOT_TEMP,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                    data=bin5_analysis_dat, REML=TRUE)
saveRDS(droplinbin5REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/null_bin5_model180REML.RDS")



linintbin5 <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin5_analysis_dat, REML=FALSE)
saveRDS(linintbin5, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin5_model180.RDS")
linintbin5 <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin5_model.RDS")
gam.check(linintbin5[[2]]) 
summary(linintbin5[[1]]) #   did not converge
summary(linintbin5[[2]]) #rsq 0.
linint5_aic <- AIC(linintbin5$lme)

linintbin5REML <- gamm(log_sum_WGTCPUE_LEN ~s(BOT_DEPTH) +
                     BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin5_analysis_dat, REML=TRUE)
saveRDS(linintbin5REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smoothint_bin5_model180REML.RDS")
summary(linintbin5REML$gam)



tv_bin5_ML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                     t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                   correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                   data=bin5_analysis_dat, REML=FALSE)
saveRDS(tv_bin5_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin5_sebs_model180ML.RDS")
summary(tv_bin5_ML$gam) #R2 = 0.21
gam.check(tv_bin5_ML[[2]]) #not as bad
summary(tv_bin5_ML[[1]]) #  
summary(tv_bin5_ML[[2]]) #
smtv5_aic <- AIC(tv_bin5_ML$lme)
plot(tv_bin5_ML$gam)

tv_bin5_REML <- gamm(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                       t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                     correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), 
                     data=bin5_analysis_dat, REML=TRUE)
saveRDS(tv_bin5_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time-varying_bin5_sebs_model180REML.RDS")
summary(tv_bin5_REML$gam) #R2 = 0.21
plot(tv_bin5_REML$gam)
draw(tv_bin5_REML$gam)

#delta aic bin 5

#REDO HERE
print(c(null5_aic, smooth5_aic, linint5_aic, smint5_aic, smtv5_aic))
AIC(smint_bin5_ML$lme, droplinbin5$lme, linintbin5$lme, dropbin5$lme, tv_bin5_ML$lme) #both smoothers and tv better, close but
#the smooth w/o int prob better b/c lower # parameters

#smooth int is lowest
delta_smoothint5 <- 0
delta_smooth5 <- smooth5_aic - smint5_aic
 delta_linint5 <- linint5_aic - smint5_aic
 delta_tv5 <- smtv5_aic - smint5_aic



#bin 5 akaike weights

 sumbin5aic <- sum(exp(-0.5*delta_smoothint5), exp(-0.5*delta_smooth5), exp(-0.5*delta_linint5), exp(-0.5*delta_tv5)) 
# 
 aw_smint5 <- exp(-0.5*delta_smoothint5)/sumbin5aic
 aw_smooth5 <- exp(-0.5*delta_smooth5)/sumbin5aic
 aw_linint5 <- exp(-0.5*delta_linint5)/sumbin5aic
 aw_tv5 <- exp(-0.5*delta_tv5)/sumbin5aic






#plot for manu

par(mfrow=c(2,3), mai=c(0.3,0.4,0.3,0.1)) 

visreg(modbin1$gam, "bottemp_anom", by="period", data=bin1_analysis_dat,
       overlay=TRUE, partial=FALSE, rug=FALSE, ylim=c(-0.5,5), legend=FALSE, xlab="",
       line=list(col=(c("#0083c9", "red"))),
       fill=list(col=(c("#0083c980", "#FF4E3780"))))

visreg(modbin2$gam, "bottemp_anom", by="period", data=bin2_analysis_dat,
       overlay=TRUE, partial=FALSE, rug=FALSE, ylim=c(-0.5,5), legend=FALSE, xlab="",
       line=list(col=(c("#0083c9", "red"))),
       fill=list(col=(c("#0083c980", "#FF4E3780"))))

visreg(modbin3$gam, "bottemp_anom", by="period", data=bin3_analysis_dat,
       overlay=TRUE, partial=FALSE, rug=FALSE, ylim=c(-0.5,5), legend=FALSE, xlab="",
       line=list(col=(c("#0083c9", "red"))),
       fill=list(col=(c("#0083c980", "#FF4E3780"))))

visreg(modbin4$gam, "bottemp_anom", by="period", data=bin4_analysis_dat,
       overlay=TRUE, partial=FALSE, rug=FALSE, ylim=c(-0.5,5), legend=FALSE, xlab="",
       line=list(col=(c("#0083c9", "red"))),
       fill=list(col=(c("#0083c980", "#FF4E3780"))))

visreg(modbin5$gam, "bottemp_anom", by="period", data=bin5_analysis_dat,
       overlay=TRUE, partial=FALSE, rug=FALSE, ylim=c(-0.5,5), legend=FALSE, xlab="",
       line=list(col=(c("#0083c9", "red"))),
       fill=list(col=(c("#0083c980", "#FF4E3780"))))

# visreg(modbin5$gam, "bottemp_anom", by="period", data=bin5_analysis_dat,
#        overlay=TRUE, partial=FALSE, rug=FALSE, ylim=c(-0.5,5), legend=FALSE, xlab="",
#        line=list(col=(c("#0083c9", "red"))),
#        fill=list(col=(c("#0083c980", "#FF4E3780"))))



#plot again with time varying------

#

par(mfrow=c(2,3), mai=c(0.3,0.4,0.3,0.1)) 


s1 <- getViz(tv_bin1_REML$gam)
p1 <- plot(sm(s1, 2)) +
  labs(title = NULL) + #l_points() +
  scale_fill_distiller(palette = "Spectral", type = "div") +
   theme(legend.position = "none")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        #  axis.title.x = element_blank(),
        #  axis.title.y = element_blank(),
        plot.title = element_blank()) + xlab("Year") +
  ylab("Bottom temperature") #+ labs(fill = "Effect")

p1 

#replot bin 1 since smooth is better?

# plot_smooth(cmod1_noint_180$gam, view="BOT_TEMP",ylab="Partial effect on log(CPUE)", xlab="Bottom Temperature",
#             ylim=c(-1,5))
# 
# plot_smooth(ptest_REML$gam, view="BOT_TEMP", plot_all=c("period"), ylab="Partial effect on log(CPUE)", xlab="Bottom Temperature",
#             ylim=c(-1,5), col=c("blue", "red")) #finally works!

s1.1 <- getViz(dropbin1REML$gam)
p1.1 <- plot(sm(s1.1, 2)) +
  labs(title = NULL) + #l_points() +
  scale_fill_distiller(palette = "Spectral", type = "div") +
  theme(legend.position = "none")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        #  axis.title.x = element_blank(),
        #  axis.title.y = element_blank(),
        plot.title = element_blank()) + xlab("Bottom temperature") +
  ylab("Partial effect on log(CPUE)") #+ labs(fill = "Effect")

p1.1 



s2 <- getViz(tv_bin2_REML$gam)
p2 <- plot(sm(s2, 2)) +
  labs(title = NULL) + #l_points() +
  scale_fill_distiller(palette = "Spectral", type = "div") +
   theme(legend.position = "none")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        #  axis.title.x = element_blank(),
        #  axis.title.y = element_blank(),
        plot.title = element_blank()) + xlab("Year") +
  ylab("Bottom temperature") #+ labs(fill = "Effect")

p2 

s3 <- getViz(tv_bin3_REML$gam)
p3 <- plot(sm(s3, 2)) +
  labs(title = NULL) + #l_points() +
  scale_fill_distiller(palette = "Spectral", type = "div") +
   theme(legend.position = "none")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        #  axis.title.x = element_blank(),
        #  axis.title.y = element_blank(),
        plot.title = element_blank()) + xlab("Year") +
  ylab("Bottom temperature") #+ labs(fill = "Effect")

p3 

s4 <- getViz(tv_bin4_REML$gam)
p4 <- plot(sm(s4, 2)) +
  labs(title = NULL) + #l_points() +
  scale_fill_distiller(palette = "Spectral", type = "div") +
   theme(legend.position = "none")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        #  axis.title.x = element_blank(),
        #  axis.title.y = element_blank(),
        plot.title = element_blank()) + xlab("Year") +
  ylab("Bottom temperature") #+ labs(fill = "Effect")

p4 

s5 <- getViz(tv_bin5_REML$gam)
p5 <- plot(sm(s5, 2)) +
labs(title = NULL) + #l_points() +
  scale_fill_distiller(palette = "Spectral", type = "div") +
  # theme(legend.position = "none")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        #  axis.title.x = element_blank(),
        #  axis.title.y = element_blank(),
        plot.title = element_blank()) + xlab("Year") +
  ylab("Bottom temperature") + labs(fill = "Effect")

p5 



#bins modeled together---------------------------------------------------------------------

#spatial correlation doesn't work b/c same haul have distance zero,
#will try with haul in yr random instead
library(gamm4)


smint_allbins_ML <- gamm4(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        s(BOT_TEMP, period, by=bin, bs="fs"),  random=~(1|YEAR/HAUL),
                      data=binmeta_clean_anom, REML=FALSE)
#saveRDS(smint_allbins_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_allbins_sebs_model.RDS")
smint_allbins_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_allbins_sebs_model.RDS")
gam.check(smint_allbins_ML$gam ) #not great
summary(smint_allbins_ML$gam)
plot(smint_allbins_ML$gam)

sm_allbins_ML <- gamm4(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                           s(BOT_TEMP, by=bin, bs="fs"),  random=~(1|YEAR/HAUL),
                         data=binmeta_clean_anom, REML=FALSE)
#saveRDS(sm_allbins_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/sm_allbins_sebs_model.RDS")
sm_allbins_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/sm_allbins_sebs_model.RDS")
gam.check(sm_allbins_ML$gam ) #not great
summary(sm_allbins_ML$gam)
plot(sm_allbins_ML$gam)

linint_allbins_ML <- gamm4(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                        BOT_TEMP:period:bin,  random=~(1|YEAR/HAUL),
                      data=binmeta_clean_anom, REML=FALSE)
#saveRDS(linint_allbins_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/linint_allbins_sebs_model.RDS")
linint_allbins_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/linint_allbins_sebs_model.RDS")
gam.check(linint_allbins_ML$gam ) #not perfect!
summary(linint_allbins_ML$gam)
plot(linint_allbins_ML$gam)

lin_allbins_ML <- gamm4(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                            BOT_TEMP:bin,  random=~(1|YEAR/HAUL),
                          data=binmeta_clean_anom, REML=FALSE)
#saveRDS(lin_allbins_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/linint_allbins_sebs_model.RDS")
lin_allbins_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/linint_allbins_sebs_model.RDS")
gam.check(lin_allbins_ML$gam ) #not perfect!
summary(lin_allbins_ML$gam)
plot(lin_allbins_ML$gam)


AIC(smint_allbins_ML$mer, sm_allbins_ML$mer,
    linint_allbins_ML$mer, lin_allbins_ML$mer) #smooth int is far far lowest
#BUT did it work!?

#what about no bins?
smint_nobins_ML <- gamm4(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                            s(BOT_TEMP, period, bs="fs"),  random=~(1|YEAR/HAUL),
                          data=binmeta_clean_anom, REML=FALSE)
#saveRDS(smint_nobins_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smint_nobins_sebs_model.RDS")


smnoint_nobins_ML <- gamm4(log_sum_WGTCPUE_LEN ~ s(BOT_DEPTH) +
                           s(BOT_TEMP),  random=~(1|YEAR/HAUL),
                         data=binmeta_clean_anom, REML=FALSE)
#saveRDS(smnoint_nobins_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/smnoint_nobins_sebs_model.RDS")

AIC(smint_allbins_ML$mer, sm_allbins_ML$mer,
    linint_allbins_ML$mer, lin_allbins_ML$mer,
    smint_nobins_ML$mer, smnoint_nobins_ML$mer) #without bins is worse

#plot smoothers------
library(cowplot)

p1 <- draw(dropbin1$gam, select=2, gg=TRUE) + scale_y_continuous(limits = c(-3, 1)) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        title = element_blank())
p2 <- draw(dropbin2$gam, select=2, gg=TRUE) + scale_y_continuous(limits = c(-3, 1))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        title = element_blank())
p3 <- draw(dropbin3$gam, select=2, gg=TRUE) + scale_y_continuous(limits = c(-3, 1))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        title = element_blank())
p4 <- draw(dropbin4$gam, select=2, gg=TRUE) + scale_y_continuous(limits = c(-3, 1))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        title = element_blank())
p5 <- draw(dropbin5$gam, select=2, gg=TRUE) + scale_y_continuous(limits = c(-3, 1))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        title = element_blank())

plot_grid(p1, p2, p3, p4, p5,
          labels = c('0-20cm', '20-30cm', '30-40cm', '40-50cm', '50cm+'))  #are these 31, 41, or 30, 40?


