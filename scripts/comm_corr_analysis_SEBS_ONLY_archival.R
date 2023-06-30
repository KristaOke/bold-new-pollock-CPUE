#======================================================================================================================
#Community correlation analyses
#
#Krista, Jan 2021
#======================================================================================================================
#Notes:
#======================================================================================================================

#load packages
library(ggplot2)
library(tidyverse)

#load data
wd <- getwd()
northsouthdata_all <- read.csv(paste(wd,"/data/survey data/combined_cleaned_north-south_1982-2021_bot_trawl_data.csv",sep=""))

northsouthselect <- northsouthdata_all[which(northsouthdata_all$SCIENTIFIC=="Chionoecetes bairdi"|
                                               northsouthdata_all$SCIENTIFIC=="Atheresthes stomias"|
                                               northsouthdata_all$SCIENTIFIC=="Hippoglossus stenolepis"|
                                               northsouthdata_all$SCIENTIFIC=="Limanda aspera"|
                                               northsouthdata_all$SCIENTIFIC=="Lepidopsetta sp."|
                                               northsouthdata_all$SCIENTIFIC=="Chionoecetes opilio"|
                                               northsouthdata_all$SCIENTIFIC=="Gadus macrocephalus"|
                                               northsouthdata_all$SCIENTIFIC=="Hippoglossoides elassodon"|
                                               northsouthdata_all$SCIENTIFIC=="Pleuronectes quadrituberculatus"|
                                               northsouthdata_all$SCIENTIFIC=="Lepidopsetta polyxystra"|
                                               northsouthdata_all$SCIENTIFIC=="Gadus chalcogrammus"),]

wide_northsouth  <- northsouthselect[,-c(8:9, 10,12)] %>% pivot_wider(names_from=SCIENTIFIC, 
                                        values_from=logCPUE)
#eek this causes some bad column names!
wide_northsouth <- wide_northsouth %>% rename(logCPUE_Chionoecetes_bairdi = "Chionoecetes bairdi",            
                                    logCPUE_Atheresthes_stomias = "Atheresthes stomias",             
                                    logCPUE_Hippoglossus_stenolepis = "Hippoglossus stenolepis",
                                    logCPUE_Limanda_aspera = "Limanda aspera",         
                                    logCPUE_Lepidopsetta_sp= "Lepidopsetta sp.",               
                                    logCPUE_Chionoecetes_opilio = "Chionoecetes opilio",             
                                    logCPUE_Gadus_macrocephalus = "Gadus macrocephalus",            
                                    logCPUE_Hippoglossoides_elassodon = "Hippoglossoides elassodon",  
                                    logCPUE_Pleuronectes_quadrituberculatus = "Pleuronectes quadrituberculatus", 
                                    logCPUE_Lepidopsetta_polyxystra = "Lepidopsetta polyxystra", 
                                    logCPUE_Gadus_chalcogrammus = "Gadus chalcogrammus"
                                    ) 




wide_northsouth$period<- NA
wide_northsouth$period[which(wide_northsouth$YEAR<2014)]<-"early"
wide_northsouth$period[which(wide_northsouth$YEAR>2013)]<-"late"
wide_northsouth$period <- as.factor(wide_northsouth$period)


#SHould NAs be zeros?
wide_northsouth$logCPUE_Atheresthes_stomias[is.na(wide_northsouth$logCPUE_Atheresthes_stomias)==TRUE] <- 0
wide_northsouth$logCPUE_Hippoglossus_stenolepis[is.na(wide_northsouth$logCPUE_Hippoglossus_stenolepis)==TRUE] <- 0
wide_northsouth$logCPUE_Hippoglossoides_elassodon[is.na(wide_northsouth$logCPUE_Hippoglossoides_elassodon)==TRUE] <- 0
wide_northsouth$logCPUE_Limanda_aspera[is.na(wide_northsouth$logCPUE_Limanda_aspera)==TRUE] <- 0
wide_northsouth$logCPUE_Lepidopsetta_polyxystra[is.na(wide_northsouth$logCPUE_Lepidopsetta_polyxystra)==TRUE] <- 0
wide_northsouth$logCPUE_Pleuronectes_quadrituberculatus[is.na(wide_northsouth$logCPUE_Pleuronectes_quadrituberculatus)==TRUE] <- 0
wide_northsouth$logCPUE_Gadus_macrocephalus[is.na(wide_northsouth$logCPUE_Gadus_macrocephalus)==TRUE] <- 0
wide_northsouth$logCPUE_Gadus_chalcogrammus[is.na(wide_northsouth$logCPUE_Gadus_chalcogrammus)==TRUE] <- 0
wide_northsouth$logCPUE_Chionoecetes_bairdi[is.na(wide_northsouth$logCPUE_Chionoecetes_bairdi)==TRUE] <- 0
wide_northsouth$logCPUE_Chionoecetes_opilio[is.na(wide_northsouth$logCPUE_Chionoecetes_opilio)==TRUE] <- 0
wide_northsouth$logCPUE_Lepidopsetta_sp[is.na(wide_northsouth$logCPUE_Lepidopsetta_sp)==TRUE] <- 0


yrs <- unique(wide_northsouth$YEAR)
#loop through species, loop through each year, get correlation and store

yr_vec <- vector(mode="numeric", length=length(yrs)*10)
sps_vec <- vector(mode="numeric", length=length(yrs)*10)
cor_vec <- vector(mode="numeric", length=length(yrs)*10)

#yrs <- unique(wide_northsouth$YEAR)
i <- 1
counter <- 1
for(i in 1:length(yrs)){
  temp_yr <- yrs[i]
  temp_wide <- wide_northsouth[which(wide_northsouth$YEAR==temp_yr),]
  #tempmat <- temp_wide[,c(38:48)]
  #tempmat <- temp_wide[,c(38:47)]
  tempmat <- temp_wide[,c(18:28)]
  
  tempcor <- cor(tempmat, method = "pearson", use = "complete.obs")
  #corrplot(tempcor)
  
  #save output
  for(l in 1:10){
    yr_vec[counter] <- temp_yr      
    sps_vec[counter] <- rownames(tempcor)[l]
    cor_vec[counter] <- tempcor[l,8]
    
    counter <- counter + 1
  }
}

cor_output <- data.frame(yr_vec, sps_vec, cor_vec)

cor_output$period <- NA
cor_output$period <- as.character(cor_output$period)
cor_output$period[which(cor_output$yr_vec<2014)]<-"early"
cor_output$period[which(cor_output$yr_vec>2013)]<- "late"
cor_output$period <- as.factor(cor_output$period)

write.csv(cor_output, file=paste(wd,"/data/community-correlations-v2_sebs_nebs.csv", sep=""))





analysis_cor <- cor_output[which(cor_output$sps_vec!="logCPUE_Gadus_chalcogrammus"),]

plot_cor <- analysis_cor
plot_cor$Species <- NA
plot_cor$Species <- as.character(plot_cor$Species)

plot_cor$Species[which(plot_cor$sps_vec=="logCPUE_Atheresthes_stomias")] <- "Arrowtooth flounder"
plot_cor$Species[which(plot_cor$sps_vec=="logCPUE_Hippoglossus_stenolepis")] <- "Pacific halibut"
plot_cor$Species[which(plot_cor$sps_vec=="logCPUE_Hippoglossoides_elassodon")] <- "Flathead sole"
plot_cor$Species[which(plot_cor$sps_vec=="logCPUE_Limanda_aspera")] <- "Yellowfin sole"
plot_cor$Species[which(plot_cor$sps_vec=="logCPUE_Lepidopsetta_polyxystra")] <- "Northern rock sole"
plot_cor$Species[which(plot_cor$sps_vec=="logCPUE_Pleuronectes_quadrituberculatus")] <- "Alaska plaice"
plot_cor$Species[which(plot_cor$sps_vec=="logCPUE_Gadus_macrocephalus")] <- "Pacific cod"
plot_cor$Species[which(plot_cor$sps_vec=="logCPUE_Chionoecetes_bairdi")] <- "Tanner (bairdi) crab"
plot_cor$Species[which(plot_cor$sps_vec=="logCPUE_Chionoecetes_opilio")] <- "Snow (opilio) crab"


ggplot(plot_cor, aes(period, cor_vec)) + geom_boxplot() + facet_wrap(~Species) + 
  geom_hline(yintercept = 0, col="red") + ylab("Correlation") + xlab("Period")

cormod <- lm(cor_vec ~ period, data = analysis_cor)
anova(cormod)

#need to do this by sps
# t.test() function makes the assumption that the variances of the two groups of samples, 
# being compared, are different. Therefore, Welch t-test is performed by default.

t_As <- t.test(cor_vec ~ period, data = analysis_cor[which(analysis_cor$sps_vec=="logCPUE_Atheresthes_stomias"),])
t_As #not sig p-value = 0.1615
p_As <- t_As$p.value

t_Hs <- t.test(cor_vec ~ period, data = analysis_cor[which(analysis_cor$sps_vec=="logCPUE_Hippoglossus_stenolepis"),])
t_Hs #p=0.008242
p_Hs <- t_Hs$p.value

t_He <- t.test(cor_vec ~ period, data = analysis_cor[which(analysis_cor$sps_vec=="logCPUE_Hippoglossoides_elassodon"),])
t_He #p=0.0.03729
p_He <- t_He$p.value

t_La <- t.test(cor_vec ~ period, data = analysis_cor[which(analysis_cor$sps_vec=="logCPUE_Limanda_aspera"),])
t_La #not sig p-value = 0.1952
p_La <- t_La$p.value

t_Lp <- t.test(cor_vec ~ period, data = analysis_cor[which(analysis_cor$sps_vec=="logCPUE_Lepidopsetta_polyxystra"),])
t_Lp #not sig p-value = 0.7689
p_Lp <- t_Lp$p.value

t_Pq <- t.test(cor_vec ~ period, data = analysis_cor[which(analysis_cor$sps_vec=="logCPUE_Pleuronectes_quadrituberculatus"),])
t_Pq #p=0.0001435
p_Pq <- t_Pq$p.value

t_Gm <- t.test(cor_vec ~ period, data = analysis_cor[which(analysis_cor$sps_vec=="logCPUE_Gadus_macrocephalus"),])
t_Gm #NOT SIG p=0.1561
p_Gm <- t_Gm$p.value

t_Cb <- t.test(cor_vec ~ period, data = analysis_cor[which(analysis_cor$sps_vec=="logCPUE_Chionoecetes_bairdi"),])
t_Cb #sig p=0.01666
p_Cb <- t_Cb$p.value

t_Co <- t.test(cor_vec ~ period, data = analysis_cor[which(analysis_cor$sps_vec=="logCPUE_Chionoecetes_opilio"),])
t_Co #not sig p-value = 0.8072
p_Co <- t_Co$p.value

ggplot(cor_output, aes(yr_vec, cor_vec, colour=sps_vec)) + geom_point() + geom_smooth() + facet_wrap(~sps_vec) + 
  geom_hline(yintercept = 0, col="red")


#FDR control=================

#following Verhoeven et al 2005 Oikos to control false discovery rate

m<- length(unique(plot_cor$Species))
alpha <- 0.05

pvec <- as.vector(c(p_As, p_Hs, p_He, p_La, p_Lp, p_Pq, p_Gm, p_Cb, p_Co))

pvec <- sort(pvec, decreasing=FALSE)

i<-1
for(i in 1:m){
  temp_p <- pvec[i]
  print(temp_p)
  temp_value <- (alpha/m)*i
  print(temp_value)
  print(temp_p < temp_value)
}
#NOW FIRST THREE remain significant, p_Cb, p_Hs, p_Pq
#Chionoecetes_bairdi
#Pleuronectes_quadrituberculatus
#Hippoglossus_stenolepis



#maps=====

library("rnaturalearth")
library("rnaturalearthdata")
library( "ggspatial" )
library("sf")



world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Atheresthes_stomias), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  


ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Hippoglossus_stenolepis), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Hippoglossoides_elassodon), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Limanda_aspera), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Lepidopsetta_polyxystra), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  


ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Pleuronectes_quadrituberculatus), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  


ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Gadus_macrocephalus), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  


ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Gadus_chalcogrammus), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  


ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Chionoecetes_bairdi), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  


ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Chionoecetes_opilio), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  



ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Lepidopsetta_sp), data=wide_northsouth) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period)  



#repeat SEBS only===============================================================

wide_sebs <- wide_northsouth[which(wide_northsouth$STRATUM!="70" &
                                     wide_northsouth$STRATUM!="71" &
                                     wide_northsouth$STRATUM!="81" &
                                     wide_northsouth$STRATUM!="0"),]


yrs <- unique(wide_sebs$YEAR)
#loop through species, loop through each year, get correlation and store

yr_vec <- vector(mode="numeric", length=length(yrs)*10)
sps_vec <- vector(mode="numeric", length=length(yrs)*10)
cor_vec <- vector(mode="numeric", length=length(yrs)*10)

#yrs <- unique(wide_northsouth$YEAR)
i <- 1
counter <- 1
for(i in 1:length(yrs)){
  temp_yr <- yrs[i]
  temp_wide <- wide_sebs[which(wide_sebs$YEAR==temp_yr),]
  #tempmat <- temp_wide[,c(38:48)]
  #tempmat <- temp_wide[,c(38:47)]
  tempmat <- temp_wide[,c(18:28)]
  
  tempcor <- cor(tempmat, method = "pearson", use = "complete.obs")
  #corrplot(tempcor)
  
  #save output
  for(l in 1:10){
    yr_vec[counter] <- temp_yr      
    sps_vec[counter] <- rownames(tempcor)[l]
    cor_vec[counter] <- tempcor[l,8]
    
    counter <- counter + 1
  }
}

cor_output_sebs <- data.frame(yr_vec, sps_vec, cor_vec)

cor_output_sebs$period <- NA
cor_output_sebs$period <- as.character(cor_output_sebs$period)
cor_output_sebs$period[which(cor_output_sebs$yr_vec<2014)]<-"early"
cor_output_sebs$period[which(cor_output_sebs$yr_vec>2013)]<- "late"
cor_output_sebs$period <- as.factor(cor_output_sebs$period)

#write.csv(cor_output_sebs, file=paste(wd,"/data/community-correlations-v2_sebs.csv", sep=""))

cor_output_sebs <- read.csv(file=paste(wd,"/data/community-correlations-v2_sebs.csv", sep=""))




analysis_cor_sebs <- cor_output_sebs[which(cor_output_sebs$sps_vec!="logCPUE_Gadus_chalcogrammus"),]

plot_cor_sebs <- analysis_cor_sebs
plot_cor_sebs$Species <- NA
plot_cor_sebs$Species <- as.character(plot_cor_sebs$Species)

plot_cor_sebs$Species[which(plot_cor_sebs$sps_vec=="logCPUE_Atheresthes_stomias")] <- "Arrowtooth flounder"
plot_cor_sebs$Species[which(plot_cor_sebs$sps_vec=="logCPUE_Hippoglossus_stenolepis")] <- "Pacific halibut"
plot_cor_sebs$Species[which(plot_cor_sebs$sps_vec=="logCPUE_Hippoglossoides_elassodon")] <- "Flathead sole"
plot_cor_sebs$Species[which(plot_cor_sebs$sps_vec=="logCPUE_Limanda_aspera")] <- "Yellowfin sole"
plot_cor_sebs$Species[which(plot_cor_sebs$sps_vec=="logCPUE_Lepidopsetta_polyxystra")] <- "Northern rock sole"
plot_cor_sebs$Species[which(plot_cor_sebs$sps_vec=="logCPUE_Pleuronectes_quadrituberculatus")] <- "Alaska plaice"
plot_cor_sebs$Species[which(plot_cor_sebs$sps_vec=="logCPUE_Gadus_macrocephalus")] <- "Pacific cod"
plot_cor_sebs$Species[which(plot_cor_sebs$sps_vec=="logCPUE_Chionoecetes_bairdi")] <- "Tanner (bairdi) crab"
plot_cor_sebs$Species[which(plot_cor_sebs$sps_vec=="logCPUE_Chionoecetes_opilio")] <- "Snow (opilio) crab"


ggplot(plot_cor_sebs, aes(period, cor_vec)) + geom_boxplot() + facet_wrap(~Species) + 
  geom_hline(yintercept = 0, col="red") + ylab("Correlation") + xlab("Period") + theme_bw()


#need to do this by sps
# t.test() function makes the assumption that the variances of the two groups of samples, 
# being compared, are different. Therefore, Welch t-test is performed by default.

t_As_sebs <- t.test(cor_vec ~ period, data = analysis_cor_sebs[which(analysis_cor_sebs$sps_vec=="logCPUE_Atheresthes_stomias"),])
t_As_sebs #not sig p-value = 0.1202
p_As_sebs <- t_As_sebs$p.value

t_Hs_sebs <- t.test(cor_vec ~ period, data = analysis_cor_sebs[which(analysis_cor_sebs$sps_vec=="logCPUE_Hippoglossus_stenolepis"),])
t_Hs_sebs #p=0.01625
p_Hs_sebs <- t_Hs_sebs$p.value

t_He_sebs <- t.test(cor_vec ~ period, data = analysis_cor_sebs[which(analysis_cor_sebs$sps_vec=="logCPUE_Hippoglossoides_elassodon"),])
t_He_sebs #p=0.0.03019
p_He_sebs <- t_He_sebs$p.value

t_La_sebs <- t.test(cor_vec ~ period, data = analysis_cor_sebs[which(analysis_cor_sebs$sps_vec=="logCPUE_Limanda_aspera"),])
t_La_sebs #not sig p-value = 0.1768
p_La_sebs <- t_La_sebs$p.value

t_Lp_sebs <- t.test(cor_vec ~ period, data = analysis_cor_sebs[which(analysis_cor_sebs$sps_vec=="logCPUE_Lepidopsetta_polyxystra"),])
t_Lp_sebs #not sig p-value = 0.2923
p_Lp_sebs <- t_Lp_sebs$p.value

t_Pq_sebs <- t.test(cor_vec ~ period, data = analysis_cor_sebs[which(analysis_cor_sebs$sps_vec=="logCPUE_Pleuronectes_quadrituberculatus"),])
t_Pq_sebs #p=0.0007213
p_Pq_sebs <- t_Pq_sebs$p.value

t_Gm_sebs <- t.test(cor_vec ~ period, data = analysis_cor_sebs[which(analysis_cor_sebs$sps_vec=="logCPUE_Gadus_macrocephalus"),])
t_Gm_sebs # p=0.03339
p_Gm_sebs <- t_Gm_sebs$p.value

t_Cb_sebs <- t.test(cor_vec ~ period, data = analysis_cor_sebs[which(analysis_cor_sebs$sps_vec=="logCPUE_Chionoecetes_bairdi"),])
t_Cb_sebs #sig p=0.01442
p_Cb_sebs <- t_Cb_sebs$p.value

t_Co_sebs <- t.test(cor_vec ~ period, data = analysis_cor_sebs[which(analysis_cor_sebs$sps_vec=="logCPUE_Chionoecetes_opilio"),])
t_Co_sebs #not sig p-value = 0.8586
p_Co_sebs <- t_Co_sebs$p.value

ggplot(cor_output_sebs, aes(yr_vec, cor_vec, colour=sps_vec)) + geom_point() + geom_smooth() + facet_wrap(~sps_vec) + 
  geom_hline(yintercept = 0, col="red")


#FDR control=================

#following Verhoeven et al 2005 Oikos to control false discovery rate

m<- length(unique(plot_cor_sebs$Species))
alpha <- 0.05

pvec <- as.vector(c(p_As_sebs, p_Hs_sebs, p_He_sebs, p_La_sebs, p_Lp_sebs, p_Pq_sebs, p_Gm_sebs, p_Cb_sebs, p_Co_sebs))

pvec <- sort(pvec, decreasing=FALSE)

i<-1
for(i in 1:m){
  temp_p <- pvec[i]
  print(temp_p)
  temp_value <- (alpha/m)*i
  print(temp_value)
  print(temp_p < temp_value)
}
#NOW FIRST THREE-FOUR remain significant, 







#repeat NEBS only===============================================================

wide_nebs <- wide_northsouth[which(wide_northsouth$STRATUM=="70" |
                                     wide_northsouth$STRATUM=="71" |
                                     wide_northsouth$STRATUM=="81"),]


yrs <- unique(wide_nebs$YEAR)
#loop through species, loop through each year, get correlation and store

yr_vec <- vector(mode="numeric", length=length(yrs)*10)
sps_vec <- vector(mode="numeric", length=length(yrs)*10)
cor_vec <- vector(mode="numeric", length=length(yrs)*10)

#yrs <- unique(wide_northsouth$YEAR)
i <- 1
counter <- 1
for(i in 1:length(yrs)){
  temp_yr <- yrs[i]
  temp_wide <- wide_nebs[which(wide_nebs$YEAR==temp_yr),]
  #tempmat <- temp_wide[,c(38:48)]
  #tempmat <- temp_wide[,c(38:47)]
  tempmat <- temp_wide[,c(18:28)]
  
  tempcor <- cor(tempmat, method = "pearson", use = "complete.obs")
  #corrplot(tempcor)
  
  #save output
  for(l in 1:10){
    yr_vec[counter] <- temp_yr      
    sps_vec[counter] <- rownames(tempcor)[l]
    cor_vec[counter] <- tempcor[l,8]
    
    counter <- counter + 1
  }
}

cor_output_nebs <- data.frame(yr_vec, sps_vec, cor_vec)

cor_output_nebs$period <- NA
cor_output_nebs$period <- as.character(cor_output_nebs$period)
cor_output_nebs$period[which(cor_output_nebs$yr_vec<2014)]<-"early"
cor_output_nebs$period[which(cor_output_nebs$yr_vec>2013)]<- "late"
cor_output_nebs$period <- as.factor(cor_output_nebs$period)

write.csv(cor_output_nebs, file=paste(wd,"/data/community-correlations-v2_nebs.csv", sep=""))





analysis_cor_nebs <- cor_output_nebs[which(cor_output_nebs$sps_vec!="logCPUE_Gadus_chalcogrammus"),]

plot_cor_nebs <- analysis_cor_nebs
plot_cor_nebs$Species <- NA
plot_cor_nebs$Species <- as.character(plot_cor_nebs$Species)

plot_cor_nebs$Species[which(plot_cor_nebs$sps_vec=="logCPUE_Atheresthes_stomias")] <- "Arrowtooth flounder"
plot_cor_nebs$Species[which(plot_cor_nebs$sps_vec=="logCPUE_Hippoglossus_stenolepis")] <- "Pacific halibut"
plot_cor_nebs$Species[which(plot_cor_nebs$sps_vec=="logCPUE_Hippoglossoides_elassodon")] <- "Flathead sole"
plot_cor_nebs$Species[which(plot_cor_nebs$sps_vec=="logCPUE_Limanda_aspera")] <- "Yellowfin sole"
plot_cor_nebs$Species[which(plot_cor_nebs$sps_vec=="logCPUE_Lepidopsetta_polyxystra")] <- "Northern rock sole"
plot_cor_nebs$Species[which(plot_cor_nebs$sps_vec=="logCPUE_Pleuronectes_quadrituberculatus")] <- "Alaska plaice"
plot_cor_nebs$Species[which(plot_cor_nebs$sps_vec=="logCPUE_Gadus_macrocephalus")] <- "Pacific cod"
plot_cor_nebs$Species[which(plot_cor_nebs$sps_vec=="logCPUE_Chionoecetes_bairdi")] <- "Tanner (bairdi) crab"
plot_cor_nebs$Species[which(plot_cor_nebs$sps_vec=="logCPUE_Chionoecetes_opilio")] <- "Snow (opilio) crab"


ggplot(plot_cor_nebs, aes(period, cor_vec)) + geom_boxplot() + facet_wrap(~Species) + 
  geom_hline(yintercept = 0, col="red") + ylab("Correlation") + xlab("Period")


#need to do this by sps
# t.test() function makes the assumption that the variances of the two groups of samples, 
# being compared, are different. Therefore, Welch t-test is performed by default.

# t_As_nebs <- t.test(cor_vec ~ period, data = analysis_cor_nebs[which(analysis_cor_nebs$sps_vec=="logCPUE_Atheresthes_stomias"),])
# t_As_nebs #not sig p-value = 0.1202
# p_As_nebs <- t_As_nebs$p.value

t_Hs_nebs <- t.test(cor_vec ~ period, data = analysis_cor_nebs[which(analysis_cor_nebs$sps_vec=="logCPUE_Hippoglossus_stenolepis"),])
t_Hs_nebs #p=0.4345
p_Hs_nebs <- t_Hs_nebs$p.value

# t_He_nebs <- t.test(cor_vec ~ period, data = analysis_cor_nebs[which(analysis_cor_nebs$sps_vec=="logCPUE_Hippoglossoides_elassodon"),])
# t_He_nebs #p=0.0.03019
# p_He_nebs <- t_He_nebs$p.value

t_La_nebs <- t.test(cor_vec ~ period, data = analysis_cor_nebs[which(analysis_cor_nebs$sps_vec=="logCPUE_Limanda_aspera"),])
t_La_nebs #not sig p-value = 0.0007298
p_La_nebs <- t_La_nebs$p.value

# t_Lp_nebs <- t.test(cor_vec ~ period, data = analysis_cor_nebs[which(analysis_cor_nebs$sps_vec=="logCPUE_Lepidopsetta_polyxystra"),])
# t_Lp_nebs #not sig p-value = 0.2923
# p_Lp_nebs <- t_Lp_nebs$p.value

t_Pq_nebs <- t.test(cor_vec ~ period, data = analysis_cor_nebs[which(analysis_cor_nebs$sps_vec=="logCPUE_Pleuronectes_quadrituberculatus"),])
t_Pq_nebs #p=0.3113
p_Pq_nebs <- t_Pq_nebs$p.value

t_Gm_nebs <- t.test(cor_vec ~ period, data = analysis_cor_nebs[which(analysis_cor_nebs$sps_vec=="logCPUE_Gadus_macrocephalus"),])
t_Gm_nebs # p=0.4459
p_Gm_nebs <- t_Gm_nebs$p.value

t_Cb_nebs <- t.test(cor_vec ~ period, data = analysis_cor_nebs[which(analysis_cor_nebs$sps_vec=="logCPUE_Chionoecetes_bairdi"),])
t_Cb_nebs #sig p=0.6641
p_Cb_nebs <- t_Cb_nebs$p.value

t_Co_nebs <- t.test(cor_vec ~ period, data = analysis_cor_nebs[which(analysis_cor_nebs$sps_vec=="logCPUE_Chionoecetes_opilio"),])
t_Co_nebs #not sig p-value = 0.06053
p_Co_nebs <- t_Co_nebs$p.value

ggplot(cor_output_nebs, aes(yr_vec, cor_vec, colour=sps_vec)) + geom_point() + geom_smooth() + facet_wrap(~sps_vec) + 
  geom_hline(yintercept = 0, col="red")


#FDR control=================

#following Verhoeven et al 2005 Oikos to control false discovery rate

m<- length(unique(plot_cor_nebs$Species))
alpha <- 0.05

pvec <- as.vector(c(p_Hs_nebs,  p_La_nebs,  p_Pq_nebs, p_Gm_nebs, p_Cb_nebs, p_Co_nebs)) #works for fewer species

pvec <- sort(pvec, decreasing=FALSE)

i<-1
for(i in 1:m){
  temp_p <- pvec[i]
  print(temp_p)
  temp_value <- (alpha/m)*i
  print(temp_value)
  print(temp_p < temp_value)
}
#just one?





