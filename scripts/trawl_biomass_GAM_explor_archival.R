#============================================================================================================================================
# GAM EXPLORATION - survey abundance data

#Created by Krista, May 2020
#============================================================================================================================================
#Notes:
#============================================================================================================================================
library(ggplot2)
library(tidyverse)
library(corrplot)
library(mgcv)
library(nlme)
library(mgcViz)
library(tidyverse)
library(mapproj)
library(visreg)



sel.trawl.dat <- read.csv("data/select_trawl_dat.csv", stringsAsFactors = FALSE, row.names = 1)

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

#convert latitude & longitude to x/y coordinates that reflect 
# actual distances using the equal-distance Albers projection:
# (I chose the parameters based on a 'rule of thumb' in the help file)
x <- mapproject(sel.trawl.dat$LONGITUDE, sel.trawl.dat$LATITUDE, "albers", param=c(55.9, 60.8))
sel.trawl.dat$long_albers <- x$x
sel.trawl.dat$lat_albers <- x$y

#want to remove SID & common which repeat data from scientific (MOST of the time!) and cause issues later
#doing this in a silly brute force way because I want to avoid indexing by col #s for better stability if data is read in repeatedly



sel.trawl.dat <- sel.trawl.dat[,c("LATITUDE", "LONGITUDE","STATION","STRATUM","YEAR","DATETIME","WTCPUE","NUMCPUE",     
                          "SCIENTIFIC","BOT_DEPTH","BOT_TEMP","SURF_TEMP","VESSEL","CRUISE","HAUL","YEAR_factor",
                          "logCPUE","long_albers","lat_albers" )]

#select early data======================================================

#we only want data pre2014 in the early period
early_dat <- sel.trawl.dat[which(sel.trawl.dat$YEAR<2014),]




#widen dataframe to have columns for each sps

# early_wide <- early_dat %>% pivot_wider(names_from=c(LONGITUDE, LATITUDE, STATION, STRATUM, YEAR, 
#                                                      DATETIME, VESSEL, CRUISE, HAUL), 
#                                         values_from=c(WTCPUE, NUMCPUE, logCPUE))

early_wide <- early_dat %>% pivot_wider(names_from=SCIENTIFIC, 
                                        values_from=c(WTCPUE, NUMCPUE, logCPUE))
#eek this causes some bad column names!
early_wide <- early_wide %>% rename(WTCPUE_Chionoecetes_bairdi = "WTCPUE_Chionoecetes bairdi",            
  WTCPUE_Atheresthes_stomia = "WTCPUE_Atheresthes stomias",             
  WTCPUE_Hippoglossus_stenolepis = "WTCPUE_Hippoglossus stenolepis",
  WTCPUE_Limanda_aspera = "WTCPUE_Limanda aspera",         
  WTCPUE_Lepidopsetta_sp= "WTCPUE_Lepidopsetta sp.",               
  WTCPUE_Chionoecetes_opilio = "WTCPUE_Chionoecetes opilio",             
  WTCPUE_Gadus_macrocephalus = "WTCPUE_Gadus macrocephalus",            
  WTCPUE_Hippoglossoides_elassodon = "WTCPUE_Hippoglossoides elassodon",  
  WTCPUE_Pleuronectes_quadrituberculatus = "WTCPUE_Pleuronectes quadrituberculatus", 
  WTCPUE_Lepidopsetta_polyxystra = "WTCPUE_Lepidopsetta polyxystra", 
  WTCPUE_Gadus_chalcogrammus = "WTCPUE_Gadus chalcogrammus", 
  NUMCPUE_Gadus_chalcogrammus = "NUMCPUE_Gadus chalcogrammus",        
  NUMCPUE_Chionoecetes_bairdi = "NUMCPUE_Chionoecetes bairdi",            
  NUMCPUE_Atheresthes_stomias = "NUMCPUE_Atheresthes stomias",           
  NUMCPUE_Hippoglossus_stenolepis = "NUMCPUE_Hippoglossus stenolepis",        
  NUMCPUE_Limanda_aspera = "NUMCPUE_Limanda aspera",                
  NUMCPUE_Lepidopsetta_sp = "NUMCPUE_Lepidopsetta sp.",               
  NUMCPUE_Chionoecetes_opilio = "NUMCPUE_Chionoecetes opilio",           
  NUMCPUE_Gadus_macrocephalus = "NUMCPUE_Gadus macrocephalus",            
  NUMCPUE_Hippoglossoides_elassodon = "NUMCPUE_Hippoglossoides elassodon",     
  NUMCPUE_Pleuronectes_quadrituberculatus = "NUMCPUE_Pleuronectes quadrituberculatus",
  NUMCPUE_Lepidopsetta_polyxystra = "NUMCPUE_Lepidopsetta polyxystra",
  logCPUE_Gadus_chalcogrammus = "logCPUE_Gadus chalcogrammus",        
  logCPUE_Chionoecetes_bairdi = "logCPUE_Chionoecetes bairdi",           
  logCPUE_Atheresthes_stomias = "logCPUE_Atheresthes stomias",            
  logCPUE_Hippoglossus_stenolepis = "logCPUE_Hippoglossus stenolepis", 
  logCPUE_Limanda_aspera = "logCPUE_Limanda aspera",        
  logCPUE_Lepidopsetta_sp = "logCPUE_Lepidopsetta sp.",              
  logCPUE_Chionoecetes_opilio = "logCPUE_Chionoecetes opilio",            
  logCPUE_Gadus_macrocephalus = "logCPUE_Gadus macrocephalus",           
  logCPUE_Hippoglossoides_elassodon = "logCPUE_Hippoglossoides elassodon",      
  logCPUE_Pleuronectes_quadrituberculatus = "logCPUE_Pleuronectes quadrituberculatus" ,
  logCPUE_Lepidopsetta_polyxystra = "logCPUE_Lepidopsetta polyxystra") 




#exclusion criteria===============

station_summary <- early_wide %>% group_by(STATION, LONGITUDE, LATITUDE) %>%
  summarize(n_yrs=n())

station_summary2 <- early_wide %>% group_by(STATION) %>%
  summarize(n_yrs=n())

joinearly <- left_join(early_wide, station_summary2)

analysis_dat <- joinearly[which(joinearly$n_yrs>5),]


#select full data=============================================

#full_dat <- sel.trawl.dat[,-c(9,11)] #drop columns that repeat info
full_dat <- sel.trawl.dat

full_dat <- full_dat[!duplicated(full_dat),] #there are 4 duplicate rows form 2016

#widen dataframe to have columns for each sps


full_wide <- full_dat %>% pivot_wider(names_from=SCIENTIFIC, 
                                        values_from=c(WTCPUE, NUMCPUE, logCPUE),
                                      values_fill=list(WTCPUE=0, NUMCPUE=0, logCPUE=0))


# full_wide <- full_dat %>% pivot_wider(names_from=SCIENTIFIC, 
#                                      # id_cols=c(CRUISE,HAUL,SID),
#                                       values_from=c(logCPUE),
#                                     #  values_fn = list(logCPUE = length),
#                                       values_fill=list(logCPUE=0)) 


#eek this causes some bad column names!
full_wide <- full_wide %>% rename(WTCPUE_Chionoecetes_bairdi = "WTCPUE_Chionoecetes bairdi",            
                                    WTCPUE_Atheresthes_stomia = "WTCPUE_Atheresthes stomias",             
                                    WTCPUE_Hippoglossus_stenolepis = "WTCPUE_Hippoglossus stenolepis",
                                    WTCPUE_Limanda_aspera = "WTCPUE_Limanda aspera",         
                                    WTCPUE_Lepidopsetta_sp= "WTCPUE_Lepidopsetta sp.",               
                                    WTCPUE_Chionoecetes_opilio = "WTCPUE_Chionoecetes opilio",             
                                    WTCPUE_Gadus_macrocephalus = "WTCPUE_Gadus macrocephalus",            
                                    WTCPUE_Hippoglossoides_elassodon = "WTCPUE_Hippoglossoides elassodon",  
                                    WTCPUE_Pleuronectes_quadrituberculatus = "WTCPUE_Pleuronectes quadrituberculatus", 
                                    WTCPUE_Lepidopsetta_polyxystra = "WTCPUE_Lepidopsetta polyxystra", 
                                    WTCPUE_Gadus_chalcogrammus = "WTCPUE_Gadus chalcogrammus", 
                                    NUMCPUE_Gadus_chalcogrammus = "NUMCPUE_Gadus chalcogrammus",        
                                    NUMCPUE_Chionoecetes_bairdi = "NUMCPUE_Chionoecetes bairdi",            
                                    NUMCPUE_Atheresthes_stomias = "NUMCPUE_Atheresthes stomias",           
                                    NUMCPUE_Hippoglossus_stenolepis = "NUMCPUE_Hippoglossus stenolepis",        
                                    NUMCPUE_Limanda_aspera = "NUMCPUE_Limanda aspera",                
                                    NUMCPUE_Lepidopsetta_sp = "NUMCPUE_Lepidopsetta sp.",               
                                    NUMCPUE_Chionoecetes_opilio = "NUMCPUE_Chionoecetes opilio",           
                                    NUMCPUE_Gadus_macrocephalus = "NUMCPUE_Gadus macrocephalus",            
                                    NUMCPUE_Hippoglossoides_elassodon = "NUMCPUE_Hippoglossoides elassodon",     
                                    NUMCPUE_Pleuronectes_quadrituberculatus = "NUMCPUE_Pleuronectes quadrituberculatus",
                                    NUMCPUE_Lepidopsetta_polyxystra = "NUMCPUE_Lepidopsetta polyxystra",
                                    logCPUE_Gadus_chalcogrammus = "logCPUE_Gadus chalcogrammus",        
                                    logCPUE_Chionoecetes_bairdi = "logCPUE_Chionoecetes bairdi",           
                                    logCPUE_Atheresthes_stomias = "logCPUE_Atheresthes stomias",            
                                    logCPUE_Hippoglossus_stenolepis = "logCPUE_Hippoglossus stenolepis", 
                                    logCPUE_Limanda_aspera = "logCPUE_Limanda aspera",        
                                    logCPUE_Lepidopsetta_sp = "logCPUE_Lepidopsetta sp.",              
                                    logCPUE_Chionoecetes_opilio = "logCPUE_Chionoecetes opilio",            
                                    logCPUE_Gadus_macrocephalus = "logCPUE_Gadus macrocephalus",           
                                    logCPUE_Hippoglossoides_elassodon = "logCPUE_Hippoglossoides elassodon",      
                                    logCPUE_Pleuronectes_quadrituberculatus = "logCPUE_Pleuronectes quadrituberculatus" ,
                                    logCPUE_Lepidopsetta_polyxystra = "logCPUE_Lepidopsetta polyxystra") 


write.csv(full_wide, file=paste(wd,"/data/full_wide_comm_dat.csv", sep=""))


#exclusion criteria===============

station_summaryF <- full_wide %>% group_by(STATION, LONGITUDE, LATITUDE) %>%
  summarize(n_yrs=n())

station_summaryF2 <- full_wide %>% group_by(STATION) %>%
  summarize(n_yrs=n())

joinfull <- left_join(full_wide, station_summaryF2)

#analysis_dat <- joinearly[which(joinearly$n_yrs>5),]

analysis_dat <- joinfull

#plot CPUE maps======

t1 <- ggplot(early_wide, aes(LONGITUDE, LATITUDE, col=logCPUE_Gadus_chalcogrammus))
t1 + geom_point() + facet_wrap(~YEAR) +
  scale_colour_gradient2(low="blue", high="red", guide="colorbar")

t2 <- ggplot(early_wide, aes(LONGITUDE, LATITUDE, fill=logCPUE_Gadus_chalcogrammus))
t2 + geom_tile() # + facet_wrap(~YEAR)

t3 <- ggplot(early_wide, aes(STATION))
t3 + geom_histogram(stat="count")




#plot covars======================================================================================

plot(early_wide$YEAR, early_wide$BOT_TEMP)
plot(early_wide$YEAR, early_wide$SURF_TEMP)

plot( early_wide$BOT_TEMP, early_wide$logCPUE_Gadus_chalcogrammus)
plot( early_wide$BOT_DEPTH, early_wide$logCPUE_Gadus_chalcogrammus)
plot( early_wide$SURF_TEMP, early_wide$logCPUE_Gadus_chalcogrammus)


ggplot(early_wide, aes(YEAR, BOT_TEMP)) + geom_point()




#NMDS=======================

#install.packages("vegan")
library(vegan)
set.seed(2)

early_comm_mat <- as.matrix(early_wide[,c(38:48)])
#going to replace NAs with zeros, should discuss whether this is a good choice or not
early_comm_mat[is.na(early_comm_mat)==TRUE]<-0
unique(is.finite(early_comm_mat)) #cool no NAs

short_early_mat <- early_comm_mat[1:1000,]
short.rel <- decostand(short_early_mat, method = "total")

short_NMDS <- metaMDS(short.rel, # Our community-by-species matrix
                      distance = "bray",
                             k=2, # The number of reduced dimensions
                             na.rm=TRUE) 
summary(short_NMDS)
short_NMDS$stress
shortptsall<-scores(short_NMDS)
scores(shortptsall, display="species")
ordiplot(short_NMDS)
goodness(short_NMDS) # Produces a results of test statistics for goodness of fit for each point
stressplot(short_NMDS) # Produces a Shepards diagram

#some do distance first, try that way
short_distmat <- 
  vegdist(short.rel, method = "bray")
short_DistanceMatrix <- 
  as.matrix(short_distmat, labels = T)
short_NMDS2 <- metaMDS(short_DistanceMatrix, # Our community-by-species matrix
                      distance = "bray",
                      k=2, # The number of reduced dimensions
                      na.rm=TRUE) 
summary(short_NMDS2)
short_NMDS2$stress
shortptsall2<-scores(short_NMDS2)
scores(shortptsall2, display="species")
ordiplot(short_NMDS2)
goodness(short_NMDS2) # Produces a results of test statistics for goodness of fit for each point
stressplot(short_NMDS2)


lessshort_early_mat <- early_comm_mat[1:5000,]
lessshort.rel <- decostand(lessshort_early_mat, method = "total")

lessshort_NMDS <- metaMDS(lessshort.rel, # Our community-by-species matrix
                      distance = "bray",
                      k=2, # The number of reduced dimensions
                      na.rm=TRUE) 
summary(lessshort_NMDS)
lessshort_NMDS$stress
lessshortptsall<-scores(lessshort_NMDS)
scores(lessshortptsall, display="species")
ordiplot(lessshort_NMDS)




#CAUTION SUPER SLOW
exampleearly_NMDS <- metaMDS(early_comm_mat, # Our community-by-species matrix
                     k=2, # The number of reduced dimensions
                     na.rm=TRUE,
                     trymax=50,
                     maxit=50) 


exampleearly_NMDSk3 <- metaMDS(early_comm_mat, # Our community-by-species matrix
                             k=3, # The number of reduced dimensions
                             na.rm=TRUE,
                             trymax=20,
                             maxit=20) 
#k=3 leads to far worse stress on most runs but lower overall? Still no convergence

#should add autotranform=FALSE to all these

exampleearly_NMDSk3_2 <- metaMDS(early_comm_mat, # Our community-by-species matrix
                               k=3, # The number of reduced dimensions
                               na.rm=TRUE,
                               trymax=40,
                               maxit=40) 

exampleearly_NMDSk3_3 <- metaMDS(early_comm_mat, # Our community-by-species matrix
                                 k=3, # The number of reduced dimensions
                                 na.rm=TRUE,
                                 trymax=60,
                                 maxit=60) 


#this one runs out of memory
exampleearly_NMDSk3_A <- metaMDS(early_comm_mat, # Our community-by-species matrix
                                 k=3, # The number of reduced dimensions
                                 na.rm=TRUE, noshare=0.2,
                                 trymax=20,
                                 maxit=20) 

exampleearly_NMDSk3_B <- metaMDS(early_comm_mat, # Our community-by-species matrix
                                 k=3, # The number of reduced dimensions
                                 na.rm=TRUE, autotransform=FALSE,
                                 trymax=3,
                                 maxit=3) 



#OK what about MDS
d <- dist(early_comm_mat ) # euclidean distances between the rows
fit1 <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim 
#IMMEDIATELY CRASHED LAPTOP
#second attempt about A half hour to run
fit1 # view results
saveRDS(fit1, file="scripts/MDS_early_output.RDS")

rownames(early_comm_mat) <- paste(early_wide$CRUISE, early_wide$HAUL, sep="-")

x1 <- fit1$points[,1]
y1 <- fit1$points[,2]
plot(x1, y1,  xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS", type="n")
text(x1, y1, labels = row.names(early_comm_mat), cex=.7)
#perhaps paste rownames split into year and haul and colour by year
mat2plot <- as.data.frame(early_comm_mat)
mat2plot$rowname <- paste(rownames(mat2plot))
mat2plot$year <- str_sub(mat2plot$rowname, 2, 5)

text(x1, y1, labels = row.names(mat2plot), cex=.7, col=mat2plot$year)

fitdf <- as.data.frame(fit1$points)

ggplot(fitdf, aes(V1, V2)) + geom_point() 

#try with short_early_mat
d2 <- dist(short_early_mat) # euclidean distances between the rows
fit1 <- cmdscale(d2,eig=TRUE, k=2) # k is the number of dim
fit1 # view results

short_early_mat3 <- early_comm_mat[1:3000,]
d3 <- dist(short_early_mat3) # euclidean distances between the rows
fit3 <- cmdscale(d3,eig=TRUE, k=2) # k is the number of dim
fit3 # view results - runs fine and under a minute!

short_early_mat5 <- early_comm_mat[1:5000,]
d5 <- dist(short_early_mat5) # euclidean distances between the rows
fit5 <- cmdscale(d5,eig=TRUE, k=2) # k is the number of dim, 3 minute run time
fit5

#is the problem at the end? Running in R says there are NAs in the matrix
11833
short_end <- early_comm_mat[5000:11833,]
de <- dist(short_end) # euclidean distances between the rows
fite <- cmdscale(de,eig=TRUE, k=2) # k is the number of dim, 3 minute run time
fite


#NMDS late only====

late_only <- full_wide[which(full_wide$YEAR>2016),]
names(late_only)
#late_mat <- as.matrix(late_only[,c(38:48)])
#late_mat <- as.matrix(late_only[,c(38:47)]) #Lepidopsetta_sp is all zeros in late period!
#also try without crabs which have a very high number of zeros and few non zeros
late_mat <- as.matrix(late_only[,c(38:45)]) 

late_NMDS <- metaMDS(late_mat, # Our community-by-species matrix
                     distance = "bray",
                             k=3, # The number of reduced dimensions
                             na.rm=TRUE,
                             noshare=TRUE,
                             trymax=200,
                             maxit=200)#,
                    # previous.best = late_NMDS) 

summary(late_NMDS)
late_NMDS$stress
lateptsall<-scores(late_NMDS)
scores(lateptsall, display="species")
ordiplot(late_NMDS)

ordiplot(late_NMDS,type="n")
orditorp(late_NMDS,display="sites",cex=1.25,air=0.01)
orditorp(late_NMDS,display="species",col="red",air=0.01)


late_highk <- metaMDS(late_mat, # Our community-by-species matrix
                     distance = "bray",
                     k=6, # The number of reduced dimensions
                     na.rm=TRUE,
                     noshare=TRUE,
                     trymax=50,
                     maxit=50) 

ordiplot(late_highk,type="n")
orditorp(late_highk,display="sites",cex=1.25,air=0.01)
orditorp(late_highk,display="species",col="red",air=0.01)

#try with only stations with pollock
pol_only <- late_only[which(late_only$logCPUE_Gadus_chalcogrammus>0),c(38:45)]
pol_mat <- as.matrix(pol_only) 

pol_highk <- metaMDS(pol_mat, # Our community-by-species matrix
                      distance = "bray",
                      k=6, # The number of reduced dimensions
                      na.rm=TRUE,
                      noshare=TRUE,
                      trymax=100,
                      maxit=100) 

ordiplot(late_highk,type="n")
orditorp(late_highk,display="sites",cex=1.25,air=0.01)
orditorp(late_highk,display="species",col="red",air=0.01)


#MDS with full data====

# full_wide <- as.data.frame(full_wide)
# summary(full_wide) 

full_comm_mat <- as.matrix(full_wide[,c(38:48)])

wd <- getwd()
write.csv(full_comm_mat, file=paste(wd,"/data/community_data_matrix.csv", sep=""))


#write to csv to run in R itself
#write.csv(full_comm_mat, file="/Users/krista/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/data/full_comm_mat.csv")

# dfull <- dist(full_comm_mat) # euclidean distances between the rows
# fitfull <- cmdscale(dfull,eig=TRUE, k=2) # k is the number of dim 
# #start1:13
# fitfull # view results
# saveRDS(fitfull, file="scripts/MDS_full_output.RDS")

#Mike ran these! Read in the RDS instead. See 'MDS_on_cluster.R' for running
fitfull <- readRDS(file="scripts/MDS_full_output.RDS")

rownames(early_comm_mat) <- paste(early_wide$CRUISE, early_wide$HAUL, sep="-")

x1 <- fitfull$points[,1]
y1 <- fitfull$points[,2]
plot(x1, y1,  xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS", type="n")
text(x1, y1, labels = row.names(full_comm_mat), cex=.7)
#perhaps paste rownames split into year and haul and colour by year
mat2plot <- as.data.frame(full_comm_mat)
mat2plot$rowname <- paste(rownames(mat2plot))
mat2plot$year <- str_sub(mat2plot$rowname, 2, 5)

text(x1, y1, labels = row.names(mat2plot), cex=.7, col=mat2plot$year)

text(x1, y1, labels = row.names(full_comm_mat), cex=.7, col=full_wide$YEAR)

text(x1, y1, labels = full_wide$YEAR, cex=.7, col=full_wide$YEAR)

full_wide$period_num <- NA
full_wide$period_num[which(full_wide$YEAR<2014)]<-1
full_wide$period_num[which(full_wide$YEAR>2013)]<-2
full_wide$period <- as.factor(full_wide$period)

text(x1, y1, labels = full_wide$YEAR, cex=.7, col=full_wide$period_num)
text(x1, y1, labels = full_wide$YEAR, cex=.7, col=full_wide$STRATUM)
text(x1, y1, labels = full_wide$STRATUM, cex=.7, col=full_wide$YEAR)
text(x1, y1, labels = full_wide$STRATUM, cex=.7, col=full_wide$STRATUM)





cov(short_early_mat)
cov(full_comm_mat)

pcatest <- prcomp(short_early_mat)



#plot climate variables========
# load environmental data
climdat <- read.csv("data/climate data.csv")

wide_join <- left_join(early_wide, climdat, by=c("YEAR" = "year"))

p2 <- ggplot(wide_join, aes(summer.cold.pool.extent, logCPUE_Gadus_chalcogrammus))
p2 + geom_smooth() + geom_point()

p3 <- ggplot(wide_join, aes(AO.jfm, logCPUE_Gadus_chalcogrammus))
p3 + geom_smooth() + geom_point()

p4 <- ggplot(wide_join, aes(south.sst.ndjfm, logCPUE_Gadus_chalcogrammus))
p4 + geom_smooth() + geom_point()

p5 <- ggplot(wide_join, aes(south.sst.amj, logCPUE_Gadus_chalcogrammus))
p5 + geom_smooth() + geom_point()







#cor models w both periods========================================================================================

#both periods inclusion criteria

both_dat <- sel.trawl.dat
#both_dat <- both_dat[,-c(9,11)] #drop columns that repeat info
both_dat <- both_dat[!duplicated(both_dat),] #few duplicated rows, remove them

#widen dataframe to have columns for each sps

both_wide <-both_dat %>% pivot_wider(names_from=SCIENTIFIC, 
                                        values_from=c(WTCPUE, NUMCPUE, logCPUE)) #
#eek this causes some bad column names!
both_wide <- both_wide %>% rename(WTCPUE_Chionoecetes_bairdi = "WTCPUE_Chionoecetes bairdi",            
                                    WTCPUE_Atheresthes_stomia = "WTCPUE_Atheresthes stomias",             
                                    WTCPUE_Hippoglossus_stenolepis = "WTCPUE_Hippoglossus stenolepis",
                                    WTCPUE_Limanda_aspera = "WTCPUE_Limanda aspera",         
                                    WTCPUE_Lepidopsetta_sp= "WTCPUE_Lepidopsetta sp.",               
                                    WTCPUE_Chionoecetes_opilio = "WTCPUE_Chionoecetes opilio",             
                                    WTCPUE_Gadus_macrocephalus = "WTCPUE_Gadus macrocephalus",            
                                    WTCPUE_Hippoglossoides_elassodon = "WTCPUE_Hippoglossoides elassodon",  
                                    WTCPUE_Pleuronectes_quadrituberculatus = "WTCPUE_Pleuronectes quadrituberculatus", 
                                    WTCPUE_Lepidopsetta_polyxystra = "WTCPUE_Lepidopsetta polyxystra", 
                                    WTCPUE_Gadus_chalcogrammus = "WTCPUE_Gadus chalcogrammus", 
                                    NUMCPUE_Gadus_chalcogrammus = "NUMCPUE_Gadus chalcogrammus",        
                                    NUMCPUE_Chionoecetes_bairdi = "NUMCPUE_Chionoecetes bairdi",            
                                    NUMCPUE_Atheresthes_stomias = "NUMCPUE_Atheresthes stomias",           
                                    NUMCPUE_Hippoglossus_stenolepis = "NUMCPUE_Hippoglossus stenolepis",        
                                    NUMCPUE_Limanda_aspera = "NUMCPUE_Limanda aspera",                
                                    NUMCPUE_Lepidopsetta_sp = "NUMCPUE_Lepidopsetta sp.",               
                                    NUMCPUE_Chionoecetes_opilio = "NUMCPUE_Chionoecetes opilio",           
                                    NUMCPUE_Gadus_macrocephalus = "NUMCPUE_Gadus macrocephalus",            
                                    NUMCPUE_Hippoglossoides_elassodon = "NUMCPUE_Hippoglossoides elassodon",     
                                    NUMCPUE_Pleuronectes_quadrituberculatus = "NUMCPUE_Pleuronectes quadrituberculatus",
                                    NUMCPUE_Lepidopsetta_polyxystra = "NUMCPUE_Lepidopsetta polyxystra",
                                    logCPUE_Gadus_chalcogrammus = "logCPUE_Gadus chalcogrammus",        
                                    logCPUE_Chionoecetes_bairdi = "logCPUE_Chionoecetes bairdi",           
                                    logCPUE_Atheresthes_stomias = "logCPUE_Atheresthes stomias",            
                                    logCPUE_Hippoglossus_stenolepis = "logCPUE_Hippoglossus stenolepis", 
                                    logCPUE_Limanda_aspera = "logCPUE_Limanda aspera",        
                                    logCPUE_Lepidopsetta_sp = "logCPUE_Lepidopsetta sp.",              
                                    logCPUE_Chionoecetes_opilio = "logCPUE_Chionoecetes opilio",            
                                    logCPUE_Gadus_macrocephalus = "logCPUE_Gadus macrocephalus",           
                                    logCPUE_Hippoglossoides_elassodon = "logCPUE_Hippoglossoides elassodon",      
                                    logCPUE_Pleuronectes_quadrituberculatus = "logCPUE_Pleuronectes quadrituberculatus" ,
                                    logCPUE_Lepidopsetta_polyxystra = "logCPUE_Lepidopsetta polyxystra") 

station_summary_both <- both_wide %>% group_by(STATION, LONGITUDE, LATITUDE) %>%
  summarize(n_yrs=n())

station_summary2_both <- both_wide %>% group_by(STATION) %>%
  summarize(n_yrs=n())

joinboth <- left_join(both_wide, station_summary2_both)

#periods_dat <- joinboth[which(joinboth$n_yrs>5),]

periods_dat <- joinboth

periods_dat$period <- NA

periods_dat$period[which(periods_dat$YEAR<2014)] <- "early"
periods_dat$period[which(periods_dat$YEAR>2013)] <- "late"

#deal w missing data=================

table(periods_dat$YEAR[which(is.na(periods_dat$BOT_TEMP)==TRUE)])
#many years have some NAs in bottom temp, a LOT in 1994

table(periods_dat$YEAR[which(is.na(periods_dat$BOT_TEMP)==FALSE)])
#but far far more with data, even in 94

#OK I should maybe loop through but I want to check these by eye so I'll do them by hand

tempmod82 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=8), data=periods_dat[which(periods_dat$YEAR==1982),] )
summary(tempmod82)
gam.check(tempmod82)
plot(tempmod82)

pre82 <- predict.gam(tempmod82, type="response")
df82 <- periods_dat[which(periods_dat$YEAR==1982),] 

missing82 <- df82[which(is.na(df82$BOT_TEMP)==TRUE),]
pred82 <- predict.gam(tempmod82, newdata=missing82, type="response")
missing82$BOT_TEMP <- pred82

tempmod83 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=8), data=periods_dat[which(periods_dat$YEAR==1983),] )
summary(tempmod83)
gam.check(tempmod83)
plot(tempmod83)

pre83 <- predict.gam(tempmod83, type="response")
df83 <- periods_dat[which(periods_dat$YEAR==1983),] 

missing83 <- df83[which(is.na(df83$BOT_TEMP)==TRUE),]
pred83 <- predict.gam(tempmod83, newdata=missing83, type="response")
missing83$BOT_TEMP <- pred83


tempmod84 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=8), data=periods_dat[which(periods_dat$YEAR==1984),] )
summary(tempmod84)
gam.check(tempmod84)
plot(tempmod84)

pre84 <- predict.gam(tempmod84, type="response")
df84 <- periods_dat[which(periods_dat$YEAR==1984),] 

missing84 <- df84[which(is.na(df84$BOT_TEMP)==TRUE),]
pred84 <- predict.gam(tempmod84, newdata=missing84, type="response")
missing84$BOT_TEMP <- pred84


tempmod85 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1985),] )
summary(tempmod85)
gam.check(tempmod85)
plot(tempmod85)

pre85 <- predict.gam(tempmod85, type="response")
df85 <- periods_dat[which(periods_dat$YEAR==1985),] 

missing85 <- df85[which(is.na(df85$BOT_TEMP)==TRUE),]
pred85 <- predict.gam(tempmod85, newdata=missing85, type="response")
missing85$BOT_TEMP <- pred85


tempmod86 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1986),] )
summary(tempmod86)
gam.check(tempmod86)
plot(tempmod86)

pre86 <- predict.gam(tempmod86, type="response")
df86 <- periods_dat[which(periods_dat$YEAR==1986),] 

missing86 <- df86[which(is.na(df86$BOT_TEMP)==TRUE),]
pred86 <- predict.gam(tempmod86, newdata=missing86, type="response")
missing86$BOT_TEMP <- pred86


tempmod87 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=8), data=periods_dat[which(periods_dat$YEAR==1987),] )
summary(tempmod87)
gam.check(tempmod87)
plot(tempmod87)

pre87 <- predict.gam(tempmod87, type="response")
df87 <- periods_dat[which(periods_dat$YEAR==1987),] 

missing87 <- df87[which(is.na(df87$BOT_TEMP)==TRUE),]
pred87 <- predict.gam(tempmod87, newdata=missing87, type="response")
missing87$BOT_TEMP <- pred87



tempmod88 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=8), data=periods_dat[which(periods_dat$YEAR==1988),] )
summary(tempmod88)
gam.check(tempmod88)
plot(tempmod88)

pre88 <- predict.gam(tempmod88, type="response")
df88 <- periods_dat[which(periods_dat$YEAR==1988),] 

missing88 <- df88[which(is.na(df88$BOT_TEMP)==TRUE),]
pred88 <- predict.gam(tempmod88, newdata=missing88, type="response")
missing88$BOT_TEMP <- pred88



tempmod89 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1989),] )
summary(tempmod89)
gam.check(tempmod89)
plot(tempmod89)

pre89 <- predict.gam(tempmod89, type="response")
df89 <- periods_dat[which(periods_dat$YEAR==1989),] 

missing89 <- df89[which(is.na(df89$BOT_TEMP)==TRUE),]
pred89 <- predict.gam(tempmod89, newdata=missing89, type="response")
missing89$BOT_TEMP <- pred89



tempmod90 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1990),] )
summary(tempmod90)
gam.check(tempmod90)
plot(tempmod90)

pre90 <- predict.gam(tempmod90, type="response")
df90 <- periods_dat[which(periods_dat$YEAR==1990),] 

missing90 <- df90[which(is.na(df90$BOT_TEMP)==TRUE),]
pred90 <- predict.gam(tempmod90, newdata=missing90, type="response")
missing90$BOT_TEMP <- pred90



tempmod91 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1991),]  )
summary(tempmod91)
gam.check(tempmod91)
plot(tempmod91)

pre91 <- predict.gam(tempmod91, type="response")
df91 <- periods_dat[which(periods_dat$YEAR==1991),] 

missing91 <- df91[which(is.na(df91$BOT_TEMP)==TRUE),]
pred91 <- predict.gam(tempmod91, newdata=missing91, type="response")
missing91$BOT_TEMP <- pred91



tempmod92 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1992),]  )
summary(tempmod92)
gam.check(tempmod92)
plot(tempmod92)

pre92 <- predict.gam(tempmod92, type="response")
df92 <- periods_dat[which(periods_dat$YEAR==1992),] 

missing92 <- df92[which(is.na(df92$BOT_TEMP)==TRUE),]
pred92 <- predict.gam(tempmod92, newdata=missing92, type="response")
missing92$BOT_TEMP <- pred92



tempmod93 <- gam(BOT_TEMP ~ s(BOT_DEPTH, k=17) + ti(long_albers, lat_albers, k=8), data=periods_dat[which(periods_dat$YEAR==1993),]  )
summary(tempmod93)
gam.check(tempmod93) #seems awfully high
plot(tempmod93) 

pre93 <- predict.gam(tempmod93, type="response")
df93 <- periods_dat[which(periods_dat$YEAR==1993),] 

missing93 <- df93[which(is.na(df93$BOT_TEMP)==TRUE),]
pred93 <- predict.gam(tempmod93, newdata=missing93, type="response")
missing93$BOT_TEMP <- pred93



tempmod94 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=7), data=periods_dat[which(periods_dat$YEAR==1994),]  )
summary(tempmod94)
gam.check(tempmod94)
plot(tempmod94)

pre94 <- predict.gam(tempmod94, type="response")
df94 <- periods_dat[which(periods_dat$YEAR==1994),] 

missing94 <- df94[which(is.na(df94$BOT_TEMP)==TRUE),]
pred94 <- predict.gam(tempmod94, newdata=missing94, type="response")
missing94$BOT_TEMP <- pred94



tempmod95 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1995),]  )
summary(tempmod95)
gam.check(tempmod95)
plot(tempmod95)

pre95 <- predict.gam(tempmod95, type="response")
df95 <- periods_dat[which(periods_dat$YEAR==1995),] 

missing95 <- df95[which(is.na(df95$BOT_TEMP)==TRUE),]
pred95 <- predict.gam(tempmod95, newdata=missing95, type="response")
missing95$BOT_TEMP <- pred95



tempmod96 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1996),]  )
summary(tempmod96)
gam.check(tempmod96)
plot(tempmod96)

pre96 <- predict.gam(tempmod96, type="response")
df96 <- periods_dat[which(periods_dat$YEAR==1996),] 

missing96 <- df96[which(is.na(df96$BOT_TEMP)==TRUE),]
pred96 <- predict.gam(tempmod96, newdata=missing96, type="response")
missing96$BOT_TEMP <- pred96



tempmod97 <- gam(BOT_TEMP ~ s(BOT_DEPTH, k=8) + ti(long_albers, lat_albers, k=11), data=periods_dat[which(periods_dat$YEAR==1997),]  )
summary(tempmod97)
gam.check(tempmod97)
plot(tempmod97)

pre97 <- predict.gam(tempmod97, type="response")
df97 <- periods_dat[which(periods_dat$YEAR==1997),] 

missing97 <- df97[which(is.na(df97$BOT_TEMP)==TRUE),]
pred97 <- predict.gam(tempmod97, newdata=missing97, type="response")
missing97$BOT_TEMP <- pred97



tempmod98 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1998),]  )
summary(tempmod98)
gam.check(tempmod98)
plot(tempmod98)

pre98 <- predict.gam(tempmod98, type="response")
df98 <- periods_dat[which(periods_dat$YEAR==1998),] 

missing98 <- df98[which(is.na(df98$BOT_TEMP)==TRUE),]
pred98 <- predict.gam(tempmod98, newdata=missing98, type="response")
missing98$BOT_TEMP <- pred98



tempmod99 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==1999),]  )
summary(tempmod99)
gam.check(tempmod99)
plot(tempmod99)

pre99 <- predict.gam(tempmod99, type="response")
df99 <- periods_dat[which(periods_dat$YEAR==1999),] 

missing99 <- df99[which(is.na(df99$BOT_TEMP)==TRUE),]
pred99 <- predict.gam(tempmod99, newdata=missing99, type="response")
missing99$BOT_TEMP <- pred99



tempmod00 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=7), data=periods_dat[which(periods_dat$YEAR==2000),]  )
summary(tempmod00)
gam.check(tempmod00)
plot(tempmod00)

pre00 <- predict.gam(tempmod00, type="response")
df00 <- periods_dat[which(periods_dat$YEAR==2000),] 

missing00 <- df00[which(is.na(df00$BOT_TEMP)==TRUE),]
pred00 <- predict.gam(tempmod00, newdata=missing00, type="response")
missing00$BOT_TEMP <- pred00



tempmod01 <- gam(BOT_TEMP ~ s(BOT_DEPTH, k=16) + ti(long_albers, lat_albers, k=11), data=periods_dat[which(periods_dat$YEAR==2001),]  )
summary(tempmod01)
gam.check(tempmod01) #still a little low but close enough
plot(tempmod01)

pre01 <- predict.gam(tempmod01, type="response")
df01 <- periods_dat[which(periods_dat$YEAR==2001),] 

missing01 <- df01[which(is.na(df01$BOT_TEMP)==TRUE),]
pred01 <- predict.gam(tempmod01, newdata=missing01, type="response")
missing01$BOT_TEMP <- pred01



tempmod02 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=8), data=periods_dat[which(periods_dat$YEAR==2002),]  )
summary(tempmod02)
gam.check(tempmod02)
plot(tempmod02)

pre02 <- predict.gam(tempmod02, type="response")
df02 <- periods_dat[which(periods_dat$YEAR==2002),] 

missing02 <- df02[which(is.na(df02$BOT_TEMP)==TRUE),]
pred02 <- predict.gam(tempmod02, newdata=missing02, type="response")
missing02$BOT_TEMP <- pred02



tempmod03 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=10), data=periods_dat[which(periods_dat$YEAR==2003),]  )
summary(tempmod03)
gam.check(tempmod03)
plot(tempmod03)

pre03 <- predict.gam(tempmod03, type="response")
df03 <- periods_dat[which(periods_dat$YEAR==2003),] 

missing03 <- df03[which(is.na(df03$BOT_TEMP)==TRUE),]
pred03 <- predict.gam(tempmod03, newdata=missing03, type="response")
missing03$BOT_TEMP <- pred03



tempmod04 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2004),]  )
summary(tempmod04)
gam.check(tempmod04)
plot(tempmod04)

pre04 <- predict.gam(tempmod04, type="response")
df04 <- periods_dat[which(periods_dat$YEAR==2004),] 

missing04 <- df04[which(is.na(df04$BOT_TEMP)==TRUE),]
pred04 <- predict.gam(tempmod04, newdata=missing04, type="response")
missing04$BOT_TEMP <- pred04



tempmod05 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2005),]  )
summary(tempmod05)
gam.check(tempmod05)
plot(tempmod05)

pre05 <- predict.gam(tempmod05, type="response")
df05 <- periods_dat[which(periods_dat$YEAR==2005),] 

missing05 <- df05[which(is.na(df05$BOT_TEMP)==TRUE),]
pred05 <- predict.gam(tempmod05, newdata=missing05, type="response")
missing05$BOT_TEMP <- pred05



tempmod06 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=10), data=periods_dat[which(periods_dat$YEAR==2006),]  )
summary(tempmod06)
gam.check(tempmod06)
plot(tempmod06)

pre06 <- predict.gam(tempmod06, type="response")
df06 <- periods_dat[which(periods_dat$YEAR==2006),] 

missing06 <- df06[which(is.na(df06$BOT_TEMP)==TRUE),]
pred06 <- predict.gam(tempmod06, newdata=missing06, type="response")
missing06$BOT_TEMP <- pred06



tempmod07 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2007),]  )
summary(tempmod07)
gam.check(tempmod07)
plot(tempmod07)

pre07 <- predict.gam(tempmod07, type="response")
df07 <- periods_dat[which(periods_dat$YEAR==2007),] 

missing07 <- df07[which(is.na(df07$BOT_TEMP)==TRUE),]
pred07 <- predict.gam(tempmod07, newdata=missing07, type="response")
missing07$BOT_TEMP <- pred07



tempmod08 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2008),]  )
summary(tempmod08)
gam.check(tempmod08)
plot(tempmod08)

pre08 <- predict.gam(tempmod08, type="response")
df08 <- periods_dat[which(periods_dat$YEAR==2008),] 

missing08 <- df08[which(is.na(df08$BOT_TEMP)==TRUE),]
pred08 <- predict.gam(tempmod08, newdata=missing08, type="response")
missing08$BOT_TEMP <- pred08



tempmod09 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=10), data=periods_dat[which(periods_dat$YEAR==2009),]  )
summary(tempmod09)
gam.check(tempmod09)
plot(tempmod09)

pre09 <- predict.gam(tempmod09, type="response")
df09 <- periods_dat[which(periods_dat$YEAR==2009),] 

missing09 <- df09[which(is.na(df09$BOT_TEMP)==TRUE),]
pred09 <- predict.gam(tempmod09, newdata=missing09, type="response")
missing09$BOT_TEMP <- pred09



tempmod10 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=10), data=periods_dat[which(periods_dat$YEAR==2010),]  )
summary(tempmod10)
gam.check(tempmod10)
plot(tempmod10)

pre10 <- predict.gam(tempmod10, type="response")
df10 <- periods_dat[which(periods_dat$YEAR==2010),] 

missing10 <- df10[which(is.na(df10$BOT_TEMP)==TRUE),]
pred10 <- predict.gam(tempmod10, newdata=missing10, type="response")
missing10$BOT_TEMP <- pred10



tempmod11 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=10), data=periods_dat[which(periods_dat$YEAR==2011),] )
summary(tempmod11)
gam.check(tempmod11)
plot(tempmod11)

pre11 <- predict.gam(tempmod11, type="response")
df11 <- periods_dat[which(periods_dat$YEAR==2011),] 

missing11 <- df11[which(is.na(df11$BOT_TEMP)==TRUE),]
pred11 <- predict.gam(tempmod11, newdata=missing11, type="response")
missing11$BOT_TEMP <- pred11



tempmod12 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=11), data=periods_dat[which(periods_dat$YEAR==2012),] )
summary(tempmod12)
gam.check(tempmod12)
plot(tempmod12)

pre12 <- predict.gam(tempmod12, type="response")
df12 <- periods_dat[which(periods_dat$YEAR==2012),] 

missing12 <- df12[which(is.na(df12$BOT_TEMP)==TRUE),]
pred12 <- predict.gam(tempmod12, newdata=missing12, type="response")
missing12$BOT_TEMP <- pred12



tempmod13 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2013),] )
summary(tempmod13)
gam.check(tempmod13)
plot(tempmod13)

pre13 <- predict.gam(tempmod13, type="response")
df13 <- periods_dat[which(periods_dat$YEAR==2013),] 

missing13 <- df13[which(is.na(df13$BOT_TEMP)==TRUE),]
pred13 <- predict.gam(tempmod13, newdata=missing13, type="response")
missing13$BOT_TEMP <- pred13



tempmod14 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2014),] )
summary(tempmod14)
gam.check(tempmod14)
plot(tempmod14)

pre14 <- predict.gam(tempmod14, type="response")
df14 <- periods_dat[which(periods_dat$YEAR==2014),] 

missing14 <- df14[which(is.na(df14$BOT_TEMP)==TRUE),]
pred14 <- predict.gam(tempmod14, newdata=missing14, type="response")
missing14$BOT_TEMP <- pred14



tempmod15 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=10), data=periods_dat[which(periods_dat$YEAR==2015),] )
summary(tempmod15)
gam.check(tempmod15)
plot(tempmod15)

pre15 <- predict.gam(tempmod15, type="response")
df15 <- periods_dat[which(periods_dat$YEAR==2015),] 

missing15 <- df15[which(is.na(df15$BOT_TEMP)==TRUE),]
pred15 <- predict.gam(tempmod15, newdata=missing15, type="response")
missing15$BOT_TEMP <- pred15



tempmod16 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2016),] )
summary(tempmod16)
gam.check(tempmod16)
plot(tempmod16)

pre16 <- predict.gam(tempmod16, type="response")
df16 <- periods_dat[which(periods_dat$YEAR==2016),] 

missing16 <- df16[which(is.na(df16$BOT_TEMP)==TRUE),]
pred16 <- predict.gam(tempmod16, newdata=missing16, type="response")
missing16$BOT_TEMP <- pred16



tempmod17 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2017),] )
summary(tempmod17)
gam.check(tempmod17)
plot(tempmod17)

pre17 <- predict.gam(tempmod17, type="response")
df17 <- periods_dat[which(periods_dat$YEAR==2017),] 

missing17 <- df17[which(is.na(df17$BOT_TEMP)==TRUE),]
pred17 <- predict.gam(tempmod17, newdata=missing17, type="response")
missing17$BOT_TEMP <- pred17



tempmod18 <- gam(BOT_TEMP ~ s(BOT_DEPTH, k=16) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2018),] )
summary(tempmod18)
gam.check(tempmod18) #still a little low on depth but close enough
plot(tempmod18)

pre18 <- predict.gam(tempmod18, type="response")
df18 <- periods_dat[which(periods_dat$YEAR==2018),] 

missing18 <- df18[which(is.na(df18$BOT_TEMP)==TRUE),]
pred18 <- predict.gam(tempmod18, newdata=missing18, type="response")
missing18$BOT_TEMP <- pred18



tempmod19 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers, k=9), data=periods_dat[which(periods_dat$YEAR==2019),] )
summary(tempmod19)
gam.check(tempmod19)
plot(tempmod19)

pre19 <- predict.gam(tempmod19, type="response")
df19 <- periods_dat[which(periods_dat$YEAR==2019),] 

missing19 <- df19[which(is.na(df19$BOT_TEMP)==TRUE),]
pred19 <- predict.gam(tempmod19, newdata=missing19, type="response")
missing19$BOT_TEMP <- pred19



#now bind together

missingall <- rbind(missing82, missing83, missing84, missing85, missing86, missing87, missing88, missing89,
    missing90, missing91, missing92, missing93, missing94, missing95,
    missing96, missing97, missing98, missing99, missing00,
    missing01, missing02, missing03, missing04, missing05,
    missing06, missing07, missing08, missing09, missing10, 
    missing11, missing12, missing13, missing14, missing15,
    missing16, missing17, missing18, missing19)


nona <- periods_dat[which(is.na(periods_dat$BOT_TEMP)==FALSE),]

filled <- rbind(nona, missingall)


#calculate anomalies and means=====================================================================

#annual 'global' (w/in dataset) mean
btempmeans <- filled %>% group_by(STATION) %>% summarize(mean_station_bottemp=mean(BOT_TEMP))
filled <- left_join(filled, btempmeans)

filled$bottemp_anom <- filled$BOT_TEMP - filled$mean_station_bottemp

periods_analysis_dat <- filled

#exclusion criteria again===============

#no longer doign this because it would drop NBS stations
# 
# station_summary3 <- periods_analysis_dat %>% group_by(STATION) %>%
#   summarize(n_yrs=n())
# 
# join2 <- left_join(periods_analysis_dat, station_summary3)
# 
# periods_analysis_dat <- join2[which(join2$n_yrs>5),]

#let's save the processed data
wd <- getwd()
write.csv(periods_analysis_dat, file=paste(wd,"/data/processed_periods_analysis_data.csv", sep=""))

#===***===***===

#Maps=====

library("rnaturalearth")
library("rnaturalearthdata")
library( "ggspatial" )
library("sf")

world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 63), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                        style = north_arrow_fancy_orienteering) +  
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Gadus_chalcogrammus)), data=periods_analysis_dat) +   
  scale_colour_gradient2(low="blue", high="red", guide="colorbar") + facet_wrap(~period) 

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 63), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=logCPUE_Gadus_chalcogrammus), data=periods_analysis_dat) + 
  facet_wrap(~period) +
  scale_colour_gradient2(low="blue", high="red")

#also depth/temp

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 63), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=BOT_DEPTH), data=periods_analysis_dat) + 
  scale_colour_gradient2(low="blue", high="red")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 63), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=periods_analysis_dat) + 
  scale_colour_gradient2(low="blue", high="red")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-180, -155), ylim = c(53, 63), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp, size=BOT_DEPTH), data=periods_analysis_dat) + 
  scale_colour_gradient2(low="blue", high="red")





