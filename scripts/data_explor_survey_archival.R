#============================================================================================================================================
# DATA EXPLORATION - survey abundance data

#Created by Krista, May 2020
#============================================================================================================================================
#Notes:
#============================================================================================================================================
library(ggplot2)
library(tidyverse)
library(corrplot)

#for now using trawl.data from 'trawl data processing.R' just to get a feel for the data
summary(trawl.data)


unique(trawl.data$COMMON)
unique(trawl.data$SCIENTIFIC) #WOW a lot of sps!!

pol.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Gadus chalcogrammus"),]
pol.trawl.dat$logCPUE <- log(pol.trawl.dat$WTCPUE + 1)

ggplot(pol.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")

table(pol.trawl.dat$STATION, pol.trawl.dat$YEAR)
length(unique(pol.trawl.dat$STATION))

ggplot(pol.trawl.dat, aes(YEAR, WTCPUE, colour=STATION)) + geom_point() + theme(legend.position = "none")

ggplot(pol.trawl.dat, aes(YEAR, WTCPUE, colour=STATION)) + geom_line() + theme(legend.position = "none")

ggplot(pol.trawl.dat, aes(YEAR, WTCPUE, colour=as.factor(STRATUM))) + geom_line() + theme(legend.position = "none")

ggplot(pol.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)





ggplot(pol.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(pol.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

pol.trawl.dat$fourthrootCPUE <- pol.trawl.dat$WTCPUE^(1/4)

ggplot(pol.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail

dat1 <- as.matrix(pol.trawl.dat[,c(7:8,12:19)])
cor1 <- cor(dat1, use="complete.obs")
corrplot(cor1)


#other sps====================================================================================

#Arctic plaice====
arcp.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Pleuronectes quadrituberculatus"),]

ggplot(arcp.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")


ggplot(arcp.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)


arcp.trawl.dat$logCPUE <- log(arcp.trawl.dat$WTCPUE + 1)

ggplot(arcp.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(arcp.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

arcp.trawl.dat$fourthrootCPUE <- arcp.trawl.dat$WTCPUE^(1/4)

ggplot(arcp.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail



#Snow crab====
snoc.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Chionoecetes opilio"),]

ggplot(snoc.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")


ggplot(snoc.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)


snoc.trawl.dat$logCPUE <- log(snoc.trawl.dat$WTCPUE + 1)

ggplot(snoc.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(snoc.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

snoc.trawl.dat$fourthrootCPUE <- snoc.trawl.dat$WTCPUE^(1/4)

ggplot(snoc.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail


#Arrowtooth====
#TWO species names!
# Atheresthes stomias --- I think it's this one
# Atheresthes evermanni
arot.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Atheresthes stomias"),]

ggplot(arot.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")


ggplot(arot.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)


arot.trawl.dat$logCPUE <- log(arot.trawl.dat$WTCPUE + 1)

ggplot(arot.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(arot.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

arot.trawl.dat$fourthrootCPUE <- arot.trawl.dat$WTCPUE^(1/4)

ggplot(arot.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail



#Flathead sole====
flhs.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Hippoglossoides elassodon"),]

ggplot(flhs.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")


ggplot(flhs.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)


flhs.trawl.dat$logCPUE <- log(flhs.trawl.dat$WTCPUE + 1)

ggplot(flhs.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(flhs.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

flhs.trawl.dat$fourthrootCPUE <- flhs.trawl.dat$WTCPUE^(1/4)

ggplot(flhs.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail




#P cod====
pcod.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Gadus macrocephalus"),]

ggplot(pcod.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")


ggplot(pcod.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)


pcod.trawl.dat$logCPUE <- log(pcod.trawl.dat$WTCPUE + 1)

ggplot(pcod.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(pcod.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

pcod.trawl.dat$fourthrootCPUE <- pcod.trawl.dat$WTCPUE^(1/4)

ggplot(pcod.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail



#Rock sole====
#THREE species names!
# Lepidopsetta polyxystra ----- Northern rock sole, newer
# Lepidopsetta bilineata ---- I think it's this one but almost no data!
# Lepidopsetta sp. ---- older
rcks.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Lepidopsetta polyxystra"|
                                     trawl.data$SCIENTIFIC=="Lepidopsetta sp." ),]

ggplot(rcks.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")


ggplot(rcks.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)


rcks.trawl.dat$logCPUE <- log(rcks.trawl.dat$WTCPUE + 1)

ggplot(rcks.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(rcks.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

rcks.trawl.dat$fourthrootCPUE <- rcks.trawl.dat$WTCPUE^(1/4)

ggplot(rcks.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail





#P cod====
yfns.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Limanda aspera"),]

ggplot(yfns.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")


ggplot(yfns.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)


yfns.trawl.dat$logCPUE <- log(yfns.trawl.dat$WTCPUE + 1)

ggplot(yfns.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(yfns.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

yfns.trawl.dat$fourthrootCPUE <- yfns.trawl.dat$WTCPUE^(1/4)

ggplot(yfns.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail



#Halibut====
halb.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Hippoglossus stenolepis"),]

ggplot(halb.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")


ggplot(halb.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)


halb.trawl.dat$logCPUE <- log(halb.trawl.dat$WTCPUE + 1)

ggplot(halb.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(halb.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

halb.trawl.dat$fourthrootCPUE <- halb.trawl.dat$WTCPUE^(1/4)

ggplot(halb.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail




#Tanner crab====
tnrc.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Chionoecetes bairdi"),]

ggplot(tnrc.trawl.dat, aes(WTCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free")


ggplot(tnrc.trawl.dat, aes(YEAR, WTCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)


tnrc.trawl.dat$logCPUE <- log(tnrc.trawl.dat$WTCPUE + 1)

ggplot(tnrc.trawl.dat, aes(YEAR, logCPUE)) + geom_point() + theme(legend.position = "none") +
  geom_smooth() + facet_wrap(~STRATUM)

ggplot(tnrc.trawl.dat, aes(logCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #little better to my eye

tnrc.trawl.dat$fourthrootCPUE <- tnrc.trawl.dat$WTCPUE^(1/4)

ggplot(tnrc.trawl.dat, aes(fourthrootCPUE)) + geom_histogram() + facet_wrap(~YEAR, scales="free") #long tail




#gather up all sps=============


sel.trawl.dat <- trawl.data[which(trawl.data$SCIENTIFIC=="Gadus chalcogrammus"|trawl.data$SCIENTIFIC=="Pleuronectes quadrituberculatus"|
                                  trawl.data$SCIENTIFIC=="Chionoecetes opilio"|trawl.data$SCIENTIFIC=="Atheresthes stomias"|
                                    trawl.data$SCIENTIFIC=="Hippoglossoides elassodon"|trawl.data$SCIENTIFIC=="Gadus macrocephalus"|
                                    trawl.data$SCIENTIFIC=="Lepidopsetta polyxystra"|trawl.data$SCIENTIFIC=="Lepidopsetta sp."|
                                   trawl.data$SCIENTIFIC=="Limanda aspera"|trawl.data$SCIENTIFIC=="Hippoglossus stenolepis"|
                                    trawl.data$SCIENTIFIC=="Chionoecetes bairdi"),]

sel.trawl.dat$logCPUE <- log(sel.trawl.dat$WTCPUE + 1)

#write.csv(sel.trawl.dat, "data/select_trawl_dat.csv")





