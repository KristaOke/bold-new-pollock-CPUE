#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#New simplified model structure on new dataset

#Krista, May 2022
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Notes:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(mgcv)
library(gamm4)
library(visreg)
library(gratia)
library("rnaturalearth")
library("rnaturalearthdata")
library( "ggspatial" )
library("sf")

#load data
wd <- getwd()
northsouthdata_all <- read.csv(paste(wd,"/data/survey data/combined_cleaned_north-south_1982-2021_bot_trawl_data.csv",sep=""))

newsebs <- northsouthdata_all[which(northsouthdata_all$STRATUM!="70" &
                                      northsouthdata_all$STRATUM!="71" &
                                      northsouthdata_all$STRATUM!="81" &
                                      northsouthdata_all$STRATUM!="0" ),]

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-178, -155), ylim = c(53, 65), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  # geom_point(aes(LONGITUDE, LATITUDE, colour=mean_station_bottemp), data=all_analysis_dat) +   
  # scale_colour_gradient2(low="blue", high="red", guide="colorbar") + 
  geom_point(aes(LONGITUDE, LATITUDE, 
                 col=as.factor(STRATUM)), data=newsebs) + theme_bw() 

#limit to pollock data
sebs_pollock <- newsebs[which(newsebs$SID=="21740"),]

#add period
sebs_pollock$period <- NA
sebs_pollock$period[which(sebs_pollock$YEAR<2014)] <- "early"
sebs_pollock$period[which(sebs_pollock$YEAR>2013)] <- "late"
sebs_pollock$period <- as.factor(sebs_pollock$period)

#fill in missing temps-----

tempmod82 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1982),] )
summary(tempmod82)
gam.check(tempmod82)
plot(tempmod82)

pre82 <- predict.gam(tempmod82, type="response")
df82 <- sebs_pollock[which(sebs_pollock$YEAR==1982),] 

missing82 <- df82[which(is.na(df82$BOT_TEMP)==TRUE),]
pred82 <- predict.gam(tempmod82, newdata=missing82, type="response")
missing82$BOT_TEMP <- pred82

tempmod83 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1983),] )
summary(tempmod83)
gam.check(tempmod83)
plot(tempmod83)

pre83 <- predict.gam(tempmod83, type="response")
df83 <- sebs_pollock[which(sebs_pollock$YEAR==1983),] 

missing83 <- df83[which(is.na(df83$BOT_TEMP)==TRUE),]
pred83 <- predict.gam(tempmod83, newdata=missing83, type="response")
missing83$BOT_TEMP <- pred83


tempmod84 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1984),] )
summary(tempmod84)
gam.check(tempmod84)
plot(tempmod84)

pre84 <- predict.gam(tempmod84, type="response")
df84 <- sebs_pollock[which(sebs_pollock$YEAR==1984),] 

missing84 <- df84[which(is.na(df84$BOT_TEMP)==TRUE),]
pred84 <- predict.gam(tempmod84, newdata=missing84, type="response")
missing84$BOT_TEMP <- pred84


tempmod85 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1985),] )
summary(tempmod85)
gam.check(tempmod85)
plot(tempmod85)

pre85 <- predict.gam(tempmod85, type="response")
df85 <- sebs_pollock[which(sebs_pollock$YEAR==1985),] 

missing85 <- df85[which(is.na(df85$BOT_TEMP)==TRUE),]
pred85 <- predict.gam(tempmod85, newdata=missing85, type="response")
missing85$BOT_TEMP <- pred85


tempmod86 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1986),] )
summary(tempmod86)
gam.check(tempmod86)
plot(tempmod86)

pre86 <- predict.gam(tempmod86, type="response")
df86 <- sebs_pollock[which(sebs_pollock$YEAR==1986),] 

missing86 <- df86[which(is.na(df86$BOT_TEMP)==TRUE),]
pred86 <- predict.gam(tempmod86, newdata=missing86, type="response")
missing86$BOT_TEMP <- pred86


tempmod87 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1987),] )
summary(tempmod87)
gam.check(tempmod87)
plot(tempmod87)

pre87 <- predict.gam(tempmod87, type="response")
df87 <- sebs_pollock[which(sebs_pollock$YEAR==1987),] 

missing87 <- df87[which(is.na(df87$BOT_TEMP)==TRUE),]
pred87 <- predict.gam(tempmod87, newdata=missing87, type="response")
missing87$BOT_TEMP <- pred87



tempmod88 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1988),] )
summary(tempmod88)
gam.check(tempmod88)
plot(tempmod88)

pre88 <- predict.gam(tempmod88, type="response")
df88 <- sebs_pollock[which(sebs_pollock$YEAR==1988),] 

missing88 <- df88[which(is.na(df88$BOT_TEMP)==TRUE),]
pred88 <- predict.gam(tempmod88, newdata=missing88, type="response")
missing88$BOT_TEMP <- pred88



tempmod89 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1989),] )
summary(tempmod89)
gam.check(tempmod89)
plot(tempmod89)

pre89 <- predict.gam(tempmod89, type="response")
df89 <- sebs_pollock[which(sebs_pollock$YEAR==1989),] 

missing89 <- df89[which(is.na(df89$BOT_TEMP)==TRUE),]
pred89 <- predict.gam(tempmod89, newdata=missing89, type="response")
missing89$BOT_TEMP <- pred89



tempmod90 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1990),] )
summary(tempmod90)
gam.check(tempmod90)
plot(tempmod90)

pre90 <- predict.gam(tempmod90, type="response")
df90 <- sebs_pollock[which(sebs_pollock$YEAR==1990),] 

missing90 <- df90[which(is.na(df90$BOT_TEMP)==TRUE),]
pred90 <- predict.gam(tempmod90, newdata=missing90, type="response")
missing90$BOT_TEMP <- pred90



tempmod91 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1991),]  )
summary(tempmod91)
gam.check(tempmod91)
plot(tempmod91)

pre91 <- predict.gam(tempmod91, type="response")
df91 <- sebs_pollock[which(sebs_pollock$YEAR==1991),] 

missing91 <- df91[which(is.na(df91$BOT_TEMP)==TRUE),]
pred91 <- predict.gam(tempmod91, newdata=missing91, type="response")
missing91$BOT_TEMP <- pred91



tempmod92 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1992),]  )
summary(tempmod92)
gam.check(tempmod92)
plot(tempmod92)

pre92 <- predict.gam(tempmod92, type="response")
df92 <- sebs_pollock[which(sebs_pollock$YEAR==1992),] 

missing92 <- df92[which(is.na(df92$BOT_TEMP)==TRUE),]
pred92 <- predict.gam(tempmod92, newdata=missing92, type="response")
missing92$BOT_TEMP <- pred92



tempmod93 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1993),]  )
summary(tempmod93)
gam.check(tempmod93) #seems awfully high
plot(tempmod93) 

pre93 <- predict.gam(tempmod93, type="response")
df93 <- sebs_pollock[which(sebs_pollock$YEAR==1993),] 

missing93 <- df93[which(is.na(df93$BOT_TEMP)==TRUE),]
pred93 <- predict.gam(tempmod93, newdata=missing93, type="response")
missing93$BOT_TEMP <- pred93



tempmod94 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1994),]  )
summary(tempmod94)
gam.check(tempmod94)
plot(tempmod94)

pre94 <- predict.gam(tempmod94, type="response")
df94 <- sebs_pollock[which(sebs_pollock$YEAR==1994),] 

missing94 <- df94[which(is.na(df94$BOT_TEMP)==TRUE),]
pred94 <- predict.gam(tempmod94, newdata=missing94, type="response")
missing94$BOT_TEMP <- pred94



tempmod95 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1995),]  )
summary(tempmod95)
gam.check(tempmod95)
plot(tempmod95)

pre95 <- predict.gam(tempmod95, type="response")
df95 <- sebs_pollock[which(sebs_pollock$YEAR==1995),] 

missing95 <- df95[which(is.na(df95$BOT_TEMP)==TRUE),]
pred95 <- predict.gam(tempmod95, newdata=missing95, type="response")
missing95$BOT_TEMP <- pred95



tempmod96 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1996),]  )
summary(tempmod96)
gam.check(tempmod96)
plot(tempmod96)

pre96 <- predict.gam(tempmod96, type="response")
df96 <- sebs_pollock[which(sebs_pollock$YEAR==1996),] 

missing96 <- df96[which(is.na(df96$BOT_TEMP)==TRUE),]
pred96 <- predict.gam(tempmod96, newdata=missing96, type="response")
missing96$BOT_TEMP <- pred96



tempmod97 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1997),]  )
summary(tempmod97)
gam.check(tempmod97)
plot(tempmod97)

pre97 <- predict.gam(tempmod97, type="response")
df97 <- sebs_pollock[which(sebs_pollock$YEAR==1997),] 

missing97 <- df97[which(is.na(df97$BOT_TEMP)==TRUE),]
pred97 <- predict.gam(tempmod97, newdata=missing97, type="response")
missing97$BOT_TEMP <- pred97



tempmod98 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1998),]  )
summary(tempmod98)
gam.check(tempmod98)
plot(tempmod98)

pre98 <- predict.gam(tempmod98, type="response")
df98 <- sebs_pollock[which(sebs_pollock$YEAR==1998),] 

missing98 <- df98[which(is.na(df98$BOT_TEMP)==TRUE),]
pred98 <- predict.gam(tempmod98, newdata=missing98, type="response")
missing98$BOT_TEMP <- pred98



tempmod99 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==1999),]  )
summary(tempmod99)
gam.check(tempmod99)
plot(tempmod99)

pre99 <- predict.gam(tempmod99, type="response")
df99 <- sebs_pollock[which(sebs_pollock$YEAR==1999),] 

missing99 <- df99[which(is.na(df99$BOT_TEMP)==TRUE),]
pred99 <- predict.gam(tempmod99, newdata=missing99, type="response")
missing99$BOT_TEMP <- pred99



tempmod00 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2000),]  )
summary(tempmod00)
gam.check(tempmod00)
plot(tempmod00)

pre00 <- predict.gam(tempmod00, type="response")
df00 <- sebs_pollock[which(sebs_pollock$YEAR==2000),] 

missing00 <- df00[which(is.na(df00$BOT_TEMP)==TRUE),]
pred00 <- predict.gam(tempmod00, newdata=missing00, type="response")
missing00$BOT_TEMP <- pred00



tempmod01 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2001),]  )
summary(tempmod01)
gam.check(tempmod01) #
plot(tempmod01)

pre01 <- predict.gam(tempmod01, type="response")
df01 <- sebs_pollock[which(sebs_pollock$YEAR==2001),] 

missing01 <- df01[which(is.na(df01$BOT_TEMP)==TRUE),]
pred01 <- predict.gam(tempmod01, newdata=missing01, type="response")
missing01$BOT_TEMP <- pred01



tempmod02 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2002),]  )
summary(tempmod02)
gam.check(tempmod02)
plot(tempmod02)

pre02 <- predict.gam(tempmod02, type="response")
df02 <- sebs_pollock[which(sebs_pollock$YEAR==2002),] 

missing02 <- df02[which(is.na(df02$BOT_TEMP)==TRUE),]
pred02 <- predict.gam(tempmod02, newdata=missing02, type="response")
missing02$BOT_TEMP <- pred02



tempmod03 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2003),]  )
summary(tempmod03)
gam.check(tempmod03)
plot(tempmod03)

pre03 <- predict.gam(tempmod03, type="response")
df03 <- sebs_pollock[which(sebs_pollock$YEAR==2003),] 

missing03 <- df03[which(is.na(df03$BOT_TEMP)==TRUE),]
pred03 <- predict.gam(tempmod03, newdata=missing03, type="response")
missing03$BOT_TEMP <- pred03



tempmod04 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2004),]  )
summary(tempmod04)
gam.check(tempmod04)
plot(tempmod04)

pre04 <- predict.gam(tempmod04, type="response")
df04 <- sebs_pollock[which(sebs_pollock$YEAR==2004),] 

missing04 <- df04[which(is.na(df04$BOT_TEMP)==TRUE),]
pred04 <- predict.gam(tempmod04, newdata=missing04, type="response")
missing04$BOT_TEMP <- pred04



tempmod05 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2005),]  )
summary(tempmod05)
gam.check(tempmod05)
plot(tempmod05)

pre05 <- predict.gam(tempmod05, type="response")
df05 <- sebs_pollock[which(sebs_pollock$YEAR==2005),] 

missing05 <- df05[which(is.na(df05$BOT_TEMP)==TRUE),]
pred05 <- predict.gam(tempmod05, newdata=missing05, type="response")
missing05$BOT_TEMP <- pred05



tempmod06 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2006),]  )
summary(tempmod06)
gam.check(tempmod06)
plot(tempmod06)

pre06 <- predict.gam(tempmod06, type="response")
df06 <- sebs_pollock[which(sebs_pollock$YEAR==2006),] 

missing06 <- df06[which(is.na(df06$BOT_TEMP)==TRUE),]
pred06 <- predict.gam(tempmod06, newdata=missing06, type="response")
missing06$BOT_TEMP <- pred06



tempmod07 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2007),]  )
summary(tempmod07)
gam.check(tempmod07)
plot(tempmod07)

pre07 <- predict.gam(tempmod07, type="response")
df07 <- sebs_pollock[which(sebs_pollock$YEAR==2007),] 

missing07 <- df07[which(is.na(df07$BOT_TEMP)==TRUE),]
pred07 <- predict.gam(tempmod07, newdata=missing07, type="response")
missing07$BOT_TEMP <- pred07



tempmod08 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2008),]  )
summary(tempmod08)
gam.check(tempmod08)
plot(tempmod08)

pre08 <- predict.gam(tempmod08, type="response")
df08 <- sebs_pollock[which(sebs_pollock$YEAR==2008),] 

missing08 <- df08[which(is.na(df08$BOT_TEMP)==TRUE),]
pred08 <- predict.gam(tempmod08, newdata=missing08, type="response")
missing08$BOT_TEMP <- pred08



tempmod09 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2009),]  )
summary(tempmod09)
gam.check(tempmod09)
plot(tempmod09)

pre09 <- predict.gam(tempmod09, type="response")
df09 <- sebs_pollock[which(sebs_pollock$YEAR==2009),] 

missing09 <- df09[which(is.na(df09$BOT_TEMP)==TRUE),]
pred09 <- predict.gam(tempmod09, newdata=missing09, type="response")
missing09$BOT_TEMP <- pred09



tempmod10 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2010),]  )
summary(tempmod10)
gam.check(tempmod10)
plot(tempmod10)

pre10 <- predict.gam(tempmod10, type="response")
df10 <- sebs_pollock[which(sebs_pollock$YEAR==2010),] 

missing10 <- df10[which(is.na(df10$BOT_TEMP)==TRUE),]
pred10 <- predict.gam(tempmod10, newdata=missing10, type="response")
missing10$BOT_TEMP <- pred10



tempmod11 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2011),] )
summary(tempmod11)
gam.check(tempmod11)
plot(tempmod11)

pre11 <- predict.gam(tempmod11, type="response")
df11 <- sebs_pollock[which(sebs_pollock$YEAR==2011),] 

missing11 <- df11[which(is.na(df11$BOT_TEMP)==TRUE),]
pred11 <- predict.gam(tempmod11, newdata=missing11, type="response")
missing11$BOT_TEMP <- pred11



tempmod12 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2012),] )
summary(tempmod12)
gam.check(tempmod12)
plot(tempmod12)

pre12 <- predict.gam(tempmod12, type="response")
df12 <- sebs_pollock[which(sebs_pollock$YEAR==2012),] 

missing12 <- df12[which(is.na(df12$BOT_TEMP)==TRUE),]
pred12 <- predict.gam(tempmod12, newdata=missing12, type="response")
missing12$BOT_TEMP <- pred12



tempmod13 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2013),] )
summary(tempmod13)
gam.check(tempmod13)
plot(tempmod13)

pre13 <- predict.gam(tempmod13, type="response")
df13 <- sebs_pollock[which(sebs_pollock$YEAR==2013),] 

missing13 <- df13[which(is.na(df13$BOT_TEMP)==TRUE),]
pred13 <- predict.gam(tempmod13, newdata=missing13, type="response")
missing13$BOT_TEMP <- pred13



tempmod14 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2014),] )
summary(tempmod14)
gam.check(tempmod14)
plot(tempmod14)

pre14 <- predict.gam(tempmod14, type="response")
df14 <- sebs_pollock[which(sebs_pollock$YEAR==2014),] 

missing14 <- df14[which(is.na(df14$BOT_TEMP)==TRUE),]
pred14 <- predict.gam(tempmod14, newdata=missing14, type="response")
missing14$BOT_TEMP <- pred14



tempmod15 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2015),] )
summary(tempmod15)
gam.check(tempmod15)
plot(tempmod15)

pre15 <- predict.gam(tempmod15, type="response")
df15 <- sebs_pollock[which(sebs_pollock$YEAR==2015),] 

missing15 <- df15[which(is.na(df15$BOT_TEMP)==TRUE),]
pred15 <- predict.gam(tempmod15, newdata=missing15, type="response")
missing15$BOT_TEMP <- pred15



tempmod16 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2016),] )
summary(tempmod16)
gam.check(tempmod16)
plot(tempmod16)

pre16 <- predict.gam(tempmod16, type="response")
df16 <- sebs_pollock[which(sebs_pollock$YEAR==2016),] 

missing16 <- df16[which(is.na(df16$BOT_TEMP)==TRUE),]
pred16 <- predict.gam(tempmod16, newdata=missing16, type="response")
missing16$BOT_TEMP <- pred16



tempmod17 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2017),] )
summary(tempmod17)
gam.check(tempmod17)
plot(tempmod17)

pre17 <- predict.gam(tempmod17, type="response")
df17 <- sebs_pollock[which(sebs_pollock$YEAR==2017),] 

missing17 <- df17[which(is.na(df17$BOT_TEMP)==TRUE),]
pred17 <- predict.gam(tempmod17, newdata=missing17, type="response")
missing17$BOT_TEMP <- pred17



tempmod18 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2018),] )
summary(tempmod18)
gam.check(tempmod18) #
plot(tempmod18)

pre18 <- predict.gam(tempmod18, type="response")
df18 <- sebs_pollock[which(sebs_pollock$YEAR==2018),] 

missing18 <- df18[which(is.na(df18$BOT_TEMP)==TRUE),]
pred18 <- predict.gam(tempmod18, newdata=missing18, type="response")
missing18$BOT_TEMP <- pred18



tempmod19 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2019),] )
summary(tempmod19)
gam.check(tempmod19)
plot(tempmod19)

pre19 <- predict.gam(tempmod19, type="response")
df19 <- sebs_pollock[which(sebs_pollock$YEAR==2019),] 

missing19 <- df19[which(is.na(df19$BOT_TEMP)==TRUE),]
pred19 <- predict.gam(tempmod19, newdata=missing19, type="response")
missing19$BOT_TEMP <- pred19



tempmod21 <- gam(BOT_TEMP ~ s(BOT_DEPTH) + ti(long_albers, lat_albers), data=sebs_pollock[which(sebs_pollock$YEAR==2021),] )
summary(tempmod21)
gam.check(tempmod21)
plot(tempmod21)

pre21 <- predict.gam(tempmod21, type="response")
df21 <- sebs_pollock[which(sebs_pollock$YEAR==2021),] 

missing21 <- df21[which(is.na(df21$BOT_TEMP)==TRUE),]
pred21 <- predict.gam(tempmod21, newdata=missing21, type="response")
missing21$BOT_TEMP <- pred21

#now bind together

missingall <- rbind(missing82, missing83, missing84, missing85, missing86, missing87, missing88, missing89,
                    missing90, missing91, missing92, missing93, missing94, missing95,
                    missing96, missing97, missing98, missing99, missing00,
                    missing01, missing02, missing03, missing04, missing05,
                    missing06, missing07, missing08, missing09, missing10, 
                    missing11, missing12, missing13, missing14, missing15,
                    missing16, missing17, missing18, missing19, missing21)


nona <- sebs_pollock[which(is.na(sebs_pollock$BOT_TEMP)==FALSE),]

sebs_model_dat <- rbind(nona, missingall)

wd <- getwd()
write_csv(sebs_model_dat, file=paste(wd,"/data/survey data/south_1982-2021_bot_trawl_data_for_sebs_models.csv", sep=""))



#models------

sebs_model_dat <- read.csv( file=paste(wd,"/data/survey data/south_1982-2021_bot_trawl_data_for_sebs_models.csv", sep=""))


startmod <- gamm4(logCPUE ~  s(BOT_DEPTH) +
               s(BOT_TEMP, by=period, bs="fs"),  random=~(1|YEAR/HAUL), 
             data=sebs_model_dat)
gam.check(startmod$gam)
summary(startmod[[1]]) 
summary(startmod[[2]]) 
anova(startmod[[2]])
plot(startmod[[2]])
#is there any spatial correlation though?

c_dat <- sebs_model_dat[which(is.finite(sebs_model_dat$BOT_DEPTH)==TRUE&
                           is.finite(sebs_model_dat$BOT_TEMP)==TRUE &
                           is.finite(sebs_model_dat$logCPUE)==TRUE ),]
c_dat$r <- NA
c_dat$r <- resid(startmod[[2]])   # Extract residuals
sp:::bubble(c_dat[which(c_dat$YEAR==2010),], zcol="r")
ggplot(c_dat, aes(long_albers, lat_albers, col=r)) + geom_point() +
  scale_colour_gradient(high="red", low="green") + facet_wrap(~YEAR)
#looks spatially correlated

gam.check(startmod[[2]]) 

cmod1 <- gamm(logCPUE ~  s(BOT_DEPTH) +
                    s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
               correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                  data=sebs_model_dat)
saveRDS(cmod1, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int.RDS")

summary(cmod1$gam)
summary(cmod1$lme)

c2_dat <- sebs_model_dat[which(is.finite(sebs_model_dat$BOT_DEPTH)==TRUE&
                                is.finite(sebs_model_dat$BOT_TEMP)==TRUE &
                                is.finite(sebs_model_dat$logCPUE)==TRUE ),]
c2_dat$r <- NA
c2_dat$r <- resid(cmod1[[2]])   # Extract residuals
ggplot(c2_dat, aes(long_albers, lat_albers, col=r)) + geom_point() +
  scale_colour_gradient(high="red", low="green") + facet_wrap(~YEAR)

cmod1S <- gamm(logCPUE ~  s(BOT_DEPTH) +
                s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
              correlation = corSpher(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
              data=sebs_model_dat)

cmod1G <- gamm(logCPUE ~  s(BOT_DEPTH) +
                 s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
               correlation = corGaus(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
               data=sebs_model_dat)

cmod1R <- gamm(logCPUE ~  s(BOT_DEPTH) +
                 s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
               correlation = corRatio(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
               data=sebs_model_dat)


AIC(cmod1R$lme, cmod1G$lme, cmod1S$lme, cmod1$lme) #exp is best

#check nugget
cmod1_nonug <- gamm(logCPUE ~  s(BOT_DEPTH) +
                s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
              correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=FALSE), method="REML",
              data=sebs_model_dat)

AIC(cmod1$lme, cmod1_nonug$lme) #better w nugget

clin_int <- gamm(logCPUE ~  s(BOT_DEPTH) +
                      BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                    data=sebs_model_dat)
saveRDS(clin_int, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int.RDS")


cmod1_noint <- gamm(logCPUE ~  s(BOT_DEPTH) +
                      s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                    data=sebs_model_dat)
saveRDS(cmod1_noint, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int.RDS")
cmod1_noint <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int.RDS")

#again w ML

cmod1_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
              correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
              data=sebs_model_dat)
#saveRDS(cmod1_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_ML.RDS")
cmod1_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_ML.RDS")

#
clin_int_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                   BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                 data=sebs_model_dat)
#saveRDS(clin_int_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_ML.RDS")
clin_int_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_ML.RDS")



cmod1_noint_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                      s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                    data=sebs_model_dat)
#saveRDS(cmod1_noint_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_ML.RDS")
cmod1_noint_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_ML.RDS")


AIC(clin_int_ML$lme, cmod1_ML$lme, cmod1_noint_ML$lme) #no int is best aic, smooth int close


#add temp as linear for null

cmod1_lin_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                         BOT_TEMP,  random=list(YEAR_factor=~1), 
                       correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                       data=sebs_model_dat)
#saveRDS(cmod1_lin_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lineartemp_ML.RDS")
cmod1_lin_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lineartemp_ML.RDS")


AIC(clin_int_ML$lme, cmod1_ML$lme, cmod1_noint_ML$lme, cmod1_lin_ML$lme) #smooths still best, this one worse than lin int

#look at best mod reml fit
gam.check(cmod1_noint$gam)
summary(cmod1_noint$lme)
summary(cmod1_noint$gam)
plot(cmod1_noint$gam)
 v <- getViz(cmod1_noint$gam)
 plot(sm(v, 2)) + l_fitLine(alpha = 0.6)

draw(cmod1_noint$gam, select=2) + theme_bw() + xlab("Bottom temperature (°C)") +
  ylab("Partial effect on log(CPUE)") + ggtitle("")


#repeat without deepest stations--------------------------

#only stations up to depths 180
#look at those deeper than 175
ggplot(sebs_model_dat[which(sebs_model_dat$BOT_DEPTH>175),], aes(STATION, BOT_DEPTH)) + geom_point()

table(sebs_model_dat$STATION[which(sebs_model_dat$BOT_DEPTH>175)])

#each of these stations have many more observations that are shallower, I will remove just the 
#rows with depths greater than 180 instead of those stations

sebs_shallower_180 <- sebs_model_dat[which(sebs_model_dat$BOT_DEPTH<180),]

#w ML

cmod1_180_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                   s(BOT_TEMP, by=as.factor(period), bs="fs"),  random=list(YEAR_factor=~1), 
                 correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                 data=sebs_shallower_180)
#saveRDS(cmod1_180_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_ML_lessthan180.RDS")
cmod1_180_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_ML_lessthan180.RDS")
aic_smint180 <- AIC(cmod1_180_ML$lme)

#
clin_int_180_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                      BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                    data=sebs_shallower_180)
#saveRDS(clin_int_180_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_ML_lessthan180.RDS")
clin_int_180_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_ML_lessthan180.RDS")
aic_linint180 <- AIC(clin_int_180_ML$lme)


cmod1_noint_180_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                         s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                       correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                       data=sebs_shallower_180)
saveRDS(cmod1_noint_180_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_ML_lessthan180.RDS")
cmod1_noint_180_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_ML_lessthan180.RDS")
aic_sm180 <- AIC(cmod1_noint_180_ML$lme)

AIC(cmod1_noint_180_ML$lme, clin_int_180_ML$lme, cmod1_180_ML$lme) #smooth no interaction still best




#as reml

cmod1_noint_180 <- gamm(logCPUE ~  s(BOT_DEPTH) +
                             s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                           correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                           data=sebs_shallower_180)
saveRDS(cmod1_noint_180, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_lessthan180.RDS")
cmod1_noint_180 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_lessthan180.RDS")


gam.check(cmod1_noint_180$gam)
summary(cmod1_noint_180$lme)
summary(cmod1_noint_180$gam)
plot(cmod1_noint_180$gam)
v <- getViz(cmod1_noint_180$gam)
plot(sm(v, 2)) + l_fitLine(alpha = 0.6)

draw(cmod1_noint_180$gam, select=2) + theme_bw() + xlab("Bottom temperature (°C)") +
  ylab("Partial effect on log(CPUE)") + ggtitle("")

draw(cmod1_noint_180$gam, select=1) + theme_bw() + xlab("Bottom depth (m)") +
  ylab("Partial effect on log(CPUE)") + ggtitle("")


#for comparision
cmod1_180_REML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                       s(BOT_TEMP, by=as.factor(period), bs="fs"),  random=list(YEAR_factor=~1), 
                     correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                     data=sebs_shallower_180)
#saveRDS(cmod1_180_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_REML.RDS")
cmod1_180_REML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_REML_lessthan180.RDS")
summary(cmod1_180_REML$gam)
anova(cmod1_180_REML$gam)

clin_int_180_REML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                          BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                        data=sebs_shallower_180)
#saveRDS(clin_int_180_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_REML_lessthan180.RDS")
clin_int_180_REML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_REML_lessthan180.RDS")
summary(clin_int_180_REML$gam)
anova(clin_int_180_REML$gam)

library(mgcViz)

v <- getViz(cmod1_180_REML$gam)

plot(sm(v, 2)) + l_fitLine(alpha = 0.6)


#are these not plotting the interaction because of the as.factor?
sebs_shallower_180$period <- as.factor(sebs_shallower_180$period)

ptest_REML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                         s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                       correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                       data=sebs_shallower_180)
#saveRDS(ptest_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/mod_test_period_factor.RDS")
#ptest_REML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/mod_test_period_factor.RDS")

vt <- getViz(ptest_REML$gam)

plot(sm(vt, 2)) + l_fitLine(alpha = 0.6)

library(itsadug)

par(mfrow=c(1,2))

plot_smooth(cmod1_noint_180$gam, view="BOT_TEMP",ylab="Partial effect on log(CPUE)", xlab="Bottom Temperature",
            ylim=c(-1,5))

plot_smooth(ptest_REML$gam, view="BOT_TEMP", plot_all=c("period"), ylab="Partial effect on log(CPUE)", xlab="Bottom Temperature",
            ylim=c(-1,5), col=c("blue", "red")) #finally works!


#akaike weights====

#now after time varying
# 
# AIC(cmod1_noint_180_ML$lme, clin_int_180_ML$lme, cmod1_180_ML$lme) #smooth no interaction still best
# 
# #smooth is lowest
# delta_smooth <-  0
# delta_smint <- aic_smint180 - aic_sm180
# delta_linint <- aic_linint180 - aic_sm180
# 
# 
# # akaike weights
# sumbinaic <- sum(exp(-0.5*delta_smooth), exp(-0.5*delta_smint), exp(-0.5*delta_linint)) 
# # 
# aw_smoothint <- exp(-0.5*delta_smint)/sumbinaic
# aw_smooth <- exp(-0.5*delta_smooth)/sumbinaic
# aw_linint <- exp(-0.5*delta_linint)/sumbinaic

#time varying review request=================================


#model suggested by reviewer

tvmod1_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                             s(YEAR, by=BOT_TEMP),  random=list(YEAR_factor=~1), 
                           correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                           data=sebs_shallower_180)
saveRDS(tvmod1_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_ML_lessthan180.RDS")
tvmod1_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_ML_lessthan180.RDS")
aic_tv180 <- AIC(tvmod1_ML$lme)
gam.check(tvmod1_ML$gam)
summary(tvmod1_ML$gam)
anova(tvmod1_ML$gam)
plot(tvmod1_ML$gam)

tvmod2_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                    t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                  correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                  data=sebs_shallower_180)
saveRDS(tvmod2_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_ML_lessthan180.RDS")
tvmod2_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_ML_lessthan180.RDS")
aic_tv2180 <- AIC(tvmod2_ML$lme)
gam.check(tvmod2_ML$gam)
summary(tvmod2_ML$gam)
anova(tvmod2_ML$gam)
plot(tvmod2_ML$gam)

AIC(cmod1_noint_180_ML$lme, clin_int_180_ML$lme, cmod1_180_ML$lme, tvmod2_ML$lme)
anova(cmod1_noint_180_ML$lme, clin_int_180_ML$lme, cmod1_180_ML$lme, tvmod2_ML$lme)

tvmod2_REML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                    t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                  correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                  data=sebs_shallower_180)
saveRDS(tvmod2_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_REML_lessthan180.RDS")
tvmod2_REML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_REML_lessthan180.RDS")
draw(tvmod2_REML$gam, select=2, dist=0.05)
summary(tvmod2_REML$gam)

#ONCE MORE w te not t2
tvte_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                    te(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                  correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                  data=sebs_shallower_180)
saveRDS(tvte_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_varTE_ML_lessthan180.RDS")
tvte_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_varTE_ML_lessthan180.RDS")
aic_tvte180 <- AIC(tvte_ML$lme)
gam.check(tvte_ML$gam)
summary(tvte_ML$gam)
anova(tvte_ML$gam)
plot(tvte_ML$gam)

AIC(cmod1_noint_180_ML$lme, clin_int_180_ML$lme, cmod1_180_ML$lme, tvmod2_ML$lme, tvte_ML$lme)
anova(cmod1_noint_180_ML$lme, clin_int_180_ML$lme, cmod1_180_ML$lme, tvmod2_ML$lme)

# tvmod2_REML <- gamm(logCPUE ~  s(BOT_DEPTH) +
#                       t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
#                     correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
#                     data=sebs_shallower_180)
# saveRDS(tvmod2_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_REML_lessthan180.RDS")
# tvmod2_REML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_REML_lessthan180.RDS")
# draw(tvmod2_REML$gam, select=2, dist=0.05)
# summary(tvmod2_REML$gam)




#akaike weights w time varying====
#tv is lowest
delta_tv <-  0
delta_sm <- aic_sm180 - aic_tv2180
delta_smint <- aic_smint180 - aic_tv2180
delta_linint <- aic_linint180 - aic_tv2180



# akaike weights
sumbinaic <- sum(exp(-0.5*delta_smooth), exp(-0.5*delta_smint), exp(-0.5*delta_linint), exp(-0.5*delta_tv)) 
# 
aw_smoothint <- exp(-0.5*delta_smint)/sumbinaic
aw_smooth <- exp(-0.5*delta_smooth)/sumbinaic
aw_linint <- exp(-0.5*delta_linint)/sumbinaic
aw_tv <- exp(-0.5*delta_tv)/sumbinaic


#plot time varying=======================================================================

b1 <- getViz(tvmod2_REML$gam)
pb <- plot(sm(b1, 2)) +
  labs(title = NULL) + #l_points() +
  scale_fill_distiller(palette = "Spectral", type = "div") +
  # theme(legend.position = "none")+ 
  theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        #  axis.title.x = element_blank(),
        #  axis.title.y = element_blank(),
        plot.title = element_blank()) + xlab("Year") +
  ylab("Bottom temperature") + labs(fill = "Effect")

pb 

fvisgam(tvmod2_REML$gam, view=c("YEAR", "BOT_TEMP"))

#plot depth effect========

draw(tvmod2_REML$gam, select=1) + theme_bw() + xlab("Bottom depth (m)") +
  ylab("Partial effect on log(CPUE)") + ggtitle("")


#only overlap============================================================================

#get percentiles

min(sebs_shallower_180$BOT_TEMP[which(sebs_shallower_180$period=="early")])  
min(sebs_shallower_180$BOT_TEMP[which(sebs_shallower_180$period=="late")])  

max(sebs_shallower_180$BOT_TEMP[which(sebs_shallower_180$period=="early")])  
max(sebs_shallower_180$BOT_TEMP[which(sebs_shallower_180$period=="late")])  

ggplot(sebs_shallower_180, aes(period, BOT_TEMP)) + geom_boxplot() + geom_point(position="jitter")

ggplot(sebs_shallower_180, aes(YEAR, BOT_TEMP)) + geom_point()

tempsummary <- sebs_shallower_180 %>% group_by(period) %>%
  summarize(mean_Btemp=mean(BOT_TEMP, na.rm=TRUE),
            q_05=quantile(BOT_TEMP, 0.05, na.rm=TRUE),
            q_95=quantile(BOT_TEMP, 0.95, na.rm=TRUE),
            q_02=quantile(BOT_TEMP, 0.02, na.rm=TRUE),
            q_98=quantile(BOT_TEMP, 0.98, na.rm=TRUE))



#analyze only data that falls within the 5th and 95th quantile of both 
overlapperiods <- sebs_shallower_180[which(sebs_shallower_180$BOT_TEMP<4.8 & sebs_shallower_180$BOT_TEMP>0),]
broadoverlapperiods <- sebs_shallower_180[which(sebs_shallower_180$BOT_TEMP<5.6 & sebs_shallower_180$BOT_TEMP>-1),]

#smaller overlap

ovrlp_noint_180 <- gamm(logCPUE ~  s(BOT_DEPTH) +
                          s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                        data=overlapperiods)
saveRDS(ovrlp_noint_180, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_no_int_lessthan180.RDS")
ovrlp_noint_180 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_no_int_lessthan180.RDS")

ovrlp_smint_180 <- gamm(logCPUE ~  s(BOT_DEPTH) +
                          s(BOT_TEMP, by=as.factor(period), bs="fs"),  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                        data=overlapperiods)
saveRDS(ovrlp_smint_180, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180.RDS")
ovrlp_smint_180 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180.RDS")

ovrlp_linint_180 <- gamm(logCPUE ~  s(BOT_DEPTH) +
                          BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                        data=overlapperiods)
saveRDS(ovrlp_linint_180, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180.RDS")
ovrlp_linint_180 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180.RDS")

AIC(ovrlp_linint_180$lme, ovrlp_smint_180$lme, ovrlp_noint_180$lme)

plot(ovrlp_smint_180$gam)
plot(ovrlp_noint_180$gam)


ovrlap_tv_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                    t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                  correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                  data=overlapperiods)
saveRDS(ovrlap_tv_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_onlyoverlap_ML_lessthan180.RDS")
ovrlap_tv_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_onlyoverlap_ML_lessthan180.RDS")
aic_overlaptv2180 <- AIC(ovrlap_tv_ML$lme)
gam.check(ovrlap_tv_ML$gam)
summary(ovrlap_tv_ML$gam)
anova(ovrlap_tv_ML$gam)
plot(ovrlap_tv_ML$gam)

ovrlap_tv_REML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                      t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                    data=overlapperiods)
saveRDS(ovrlap_tv_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_onlyoverlap_REML_lessthan180.RDS")
ovrlap_tv_REML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_onlyoverlap_REML_lessthan180.RDS")
draw(ovrlap_tv_REML$gam, select=2, dist=0.05)
summary(ovrlap_tv_REML$gam)

#broader overlap

bigovrlp_noint_180 <- gamm(logCPUE ~  s(BOT_DEPTH) +
                          s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                        data=broadoverlapperiods)
saveRDS(bigovrlp_noint_180, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/big-overlap_no_int_lessthan180.RDS")
bigovrlp_noint_180 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/big-overlap_no_int_lessthan180.RDS")

bigovrlp_smint_180 <- gamm(logCPUE ~  s(BOT_DEPTH) +
                          s(BOT_TEMP, by=as.factor(period), bs="fs"),  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                        data=broadoverlapperiods)
saveRDS(bigovrlp_smint_180, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180.RDS")
bigovrlp_smint_180 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180.RDS")

bigovrlp_linint_180 <- gamm(logCPUE ~  s(BOT_DEPTH) +
                           BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                         correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                         data=broadoverlapperiods)
saveRDS(bigovrlp_linint_180, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180.RDS")
bigovrlp_linint_180 <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180.RDS")

AIC(bigovrlp_linint_180$lme, bigovrlp_smint_180$lme, bigovrlp_noint_180$lme)

plot(bigovrlp_smint_180$gam)
plot(bigovrlp_noint_180$gam)



bigovrlap_tv_ML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                       t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                     correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                     data=broadoverlapperiods)
saveRDS(bigovrlap_tv_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_bigoverlap_ML_lessthan180.RDS")
bigovrlap_tv_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_bigoverlap_ML_lessthan180.RDS")
aic_bigoverlaptv2180 <- AIC(bigovrlap_tv_ML$lme)
gam.check(bigovrlap_tv_ML$gam)
summary(bigovrlap_tv_ML$gam)
anova(bigovrlap_tv_ML$gam)
plot(bigovrlap_tv_ML$gam)

bigovrlap_tv_REML <- gamm(logCPUE ~  s(BOT_DEPTH) +
                         t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                       correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                       data=broadoverlapperiods)
saveRDS(bigovrlap_tv_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_bigoverlap_REML_lessthan180.RDS")
bigovrlap_tv_REML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_bigoverlap_REML_lessthan180.RDS")
draw(bigovrlap_tv_REML$gam, select=2, dist=0.05)
summary(bigovrlap_tv_REML$gam)


#refit w reml

ovrlp_noint_reml <- gamm(logCPUE ~  s(BOT_DEPTH) +
                          s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                        data=overlapperiods)
saveRDS(ovrlp_noint_reml, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_no_int_lessthan180_reml.RDS")
ovrlp_noint_reml <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_no_int_lessthan180_reml.RDS")

overlapperiods$period <- as.factor(overlapperiods$period)

ovrlp_smint_reml <- gamm(logCPUE ~  s(BOT_DEPTH) +
                          s(BOT_TEMP, by=period, bs="fs"),  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                        data=overlapperiods)
saveRDS(ovrlp_smint_reml, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180_reml.RDS")
ovrlp_smint_reml <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/overlap_smooth_int_lessthan180_reml.RDS")



par(mfrow=c(1,2))

plot_smooth(ovrlp_noint_reml$gam, view="BOT_TEMP",ylab="Partial effect on log(CPUE)", xlab="Bottom Temperature",
            ylim=c(-1,5))

plot_smooth(ovrlp_smint_reml$gam, view="BOT_TEMP", plot_all=c("period"), ylab="Partial effect on log(CPUE)", xlab="Bottom Temperature",
            ylim=c(-1,5), col=c("blue", "red")) #finally works!


AIC(ovrlp_smint_reml$lme, ovrlp_noint_reml$lme)



#repeat without log(CPUE+1)============================

#add a new column responding to reviewer comment

sebs_shallower_180$logCPUE_minconstant <- log(sebs_shallower_180$WTCPUE + 0.0000001)

#w ML

cmod1_log_ML <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                       s(BOT_TEMP, by=as.factor(period), bs="fs"),  random=list(YEAR_factor=~1), 
                     correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                     data=sebs_shallower_180)
#saveRDS(cmod1_log_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_ML_log.RDS")
cmod1_log_ML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_ML_log.RDS")
aic_smintlog <- AIC(cmod1_log_ML$lme)

#
clin_int_log_ML <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                          BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                        data=sebs_shallower_180)
#saveRDS(clin_int_log_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_ML_log.RDS")
clin_int_log_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_ML_log.RDS")
aic_linintlog <- AIC(clin_int_log_ML$lme)


cmod1_noint_log_ML <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                             s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                           correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                           data=sebs_shallower_180)
saveRDS(cmod1_noint_log_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_ML_log.RDS")
cmod1_noint_log_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_ML_log.RDS")
aic_smlog <- AIC(cmod1_noint_log_ML$lme)

AIC(cmod1_noint_log_ML$lme, clin_int_log_ML$lme, cmod1_log_ML$lme) #




#as reml

cmod1_noint_log <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                          s(BOT_TEMP),  random=list(YEAR_factor=~1), 
                        correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                        data=sebs_shallower_180)
saveRDS(cmod1_noint_log, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_log.RDS")
cmod1_noint_log <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_no_int_log.RDS")


gam.check(cmod1_noint_log$gam)
summary(cmod1_noint_log$lme)
summary(cmod1_noint_log$gam)
plot(cmod1_noint_log$gam)
v <- getViz(cmod1_noint_log$gam)
plot(sm(v, 2)) + l_fitLine(alpha = 0.6)

draw(cmod1_noint_log$gam, select=2) + theme_bw() + xlab("Bottom temperature (°C)") +
  ylab("Partial effect on log(CPUE)") + ggtitle("")

draw(cmod1_noint_log$gam, select=1) + theme_bw() + xlab("Bottom depth (m)") +
  ylab("Partial effect on log(CPUE)") + ggtitle("")


#for comparision
cmod1_log_REML <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                         s(BOT_TEMP, by=as.factor(period), bs="fs"),  random=list(YEAR_factor=~1), 
                       correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                       data=sebs_shallower_180)
#saveRDS(cmod1_log_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_REML_log.RDS")
cmod1_log_REML <- readRDS(file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_corexp_sm_int_REML_log.RDS")
summary(cmod1_log_REML$gam)
anova(cmod1_log_REML$gam)

clin_int_log_REML <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                            BOT_TEMP:period,  random=list(YEAR_factor=~1), 
                          correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                          data=sebs_shallower_180)
#saveRDS(clin_int_log_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_REML_log.RDS")
clin_int_log_REML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/new_lin_int_REML_log.RDS")
summary(clin_int_log_REML$gam)
anova(clin_int_log_REML$gam)

library(mgcViz)

v <- getViz(cmod1_log_REML$gam)

plot(sm(v, 2)) + l_fitLine(alpha = 0.6)



#model suggested by reviewer

tvmod1log_ML <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                    s(YEAR, by=BOT_TEMP),  random=list(YEAR_factor=~1), 
                  correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                  data=sebs_shallower_180)
saveRDS(tvmod1log_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_ML_log.RDS")
tvmod1log_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var_ML_log.RDS")
aic_tvlog <- AIC(tvmod1log_ML$lme)
gam.check(tvmod1log_ML$gam)
summary(tvmod1log_ML$gam)
anova(tvmod1log_ML$gam)
plot(tvmod1log_ML$gam)

tvmod2log_ML <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                    t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                  correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                  data=sebs_shallower_180)
saveRDS(tvmod2log_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_ML_log.RDS")
tvmod2log_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_ML_log.RDS")
aic_tv2log <- AIC(tvmod2log_ML$lme)
gam.check(tvmod2log_ML$gam)
summary(tvmod2log_ML$gam)
anova(tvmod2log_ML$gam)
plot(tvmod2log_ML$gam)

AIC(cmod1_noint_log_ML$lme, clin_int_log_ML$lme, cmod1_log_ML$lme, tvmod2log_ML$lme)
anova(cmod1_noint_log_ML$lme, clin_int_log_ML$lme, cmod1_log_ML$lme, tvmod2log_ML$lme)

tvmod2log_REML <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                      t2(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                    correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="REML",
                    data=sebs_shallower_180)
saveRDS(tvmod2log_REML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_REML_log.RDS")
tvmod2log_REML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_var2_REML_log.RDS")
draw(tvmod2log_REML$gam, select=2, dist=0.05)
summary(tvmod2log_REML$gam)

#ONCE MORE w te not t2
tvtelog_ML <- gamm(logCPUE_minconstant ~  s(BOT_DEPTH) +
                  te(YEAR, BOT_TEMP),  random=list(YEAR_factor=~1), 
                correlation = corExp(form=~ long_albers + lat_albers|YEAR_factor, nugget=TRUE), method="ML",
                data=sebs_shallower_180)
saveRDS(tvtelog_ML, file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_varTE_ML_log.RDS")
tvtelog_ML <- readRDS( file="~/Dropbox/Work folder/Pollock Analyses/bold-new-pollock/scripts/time_varTE_ML_log.RDS")
aic_tvtelog <- AIC(tvtelog_ML$lme)
gam.check(tvtelog_ML$gam)
summary(tvtelog_ML$gam)
anova(tvtelog_ML$gam)
plot(tvtelog_ML$gam)

AIC(cmod1_noint_log_ML$lme, clin_int_log_ML$lme, cmod1_log_ML$lme, tvmod2log_ML$lme, tvtelog_ML$lme)
anova(cmod1_noint_log_ML$lme, clin_int_log_ML$lme, cmod1_log_ML$lme, tvmod2log_ML$lme)


#akaike weights 
#tv is lowest
delta_tvlog <-  0
delta_smlog <- aic_smlog - aic_tv2log
delta_smintlog <- aic_smintlog - aic_tv2log
delta_linintlog <- aic_linintlog - aic_tv2log


#UPDATE BELOW
# akaike weights
sumbinaic <- sum(exp(-0.5*delta_smooth), exp(-0.5*delta_smint), exp(-0.5*delta_linint), exp(-0.5*delta_tv)) 
# 
aw_smoothint <- exp(-0.5*delta_smint)/sumbinaic
aw_smooth <- exp(-0.5*delta_smooth)/sumbinaic
aw_linint <- exp(-0.5*delta_linint)/sumbinaic
aw_tv <- exp(-0.5*delta_tv)/sumbinaic





