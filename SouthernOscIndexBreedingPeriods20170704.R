##########################################
############ step 1. get set up to work ########
##########################################

rm(list = ls())
library(data.table)
library(raster)
library(circular)
library(gtools)
library(Hmisc)
library(plyr)
library(Hmisc)
library(multcomp)
library(gplots)
library(reshape2)
library(reporttools)
#figure
library(lme4)
#library(afex)
library(ENmisc)
library(arm)
library(car)
library(heavy)
calcskew <- function(circulardates){
  datesC<-circular( circulardates,modulo ="2pi", units="radians", rotation="counter")
  Rbar<-rho.circular(datesC)#average clustering 0 is uncluseted 1 is all at same location
  V<-1-Rbar#sample circular variance
  t2t<- trigonometric.moment(datesC, p=2, center=TRUE)
  bbar2 <- t2t$sin
  skewness <- bbar2/(V**(3/2)) #skewness 
  return(round(skewness,2))
}



##########################################
############ step 3. set up datafiles ########
##########################################

obs.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
# Breeding quantiles
dat<-fread(paste0(obs.dir,'PointOfLayDayOfYear2016-09-20.csv'))[,c( "Scientific.Name",
                                                                    "lat",
                                                                    "lon",
                                                                    "sourceName",
                                                                    "type",
                                                                    "DOY_PL",
                                                                    "year")]

dat$DOY_PL[dat$DOY_PL==366]<-365
dat$month<-month(as.Date(as.character(dat$DOY_PL),format="%j"))



# traits
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2016-10-07.csv')
taxon<-inc[c("genus","Clem.Family","ScientificName","Clements.Scientific.name","sort.v2016","family","English.name","Order")]
dat<-merge(dat,taxon,all.x=TRUE,by.x="Scientific.Name",by.y="ScientificName")
#add koeppen zones
koeppen<-raster(paste0('/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'),
                proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
#combine equatorial and tropical: 41 Equatorial, 35 Tropical, extract koeppen zone
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
koeppen<- reclassify(koeppen, rclmat)
dat$koeppen<-raster::extract(koeppen,cbind(dat$lon, dat$lat))


#convert data to radians to make circular vector
dat$Radians <-(dat$DOY_PL/365*360)*pi / 180
#Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,3-Grassland,3-Temperate
TEMP<-subset(dat, koeppen == 3) #temperate region Breeding Period
TEMP$month<-formatC(TEMP$month,width=2,flag='0') #format months 1 becomes 01
TEMP$date<-paste0(TEMP$year,TEMP$month) #make date string that matches NOAA SOI
TEMP<-subset(TEMP, date >=190001) #limit time to same as NOAA data
#read in BOM ENSO events
#BOM
SOI<-read.csv("/Users/daisy/Google Drive/PhD/Oscilations/Data/BOM_SOI.csv")
SOI$month2<-formatC(SOI$monthNumeric,width=2,flag='0') #format months 1 becomes 01
SOI$date<-paste0(SOI$Year,SOI$month2) #make date string that matches NOAA SOI
SOI<-SOI[,c("date","BreedingENSO","Year","monthNumeric","month2")]
# #NOAA - no longer using
# NOAA<-read.csv("/Users/daisy/Google Drive/PhD/Data/SOI/NOAA_OceanicNinoIndex5MonthPhase.csv")
# colnames(NOAA)<-c("Year",1:12)
# NOAA<-melt(NOAA, id="Year")
# NOAA$month<-formatC(NOAA$variable,width=2,flag='0') #format months 1 becomes 01
# NOAA$date<-paste0(NOAA$Year,NOAA$month)
# NOAA<-NOAA[,c("value","date")]
# colnames(NOAA)<-c("NOAA","date")
# SOI<-merge(SOI,NOAA,all.x=TRUE,by="date")


TEMP<-merge(TEMP,SOI, by.x="date",by.y="date",all.x=TRUE)
#TEMP<-subset(TEMP,!is.na(TEMP$SOI_sustained5months))


# #See if there is bias in the months of data
#
# el<-subset(SOI,NOAA=='E')
# la<-subset(SOI,NOAA=='L')
# nu<-subset(SOI,NOAA=='N')
# par(mfrow = c(1,3))
# hist(as.numeric(el$month2), main = "El Nino",ylim=c(0,45))
# hist(as.numeric(la$month2), main = "La Nina",ylim=c(0,45))
# hist(as.numeric(nu$month2), main = "Neutral",ylim=c(0,45))
#

##########################################
############ step 4. get species with more than 100 obs in each phase ########
##########################################

#get el Nino obs and species with more than 100 observations
eN<-subset(TEMP,BreedingENSO=='NINO')
spSumeN<-as.data.frame(table(eN$Scientific.Name))
sp_eN<-as.vector(droplevels(subset(spSumeN,Freq>=100)$Var1))
lN<-subset(TEMP,BreedingENSO=='NINA')
spSumlN<-as.data.frame(table(lN$Scientific.Name))
sp_lN<-as.vector(droplevels(subset(spSumlN,Freq>=100)$Var1))
# nN<-subset(TEMP,NOAA=="N")
# spSumnN<-as.data.frame(table(nN$Scientific.Name))
# sp_nN<-as.vector(droplevels(subset(spSumnN,Freq>=100)$Var1))
#name of species in both, limit data to these species
fin_sp<-intersect(sp_eN,sp_lN)
#fin_sp<-intersect(fin_sp,sp_nN)#neutral
eN<-droplevels(eN[eN$Scientific.Name %in% fin_sp,])
eN$SOI<-as.character(eN$SOI)
eN$SOI<-"elNino"
lN<-lN[lN$Scientific.Name %in% fin_sp,]
lN$SOI<-as.character(lN$SOI)
lN$SOI<-"LaNina"
# nN<-nN[nN$Scientific.Name %in% fin_sp,]
# nN$SOI<-as.character(nN$SOI)
# nN$SOI<-"Neutral"
# #make histograms of first egg dates
# par(mfrow = c(1,2))
# hist(lN$DOY_PL,main="NINA",xlab="",ylab = "",xlim = c(1,365),freq=F,breaks=seq(1,365,length = 12),ylim=c(0,0.009))
# hist(eN$DOY_PL,main="NINO",xlab="",ylab = "",xlim = c(1,365),freq=F,breaks=seq(1,365,length = 12),ylim=c(0,0.009))
# #hist(nN$DOY_PL,main="NEUT",xlab="",ylab = "",xlim = c(1,365),freq=TRUE,breaks=seq(1,365,length = 12))


########################################
###########Step 4b. remove 'unknown' breeding observations when there are more than 100 from the other catagories
#######################################


#get NINO species with more than 100 Egg or young or multi observations
  datAcc<-subset(eN, type == "solitary-egg" | type == "multi"| type == "young")
  datAcc_sp<-as.data.frame(table(datAcc$Scientific.Name))
  dat_eN_Acc_sp<-as.vector(droplevels(subset(datAcc_sp,Freq>=100)$Var1))
#get NINA species with more than 100 Egg or young or multi observations
  dat_LA_Acc<-subset(lN, type == "solitary-egg" | type == "multi"| type == "young")
  dat_LA_Acc_sp<-as.data.frame(table(dat_LA_Acc$Scientific.Name))
  dat_LA_Acc_sp<-as.vector(droplevels(subset(dat_LA_Acc_sp,Freq>=100)$Var1))
#species in both la nina and el nino
  fin_sp_Acc<-intersect(dat_eN_Acc_sp,dat_LA_Acc_sp) #25 species in total
#get the high accuracy data
  dat_eN_Acc<-datAcc[datAcc$Scientific.Name %in% fin_sp_Acc,]
  dat_LA_Acc<-dat_LA_Acc[dat_LA_Acc$Scientific.Name %in% fin_sp_Acc,]
# #high accuracy neutral
#   dat_nN_Acc<-subset(nN, type == "solitary-egg" | type == "multi"| type == "young")
#   dat_Nu_Acc<-dat_nN_Acc[dat_nN_Acc$Scientific.Name %in% fin_sp_Acc,]
#put the data back together
  high<-rbind(dat_LA_Acc,dat_eN_Acc)#,dat_Nu_Acc)
  high$Acc<-"high"
#get the data for the rest of the species
  dat_eN_LOW<-eN[eN$Scientific.Name %nin% fin_sp_Acc,]
    length(unique(dat_eN_LOW$Scientific.Name))
  dat_LA_LOW<-lN[lN$Scientific.Name %nin% fin_sp_Acc,]
    length(unique(dat_LA_LOW$Scientific.Name))
  # dat_Nu_LOW<-nN[nN$Scientific.Name %nin% fin_sp_Acc,]
  #   length(unique(dat_Nu_LOW$Scientific.Name))
  low<-rbind(dat_eN_LOW,dat_LA_LOW)#,dat_Nu_LOW)
  low$Acc<-"low"

  #put data back together
  PhaseDat<-droplevels(rbind(high,low))
  
##########
elevation<-raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Elevation/eMAST_ANUClimate_fx_elev_v1m0.nc")
distToCoast<-raster("/Users/daisy/Google Drive/PhD/Data/Spatial/DistanceToCoast/eMAST_ANUClimate_fx_dist_v1m0.nc")
  
PhaseDat$elevation<-raster::extract(elevation,cbind(PhaseDat$lon, PhaseDat$lat))
PhaseDat$distToCoast<-raster::extract(distToCoast,cbind(PhaseDat$lon, PhaseDat$lat))

  ##########################################
############ step 5. Calc breeding quantiles during differnt phases ########
##########################################

#combine data
SpeciesSummary<-list()
state<-c("elNino","LaNina")#,"Neutral")
alldat90<-list()
for (i in 1:length(fin_sp)){
  phaseSummary<-list()
  phasealldat90<-list()
  for (ii in 1:length(state)){
    #get species data in specific phase  
    spdat<-subset(PhaseDat,Scientific.Name==fin_sp[i] & SOI == state[ii])#get species data
    try(if(nrow(spdat) < 100) stop("not enough observations"))
    NoObs<-as.numeric(nrow(spdat))#get observation count
    CommonName<-as.character(subset(inc,ScientificName==fin_sp[i],"Common.me",drop=TRUE))
    Order<-as.character(subset(inc,ScientificName==fin_sp[i],order,drop=TRUE))
    family<-as.character(subset(inc,ScientificName==fin_sp[i],family,drop=TRUE))
    #clades<-as.character(subset(inc,ScientificName==fin_sp[i],"Patch.clade",drop=TRUE))
    obsRadians<-spdat$Radians
    #date of quantiles
    quant<-quantile.circular(obsRadians,c(0,.05,.5,.95),type=8)
    Per0<-round((quant[[1]]*180/pi)/360*365)
    Per5<-round((quant[[2]]*180/pi)/360*365)
    Per50<-round((quant[[3]]*180/pi)/360*365)
    Per95<-round((quant[[4]]*180/pi)/360*365)
    #Breeding Period Length
    StartDate<-as.numeric(Per5)
    EndDate<-as.numeric(Per95)
    BPL<-ifelse (EndDate < StartDate, 365-StartDate+EndDate,EndDate - StartDate )
    ###90% of data, remove obs outside of peak breeding
    spdat90<-if (quant[[2]]<quant[[4]]){
      subset(spdat, Radians >=quant[[2]] & Radians <= quant[[4]])
    }else rbind(subset(spdat, Radians >=quant[[2]]),subset(spdat, Radians <=quant[[4]]))
    #calculate skewness and clumpiness
    obsRadians90C<-circular( spdat90$Radians,modulo ="2pi", units="radians", rotation="counter")
    Rbar90Per<-round(rho.circular( obsRadians90C),2)#average clustering 0 is uncluseted 1 is all at same location
    skew90<-calcskew(spdat90$Radians)#negative numbers skewed counter clockwise, poitive clockwise, 0 = not skewed
    #put the data together and name it
    Perdat<-data.frame(paste(state[ii]),paste(fin_sp[i]),NoObs,CommonName,Order,family,length(obsRadians),
                  Per0,Per5,Per50,Per95,BPL,Rbar90Per,skew90)
    colnames(Perdat)<-c('Phase','Species','NoObs','CommonName','Order','family','ObservationCount','Quantile0','Quantile5',
                        'Quantile50','Quantile95','BreedingPeriod','Rbar90','Skew90')
    Perdat$Accu<-unique(spdat$Acc)
    phaseSummary[[ii]]<-Perdat
    phasealldat90[[ii]]<-spdat90
  }
  SpeciesSummary[[i]]<-do.call('rbind',phaseSummary)
  alldat90[[i]]<-do.call('rbind',phasealldat90)
}

alldat90<-as.data.frame(do.call('rbind',alldat90))
allDat<-as.data.frame(do.call('rbind',SpeciesSummary))

try(if(nrow(alldat90) < nrow(PhaseDat)-nrow(PhaseDat)/10) stop("more than 10% removed"))

######################################################
########test to see if el nino and la nina circualr data has same distributions
######################################################


modDat<-list()
for(j in 1:length(fin_sp)){
  spdat<-subset(alldat90,Scientific.Name==fin_sp[j])
  summar<-subset(allDat,Species==fin_sp[j])
  elN<-subset(summar,Phase=="elNino",select=c("Species","CommonName","Accu","Quantile5","Quantile50","BreedingPeriod","Rbar90","Skew90"))
  
  laN<-subset(summar,Phase=="LaNina",select=c("Quantile5","Quantile50","BreedingPeriod","Rbar90","Skew90"))
  colnames(laN)<-c("LNQuantile5","LNQuantile50","LNBreedingPeriod","LNRbar90","LNSkew90")
  #make data circular
  data<-list(
    el=circular(subset(spdat,SOI=="elNino",select= Radians),units="radians",template="geographics"),
    la=circular(subset(spdat,SOI=="LaNina",select= Radians),units="radians",template="geographics")
    )
  

  #test for significant differences in the two distributions, 
  #Watson's two-sample U2 tests (Watson 1962), 
  #tests if two samples of circular data come from same population
      datetStat<-round(watson.two.test(data$la,data$el,alpha=0.05)$statistic,3)
      datepVal<-if (datetStat > 0.385) {
        "< 0.001"} else 
          if (datetStat > 0.268) {
        "< 0.01"} else 
          if (datetStat >= 0.187) {
        "< 0.05"} else 
        {"> 0.05"}
      
      Species<-fin_sp[j]
      ObsCount<-sum(summar$NoObs)
      
modDat[[j]]<-cbind(elN,laN,datetStat,datepVal,ObsCount)
}

modDat<-do.call("rbind",modDat)
View(modDat)

#50% of the species tested as having significantly differnt distributions

############################
####Modify date for lm
############################

ModifiedDate<-list()
for(m in 1:length(fin_sp)){
spdat<-subset(alldat90,Scientific.Name==fin_sp[m])
  earlstart<-min(subset(allDat,Species==fin_sp[m])$Quantile5)
  latestend<-max(subset(allDat,Species==fin_sp[m])$Quantile95)
  #adjust dates so those at the begining of the year come after 365

    df<-subset(spdat,DOY_PL<earlstart)#dates to adjust
    df$modifiedDOY<-df$DOY_PL+365
    df2<-subset(spdat,DOY_PL>earlstart)#dates that are fine
    df2$modifiedDOY<-df2$DOY_PL
    ModifiedDate[[m]]<-rbind(df,df2)
}

ModifiedDate<-do.call("rbind",ModifiedDate)

######################################################
#########lm for each species 
######################################################


for(j in 1:length(fin_sp)){
  spdat<-subset(ModifiedDate,Scientific.Name==fin_sp[j])
  fit<-lm(formula = modifiedDOY ~ BreedingENSO+elevation+distToCoast+lat,data = spdat)
  summary(fit)
  #fit<-lm(formula = elevation ~ BreedingENSO,data = spdat)
  #qqPlot(residuals(lmAllSp),main=("start"))# you don't need to worry about residuals because we are not using this for prediction
  fitAnoc<-Anova(fit) # anova table
  p<-round(Anova(fit)$"Pr(>F)"[1],digits = 3)#get p values from lmer
  #text(130,-40,p)
  modAnova<-Anova(fit, test.statistic="F")
  f<-modAnova$F[1]
  df<-modAnova$Df
}
  
visreg(fit,"BreedingENSO",type="conditional")  
visreg(fit,"elevation",by="BreedingENSO")







qqPlot(residuals(lmAllSp),main=("start"))
Anova(lmAllSp, test.statistic="F")

##########################################
############ step 5. Are there differences in Breeding period lengths? ########
##########################################

levels(allDat$Phase) <- abbreviate(levels(allDat$Phase))
allDat$Phase <- with(allDat, reorder(Phase,
                            BreedingPeriod, mean))
lmBP<-lmer(formula = BreedingPeriod ~ Phase+ (1+Phase|family),REML=T, data = allDat)
#anova(lmBP,lmO) - #includeing Order did not make significant difference
#set up general linear hypothesis
qqPlot(residuals(lmBP),main=("start"))
Anova(lmBP, test.statistic="F")


TukeyRegion2<- glht(lmBP, linfct=mcp(Phase="Tukey"))
summary(TukeyRegion2)
lets <- cld(TukeyRegion2)$mcletters$Letters
library(doBy)
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(BreedingPeriod ~ Phase, data=allDat,
                FUN=c(mean,SE,sd,min,max))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)



##########################################
############ step 6. Are there differences in start dates? ########

#NO SPECIES HAVE START DATES THAT CROSS OVER THE YEAR,
#DO NOT NEED TO WORRY ABOUT CIRCULAR NATURE OF DATA
##########################################

b<-subset(allDat,Phase=="LaNn",c(Quantile5,Species))
a<-subset(allDat,Phase=="elNn",c(Quantile5,Species))
colnames(a)[1]<-c("el")
colnames(b)[1]<-c("la")
c<-merge(a,b,by="Species")
c$diff<-c$el-c$la


allDat$Phase <- with(allDat, reorder(Phase,
                                     Quantile5, mean))
lmstart<-lmer(formula = Quantile5 ~ Phase + (1+Phase|family), data = allDat,REML=T)
#anova(lmstart,lmstart2)
qqPlot(residuals(lmstart),main=("start"))

Anova(lmstart, test.statistic="F")
#set up general linear hypothesis
TukeyRegion2<- glht(lmstart, linfct=mcp(Phase="Tukey"))
lets <- cld(TukeyRegion2)$mcletters$Letters
library(doBy)
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(Quantile5 ~ Phase, data=allDat,
                FUN=c(mean,SE))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
b <- with(mn, barplot2(Quantile5.mean,
                       names.arg=Phase,
                       #main=spSum[i],
                       plot.ci=TRUE, ylim=c(0,300),
                       ci.l=Quantile5.mean - Quantile5.SE,
                       ci.u=Quantile5.mean + Quantile5.SE))
text(b, 280, lets)



##########################################
############ step 7. Are there differences in end dates? ########
##########################################

#modify date
allDat$MODIFIED_Quantile95<-with(allDat,ifelse(Quantile95<Quantile5,Quantile95+365,Quantile95))
#check that dates make sense
b<-subset(allDat,Phase=="LaNn",c(MODIFIED_Quantile95,Species))
a<-subset(allDat,Phase=="elNn",c(MODIFIED_Quantile95,Species))
colnames(a)[1]<-c("el")
colnames(b)[1]<-c("la")
c<-merge(a,b,by="Species")
c$diff<-c$el-c$la #modified dates so no longer are circular in nature


allDat$Phase <- with(allDat, reorder(Phase,
                                     MODIFIED_Quantile95, mean))
lmEnd<-lmer(formula = MODIFIED_Quantile95 ~ Phase+ (1+Phase|family), data = allDat,REML = TRUE) 
lmEnd2<-lmer(formula = MODIFIED_Quantile95 ~ Phase+ (1+Phase|family) + (1+Phase|Order), data = allDat,REML=TRUE)
anova(lmEnd,lmEnd2) #lmEnd is better model
Anova(lmEnd, test.statistic="F")
qqPlot(residuals(lmEnd),main=("mod1"))

#set up general linear hypothesis
TukeyRegion2<- glht(lmEnd, linfct=mcp(Phase="Tukey"))
lets <- cld(TukeyRegion2)$mcletters$Letters
library(doBy)
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(MODIFIED_Quantile95 ~ Phase, data=allDat,
                FUN=c(mean,SE,sd,min,max))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
b <- with(mn, barplot2(MODIFIED_Quantile95.mean,
                       names.arg=Phase,
                       #main=spSum[i],
                       plot.ci=TRUE, ylim=c(0,360),
                       ci.l=MODIFIED_Quantile95.mean - MODIFIED_Quantile95.SE,
                       ci.u=MODIFIED_Quantile95.mean + MODIFIED_Quantile95.SE))
text(b, 280, lets)


##########################################
############ step 8. Are there differences in rbar? ########
##########################################

#diff in rBar, el Nino less clumped
allDat$Phase <- with(allDat, reorder(Phase,
                                     Rbar90, mean))
lmClades<-lm(formula = Rbar90 ~ Phase, data = allDat)
anova(lmClades)
#set up general linear hypothesis
TukeyRegion2<- glht(lmClades, linfct=mcp(Phase="Tukey"))
lets <- cld(TukeyRegion2)$mcletters$Letters
library(doBy)
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(Rbar90 ~ Phase, data=allDat,
                FUN=c(mean,SE))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
b <- with(mn, barplot2(Rbar90.mean,
                       names.arg=Phase,
                       #main=spSum[i],
                       plot.ci=TRUE, ylim=c(0,1),
                       ci.l=Rbar90.mean - Rbar90.SE,
                       ci.u=Rbar90.mean + Rbar90.SE))
text(b, .9, lets)

##########################################
############ step 8. Are there differences in skew? ########
##########################################

allDat$Phase <- with(allDat, reorder(Phase,
                                     Skew90, mean))
lmClades<-lm(formula = Skew90 ~ Phase, data = allDat)
anova(lmClades)
#set up general linear hypothesis
TukeyRegion2<- glht(lmClades, linfct=mcp(Phase="Tukey"))
lets <- cld(TukeyRegion2)$mcletters$Letters
library(doBy)
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(Skew90 ~ Phase, data=allDat,
                FUN=c(mean,SE))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
b <- with(mn, barplot2(Skew90.mean,
                       names.arg=Phase,
                       #main=spSum[i],
                       plot.ci=TRUE, ylim=c(-1,1),
                       ci.l=Skew90.mean - Skew90.SE,
                       ci.u=Skew90.mean + Skew90.SE))
text(b, .9, lets)


##########################################
############ step 9. how many species have significant differnces between el and la? ########
##########################################


a<-subset(allDat,Phase=="elNn",select=c(Species,NoObs,BreedingPeriod,Quantile5,Quantile95))
colnames(a)<-c("Species","NoObsEL","BreedingPeriod_elNn","Quantile5_elNn","Quantile95_elNn")
b<-subset(allDat,Phase=="LaNn",select=c(Species,NoObs,BreedingPeriod,Quantile5,Quantile95))
colnames(b)<-c("Species","NoObsLa","BreedingPeriod_LaNn","Quantile5_LaNn","Quantile95_LaNn")
supplTAB<-merge(a,b,by = "Species")

#get data and species
ELDat<-subset(alldat90,SOI=="elNino"|SOI=="LaNina")
ELsp<-as.vector(unique(ELDat$Scientific.Name))
spPval<-list()
for (i in 1:length(ELsp)){
  #modify dates
  startEL<-subset(supplTAB,Species==ELsp[i],Quantile5_elNn)[,1]
  startLA<-subset(supplTAB,Species==ELsp[i],Quantile5_LaNn)[,1]
  ELsub<-subset(ELDat,Scientific.Name==ELsp[i] & SOI=="elNino")
  ELLat<-mean(ELsub$lat)
  ELLon<-mean(ELsub$lon)
  LAsub<-subset(ELDat,Scientific.Name==ELsp[i] & SOI=="LaNina")
  LALat<-mean(LAsub$lat)
  LAlon<-mean(LAsub$lon)
  ELsub$MOD_PL<-with(ELsub,ifelse(DOY_PL<startEL,DOY_PL+365,DOY_PL))
  LAsub$MOD_PL<-with(LAsub,ifelse(DOY_PL<startLA,DOY_PL+365,DOY_PL))
  ELLA<-rbind(ELsub,LAsub)
  ELLA$SOI<-factor(ELLA$SOI)
  #get pvals
  lmDay<-lm(formula = MOD_PL ~ SOI, data = ELLA)
  lmLat<-lm(formula = lat ~ SOI, data = ELLA)
  lmLon<-lm(formula = lon ~ SOI, data = ELLA)
  aDay<-Anova(lmDay)
  aLat<-Anova(lmLat)
  aLon<-Anova(lmLon)
  pDay<-formatPval(aDay$'Pr(>F)',3,eps=0.05)[1]
  pLat<-formatPval(aLat$'Pr(>F)',3,eps=0.05)[1]
  pLon<-formatPval(aLon$'Pr(>F)',3,eps=0.05)[1]
  
  spPval[[i]]<-data.frame(t(c(ELsp[i], pDay)))
  
}
  
SIG<-do.call("rbind",spPval)
colnames(SIG)<-c("Species","p")


supplTAB<-merge(supplTAB,SIG,by="Species")


#get data ready for supplementary information, make sure names are up to date. 
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
#keep only species of interest
inc<-subset(inc, remove !=1)
taxon<-inc[c("genus","family","ScientificName","Common.me")]
s1<-merge(supplTAB,taxon,all.x=TRUE,by.x="Species",by.y="ScientificName")

#check names are okay with Australian Bird Data Version 1
birdData<-read.csv("/Users/daisy/Google Drive/PhD/birdTraits/PublishedData/Australian_Bird_Data_Version_1.csv")
#keep only full species
birdData2<-subset(birdData, is.na(X6_Subspecies_name_2))[c(1:5)]

birdData2$Species<-with(birdData2,paste(X4_Genus_name_2,X5_Species_name_2,sep=" "))
s2<-merge(s1,birdData2,all.x=TRUE,by="Species")
bad<-subset(s2,is.na(X4_Genus_name_2))
good<-subset(s2,!is.na(X4_Genus_name_2))

good<-good[,c("Species",
            "Common.me",
            "NoObsEL",
            "NoObsLa",
            "BreedingPeriod_elNn",
            "Quantile5_elNn",
            "BreedingPeriod_LaNn",
            "Quantile5_LaNn",
            "p")]

i <- sapply(good, is.factor)
good[i] <- lapply(good[i], as.character)


bad<-merge(bad,birdData2,all.x=TRUE,by.x="Common.me",by.y="X3_Taxon_common_name_2")
bad$Species<-bad$Species.y
bad<-bad[,c("Species",
            "Common.me",
            "NoObsEL",
            "NoObsLa",
            "BreedingPeriod_elNn",
            "Quantile5_elNn",
            "BreedingPeriod_LaNn",
            "Quantile5_LaNn",
            "p")]
#get rid of factors
i <- sapply(bad, is.factor)
bad[i] <- lapply(bad[i], as.character)

all<-rbind(good,bad)                                                                                                                        

write.csv(all,"/Users/daisy/Google Drive/PhD/BreedingTiming/manuscript/TablesFigures/SupplementaryData2.csv")






#write.csv(supplTAB,"/Users/daisy/Google Drive/PhD/BreedingTiming/manuscript/TablesFigures/BreedingPeriod20160105.csv")






##########################################
############ step 11. Figure - Breeding Period by Phase   ########
##########################################

# library(RColorBrewer)
# color.gen<-c(brewer.pal(n=11, name="RdYlBu"),
#              brewer.pal(n=11, name="RdBu"),
#              brewer.pal(n=12, name="Paired"),
#              brewer.pal(n=11, name="PuOr"),
#              brewer.pal(n=11, name="BrBG"))


boxSlopes <- function(condition,group,labelY,sig,sigsize,figNumber,Sp,dat,axesTF,at,labels){
  wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
              ylab=labelY, 
              main="",names=c("",""),outcol="white", 
              boxwex=0.5,border = c("black","black"),
              cex.lab=1.3, axes=axesTF)
  box()
  axis(1, at=1:2,labels=FALSE)
  mtext(expression(El~Ni*tilde(n)*o), side=1, line=1, at=1, cex=.8)
  mtext(expression(La~Ni*tilde(n)*a), side=1, line=1, at=2, cex=.8)
  if(axesTF==FALSE){
    axis(2, at=at, 
         labels=labels)
  }
  mtext(figNumber, side=3, line=1, at=1.5,  adj=8,cex=.8,font=2)
  #add line for each species
  for (i in 1:length(Sp)) {
    el<-subset(dat,Species==Sp[i] & Phase =="elNn",select = paste(condition))[,1]
    la<-subset(dat,Species  == Sp[i] & Phase =="LaNn",select = condition)[,1]
    color="grey40"
    #color=color.gen[i] #for coloured lines
    #segments(1,slopes.gen$intercept[i]+slopes.gen$slope[i],2,slopes.gen$intercept[i], col=color, lwd=0.7, lty=1)
    segments(1,el,2,la, col=color, lwd=0.7, lty=1)
  }
  #add in the box plots again so that they are not covered by the lines
  wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
  ylab=labelY, 
  main="",names=c("",""), 
  boxwex=0.5,border = c("black","black"),
  cex.lab=1.3, axes=axesTF, add = TRUE, col ="#FFFFFF00" )
  mtext(sig, side=3, line=-2.5, at=1.5,  adj=0.5,cex=sigsize)
}

#prep data
allDat$Species<-lapply(allDat$Species, as.character)
SpeciesTest<-unlist(allDat$Species)
allDat$Genus<-sapply(strsplit(SpeciesTest, " "), "[", 1)
LEallDat<-droplevels(subset(allDat,Phase=="elNn" | Phase =="LaNn"))

#reorder data so elNina is first
LEallDat$Phase <- with(LEallDat, reorder(Phase,
                                     BreedingPeriod, mean))
#set up PDF
pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/BreedingPeriodByPhase_Significant20160424.pdf",
    width = 8, height = 4)

par(mfrow=c(1,3))


#Three figures 1 row 3 columns, breeding period, start, conclusion
boxSlopes(condition="BreedingPeriod",
          group="Phase",
          labelY="ELP (day)",
          sig="*",
          sigsize=1.5,
          figNumber="(a)",
          Sp=unlist(unique(LEallDat$Species)),
          dat=LEallDat,
          axesTF=TRUE)
         
boxSlopes(condition="Quantile5",
          group="Phase",
          labelY="Start of ELP",
          sig="n.s.",
          sigsize=.9,
          figNumber="(b)",
          Sp=unlist(unique(LEallDat$Species)),
          dat=LEallDat,
          axesTF=FALSE,
          at=c(1,32,61,92,121,152,182,214,245,275,335), 
          labels=c("Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))


#modify date
LEallDat$MODIFIED_Quantile95<-with(LEallDat,ifelse(Quantile95<Quantile5,Quantile95+365,Quantile95))
boxSlopes(condition="MODIFIED_Quantile95",
          group="Phase",
          labelY="Termination of ELP",
          sig="*",
          sigsize=1.5,
          figNumber="(c)",
          Sp=unlist(unique(LEallDat$Species)),
          dat=LEallDat,
          axesTF=FALSE,
          at=c(245,275,305,335,366,397,425,456,487,518,548), 
          labels=c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","July"))

# LEallDat$MODIFIED_Quantile50<-with(LEallDat,ifelse(Quantile50<Quantile5,Quantile50+365,Quantile50))
# 
# boxSlopes(condition="MODIFIED_Quantile50",
#           group="Phase",
#           labelY="End dates",
#           sig="*",
#           sigsize=1.5,
#           figNumber="(c)",
#           Sp=unlist(unique(LEallDat$Species)),
#           dat=LEallDat,
#           axesTF=FALSE,
#           at=c(245,275,305,335,366,397,425,456,487,518,548), 
#           labels=c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","July"))



dev.off()




##########################################
############ step 12. Figure for Species on the move- Breeding Period by Phase   ########
##########################################

boxSlopes <- function(condition,group,labelY,sig,sigsize,figNumber,Sp,dat,axesTF,at,labels){
  wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
              ylab=labelY, 
              main="",names=c("",""),outcol="white", 
              boxwex=0.5,border = c("black","black"),
              cex.lab=1, axes=axesTF)
  box()
  axis(1, at=1:2,labels=FALSE)
  mtext(expression(El~Ni*tilde(n)*o), side=1, line=1, at=1, cex=1)
  mtext(expression(La~Ni*tilde(n)*a), side=1, line=1, at=2, cex=1)
  if(axesTF==FALSE){
    axis(2, at=at, 
         labels=labels)
  }
  mtext(figNumber, side=3, line=1, at=1.5,cex=1)
  #add line for each species
  
  el<-subset(dat,Phase =="elNn",select = paste(condition))[,1]
  la<-subset(dat,Phase =="LaNn",select = condition)[,1]
  absDiff<-abs(el-la)
  
  dat2<-as.data.frame(cbind(el,la, absDiff))
 dat2<-dat2[order(absDiff),] 
  
  
  for (i in 1:nrow(dat2)) {
    el<-dat2[i,1]#subset(dat,Species==Sp[i] & Phase =="elNn",select = paste(condition))[,1]
    la<-dat2[i,2]#subset(dat,Species  == Sp[i] & Phase =="LaNn",select = condition)[,1]
    
    color=ifelse(la-el>=20 | el-la>=20, "dodgerblue3", "grey60")
    #color=color.gen[i] #for coloured lines
    #segments(1,slopes.gen$intercept[i]+slopes.gen$slope[i],2,slopes.gen$intercept[i], col=color, lwd=0.7, lty=1)
    segments(1,el,2,la, col=color, lwd=0.8, lty=1)
  }
  #add in the box plots again so that they are not covered by the lines
  wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
              ylab=labelY, 
              main="",names=c("",""), 
              boxwex=0.5,border = c("black","black"),
              cex.lab=1, axes=axesTF, add = TRUE, col ="#FFFFFF00" )
  mtext(sig, side=3, line=-2.5, at=1.5,  adj=0.5,cex=sigsize)
}

#prep data
allDat$Species<-lapply(allDat$Species, as.character)
SpeciesTest<-unlist(allDat$Species)
allDat$Genus<-sapply(strsplit(SpeciesTest, " "), "[", 1)
LEallDat<-droplevels(subset(allDat,Phase=="elNn" | Phase =="LaNn"))

#reorder data so elNina is first
LEallDat$Phase <- with(LEallDat, reorder(Phase,
                                         BreedingPeriod, mean))
#set up PDF
pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/SpeciesOnTheMoveBreedingPeriodByPhase_Significant.pdf",
   width = 8, height = 4)

par(mfrow=c(1,2))


#Three figures 1 row 2 columns, breeding period, start
boxSlopes(condition="BreedingPeriod",
          group="Phase",
          labelY="Length (days)",
          sig="*",
          sigsize=1.5,
          figNumber="Egg-laying Period",
          Sp=unlist(unique(LEallDat$Species)),
          dat=LEallDat,
          axesTF=TRUE)

boxSlopes(condition="Quantile5",
          group="Phase",
          labelY="Date",
          sig="n.s.",
          sigsize=.9,
          figNumber="Start date",
          Sp=unlist(unique(LEallDat$Species)),
          dat=LEallDat,
          axesTF=FALSE,
          at=c(1,32,61,92,121,152,182,214,245,275,335), 
          labels=c("Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))





dev.off()










##########################################
############ Step 8. What are differences between La Nina and El Nino ########
##########################################
supplTAB$MODIFIED_LaQuantile95<-with(supplTAB,ifelse(Quantile95_LaNn<Quantile95_elNn & Quantile95_elNn >= 100 & Quantile95_LaNn <= 100,Quantile95_LaNn+365,Quantile95_LaNn))
supplTAB$start.diff=with(supplTAB,Quantile5_elNn-Quantile5_LaNn)
supplTAB$end.diff=with(supplTAB,Quantile95_elNn-MODIFIED_LaQuantile95)
supplTAB$period.diff=with(supplTAB,BreedingPeriod_elNn-BreedingPeriod_LaNn)


#how many sppecies with significant differences have breeding periods at least 20 days longer during La Nina?
nrow(subset(supplTAB,period.diff<=-20)) #23
nrow(subset(supplTAB,period.diff<=-10)) #35
nrow(subset(supplTAB,period.diff>=20)) #0 El Nina
nrow(subset(supplTAB,period.diff>=10))

#how many have little change
nrow(subset(supplTAB,period.diff<10 & period.diff>-10))
a<-subset(supplTAB,period.diff<10 & period.diff>-10)
#how many species with significant differences have end dates at least 20 days later during La Nina?
nrow(subset(supplTAB,end.diff<=-20))  #25
nrow(subset(supplTAB,end.diff<=-10))  #32
#how many species with significant differences have end dates at least 10 days earlier during La Nina?
nrow(subset(supplTAB,end.diff>=10 & p =="< 0.05"))    #0

#how many species with significant differences have start dates at least 20 earlier during La Nina?
nrow(subset(supplTAB,start.diff<=-20)) #1
nrow(subset(supplTAB,start.diff<=-10)) #8

#how many species with significant differences have start dates at least 20 earlier during El Nina?
nrow(subset(supplTAB,start.diff>=20)) #3
nrow(subset(supplTAB,start.diff>=10)) #11











##########################################
############ Step 8. How much earlier is the breeding period? ########
##########################################

library(visreg)
library(reporttools)
#
mod<-lmer(DOY_PL ~ SOI + (1+SOI|Scientific.Name),data=ELDat)
coef(mod)$Scientific.Name
summary(mod)$coefficients[1,1]  #predicted value outcrossers
summary(mod)$coefficients[1,1]  + summary(mod)$coefficients[2,1] #predicted value selfers 
t=list(mod)
rsquared.glmm(list(mod)) #proportion of variance explained by mating system and random effects (marginal,conditional)
coef(mod)$Scientific.Name
modANOVA<-Anova(mod) #best way to get p values from lmer

modANOVA$'Pr(>Chisq)'
#visreg(mod)

formatPval(modANOVA$'Pr(>Chisq)',3,eps=0.05)


mod=lmer(BreedingPeriod ~ Phase + (1+Phase|Clade),data=LEallDat)
summary(mod)$coefficients[1,1]  #predicted value outcrossers
summary(mod)$coefficients[1,1]  + summary(mod)$coefficients[2,1] #predicted value selfers 
t=list(mod)
rsquared.glmm(list(mod)) #proportion of variance explained by mating system and random effects (marginal,conditional)
coef(mod)$Clade

