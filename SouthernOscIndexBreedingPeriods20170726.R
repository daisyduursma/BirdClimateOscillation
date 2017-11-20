##########################################
############ step 1. get set up to work ########
##########################################

##########################################
############ step 1. get set up to work ########
##########################################

rm(list = ls())
library(ggplot2)
library(doBy)
library(visreg)
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
library(MuMIn)
#library(afex)
library(ENmisc)
#library(arm)
library(car)
#library(heavy)
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
#Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,13-Grassland,3-Temperate
TEMP<-subset(dat, koeppen == 22 |koeppen == 13) #temperate region Breeding Period
TEMP$month<-formatC(TEMP$month,width=2,flag='0') #format months 1 becomes 01
TEMP$date<-paste0(TEMP$year,TEMP$month) #make date string that matches NOAA SOI
TEMP<-subset(TEMP, date >=190001) #limit time to same as NOAA data
#read in BOM ENSO events
#BOM
SOI<-read.csv("/Users/daisy/Google Drive/PhD/ENSO/Data/BOM_SOItest.csv")
SOI$month2<-formatC(SOI$monthNumeric,width=2,flag='0') #format months 1 becomes 01
SOI$date<-paste0(SOI$Year,SOI$month2) #make date string that matches NOAA SOI
SOI<-SOI[,c("date","BreedingENSO","Year","monthNumeric","month2","X")]
#how many breeding seasons in each phase
July<-subset(SOI,month2=="07")
table(July$BreedingENSO)

TEMP<-merge(TEMP,SOI, by.x="date",by.y="date",all.x=TRUE)
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
nN<-subset(TEMP,BreedingENSO=='neutral')
spSumnN<-as.data.frame(table(nN$Scientific.Name))
sp_nN<-as.vector(droplevels(subset(spSumnN,Freq>=100)$Var1))
#name of species in both, limit data to these species
fin_sp<-intersect(sp_eN,sp_lN)
fin_sp<-intersect(fin_sp,sp_nN)#neutral
eN<-droplevels(eN[eN$Scientific.Name %in% fin_sp,])
eN$SOI<-as.character(eN$SOI)
eN$SOI<-"elNino"
lN<-lN[lN$Scientific.Name %in% fin_sp,]
lN$SOI<-as.character(lN$SOI)
lN$SOI<-"LaNina"
nN<-nN[nN$Scientific.Name %in% fin_sp,]
nN$SOI<-as.character(nN$SOI)
nN$SOI<-"Neutral"

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
#get Neutral species with more than 100 Egg or young or multi observations
dat_NU_Acc<-subset(nN, type == "solitary-egg" | type == "multi"| type == "young")
dat_NU_Acc_sp<-as.data.frame(table(dat_NU_Acc$Scientific.Name))
dat_NU_Acc_sp<-as.vector(droplevels(subset(dat_NU_Acc_sp,Freq>=100)$Var1))


#species in both la nina and el nino
fin_sp_Acc<-intersect(dat_eN_Acc_sp,dat_LA_Acc_sp) #25 species in total
fin_sp_Acc<-intersect(fin_sp_Acc,dat_NU_Acc_sp)
#get the high accuracy data
dat_eN_Acc<-datAcc[datAcc$Scientific.Name %in% fin_sp_Acc,]
dat_LA_Acc<-dat_LA_Acc[dat_LA_Acc$Scientific.Name %in% fin_sp_Acc,]
dat_NU_Acc<-dat_NU_Acc[dat_NU_Acc$Scientific.Name %in% fin_sp_Acc,]
#put the data back together
high<-rbind(dat_LA_Acc,dat_eN_Acc)#
high<-rbind(high,dat_NU_Acc)
high$Acc<-"high"
#get the data for the rest of the species
dat_eN_LOW<-eN[eN$Scientific.Name %nin% fin_sp_Acc,]
length(unique(dat_eN_LOW$Scientific.Name))
dat_LA_LOW<-lN[lN$Scientific.Name %nin% fin_sp_Acc,]
length(unique(dat_LA_LOW$Scientific.Name))
dat_Nu_LOW<-nN[nN$Scientific.Name %nin% fin_sp_Acc,]
length(unique(dat_Nu_LOW$Scientific.Name))

low<-rbind(dat_eN_LOW,dat_LA_LOW)#,
low<-rbind(low,dat_Nu_LOW)
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
state<-c("elNino","LaNina","Neutral")
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


#50% of the species tested as having significantly differnt distributions


#Of the species that were significant, 
sigdat<-subset(modDat,datepVal=="< 0.05" |datepVal=="< 0.01"|datepVal=="< 0.001")
#how many had shorter beeding periods during el nino
short<-with(sigdat,BreedingPeriod -LNBreedingPeriod)
mean(short[short < 0])#how much shorter
sd(short[short < 0])
mean(short[short > 0])#how much longer
sd(short[short >0])
#how many species had earler peaks during El nino
peak<-with(sigdat,Quantile50 -LNQuantile50)
mean(peak[peak < 0])#how much shorter
sd(peak[peak < 0])
mean(peak[peak > 0])#how much shorter
sd(peak[peak > 0])

start<-sort(with(sigdat,Quantile5 -LNQuantile5))
mean(start[start < 0])#how much shorter
sd(start[start < 0])
mean(start[start > 0])#how much shorter
sd(start[start > 0])

ccls<-sort(with(sigdat,(Quantile5+BreedingPeriod) -(LNQuantile5+LNBreedingPeriod)))
mean( ccls[ ccls < 0])#how much shorter
sd( ccls[ ccls < 0])
mean( ccls[ ccls > 0])#how much shorter
sd( ccls[ ccls > 0])









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


ModifiedDate2<-list()
for(m in 1:length(fin_sp)){
  spdat<-subset(PhaseDat,Scientific.Name==fin_sp[m])
  earlstart<-min(subset(allDat,Species==fin_sp[m])$Quantile5)
  latestend<-max(subset(allDat,Species==fin_sp[m])$Quantile95)
  #adjust dates so those at the begining of the year come after 365
  
  df<-subset(spdat,DOY_PL<earlstart)#dates to adjust
  df$modifiedDOY<-df$DOY_PL+365
  df2<-subset(spdat,DOY_PL>earlstart)#dates that are fine
  df2$modifiedDOY<-df2$DOY_PL
  ModifiedDate2[[m]]<-rbind(df,df2)
}

PhaseDat<-do.call("rbind",ModifiedDate2)

##########################################
############ Are there differences in mean during El Nino and La Nina? 
##########################################


ModifiedDate$Phase <- factor(as.character(ModifiedDate$BreedingENSO), levels = c("NINO", "NINA", "neutral"))

ModifiedDate$elevation[ModifiedDate$elevation <= 0] <- 1
ModifiedDate$elevation[is.na(ModifiedDate$elevation) ] <- 0.01
ModifiedDate$elevation<-ModifiedDate$elevation+1

ModifiedDate$distToCoast[ModifiedDate$distToCoast <= 0] <- 1
ModifiedDate$distToCoast[is.na(ModifiedDate$distToCoast)] <-1
ModifiedDate$distToCoast<-ModifiedDate$distToCoast+1

ModifiedDate$logelev <- log10(ModifiedDate$elevation)
ModifiedDate$logdist <- log10(ModifiedDate$distToCoast)


#FULL - all fixed effects and random model
mF<-lmer(modifiedDOY ~ BreedingENSO + logelev*logdist + lat+(1|Acc)+(1|Scientific.Name)+(1|Order),REML=T, data = ModifiedDate )
summary(mF)
qqPlot(residuals(mF),main=("DOY"))
Anova(mF,test.statistic="Chisq")#get p value
anova(mF)#get F statement
r.squaredGLMM(mF)


TukeyENSO<- glht(mF, linfct=mcp(BreedingENSO="Tukey"))
summary(TukeyENSO)
# summaryBy(ChangeBPEL~season, data=wide,
#           FUN=c(mean,sd)) #find mean and sd
#*model shows significant variation but not between El Nina and Neutral


#models calculated with maximum likelihood to get AIC BIC values
m0ML<-lmer(modifiedDOY ~1 + (1|Scientific.Name) + (1|Order),REML=F, data = ModifiedDate )
mFML<-lmer(modifiedDOY ~ BreedingENSO + lat+ logelev+logdist+ (1|Scientific.Name)+(1|Order),REML=F, data = ModifiedDate )
summary(mFML)
Anova(mFML,test.statistic="Chisq")


######################################
#species that meet El Nino Hypothesis
######################################

#change direction of data
wide<-reshape(allDat, idvar = c("Species","CommonName","Order","family") , timevar = "Phase", direction = "wide")
#modify end dates, account for those species that breed over turn of year
wide$modified95<-ifelse( wide$Quantile95.elNino<200, wide$Quantile95.elNino+365, wide$Quantile95.elNino)
wide$modified95LA<-ifelse( wide$Quantile95.LaNina<150, wide$Quantile95.LaNina+365, wide$Quantile95.LaNina)
wide$modified95NEU<-ifelse( wide$Quantile95.Neutral<175,wide$Quantile95.Neutral+365, wide$Quantile95.Neutral)


#We hypothesised that there would be earlier peaks and terminations 
#of the breeding period during El Nino and the breeding periods would be shorter. 
wide$earlierELstart<-with(wide,Quantile5.elNino-Quantile5.LaNina)
wide$earlierELpeak<-with(wide,Quantile50.elNino-Quantile50.LaNina)
#modifiy dates for end
wide$earlierELend<-with(wide,modified95-modified95LA)
wide$shorterBPEL<-with(wide,BreedingPeriod.elNino-BreedingPeriod.LaNina)
wide$ChangeBPEL<-with(wide,(((BreedingPeriod.elNino-BreedingPeriod.LaNina)/BreedingPeriod.elNino)*100))
#number of species that meet hypothesis
nrow(subset(wide,earlierELend < 0 & earlierELpeak < 0&  shorterBPEL <0))#La Nina


#for neutral
wide$earlierELstartNeu<-with(wide,Quantile5.elNino-Quantile5.Neutral)
wide$earlierELendNeu<-with(wide,modified95-modified95NEU)
wide$earlierELpeakNeu<-with(wide,Quantile50.elNino-Quantile50.Neutral)
wide$shorterBPELNeu<-with(wide,BreedingPeriod.elNino-BreedingPeriod.Neutral)
wide$ChangeBPELNeu<-with(wide,(((BreedingPeriod.elNino-BreedingPeriod.Neutral)/BreedingPeriod.elNino)*100))
nrow(subset(wide,earlierELendNeu < 0 & earlierELpeakNeu < 0&  shorterBPELNeu <0))
wide$ChangeBPLANeu<-with(wide,(((BreedingPeriod.LaNina-BreedingPeriod.Neutral)/BreedingPeriod.LaNina)*100))

#number of species that meet hypothesis
nrow(subset(wide,earlierELend < 0 & earlierELpeak < 0&  shorterBPEL <0))#La Nina
nrow(subset(wide,earlierELendNeu < 0 & earlierELpeakNeu < 0&  shorterBPELNeu <0))
#% species with shorter breeding periods
nrow(subset(wide, shorterBPEL <0))/60*100
nrow(subset(wide, shorterBPELNeu <0))/60*100
# xx% of species having later starts 
nrow(subset(wide, earlierELstart >0))/60*100
nrow(subset(wide, earlierELstartNeu >0))/60*100

#xx% of species having earlier peaks and
nrow(subset(wide, wide$earlierELpeak <0))/60*100
nrow(subset(wide, wide$earlierELpeakNeu <0))/60*100
            
#% of species having earlier terminations to egg-laying period.
nrow(subset(wide, wide$earlierELend <0))/60*100
nrow(subset(wide, wide$earlierELendNeu <0))/60*100

#During La Nina events, we hypothesised that favourable climatic 
#conditions should extend the breeding period  resulting in later peaks 
#and terminations (see Figure 1). 
#number of species that meet hypothesis
nrow(subset(wide,earlierELend > 0 & earlierELpeak >0 & shorterBPEL > 0 ))
nrow(subset(wide,earlierELendNeu > 0 & earlierELpeakNeu >0 & shorterBPELNeu > 0 ))






# ######################################
# #species that meet El Nino Hypothesis
# ######################################
# 
# 
# wide<-reshape(allDat, idvar = c("Species","CommonName","Order","family") , timevar = "Phase", direction = "wide")
# 
# #We hypothesised that there would be earlier peaks and terminations 
# #of the breeding period during El Nino and the breeding periods would be shorter. 
# 
# wide$earlierELstart<-with(wide,Quantile5.elNino-Quantile5.LaNina)
# wide$earlierELend<-with(wide,Quantile95.elNino-Quantile95.LaNina)
# wide$earlierELpeak<-with(wide,Quantile50.elNino-Quantile50.LaNina)
# wide$shorterBPEL<-with(wide,BreedingPeriod.elNino-BreedingPeriod.LaNina)
# wide$ChangeBPEL<-with(wide,(((BreedingPeriod.elNino-BreedingPeriod.LaNina)/BreedingPeriod.elNino)*100))
# 
# #number of species that meet hypothesis
# nrow(subset(wide,earlierELend < 0 & earlierELpeak < 0&  shorterBPEL <0))
# #% species with shorter breeding periods
# nrow(subset(wide, shorterBPEL <0))/60*100
# # xx% of species having later starts 
# nrow(subset(wide, earlierELstart >0))/60*100
# 
# #xx% of species having earlier peaks and
# nrow(subset(wide, wide$earlierELpeak <0))/60*100
# 
# #% of species having earlier terminations to egg-laying period.
# nrow(subset(wide, wide$earlierELend <0))/60*100
# 
# #During La Nina events, we hypothesised that favourable climatic 
# #conditions should extend the breeding period  resulting in later peaks 
# #and terminations (see Figure 1). 
# #number of species that meet hypothesis
# nrow(subset(wide,earlierELend > 0 & earlierELpeak >0 & shorterBPEL > 0 ))
# nrow(subset(wide,shorterBPEL < 0 ))



# 
# ##########################################
# #% change in egglaying period from neutral
# ##########################################
# 
# ELchange<-data.frame(change=with(wide,(((BreedingPeriod.elNino-BreedingPeriod.Neutral)/BreedingPeriod.elNino)*100)),
#                      SOI="NINO")
# LAchange<-data.frame(change=with(wide,(((BreedingPeriod.LaNina-BreedingPeriod.Neutral)/BreedingPeriod.LaNina)*100)),
#                      SOI="NINA")
# LAELchange<-data.frame(change=with(wide,(((BreedingPeriod.elNino-BreedingPeriod.LaNina)/BreedingPeriod.elNino)*100)),
#                      SOI="NINONINA")
# dfChange<-rbind(ELchange,LAchange)
# dfChange<-rbind(dfChange,LAELchange)
# 
# lmChange<-lm(change~SOI ,data=dfChange)
# summary(lmChange)#end data explains about 16% if variation we observe, with those concluding earlier in the year
# 
# beanplot(change~SOI,data=dfChange, ll = 0.04, main = NA,
#          col="white",method="jitter",
#          ylab=expression(paste("%",Delta)),
#          show.names=FALSE)
# 
# 
# 
# 
# 
# ##########################################
# # change start in egglaying period from neutral
# ##########################################
# ELstart<-data.frame(start=with(wide,Quantile5.Neutral-Quantile5.elNino),
#                      SOI="NINO")
# LAstart<-data.frame(start=with(wide,Quantile5.Neutral-Quantile5.LaNina),
#                      SOI="NINA")
# dfstart<-rbind(ELstart,LAstart)
# 
# lmstart<-lm(start~SOI ,data=dfstart)
# summary(lmstart)
# 
# ##########################################
# # change peak in egglaying period from neutral
# ##########################################
# ELpeak<-data.frame(peak=with(wide,Quantile50.Neutral-Quantile50.elNino),
#                     SOI="NINO")
# LApeak<-data.frame(peak=with(wide,Quantile50.Neutral-Quantile50.LaNina),
#                     SOI="NINA")
# dfpeak<-rbind(ELpeak,LApeak)
# 
# lmpeak<-lm(peak~SOI ,data=dfpeak)
# summary(lmpeak)



##########################################
# change end in egglaying period from neutral
# ##########################################
# ELend<-data.frame(end=with(wide,modified95NEU-modified95),
#                    SOI="NINO")
# LAend<-data.frame(end=with(wide,modified95NEU-modified95LA),
#                    SOI="NINA")
# dfend<-rbind(ELend,LAend)
# 
# lmend<-lm(end~SOI ,data=dfend)
# summary(lmend)


#Can we explain those species that have biggest change during El Nino
lmChange<-lm(ChangeBPEL~modified95NEU  ,data=wide)#El Nino La Nina
lmChange<-lm(ChangeBPEL~modified95  ,data=wide)#EL to Nu
summary(lmChange)#end data explains about 16% if variation we observe, with those concluding earlier in the year
Anova(lmChange)
visreg(lmChange)
qqPlot(residuals(lmChange),main=("DOY"))


lmChange<-lm(ChangeBPEL~modified95  ,data=wide)


###############################################

#No evidence that breeding periods are more skewed later in the year
#for either La Nina or El Nino breeding observations
###############################################

lm1<-lm(Skew90.LaNina~modified95NEU,data=wide)
summary(lm1)
lm2<-lm(Skew90.elNino~modified95,data=wide)
summary(lm2)

##########################################
#La Nina significnatly more skewed than either El Nino or newutral phase
#La Nina years are left skewed, and others right skewed

#########################################

skewLong<-rbind(data.frame(x=wide$Skew90.Neutral,y="Neutral"),
                rbind(data.frame(x=wide$Skew90.LaNina,y="LaNina"),
                      data.frame(x=wide$Skew90.elNino,y="ElNino")))
lm3<-lm(x~y  ,data=skewLong)
summary(lm3)
Tukeyskew<- glht(lm3, linfct=mcp(y="Tukey"))
summary(lm3)
lm4<-lm(Skew90.LaNina~Skew90.Neutral  ,data=wide)
summary(lm4)
lm5<-lm(Skew90.Neutral~Skew90.elNino  ,data=wide)
summary(lm5)
qqPlot(residuals(lm5))

summaryBy(x~y,data=skewLong)
###################################################
#seasonality of change in egg laying between El an LA
##################################################
#wide$season <- cut(wide$Quantile50.Neutral, 3, labels=c("Early", "Mid","Late"))
wide$season <- cut(wide$modified95, 3, labels=c("Early", "Mid","Late"))
fit.season<-lmer(ChangeBPEL~season+(1|Order) ,data=wide)
Anova(fit.season,test.statistic="F")
summary(fit.season)
TukeyRegion2<- glht(fit.season, linfct=mcp(season="Tukey"))
summary(TukeyRegion2)
visreg(fit.season,"season", type = "conditional")
summaryBy(ChangeBPEL~season, data=wide,
          FUN=c(mean,sd)) #find mean and sd

fit.seasonELtoNU<-lmer(ChangeBPELNeu~season+(1|Order) ,data=wide)
Anova(fit.seasonELtoNU, test.statistic="F")
#find differencs between regions
TukeyRegion2<- glht(fit.seasonELtoNU, linfct=mcp(season="Tukey"))
summary(TukeyRegion2)
summaryBy(ChangeBPELNeu~season, data=wide,
          FUN=c(mean,sd)) #find mean and sd

fit.seasonLAtoNU<-lmer(ChangeBPLANeu~season+(1|Order) ,data=wide)
Anova(fit.seasonLAtoNU, test.statistic="F")
visreg(fit.seasonLAtoNU,"season", type = "conditional")
#find differencs between regions
TukeyRegion2<- glht(fit.seasonLAtoNU, linfct=mcp(season="Tukey"))
summary(TukeyRegion2)
summaryBy(ChangeBPLANeu~season, data=wide,
          FUN=c(mean,sd)) #find mean and sd




#seasonal trend in length of egg-laying period
fit.season.length<-lm(BreedingPeriod.elNino~season ,data=wide)
Anova(fit.season.length, test.statistic="F")
#find differencs between regions
TukeyRegion2<- glht(fit.season.length, linfct=mcp(season="Tukey"))
summary(TukeyRegion2)
summaryBy(BreedingPeriod.elNino~season, data=wide,
          FUN=c(mean,sd)) #find mean and sd


#seasonal trend in length of egg-laying period
fit.season.lengthLA<-lm(BreedingPeriod.LaNina~season ,data=wide)
Anova(fit.season.lengthLA, test.statistic="F")
#find differencs between regions
TukeyRegion2<- glht(fit.season.lengthLA, linfct=mcp(season="Tukey"))
summary(TukeyRegion2)
summaryBy(BreedingPeriod.LaNina~season, data=wide,
          FUN=c(mean,sd)) #find mean and sd

#seasonal trend in length of egg-laying period
fit.season.lengthNu<-lmer(BreedingPeriod.Neutral~season+(1|Order) ,data=wide)
Anova(fit.season.lengthNu, test.statistic="F")
#find differencs between regions
TukeyRegion2<- glht(fit.season.lengthNu, linfct=mcp(season="Tukey"))
summary(TukeyRegion2)
summaryBy(BreedingPeriod.Neutral~season, data=wide,
          FUN=c(mean,sd)) #find mean and sdwq


la<-wide[c("season","BreedingPeriod.LaNina")]
colnames(la)<-c("Season","Length")
la$Season<-paste(la$Season,"1")

el<-wide[c("season","BreedingPeriod.elNino")]
colnames(el)<-c("Season","Length")
el$Season<-paste(el$Season,"2")

df1<-rbind(la,el)

df1$Season<-factor(df1$Season)
df1$Season = factor(df1$Season,levels(df1$Season)[c(2,1,6,5,4,3)])

library("vioplot")
library(beanplot)
library(doBy)

summaryBy(modified95~season,data=wide,FUN=c(min,max))

pdf(file = paste0("/Users/daisy/Google Drive/PhD/ENSO/Manuscript/GlobalChangeBiology/Figures/",as.Date(Sys.time()),".pdf"),
    width = 3.2, height = 6)
par(mfrow = c(2,1),
    mar = c(2,4,2,1))

beanplot(ChangeBPEL~season,data=wide, ll = 0.04, main = NA,
         col="white",method="jitter",
         ylab=expression(paste("%",Delta)),
         show.names=FALSE)
mtext(expression(bold(paste("(a) %",Delta," egg-laying period"))),line=-1.7,outer=TRUE,side=3,adj=0)
text(.95,45,"a",cex=.8)
text(1.95,45,"b",cex=.8)
text(2.95,45,"b",cex=.8)

beanplot(Length ~ Season, data = df1, ll = 0.04,log="",
                    ylab = "days", side = "both",
                    col = list("black", c("grey", "white")))

legend("bottomright", fill = c("black", "grey"),
           legend = c("El Nino", "La Nina"),bty="n",cex=.8)

text(.95,335,"a",cex=.8)
text(1.95,335,"a",cex=.8)
text(2.95,335,"b",cex=.8)

mtext(expression(bold(paste("(b) Length of the egg-laying period"))),
      line=-16,outer=TRUE,side=3,adj=0)



dev.off()


##########################################
##Make histograms of each species and label if it is early mid or late breeder ########
##########################################

#combine data

state<-c("elNino","LaNina","Neutral")

early<-subset(wide,season=="Early")
mid<-subset(wide,season=="Mid")
late<-subset(wide,season=="Late")

wide2<-rbind(early,rbind(mid,late))

plot_list = list()
for (i in 1:length(wide2$Species)){
  Pdat<-subset(PhaseDat,Scientific.Name == wide2$Species[i],
               select=c("modifiedDOY","SOI","English.name"))
  p<-ggplot(Pdat, aes(x = modifiedDOY, fill = SOI)) + 
    geom_density(alpha = 0.5, bw=30) +
    ggtitle(paste(unique(Pdat$English.name)),wide2$season[i])
  plot_list[[i]] = p
}

pdf("/Users/daisy/Google Drive/PhD/ENSO/Manuscript/GlobalChangeBiology/Figures/desert_KDE_Breedingplots.pdf")
for (i in 1:length(wide2$Species)){
    print(plot_list[[i]])
    
}
dev.off()


######################################################

#Are species breeding at signficantly differnt latitudes during el Nino and la Nina years

####
mEle<-lmer(logelev ~ BreedingENSO + lat+ logdist+modifiedDOY+ (1|Scientific.Name)+(1|Order),REML=F, data = ModifiedDate )
summary(mEle)
Anova(mEle,test.statistic="Chisq")
anova(mEle)

mlat<-lmer(lat ~ BreedingENSO + logelev+ logdist+modifiedDOY+ (1|Scientific.Name)+(1|Order),REML=F, data = ModifiedDate )
summary(mlat)
Anova(mlat,test.statistic="Chisq")
anova(mlat)
Tukeylat<- glht(mlat, linfct=mcp(BreedingENSO="Tukey"))
summary(Tukeylat)

summaryBy(lat ~ BreedingENSO, data=ModifiedDate,
          FUN=c(mean,sd)) #find mean and sd





mdist<-lmer(logdist ~ BreedingENSO + logelev+ lat+modifiedDOY+ (1|Scientific.Name)+(1|Order),REML=F, data = ModifiedDate )
summary(mdist)
Anova(mdist,test.statistic="Chisq")
anova(mdist)






















# 
# 
# 
# 
# #get data ready for supplementary information, make sure names are up to date. 
# inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
# #keep only species of interest
# inc<-subset(inc, remove !=1)
# taxon<-inc[c("genus","family","ScientificName","Common.me")]
# s1<-merge(supplTAB,taxon,all.x=TRUE,by.x="Species",by.y="ScientificName")
# 
# #check names are okay with Australian Bird Data Version 1
# birdData<-read.csv("/Users/daisy/Google Drive/PhD/birdTraits/PublishedData/Australian_Bird_Data_Version_1.csv")
# #keep only full species
# birdData2<-subset(birdData, is.na(X6_Subspecies_name_2))[c(1:5)]
# 
# birdData2$Species<-with(birdData2,paste(X4_Genus_name_2,X5_Species_name_2,sep=" "))
# s2<-merge(s1,birdData2,all.x=TRUE,by="Species")
# bad<-subset(s2,is.na(X4_Genus_name_2))
# good<-subset(s2,!is.na(X4_Genus_name_2))
# 
# good<-good[,c("Species",
#             "Common.me",
#             "NoObsEL",
#             "NoObsLa",
#             "BreedingPeriod_elNn",
#             "Quantile5_elNn",
#             "BreedingPeriod_LaNn",
#             "Quantile5_LaNn",
#             "p")]
# 
# i <- sapply(good, is.factor)
# good[i] <- lapply(good[i], as.character)
# 
# 
# bad<-merge(bad,birdData2,all.x=TRUE,by.x="Common.me",by.y="X3_Taxon_common_name_2")
# bad$Species<-bad$Species.y
# bad<-bad[,c("Species",
#             "Common.me",
#             "NoObsEL",
#             "NoObsLa",
#             "BreedingPeriod_elNn",
#             "Quantile5_elNn",
#             "BreedingPeriod_LaNn",
#             "Quantile5_LaNn",
#             "p")]
# #get rid of factors
# i <- sapply(bad, is.factor)
# bad[i] <- lapply(bad[i], as.character)
# 
# all<-rbind(good,bad)                                                                                                                        
# 
# write.csv(all,"/Users/daisy/Google Drive/PhD/BreedingTiming/manuscript/TablesFigures/SupplementaryData2.csv")

# 
# 
# 
# 
# 
# #write.csv(supplTAB,"/Users/daisy/Google Drive/PhD/BreedingTiming/manuscript/TablesFigures/BreedingPeriod20160105.csv")
# 
# 
# 
# 
# 
# 
# ##########################################
# ############ step 11. Figure - Breeding Period by Phase   ########
# ##########################################
# 
# # library(RColorBrewer)
# # color.gen<-c(brewer.pal(n=11, name="RdYlBu"),
# #              brewer.pal(n=11, name="RdBu"),
# #              brewer.pal(n=12, name="Paired"),
# #              brewer.pal(n=11, name="PuOr"),
# #              brewer.pal(n=11, name="BrBG"))
# 
# 
# boxSlopes <- function(condition,group,labelY,sig,sigsize,figNumber,Sp,dat,axesTF,at,labels){
#   wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
#               ylab=labelY, 
#               main="",names=c("",""),outcol="white", 
#               boxwex=0.5,border = c("black","black"),
#               cex.lab=1.3, axes=axesTF)
#   box()
#   axis(1, at=1:2,labels=FALSE)
#   mtext(expression(El~Ni*tilde(n)*o), side=1, line=1, at=1, cex=.8)
#   mtext(expression(La~Ni*tilde(n)*a), side=1, line=1, at=2, cex=.8)
#   if(axesTF==FALSE){
#     axis(2, at=at, 
#          labels=labels)
#   }
#   mtext(figNumber, side=3, line=1, at=1.5,  adj=8,cex=.8,font=2)
#   #add line for each species
#   for (i in 1:length(Sp)) {
#     el<-subset(dat,Species==Sp[i] & Phase =="elNn",select = paste(condition))[,1]
#     la<-subset(dat,Species  == Sp[i] & Phase =="LaNn",select = condition)[,1]
#     color="grey40"
#     #color=color.gen[i] #for coloured lines
#     #segments(1,slopes.gen$intercept[i]+slopes.gen$slope[i],2,slopes.gen$intercept[i], col=color, lwd=0.7, lty=1)
#     segments(1,el,2,la, col=color, lwd=0.7, lty=1)
#   }
#   #add in the box plots again so that they are not covered by the lines
#   wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
#   ylab=labelY, 
#   main="",names=c("",""), 
#   boxwex=0.5,border = c("black","black"),
#   cex.lab=1.3, axes=axesTF, add = TRUE, col ="#FFFFFF00" )
#   mtext(sig, side=3, line=-2.5, at=1.5,  adj=0.5,cex=sigsize)
# }
# 
# #prep data
# allDat$Species<-lapply(allDat$Species, as.character)
# SpeciesTest<-unlist(allDat$Species)
# allDat$Genus<-sapply(strsplit(SpeciesTest, " "), "[", 1)
# LEallDat<-droplevels(subset(allDat,Phase=="elNn" | Phase =="LaNn"))
# 
# #reorder data so elNina is first
# LEallDat$Phase <- with(LEallDat, reorder(Phase,
#                                      BreedingPeriod, mean))
# #set up PDF
# pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/BreedingPeriodByPhase_Significant20160424.pdf",
#     width = 8, height = 4)
# 
# par(mfrow=c(1,3))
# 
# 
# #Three figures 1 row 3 columns, breeding period, start, conclusion
# boxSlopes(condition="BreedingPeriod",
#           group="Phase",
#           labelY="ELP (day)",
#           sig="*",
#           sigsize=1.5,
#           figNumber="(a)",
#           Sp=unlist(unique(LEallDat$Species)),
#           dat=LEallDat,
#           axesTF=TRUE)
#          
# boxSlopes(condition="Quantile5",
#           group="Phase",
#           labelY="Start of ELP",
#           sig="n.s.",
#           sigsize=.9,
#           figNumber="(b)",
#           Sp=unlist(unique(LEallDat$Species)),
#           dat=LEallDat,
#           axesTF=FALSE,
#           at=c(1,32,61,92,121,152,182,214,245,275,335), 
#           labels=c("Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))
# 
# 
# #modify date
# LEallDat$MODIFIED_Quantile95<-with(LEallDat,ifelse(Quantile95<Quantile5,Quantile95+365,Quantile95))
# boxSlopes(condition="MODIFIED_Quantile95",
#           group="Phase",
#           labelY="Termination of ELP",
#           sig="*",
#           sigsize=1.5,
#           figNumber="(c)",
#           Sp=unlist(unique(LEallDat$Species)),
#           dat=LEallDat,
#           axesTF=FALSE,
#           at=c(245,275,305,335,366,397,425,456,487,518,548), 
#           labels=c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","July"))
# 
# # LEallDat$MODIFIED_Quantile50<-with(LEallDat,ifelse(Quantile50<Quantile5,Quantile50+365,Quantile50))
# # 
# # boxSlopes(condition="MODIFIED_Quantile50",
# #           group="Phase",
# #           labelY="End dates",
# #           sig="*",
# #           sigsize=1.5,
# #           figNumber="(c)",
# #           Sp=unlist(unique(LEallDat$Species)),
# #           dat=LEallDat,
# #           axesTF=FALSE,
# #           at=c(245,275,305,335,366,397,425,456,487,518,548), 
# #           labels=c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","July"))
# 
# 
# 
# dev.off()
# 
# 
# 
# 
# ##########################################
# ############ step 12. Figure for Species on the move- Breeding Period by Phase   ########
# ##########################################
# 
# boxSlopes <- function(condition,group,labelY,sig,sigsize,figNumber,Sp,dat,axesTF,at,labels){
#   wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
#               ylab=labelY, 
#               main="",names=c("",""),outcol="white", 
#               boxwex=0.5,border = c("black","black"),
#               cex.lab=1, axes=axesTF)
#   box()
#   axis(1, at=1:2,labels=FALSE)
#   mtext(expression(El~Ni*tilde(n)*o), side=1, line=1, at=1, cex=1)
#   mtext(expression(La~Ni*tilde(n)*a), side=1, line=1, at=2, cex=1)
#   if(axesTF==FALSE){
#     axis(2, at=at, 
#          labels=labels)
#   }
#   mtext(figNumber, side=3, line=1, at=1.5,cex=1)
#   #add line for each species
#   
#   el<-subset(dat,Phase =="elNn",select = paste(condition))[,1]
#   la<-subset(dat,Phase =="LaNn",select = condition)[,1]
#   absDiff<-abs(el-la)
#   
#   dat2<-as.data.frame(cbind(el,la, absDiff))
#  dat2<-dat2[order(absDiff),] 
#   
#   
#   for (i in 1:nrow(dat2)) {
#     el<-dat2[i,1]#subset(dat,Species==Sp[i] & Phase =="elNn",select = paste(condition))[,1]
#     la<-dat2[i,2]#subset(dat,Species  == Sp[i] & Phase =="LaNn",select = condition)[,1]
#     
#     color=ifelse(la-el>=20 | el-la>=20, "dodgerblue3", "grey60")
#     #color=color.gen[i] #for coloured lines
#     #segments(1,slopes.gen$intercept[i]+slopes.gen$slope[i],2,slopes.gen$intercept[i], col=color, lwd=0.7, lty=1)
#     segments(1,el,2,la, col=color, lwd=0.8, lty=1)
#   }
#   #add in the box plots again so that they are not covered by the lines
#   wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
#               ylab=labelY, 
#               main="",names=c("",""), 
#               boxwex=0.5,border = c("black","black"),
#               cex.lab=1, axes=axesTF, add = TRUE, col ="#FFFFFF00" )
#   mtext(sig, side=3, line=-2.5, at=1.5,  adj=0.5,cex=sigsize)
# }
# 
# #prep data
# allDat$Species<-lapply(allDat$Species, as.character)
# SpeciesTest<-unlist(allDat$Species)
# allDat$Genus<-sapply(strsplit(SpeciesTest, " "), "[", 1)
# LEallDat<-droplevels(subset(allDat,Phase=="elNn" | Phase =="LaNn"))
# 
# #reorder data so elNina is first
# LEallDat$Phase <- with(LEallDat, reorder(Phase,
#                                          BreedingPeriod, mean))
# #set up PDF
# pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/SpeciesOnTheMoveBreedingPeriodByPhase_Significant.pdf",
#    width = 8, height = 4)
# 
# par(mfrow=c(1,2))
# 
# 
# #Three figures 1 row 2 columns, breeding period, start
# boxSlopes(condition="BreedingPeriod",
#           group="Phase",
#           labelY="Length (days)",
#           sig="*",
#           sigsize=1.5,
#           figNumber="Egg-laying Period",
#           Sp=unlist(unique(LEallDat$Species)),
#           dat=LEallDat,
#           axesTF=TRUE)
# 
# boxSlopes(condition="Quantile5",
#           group="Phase",
#           labelY="Date",
#           sig="n.s.",
#           sigsize=.9,
#           figNumber="Start date",
#           Sp=unlist(unique(LEallDat$Species)),
#           dat=LEallDat,
#           axesTF=FALSE,
#           at=c(1,32,61,92,121,152,182,214,245,275,335), 
#           labels=c("Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))
# 
# 
# 
# 
# 
# dev.off()
# 
# 
# 
# 
# 





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





