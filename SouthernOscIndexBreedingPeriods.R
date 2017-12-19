

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
library(beanplot)
library(reporttools)
library(lmerTest)
library(spdep)
library(lme4)
library(MuMIn)
library(nlme)
#library(afex)
library(ENmisc)
#library(arm)
library(car)
#library(heavy)
calcskew <- function(circulardates) {
  datesC <-
    circular(
      circulardates,
      modulo = "2pi",
      units = "radians",
      rotation = "counter"
    )
  Rbar <-
    rho.circular(datesC)#average clustering 0 is uncluseted 1 is all at same location
  V <- 1 - Rbar#sample circular variance
  t2t <- trigonometric.moment(datesC, p = 2, center = TRUE)
  bbar2 <- t2t$sin
  skewness <- bbar2 / (V ** (3 / 2)) #skewness
  return(round(skewness, 2))
}



##########################################
############ step 3. set up datafiles ########
##########################################

obs.dir <-
  '/Users/daisy/GoogleDrive/PhD/Data/Observaitons/Cleaned/Breeding/'
# Breeding quantiles
dat <-
  fread(paste0(obs.dir, 'PointOfLayDayOfYear2016-09-20.csv'))[, c("Scientific.Name",
                                                                  "lat",
                                                                  "lon",
                                                                  "sourceName",
                                                                  "type",
                                                                  "DOY_PL",
                                                                  "year",
                                                                  "month")]



dat$DOY_PL[dat$DOY_PL == 366] <- 365
dat$month <- month(as.Date(as.character(dat$DOY_PL), format = "%j"))



# traits
inc <-
  read.csv(
    '/Users/daisy/GoogleDrive/PhD/BreedingTiming/tables/SpeciesOfInterest_2016-10-07.csv'
  )
taxon <-
  inc[c(
    "genus",
    "Clem.Family",
    "ScientificName",
    "Clements.Scientific.name",
    "sort.v2016",
    "family",
    "English.name",
    "Order"
  )]
dat <-
  merge(dat,
        taxon,
        all.x = TRUE,
        by.x = "Scientific.Name",
        by.y = "ScientificName")
#add koeppen zones
koeppen <-
  raster(
    paste0(
      '/Users/daisy/GoogleDrive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'
    ),
    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")
  )
#combine equatorial and tropical: 41 Equatorial, 35 Tropical, extract koeppen zone
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol = 3, byrow = TRUE)
koeppen <- reclassify(koeppen, rclmat)
dat$koeppen <- raster::extract(koeppen, cbind(dat$lon, dat$lat))


#convert data to radians to make circular vector
dat$Radians <- (dat$DOY_PL / 365 * 360) * pi / 180

##########################
#Chose the region of interest
#Analysis is identical for both regions
#Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,13-Grassland,3-Temperate

##########################

#TEMP<-subset(dat, koeppen == 22 |koeppen == 13) #arid region Breeding Period
TEMP <- subset(dat, koeppen == 3)#temperate region
TEMP$month <-
  formatC(TEMP$month, width = 2, flag = '0') #format months 1 becomes 01
TEMP$date <-
  paste0(TEMP$year, TEMP$month) #make date string that matches NOAA SOI
TEMP <- subset(TEMP, date >= 190001) #limit time to same as NOAA data
#read in BOM ENSO events
#BOM
SOI <-
  read.csv("/Users/daisy/GoogleDrive/PhD/ENSO/Data/BOM_SOItest.csv")
SOI$month2 <-
  formatC(SOI$monthNumeric, width = 2, flag = '0') #format months 1 becomes 01
SOI$date <-
  paste0(SOI$Year, SOI$month2) #make date string that matches NOAA SOI
SOI <-
  SOI[, c("date", "BreedingENSO", "Year", "monthNumeric", "month2", "X")]
#how many breeding seasons in each phase
July <- subset(SOI, month2 == "07")
table(July$BreedingENSO)

TEMP <- merge(TEMP,
              SOI,
              by.x = "date",
              by.y = "date",
              all.x = TRUE)
##########################################
############ step 4. get species with more than 100 obs in each phase ########
##########################################

#get el Nino obs and species with more than 100 observations
eN <- subset(TEMP, BreedingENSO == 'NINO')
spSumeN <- as.data.frame(table(eN$Scientific.Name))
sp_eN <- as.vector(droplevels(subset(spSumeN, Freq >= 100)$Var1))
lN <- subset(TEMP, BreedingENSO == 'NINA')
spSumlN <- as.data.frame(table(lN$Scientific.Name))
sp_lN <- as.vector(droplevels(subset(spSumlN, Freq >= 100)$Var1))
nN <- subset(TEMP, BreedingENSO == 'neutral')
spSumnN <- as.data.frame(table(nN$Scientific.Name))
sp_nN <- as.vector(droplevels(subset(spSumnN, Freq >= 100)$Var1))
#name of species in both, limit data to these species
fin_sp <- intersect(sp_eN, sp_lN)
fin_sp <- intersect(fin_sp, sp_nN)#neutral
eN <- droplevels(eN[eN$Scientific.Name %in% fin_sp, ])
eN$SOI <- as.character(eN$SOI)
eN$SOI <- "elNino"
lN <- lN[lN$Scientific.Name %in% fin_sp, ]
lN$SOI <- as.character(lN$SOI)
lN$SOI <- "LaNina"
nN <- nN[nN$Scientific.Name %in% fin_sp, ]
nN$SOI <- as.character(nN$SOI)
nN$SOI <- "Neutral"

########################################
###########Step 4b. remove 'unknown' breeding observations when there are more than 100 from the other catagories
#######################################


#get NINO species with more than 100 Egg or young or multi observations
datAcc <-
  subset(eN, type == "solitary-egg" |
           type == "multi" | type == "young")
datAcc_sp <- as.data.frame(table(datAcc$Scientific.Name))
dat_eN_Acc_sp <-
  as.vector(droplevels(subset(datAcc_sp, Freq >= 100)$Var1))
#get NINA species with more than 100 Egg or young or multi observations
dat_LA_Acc <-
  subset(lN, type == "solitary-egg" |
           type == "multi" | type == "young")
dat_LA_Acc_sp <- as.data.frame(table(dat_LA_Acc$Scientific.Name))
dat_LA_Acc_sp <-
  as.vector(droplevels(subset(dat_LA_Acc_sp, Freq >= 100)$Var1))
#get Neutral species with more than 100 Egg or young or multi observations
dat_NU_Acc <-
  subset(nN, type == "solitary-egg" |
           type == "multi" | type == "young")
dat_NU_Acc_sp <- as.data.frame(table(dat_NU_Acc$Scientific.Name))
dat_NU_Acc_sp <-
  as.vector(droplevels(subset(dat_NU_Acc_sp, Freq >= 100)$Var1))


#species in both la nina and el nino
fin_sp_Acc <-
  intersect(dat_eN_Acc_sp, dat_LA_Acc_sp) #25 species in total
fin_sp_Acc <- intersect(fin_sp_Acc, dat_NU_Acc_sp)
#get the high accuracy data
dat_eN_Acc <- datAcc[datAcc$Scientific.Name %in% fin_sp_Acc, ]
dat_LA_Acc <- dat_LA_Acc[dat_LA_Acc$Scientific.Name %in% fin_sp_Acc, ]
dat_NU_Acc <- dat_NU_Acc[dat_NU_Acc$Scientific.Name %in% fin_sp_Acc, ]
#put the data back together
high <- rbind(dat_LA_Acc, dat_eN_Acc)#
high <- rbind(high, dat_NU_Acc)
high$Acc <- "high"
length(unique(high$Scientific.Name))#number of species with only high accuracy observation
#get the data for the rest of the species
dat_eN_LOW <- eN[eN$Scientific.Name %nin% fin_sp_Acc, ]
length(unique(dat_eN_LOW$Scientific.Name))
dat_LA_LOW <- lN[lN$Scientific.Name %nin% fin_sp_Acc, ]
length(unique(dat_LA_LOW$Scientific.Name))
dat_Nu_LOW <- nN[nN$Scientific.Name %nin% fin_sp_Acc, ]
length(unique(dat_Nu_LOW$Scientific.Name))

low <- rbind(dat_eN_LOW, dat_LA_LOW)#,
low <- rbind(low, dat_Nu_LOW)
low$Acc <- "low"

#put data back together
PhaseDat <- droplevels(rbind(high, low))

##########
elevation <-
  raster(
    "/Users/daisy/GoogleDrive/PhD/Data/Spatial/Elevation/eMAST_ANUClimate_fx_elev_v1m0.nc"
  )
distToCoast <-
  raster(
    "/Users/daisy/GoogleDrive/PhD/Data/Spatial/DistanceToCoast/eMAST_ANUClimate_fx_dist_v1m0.nc"
  )

PhaseDat$elevation <-
  raster::extract(elevation, cbind(PhaseDat$lon, PhaseDat$lat))

PhaseDat$distToCoast <-
  raster::extract(distToCoast, cbind(PhaseDat$lon, PhaseDat$lat))


##########################################
############ step 5. Calc breeding quantiles during differnt phases ########
##########################################

#combine data
SpeciesSummary <- list()
state <- c("elNino", "LaNina", "Neutral")
alldat90 <- list()
for (i in 1:length(fin_sp)) {
  phaseSummary <- list()
  phasealldat90 <- list()
  for (ii in 1:length(state)) {
    #get species data in specific phase
    spdat <-
      subset(PhaseDat, Scientific.Name == fin_sp[i] &
               SOI == state[ii])#get species data
    try(if (nrow(spdat) < 100)
      stop("not enough observations"))
    NoObs <- as.numeric(nrow(spdat))#get observation count
    CommonName <-
      as.character(subset(inc, ScientificName == fin_sp[i], "Common.me", drop =
                            TRUE))
    Order <-
      as.character(subset(inc, ScientificName == fin_sp[i], order, drop = TRUE))
    family <-
      as.character(subset(inc, ScientificName == fin_sp[i], family, drop = TRUE))
    #clades<-as.character(subset(inc,ScientificName==fin_sp[i],"Patch.clade",drop=TRUE))
    obsRadians <- spdat$Radians
    #date of quantiles
    quant <- quantile.circular(obsRadians, c(0, .05, .5, .95), type = 8)
    Per0 <- round((quant[[1]] * 180 / pi) / 360 * 365)
    Per5 <- round((quant[[2]] * 180 / pi) / 360 * 365)
    Per50 <- round((quant[[3]] * 180 / pi) / 360 * 365)
    Per95 <- round((quant[[4]] * 180 / pi) / 360 * 365)
    #Breeding Period Length
    StartDate <- as.numeric(Per5)
    EndDate <- as.numeric(Per95)
    BPL <-
      ifelse (EndDate < StartDate,
              365 - StartDate + EndDate,
              EndDate - StartDate)
    ###90% of data, remove obs outside of peak breeding
    spdat90 <- if (quant[[2]] < quant[[4]]) {
      subset(spdat, Radians >= quant[[2]] & Radians <= quant[[4]])
    } else
      rbind(subset(spdat, Radians >= quant[[2]]),
            subset(spdat, Radians <= quant[[4]]))
    #calculate skewness and clumpiness
    obsRadians90C <-
      circular(
        spdat90$Radians,
        modulo = "2pi",
        units = "radians",
        rotation = "counter"
      )
    Rbar90Per <-
      round(rho.circular(obsRadians90C), 2)#average clustering 0 is uncluseted 1 is all at same location
    skew90 <-
      calcskew(spdat90$Radians)#negative numbers skewed counter clockwise, poitive clockwise, 0 = not skewed
    #put the data together and name it
    Perdat <-
      data.frame(
        paste(state[ii]),
        paste(fin_sp[i]),
        NoObs,
        CommonName,
        Order,
        family,
        length(obsRadians),
        Per0,
        Per5,
        Per50,
        Per95,
        BPL,
        Rbar90Per,
        skew90
      )
    colnames(Perdat) <-
      c(
        'Phase',
        'Species',
        'NoObs',
        'CommonName',
        'Order',
        'family',
        'ObservationCount',
        'Quantile0',
        'Quantile5',
        'Quantile50',
        'Quantile95',
        'BreedingPeriod',
        'Rbar90',
        'Skew90'
      )
    Perdat$Accu <- unique(spdat$Acc)
    phaseSummary[[ii]] <- Perdat
    phasealldat90[[ii]] <- spdat90
  }
  SpeciesSummary[[i]] <- do.call('rbind', phaseSummary)
  alldat90[[i]] <- do.call('rbind', phasealldat90)
}

alldat90 <- as.data.frame(do.call('rbind', alldat90))
allDat <- as.data.frame(do.call('rbind', SpeciesSummary))

try(if (nrow(alldat90) < nrow(PhaseDat) - nrow(PhaseDat) / 10)
  stop("more than 10% removed"))


############################
####Modify date for lm
############################

ModifiedDate <- list()
for (m in 1:length(fin_sp)) {
  spdat <- subset(alldat90, Scientific.Name == fin_sp[m])
  earlstart <- min(subset(allDat, Species == fin_sp[m])$Quantile5)
  latestend <- max(subset(allDat, Species == fin_sp[m])$Quantile95)
  #adjust dates so those at the begining of the year come after 365

  df <- subset(spdat, DOY_PL < earlstart)#dates to adjust
  df$modifiedDOY <- df$DOY_PL + 365
  df2 <- subset(spdat, DOY_PL > earlstart)#dates that are fine
  df2$modifiedDOY <- df2$DOY_PL
  ModifiedDate[[m]] <- rbind(df, df2)
}

ModifiedDate <- do.call("rbind", ModifiedDate)


ModifiedDate2 <- list()
for (m in 1:length(fin_sp)) {
  spdat <- subset(PhaseDat, Scientific.Name == fin_sp[m])
  earlstart <- min(subset(allDat, Species == fin_sp[m])$Quantile5)
  latestend <- max(subset(allDat, Species == fin_sp[m])$Quantile95)
  #adjust dates so those at the begining of the year come after 365

  df <- subset(spdat, DOY_PL < earlstart)#dates to adjust
  df$modifiedDOY <- df$DOY_PL + 365
  df2 <- subset(spdat, DOY_PL > earlstart)#dates that are fine
  df2$modifiedDOY <- df2$DOY_PL
  ModifiedDate2[[m]] <- rbind(df, df2)
}

PhaseDat <- do.call("rbind", ModifiedDate2)

##########################################
############ Are there differences in mean during El Nino and La Nina?
##########################################


ModifiedDate$Phase <-
  factor(as.character(ModifiedDate$BreedingENSO),
         levels = c("NINO", "NINA", "neutral"))

ModifiedDate$elevation[ModifiedDate$elevation <= 0] <- 1
ModifiedDate$elevation[is.na(ModifiedDate$elevation)] <- 0.01
ModifiedDate$elevation <- ModifiedDate$elevation + 1

ModifiedDate$distToCoast[ModifiedDate$distToCoast <= 0] <- 1
ModifiedDate$distToCoast[is.na(ModifiedDate$distToCoast)] <- 1
ModifiedDate$distToCoast <- ModifiedDate$distToCoast + 1

ModifiedDate$logelev <- log10(ModifiedDate$elevation)
ModifiedDate$logdist <- log10(ModifiedDate$distToCoast)



# ModifiedDate$lat<-round(ModifiedDate$lat,2)
# ModifiedDate$lon<-round(ModifiedDate$lon,2)
# ModifiedDate<-unique(ModifiedDate)


#lme linear mixed-effects regression command in the nlme R package allows the user to fit a regression model in which the outcome and the expected errors are spatially autocorrelated.

ensolm <- lme(fixed = modifiedDOY ~ BreedingENSO + logelev * logdist + lat,
              data = ModifiedDate,
              random = ~ 1 |Order/Scientific.Name)

library(visreg)

visreg(ensolm, "BreedingENSO")

ensolm2 <- lme(fixed = modifiedDOY ~ BreedingENSO + logelev * logdist + lat,
              data = ModifiedDate,
              random = ~ 1 |Scientific.Name)




summary(ensolm)
r.squaredGLMM(ensolm)

soil.gau <- update(ensolm, correlation = corGaus(1, form = ~lat+lon))
summary(soil.gau)
r.squaredGLMM(ensolm)
            

mF <-
  lmer(
    modifiedDOY ~ BreedingENSO + logelev * logdist + 
      lat + (1|Year) + (1 + BreedingENSO| Scientific.Name),
    REML = T,
    data = ModifiedDate
  )




summary(mF)
qqPlot(residuals(mF), main = ("DOY"))
Anova(mF, test.statistic = "Chisq")#get p value
r.squaredGLMM(mF)
anova(mF, test = "marginal", type = 2)#get F statement
TukeyENSO <- glht(mF, linfct = mcp(BreedingENSO = "Tukey"))
summary(TukeyENSO)

library(sjPlot)
library(sjmisc)
sjt.lmer(mF)

#get mean breeding days for differnt phase
lsmeansLT(mF)


visreg(mF,
       "modifiedDOY",
       by = "BreedingENSO",
       overlay = TRUE,
       alpha = .09)


######################################
#species that meet La Nina Hypothesis
######################################

#change direction of data
wide <-
  reshape(
    allDat,
    idvar = c("Species", "CommonName", "Order", "family") ,
    timevar = "Phase",
    direction = "wide"
  )
#modify end dates, account for those species that breed over turn of year
wide$modified95 <-
  ifelse(wide$Quantile95.elNino < 200,
         wide$Quantile95.elNino + 365,
         wide$Quantile95.elNino)
wide$modified95LA <-
  ifelse(wide$Quantile95.LaNina < 150,
         wide$Quantile95.LaNina + 365,
         wide$Quantile95.LaNina)
wide$modified95NEU <-
  ifelse(
    wide$Quantile95.Neutral < 175,
    wide$Quantile95.Neutral + 365,
    wide$Quantile95.Neutral
  )


#We hypothesised
#earlier start during La Nina
#length of egg-laying period should be longer during la nina
#peak and conclusion later in the year
wide$earlierstart <- with(wide, Quantile5.LaNina - Quantile5.elNino)
wide$longerBPEL <-
  with(wide, BreedingPeriod.LaNina - BreedingPeriod.elNino)
wide$laterpeak <- with(wide, Quantile50.LaNina - Quantile50.elNino)
wide$laterconclusion <- with(wide, modified95LA - modified95)

wide$earlierstartNEU <- with(wide, Quantile5.LaNina - Quantile5.Neutral)
wide$longerBPNEU <-
  with(wide, BreedingPeriod.LaNina - BreedingPeriod.Neutral)
wide$laterpeakNEU <- with(wide, Quantile50.LaNina - Quantile50.Neutral)
wide$laterconclusionNEU <- with(wide, modified95LA - modified95NEU)

#write.csv(wide, '/Users/daisy/GoogleDrive/PhD/ENSO/Tables/TemperateWideSummary.csv', row.names=F)

#% of species with earlier starts during La Nina
nrow(subset(wide, earlierstart < 0)) / nrow(wide) * 100
nrow(subset(wide, earlierstartNEU < 0)) / nrow(wide) * 100

#% species with longer breeding periods
nrow(subset(wide, longerBPEL > 0)) / nrow(wide) * 100
nrow(subset(wide, longerBPNEU > 0)) / nrow(wide) * 100

#xx% of species having later peaks and
nrow(subset(wide, wide$laterpeak > 0)) / nrow(wide) * 100
nrow(subset(wide, wide$laterpeakNEU < 0)) / nrow(wide) * 100

#% of species having later terminations to egg-laying period.
nrow(subset(wide, wide$laterconclusion > 0)) / nrow(wide) * 100
nrow(subset(wide, wide$laterconclusionNEU < 0)) / nrow(wide) * 100


wide$ChangeBPEL <-
  with(wide, (((BreedingPeriod.elNino - BreedingPeriod.LaNina) / BreedingPeriod.elNino
  ) * 100))
wide$ChangeBPELNeu <-
  with(wide, (((BreedingPeriod.elNino - BreedingPeriod.Neutral) / BreedingPeriod.elNino
  ) * 100))
wide$ChangeBPLANeu <-
  with(wide, (((BreedingPeriod.LaNina - BreedingPeriod.Neutral) / BreedingPeriod.LaNina
  ) * 100))
wide$ChangeBPLAEL <-
  with(wide, (((BreedingPeriod.LaNina - BreedingPeriod.elNino) / BreedingPeriod.LaNina
  ) * 100))


BL <-
  subset(
    wide,
    select = c(
      "Species",
      "CommonName",
      "BreedingPeriod.LaNina",
      "BreedingPeriod.Neutral",
      "BreedingPeriod.elNino"
    )
  )
colMeans(BL[3:5])






#####################
#Can we explain those species that have biggest change in breeding in relation to the time of year
#those species that during neutral conditions breed before the end of the year
#vs those that breed after the turn of the year
#################


#reformat data
a <- subset(wide, select = c("ChangeBPELNeu", "modified95NEU", "Order"))
colnames(a) <- c("change", "end", "Order")
a$time <- "EtoNeu"
b <- subset(wide, select = c("ChangeBPEL", "modified95NEU", "Order"))
colnames(b) <- c("change", "end", "Order")
b$time <- "EtoL"
c <- rbind(a, b)
d <- subset(wide, select = c("ChangeBPLANeu", "modified95NEU", "Order"))
colnames(d) <- c("change", "end", "Order")
d$time <- "LtoN"
c <- rbind(c, d)
c$time <- as.factor(c$time)
c$EL<-as.factor(ifelse(c$end<=365,"A","B"))
c$phaseTime<-paste0(c$time,c$EL)

changeDat<-subset(c, time == "EtoNeu" | time == "LtoN")

#model to assess if species that breed primarily in the spring months have a greater response to La Niña and El Niño events, when compared to Neutral conditions
summary(fit <- lmer(change ~ phaseTime + (1+phaseTime|Order), data = changeDat))


lsmeansLT(fit)
qqPlot(residuals(fit), main = ("DOY"))
Anova(fit, test.statistic = "Chisq")#get p value
anova(fit, test = "marginal", type = 2)#get F statement
r.squaredGLMM(fit)
lsmeansLT(fit)

###############################
#alternative analysis
#############################
#change data around to use with visreag
a <- subset(wide, select = c("ChangeBPELNeu", "modified95NEU", "Order"))
colnames(a) <- c("change", "end", "Order")
a$time <- "EtoNeu"
# b <- subset(wide, select = c("ChangeBPEL", "modified95NEU", "Order"))
# colnames(b) <- c("change", "end", "Order")
# b$time <- "EtoL"
# c <- rbind(a, b)
d <- subset(wide, select = c("ChangeBPLANeu", "modified95NEU", "Order"))
colnames(d) <- c("change", "end", "Order")
d$time <- "LtoN"
c <- rbind(a, d)
c$time <- as.factor(c$time)

# cArid<-c
# fitArid <- lmer(change ~ time * end +(1+time|Order), data = cArid)

fit <- lmer(change ~ time * end +(1+time|Order), data = c)
qqPlot(residuals(fit))
par(mfrow = c(1, 2))
plot(fit, which = c(3, 2))
par(mfrow = c(1, 1))
visreg(fit,
       "end",
       by = "time",
       overlay = TRUE,
       alpha = .095)

summary(fit)
anova(fit, test = "marginal", type = 2)#get F statement
r.squaredGLMM(fit)
Anova(fit, test.statistic = "Chisq")#get p value
r.squaredGLMM(fit)
lsmeansLT(fit)

#fitTemp<-fit

pdf(
  file = paste0(
    "/Users/daisy/GoogleDrive/PhD/ENSO/Figures/TemperateChangeELP",
    as.Date(Sys.time()),
    ".pdf"
  ),
  width = 3.3,
  height = 5
)

par(mfrow = c(2, 1),
    mar = c(5, 4, 1, 1))


visreg(
  fit,
  "end",
  by = "time",
  overlay = TRUE,
  alpha = .095,
  line = list(col = c("grey50","goldenrod1")),
  points = list(
    cex = .7,
    pch = 20,
    col = c("grey50","goldenrod1")
  ),
  fill = list(col = c( "#7F7F7F7F", "#FFC1257F")),
  ylab = "% change",
  xlab = "month",
  legend = FALSE,
  xaxt = 'n'
)
axis(
  side = 1,
  at = c(258, 289, 320, 350, 381, 413, 440, 470),
  c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
)


legend(
  360,
  -15,
  c("La Nina", "El Nino"),
  lwd = 4,
  col = rev(c("grey50","goldenrod1")),
  bty = 'n'
)


visreg(
  fitArid,
  "end",
  by = "time",
  overlay = TRUE,
  alpha = .095,
  line = list(col = c("grey50","goldenrod1")),
  points = list(
    cex = .7,
    pch = 20,
    col = c("grey50","goldenrod1")
  ),
  fill = list(col = c( "#7F7F7F7F", "#FFC1257F")),
  ylab = "% change",
  xlab = "month",
  legend = FALSE,
  xaxt = 'n'
)
axis(
  side = 1,
  at = c(258, 289, 320, 350, 381, 413, 440, 470),
  c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
)

dev.off()

# #change in breedin length compared to La nina, breeding that ends by dec 1
# dec <- subset(wide, modified95NEU < 365)
# mean(dec$ChangeBPLANeu) #La nina to neutral
# sd(dec$ChangeBPLANeu)
# mean(dec$ChangeBPELNeu) #La nina to El Nino
# sd(dec$ChangeBPELNeu)
# 
# 




