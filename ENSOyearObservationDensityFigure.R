#make histogram showing observation density by phase and years that were differnt phase

library(ggplot2)
library(tidyr)

# get data together using other script, make sure data is for both arid region and temperate region species

#PhaseDatAll - name of data


obYear<- ggplot(PhaseDatAll, aes(Year, colour = SOI)) +
  geom_density(show.legend = F, adjust = 2) +
  theme_minimal() +
  scale_color_manual(values = c(LaNina = "#FFC1257F", elNino = "#7F7F7F7F", Neutral = "#0000807F")) +
  scale_size_manual(values = c(1,1, 1))+
  
  geom_segment(aes(x = Year, y = 0, xend = Year, yend = 0.0025))

#add in lines for years of diffent phases

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
SOI<-subset(July, BreedingENSO != "")


pdf(file = paste0("/Users/daisy/GoogleDrive/PhD/ENSO/Manuscript/DiversityAndDist/Figures/ENSOPhase",as.Date(Sys.time()),".pdf"),
    width = 3.3, height = 3)

 ggplot(PhaseDatAll, aes(Year, colour = SOI)) +
  geom_density(show.legend = T, adjust = 2,size=1) +
   theme_bw() +
   theme(axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         legend.position = c(0.2,0.57)) +
  scale_color_manual(values = c(LaNina = "#FFC125", elNino = "#7F7F7F", Neutral = "#000080")) +
  
  
  geom_segment(aes(x = Year, y = -0.002, xend = Year, yend = 0),show.legend = T)



dev.off()

