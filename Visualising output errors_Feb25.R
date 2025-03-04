

rm(list=ls())
library(runjags)
library(fitdistrplus)
library(epiR) # for betabuster
library(tidyverse)
library(HDInterval)
library(dplyr)
library(modeest)
library(ggplot2)
library(ggmcmc)
library(coda)
library(grid)

df <- read.csv() #Read in the model output/ error data file.


df <- df %>% mutate(sensdiff = Actual.sensitivity.one-Actual.sensitivity.two) #Difference between the sensitivity values of the two symptoms
sediff <- abs(df$sensdiff)

df <- df %>% mutate(specdiff = Actual.specificity.one-Actual.specificity.two) #Difference between the specificity values of the two symptoms
spdiff <- abs(df$specdiff)

df <- df %>% mutate(prevdiff = High.Prevalence-Low.Prevalence) #Difference between the sensitivity values of the two symptoms




#Adding these into the dataframe
dfnew<- cbind(df,sediff)
dfnew <- cbind(dfnew,spdiff)



################################################################################################################################################

#Dividing up the dataframe based on prior knowledge#

pdf <- subset(dfnew, Model.Knowledge == "Poor") #Both poor priors

godf <- subset(dfnew, Model.Knowledge == "Good one") #One very narrow prior

gbdf <- subset(dfnew, Model.Knowledge == "Good both") #Both good priors

###############################################################################################################################################

#Now finding where the simulated sensitivity and specificity values fall on the prior distributions#

#Both poor priors first#

#Sensitivity values#
pse1 <- c()
for(i in 1:150){
  pse1[i] <- pbeta(pdf$Actual.sensitivity.one[i], 1.607, 1.7284)
}
pse2 <- c()
for(i in 1:150){
  pse2[i] <- pbeta(pdf$Actual.sensitivity.two[i], 1.607, 1.7284)
}

#Specificity values#
psp1 <- c()
for(i in 1:150){
  psp1[i] <- pbeta(pdf$Actual.specificity.one[i], 1.394, 1.201)
}
psp2 <- c()
for(i in 1:150){
  psp2[i] <- pbeta(pdf$Actual.specificity.two[i], 1.394, 1.201)
}

#Prevalence values#
pprev1 <- c()
for(i in 1:150){
  pprev1[i] <- pbeta(pdf$Low.Prevalence[i], 3, 8)
}
pprev2 <- c()
for(i in 1:150){
  pprev2[i] <- pbeta(pdf$High.Prevalence[i], 6, 4)
}

#Covariance values#
pcovp <- c()
for(i in 1:150){
  pcovp[i] <- punif(pdf$Covdp[i], pdf$Covdp.min[i], pdf$Covdp.max[i])
}
pcovn <- c()
for(i in 1:150){
  pcovn[i] <- punif(pdf$Covdn[i], pdf$Covdn.min[i], pdf$Covdn.max[i])
}


#Adding these values to the dataframe#
pdf <- cbind(pdf, pse1)
pdf <- cbind(pdf, pse2)
pdf <- cbind(pdf, psp1)
pdf <- cbind(pdf, psp2)
pdf <- cbind(pdf, pprev1)
pdf <- cbind(pdf, pprev2)
pdf <- cbind(pdf, pcovp)
pdf <- cbind(pdf, pcovn)


#Calculating the difference of these values from the median of the priors. i.e. a measure of how far they are at the tail#
for(i in 1:150){
pdf$npse1[i] <- 0.5-pdf$pse1[i]
if(pdf$npse1[i]<0){
  pdf$npse1[i] <- pdf$npse1[i]*(-1)
}
}

for(i in 1:150){
  pdf$npsp1[i] <- 0.5-pdf$psp1[i]
  if(pdf$npsp1[i]<0){
    pdf$npsp1[i] <- pdf$npsp1[i]*(-1)
  }
}


for(i in 1:150){
  pdf$npse2[i] <- 0.5-pdf$pse2[i]
  if(pdf$npse2[i]<0){
    pdf$npse2[i] <- pdf$npse2[i]*(-1)
  }
}

for(i in 1:150){
  pdf$npsp2[i] <- 0.5-pdf$psp2[i]
  if(pdf$npsp2[i]<0){
    pdf$npsp2[i] <- pdf$npsp2[i]*(-1)
  }
}

for(i in 1:150){
  pdf$npprev1[i] <- 0.5-pdf$pprev1[i]
  if(pdf$npprev1[i]<0){
    pdf$npprev1[i] <- pdf$npprev1[i]*(-1)
  }
}

for(i in 1:150){
  pdf$npprev2[i] <- 0.5-pdf$pprev2[i]
  if(pdf$npprev2[i]<0){
    pdf$npprev2[i] <- pdf$npprev2[i]*(-1)
  }
}

for(i in 1:150){
  pdf$nppcovp[i] <- 0.5-pdf$pcovp[i]
  if(pdf$nppcovp[i]<0){
    pdf$nppcovp[i] <- pdf$nppcovp[i]*(-1)
  }
}

for(i in 1:150){
  pdf$nppcovn[i] <- 0.5-pdf$pcovn[i]
  if(pdf$nppcovn[i]<0){
    pdf$nppcovn[i] <- pdf$nppcovn[i]*(-1)
  }
}

#Now calculating the maximum distance from the median in a model (i.e. a measure of how far at the tail)
for(i in 1:150){
pdf$maxp[i] <- max(c(pdf$npse1[i], pdf$npse2[i], pdf$npsp1[i], pdf$npsp2[i], pdf$npprev1[i], pdf$npprev2[i],
                   pdf$nppcovp[i], pdf$nppcovn[i]), na.rm = TRUE)}


#Now calculate the minimum difference present between respective sensitivity or specificity values#
for(i in 1:150){
  pdf$mind[i] <- min(pdf$sediff[i], pdf$spdiff[i])}


#Now calculate the largest error between the sensitivity and specificity in the model to later visualise if this has an impact#
for(i in 1:150){
  pdf$maxerr[i] <- max(pdf$Sensitivity.Error[i], pdf$Specificity.Error[i])}


#The aim is then to examine whether large error is associated with either a sensitivity/specificity value from the
#tail of a distribution, or when sensitivity/ specificity values are close together


####################################################################################################################################################################


#Now finding where the simulated sensitivity and specificity values fall on the prior distributions#

#Now for very one good prior (symptom 2)#

#Sensitivity values#
pse1 <- c()
for(i in 1:150){
  pse1[i] <- pbeta(godf$Actual.sensitivity.one[i], 1.607, 1.7284)
}
pse2 <- c()
for(i in 1:150){
  pse2[i] <- pbeta(godf$Actual.sensitivity.two[i], 10.654, 5.246)
}

#Specificity values#
psp1 <- c()
for(i in 1:150){
  psp1[i] <- pbeta(godf$Actual.specificity.one[i], 1.394, 1.201)
}
psp2 <- c()
for(i in 1:150){
  psp2[i] <- pbeta(godf$Actual.specificity.two[i], 11.884, 1.829)
}

#Prevalence values#
pprev1 <- c()
for(i in 1:150){
  pprev1[i] <- pbeta(godf$Low.Prevalence[i], 3, 8)
}
pprev2 <- c()
for(i in 1:150){
  pprev2[i] <- pbeta(godf$High.Prevalence[i], 6, 4)
}

#Covariance values#
pcovp <- c()
for(i in 1:150){
  pcovp[i] <- punif(godf$Covdp[i], godf$Covdp.min[i], godf$Covdp.max[i])
}
pcovn <- c()
for(i in 1:150){
  pcovn[i] <- punif(godf$Covdn[i], godf$Covdn.min[i], godf$Covdn.max[i])
}

#Adding these values to the dataframe#
godf <- cbind(godf, pse1)
godf <- cbind(godf, pse2)
godf <- cbind(godf, psp1)
godf <- cbind(godf, psp2)
godf <- cbind(godf, pprev1)
godf <- cbind(godf, pprev2)
godf <- cbind(godf, pcovp)
godf <- cbind(godf, pcovn)

#Calculating the difference of these values from the median of the priors. i.e. a measure of how far they are at the tail#
for(i in 1:150){
  godf$npse1[i] <- 0.5-godf$pse1[i]
  if(godf$npse1[i]<0){
    godf$npse1[i] <- godf$npse1[i]*(-1)
  }
}

for(i in 1:150){
  godf$npsp1[i] <- 0.5-godf$psp1[i]
  if(godf$npsp1[i]<0){
    godf$npsp1[i] <- godf$npsp1[i]*(-1)
  }
}


for(i in 1:150){
  godf$npse2[i] <- 0.5-godf$pse2[i]
  if(godf$npse2[i]<0){
    godf$npse2[i] <- godf$npse2[i]*(-1)
  }
}

for(i in 1:150){
  godf$npsp2[i] <- 0.5-godf$psp2[i]
  if(godf$npsp2[i]<0){
    godf$npsp2[i] <- godf$npsp2[i]*(-1)
  }
}

for(i in 1:150){
  godf$npprev1[i] <- 0.5-godf$pprev1[i]
  if(godf$npprev1[i]<0){
    godf$npprev1[i] <- godf$npprev1[i]*(-1)
  }
}

for(i in 1:150){
  godf$npprev2[i] <- 0.5-godf$pprev2[i]
  if(godf$npprev2[i]<0){
    godf$npprev2[i] <- godf$npprev2[i]*(-1)
  }
}

for(i in 1:150){
  godf$nppcovp[i] <- 0.5-godf$pcovp[i]
  if(godf$nppcovp[i]<0){
    godf$nppcovp[i] <- godf$nppcovp[i]*(-1)
  }
}

for(i in 1:150){
  godf$nppcovn[i] <- 0.5-godf$pcovn[i]
  if(godf$nppcovn[i]<0){
    godf$nppcovn[i] <- godf$nppcovn[i]*(-1)
  }
}

#Now calculating the maximum distance from the median in a model (i.e. a measure of how far at the tail)
for(i in 1:150){
  godf$maxp[i] <- max(c(godf$npse1[i], godf$npse2[i], godf$npsp1[i], godf$npsp2[i], godf$npprev1[i], godf$npprev2[i],
                      godf$nppcovp[i], godf$nppcovn[i]), na.rm = TRUE)}


#Now calculate the minimum difference present between respective sensitivity or specificity values#
for(i in 1:150){
  godf$mind[i] <- min(godf$sediff[i], godf$spdiff[i])}


#Ok now calculate the largest error between the sensitivity and specificity in the model#
for(i in 1:150){
  godf$maxerr[i] <- max(godf$Sensitivity.Error[i], godf$Specificity.Error[i])}


#The aim is then to examine whether large error is associated with either a sensitivity/specificity value from the
#tail of a distribution, or when sensitivity/ specificity values are close together



###################################################################################################################################################################


#Now finding where the simulated sensitivity and specificity values fall on the prior distributions#

#Now for very one good prior (symptom 2)#

#Sensitivity values#
pse1 <- c()
for(i in 1:150){
  pse1[i] <- pbeta(gbdf$Actual.sensitivity.one[i], 1.607, 1.7284)
}
pse2 <- c()
for(i in 1:150){
  pse2[i] <- pbeta(gbdf$Actual.sensitivity.two[i], 1.607, 1.7284)
}

#Specificity values#
psp1 <- c()
for(i in 1:150){
  psp1[i] <- pbeta(gbdf$Actual.specificity.one[i], 1.394, 1.201)
}
psp2 <- c()
for(i in 1:150){
  psp2[i] <- pbeta(gbdf$Actual.specificity.two[i], 1.394, 1.201)
}

#Prevalence values#
pprev1 <- c()
for(i in 1:150){
  pprev1[i] <- pbeta(gbdf$Low.Prevalence[i], 3, 8)
}
pprev2 <- c()
for(i in 1:150){
  pprev2[i] <- pbeta(gbdf$High.Prevalence[i], 6, 4)
}

#Covariance values#
pcovp <- c()
for(i in 1:150){
  pcovp[i] <- punif(gbdf$Covdp[i], gbdf$Covdp.min[i], gbdf$Covdp.max[i])
}
pcovn <- c()
for(i in 1:150){
  pcovn[i] <- punif(gbdf$Covdn[i], gbdf$Covdn.min[i], gbdf$Covdn.max[i])
}

#Adding these values to the dataframe#
gbdf <- cbind(gbdf, pse1)
gbdf <- cbind(gbdf, pse2)
gbdf <- cbind(gbdf, psp1)
gbdf <- cbind(gbdf, psp2)
gbdf <- cbind(gbdf, pprev1)
gbdf <- cbind(gbdf, pprev2)
gbdf <- cbind(gbdf, pcovp)
gbdf <- cbind(gbdf, pcovn)


#Calculating the difference of these values from the median of the priors. i.e. a measure of how far they are at the tail#
for(i in 1:150){
  gbdf$npse1[i] <- 0.5-gbdf$pse1[i]
  if(gbdf$npse1[i]<0){
    gbdf$npse1[i] <- gbdf$npse1[i]*(-1)
  }
}

for(i in 1:150){
  gbdf$npsp1[i] <- 0.5-gbdf$psp1[i]
  if(gbdf$npsp1[i]<0){
    gbdf$npsp1[i] <- gbdf$npsp1[i]*(-1)
  }
}


for(i in 1:150){
  gbdf$npse2[i] <- 0.5-gbdf$pse2[i]
  if(gbdf$npse2[i]<0){
    gbdf$npse2[i] <- gbdf$npse2[i]*(-1)
  }
}

for(i in 1:150){
  gbdf$npsp2[i] <- 0.5-gbdf$psp2[i]
  if(gbdf$npsp2[i]<0){
    gbdf$npsp2[i] <- gbdf$npsp2[i]*(-1)
  }
}

for(i in 1:150){
  gbdf$npprev1[i] <- 0.5-gbdf$pprev1[i]
  if(gbdf$npprev1[i]<0){
    gbdf$npprev1[i] <- gbdf$npprev1[i]*(-1)
  }
}

for(i in 1:150){
  gbdf$npprev2[i] <- 0.5-gbdf$pprev2[i]
  if(gbdf$npprev2[i]<0){
    gbdf$npprev2[i] <- gbdf$npprev2[i]*(-1)
  }
}

for(i in 1:150){
  gbdf$nppcovp[i] <- 0.5-gbdf$pcovp[i]
  if(gbdf$nppcovp[i]<0){
    gbdf$nppcovp[i] <- gbdf$nppcovp[i]*(-1)
  }
}

for(i in 1:150){
  gbdf$nppcovn[i] <- 0.5-gbdf$pcovn[i]
  if(gbdf$nppcovn[i]<0){
    gbdf$nppcovn[i] <- gbdf$nppcovn[i]*(-1)
  }
}

#Now calculating the maximum distance from the median in a model (i.e. a measure of how far at the tail)
for(i in 1:150){
  gbdf$maxp[i] <- max(c(gbdf$npse1[i], gbdf$npse2[i], gbdf$npsp1[i], gbdf$npsp2[i], gbdf$npprev1[i], gbdf$npprev2[i],
                      gbdf$nppcovp[i],gbdf$nppcovn[i]), na.rm = TRUE)}


#Now calculate the minimum difference present between respective sensitivity or specificity values#
for(i in 1:150){
  gbdf$mind[i] <- min(gbdf$sediff[i], gbdf$spdiff[i])}


#Ok now calculate the largest error between the sensitivity and specificity in the model#
for(i in 1:150){
  gbdf$maxerr[i] <- max(gbdf$Sensitivity.Error[i], gbdf$Specificity.Error[i])}

#The aim is then to examine whether large error is associated with either a sensitivity/specificity value from the
#tail of a distribution, or when sensitivity/ specificity values are close together

################## Plot ##################

#Combine all#

maindf <- rbind(gbdf,godf)
maindf <- rbind(pdf,maindf)

#maindf$errordiffb <- ifelse(maindf$mind < 0.05, "Less than 0.05", "Greater than 0.05")  #values based on looking at data - there is no apparent effect of how close sensitivity/ specificity values are together


mainplot <- ggplot(maindf, aes(x = maxp, y = maxerr)) +
  #geom_point(aes(color = errordiffb, shape = errordiffb), size = 1.2) +
  geom_point()+
  labs(color = "Min difference\nbetween different\nsymptom sens/spec\nvalues", shape = "Min difference\nbetween different\nsymptom sens/spec\nvalues") +
  #scale_color_manual(values=c('Less than 0.05'="red", 'Greater than 0.05'="blue"), name="Min difference\n between values")+
  scale_y_continuous(name = "Max median estimate error", limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1))+
  scale_x_continuous(name = "Max percentile distance\nfrom prior median", limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1))+
  theme_minimal()+
  theme(text=element_text(size=12))
mainplot



