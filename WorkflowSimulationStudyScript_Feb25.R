
#This script is for examining a workflow where 25 surveyors assess 80 trees for two symptoms with no covariance and poor
#prior information across two sites

#Adjustments required to increase number of surveyors, number of trees, include covariance,
#and run with good, or very good priors are indicated.

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


#Variables we are going to save in loop#
actsevol1 <- c()
actsevol2 <- c()
actspvol1 <- c()
actspvol2 <- c()
actualsevol1 <- c()
actualsevol2 <- c()
actualspvol1 <- c()
actualspvol2 <- c()

estsealpha1 <- c()
estsealpha1u95 <- c()
estsealpha1l95 <- c()
estsebeta1 <- c()
estsebeta1u95 <- c()
estsebeta1l95 <- c()
estspalpha1<- c()
estspalpha1u95 <- c()
estspalpha1l95 <- c()
estspbeta1 <- c()
estspbeta1u95 <- c()
estspbeta1l95 <- c()

estsealpha1psrf <- c()
estsealpha1MCofSD <- c()
estsealpha1AC50 <- c()

estsebeta1psrf <- c()
estsebeta1MCofSD <- c()
estsebeta1AC50 <- c()

estspalpha1psrf <- c()
estspalpha1MCofSD <- c()
estspalpha1AC50 <- c()

estspbeta1psrf <- c()
estspbeta1MCofSD <- c()
estspbeta1AC50 <- c()

estsealpha2 <- c()
estsealpha2u95 <- c()
estsealpha2l95 <- c()
estsebeta2 <- c()
estsebeta2u95 <- c()
estsebeta2l95 <- c()
estspalpha2<- c()
estspalpha2u95 <- c()
estspalpha2l95 <- c()
estspbeta2 <- c()
estspbeta2u95 <- c()
estspbeta2l95 <- c()

estsealpha2psrf <- c()
estsealpha2MCofSD <- c()
estsealpha2AC50 <- c()

estsebeta2psrf <- c()
estsebeta2MCofSD <- c()
estsebeta2AC50 <- c()

estspalpha2psrf <- c()
estspalpha2MCofSD <- c()
estspalpha2AC50 <- c()

estspbeta2psrf <- c()
estspbeta2MCofSD <- c()
estspbeta2AC50 <- c()


estprev1 <- c()
estprev1psrf <- c()
estprev1MCofSD <- c()
estprev1AC50 <- c()
estprev1u95 <- c()
estprev1l95 <- c()

estprev2 <- c()
estprev2psrf <- c()
estprev2MCofSD <- c()
estprev2AC50 <- c()
estprev2u95 <- c()
estprev2l95 <- c()

actualprev1 <- c()
actualprev2 <- c()

#Saving the chains from the model runs#
trialneww <- c()
result <- c()
finalmc <- list()
trialfinal <- mcmc()
trialnew <- c()


######Simulate surveyor Sensitivity and Specificity for two symptoms on a tree, and there observations from these trees (i.e. presence/ absence of a symptom) #########

for(j in 1:25){   #Adjust number of surveyors assessed here (this run is for 25)
  
site1symptom <- rbinom(80, 1, prob = 0.30) #We are selecting 80 trees a surveyor will sample from site 1 with a probability of 0.30 the tree will actually be diseased
  
site2symptom <- rbinom(80, 1, prob = 0.60) #We are selecting 80 trees a surveyor will sample from site 1 with a probability of 0.30 the tree will actually be diseased
  
#Label trees on each site#
site1tree <- as.factor(seq(1,80))
site2tree <- as.factor(seq(1,80))
site1tree <- paste("a", site1tree, sep="")
site2tree <- paste("b", site1tree, sep="")  
  
  
##############Simulating surveyor sensitivity and specificity values from the observed AOD stem bleed data (i.e. a useful symptom to assess)#############

allsite <- read.csv()#read in AOD data in AODprocesseddata.csv file

#Specificity first#
#An estimate of confidence intervals
allsite$pe <- sqrt(allsite$BleedSpecificity*(1-allsite$BleedSpecificity)/allsite$nnegbleed)
allsite$ci <- qnorm(0.975)*allsite$pe

bleedspec <- allsite[,c(12,15,16)]
bleedspec$lower <- bleedspec$BleedSpecificity-bleedspec$ci
bleedspec$upper <- bleedspec$BleedSpecificity+bleedspec$ci

#Provide some broad starting parameters to estimate beta distribution parameters for the distribution of specificity#
startest <- epi.betabuster(median(bleedspec$BleedSpecificity), 0.95, TRUE, min(bleedspec$lower), 0.95)

#Now fit this to the specificity data to estimate the beta distribution parameters for the specificity of the observed AOD stem bleed data
bleedsp <- bleedspec$BleedSpecificity
trialsp <- fitdist(bleedsp,dist="beta",method="mle",start=startest[c(1,2)])
trialsp

#Now simulate specificity values from this
simd <- rbeta(10000, trialsp$estimate[[1]], trialsp$estimate[[2]])

#and sample two specificity values (one for symptom one, another for symptom two)
simvolsp <- sample(simd,size=2,replace=FALSE)

#NOW REPEAT FOR SENSITIVITY#
rm(allsite,startest,bleedsp,simd)

allsite <- read.csv()#read in AOD data in AODprocesseddata.csv file

allsite$pe <- sqrt(allsite$BleedSensitivity*(1-allsite$BleedSensitivity)/allsite$nposbleed)

allsite$ci <- qnorm(0.975)*allsite$pe

bleedsens <- allsite[,c(9,15,16)]
bleedsens$lower <- bleedsens$BleedSensitivity-bleedsens$ci
bleedsens$upper <- bleedsens$BleedSensitivity+bleedsens$ci
min(bleedsens$lower)

startest <- epi.betabuster(median(bleedsens$BleedSensitivity), 0.95, TRUE, min(bleedsens$lower), 0.95)

bleedse <- bleedsens$BleedSensitivity
trialse <- fitdist(bleedse,dist="beta",method="mle",start=startest[c(1,2)])
trialse

simd <- rbeta(10000, trialse$estimate[[1]], trialse$estimate[[2]])

simvolse <- sample(simd,size=2,replace=FALSE)

rm(allsite,startest,bleedse,simd)


#The parameters for the distributions are:
#Sensitivity: alpha = 10.654; beta = 5.246
#Specificity: alpha = 11.884; beta = 1.829


#######Now sample the trees for two symptoms with the simulated Sensitivity and Specificity value for each symptom#######

#Please note: I call them vol1 and vol2, but these actually represent symptom 1 and symptom 2#

#For Site 1#
vol1 <- c()
for(i in 1:length(site1symptom)){if(site1symptom[i]>0){
  vol1[i] <- rbinom(1,1,simvolse[1])
}else{vol1[i] <- rbinom(1,1,1-simvolsp[1])}}


#Without covariance:
vol2 <- c()
for(i in 1:length(site1symptom)){if(site1symptom[i]>0){
  vol2[i] <- rbinom(1,1,simvolse[2])
}else{vol2[i] <- rbinom(1,1,1-simvolsp[2])}}

#Below is the code for including covariance#
#vol2 <- c()
#for(i in 1:length(site1symptom)){if(site1symptom[i]>0 & vol1[i]>0){
  #vol2[i] <- rbinom(1,1,0.85)
#}else if(site1symptom[i]>0 & vol1[i]<1){
 # vol2[i] <- rbinom(1,1,0.15)}
  #else if (site1symptom[i]<1 & vol1[i]<1){
   # vol2[i] <- rbinom(1,1,0.15)}
 # else if (site1symptom[i]<1 & vol1[i]>0){
  #  vol2[i] <- rbinom(1,1,0.85) 
  #}}

#Now calculate the actual Sensitivity and Specificity of symptom 1 and symptom 2 on Site 1#

#Symptom 2#
scoresevol2 <- c()
scorespvol2 <- c()
for(i in 1:length(vol2)){if(vol2[i]==site1symptom[i]&site1symptom[i]==1){
  scoresevol2[i] <- 1} else if(vol2[i]!=site1symptom[i]&site1symptom[i]==1){scoresevol2[i] <- 0}
  else if(vol2[i]==site1symptom[i]&site1symptom[i]==0){
    scorespvol2[i] <- 1} else if(vol2[i]!=site1symptom[i]&site1symptom[i]==0){scorespvol2[i] <- 0}
}
scoresevol2 <- na.omit(scoresevol2)
scorespvol2  <- na.omit(scorespvol2)
scoresevol2 <- rep(sum(scoresevol2)/length(scoresevol2),time=80) #Adjust value here if sampling more than 80 trees#
scorespvol2 <- rep(sum(scorespvol2)/length(scorespvol2),time=80) #Adjust value here if sampling more than 80 trees#

#Symptom 1#
scoresevol1 <- c()
scorespvol1 <- c()
for(i in 1:length(vol1)){if(vol1[i]==site1symptom[i]&site1symptom[i]==1){
  scoresevol1[i] <- 1} else if(vol1[i]!=site1symptom[i]&site1symptom[i]==1){scoresevol1[i] <- 0}
  else if(vol1[i]==site1symptom[i]&site1symptom[i]==0){
    scorespvol1[i] <- 1} else if(vol1[i]!=site1symptom[i]&site1symptom[i]==0){scorespvol1[i] <- 0}
}
scoresevol1 <- na.omit(scoresevol1)
scorespvol1  <- na.omit(scorespvol1)
scoresevol1 <- rep(sum(scoresevol1)/length(scoresevol1),time=80) #Adjust value here if sampling more than 80 trees#
scorespvol1 <- rep(sum(scorespvol1)/length(scorespvol1),time=80) #Adjust value here if sampling more than 80 trees#

#now make the dataframe for site 1
datasite1 <- data.frame(site1tree,vol1,vol2,scoresevol1,scoresevol2,scorespvol1,scorespvol2)

rm(scoresevol1, scoresevol2, scorespvol1, scorespvol2) #Delete them because saved above#

#Now repeat for Site 2#

#Symptom 1#
vol1b <- c()
for(i in 1:length(site2symptom)){if(site2symptom[i]>0){
  vol1b[i] <- rbinom(1,1,simvolse[1])
}else{vol1b[i] <- rbinom(1,1,1-simvolsp[1])}}

#Symptom 2 without covariance#
vol2b <- c()
for(i in 1:length(site2symptom)){if(site2symptom[i]>0){
 vol2b[i] <- rbinom(1,1,simvolse[2])
 }else{vol2b[i] <- rbinom(1,1,1-simvolsp[2])}}

#Symptom 2 with covariance#
#vol2b <- c()
#for(i in 1:length(site2symptom)){if(site2symptom[i]>0 & vol1b[i]>0){
#  vol2b[i] <- rbinom(1,1,0.85)
#}else if(site2symptom[i]>0 & vol1b[i]<1){
#  vol2b[i] <- rbinom(1,1,0.15)}
#  else if (site2symptom[i]<1 & vol1b[i]<1){
#    vol2b[i] <- rbinom(1,1,0.15)}
#  else if (site2symptom[i]<1 & vol1b[i]>0){
#    vol2b[i] <- rbinom(1,1,0.85) 
#  }}

#Now calculate the actual Sensitivity and Specificity of volunteer 1 and volunteer 2 on Site 2#

scoresevol2 <- c()
scorespvol2 <- c()
for(i in 1:length(vol2b)){if(vol2b[i]==site2symptom[i]&site2symptom[i]==1){
  scoresevol2[i] <- 1} else if(vol2b[i]!=site2symptom[i]&site2symptom[i]==1){scoresevol2[i] <- 0}
  else if(vol2b[i]==site2symptom[i]&site2symptom[i]==0){
    scorespvol2[i] <- 1} else if(vol2b[i]!=site2symptom[i]&site2symptom[i]==0){scorespvol2[i] <- 0}
}
scoresevol2 <- na.omit(scoresevol2)
scorespvol2  <- na.omit(scorespvol2)
scoresevol2 <- rep(sum(scoresevol2)/length(scoresevol2),time=80) #Adjust value here if sampling more than 80 trees#
scorespvol2 <- rep(sum(scorespvol2)/length(scorespvol2),time=80) #Adjust value here if sampling more than 80 trees#

#vol1b
scoresevol1 <- c()
scorespvol1 <- c()
for(i in 1:length(vol1b)){if(vol1b[i]==site2symptom[i]&site2symptom[i]==1){
  scoresevol1[i] <- 1} else if(vol1b[i]!=site2symptom[i]&site2symptom[i]==1){scoresevol1[i] <- 0}
  else if(vol1b[i]==site2symptom[i]&site2symptom[i]==0){
    scorespvol1[i] <- 1} else if(vol1b[i]!=site2symptom[i]&site2symptom[i]==0){scorespvol1[i] <- 0}
}
scoresevol1 <- na.omit(scoresevol1)
scorespvol1  <- na.omit(scorespvol1)
scoresevol1 <- rep(sum(scoresevol1)/length(scoresevol1),time=80) #Adjust value here if sampling more than 80 trees#
scorespvol1 <- rep(sum(scorespvol1)/length(scorespvol1),time=80) #Adjust value here if sampling more than 80 trees#

#now make the dataframe for site 1
datasite2 <- data.frame(site2tree,vol1b,vol2b,scoresevol1,scoresevol2,scorespvol1,scorespvol2)

rm(scoresevol1, scoresevol2, scorespvol1, scorespvol2) #Delete them because saved above#

#Calculate covariance for the simulated survey dataset#
#Bind the two sites' symptom and observation data

allvol1 <- c(vol1,vol1b)
allvol2 <- c(vol2,vol2b)
allsitesymptom <- c(site1symptom, site2symptom)


allscores <- cbind(allsitesymptom,allvol1,allvol2)
allscoresdf <- data.frame(allscores)

#For calculation of the covariance - see Branscum et al., 2005 paper and Dendukuria and Joseph, 2001 paper for equations#
#Sensitivity#
se11 <- c() #Both symptom observations are TP#
se11 <- 1:160 #Adjust value here if sampling more than 80 trees i.e. 200 for 100 trees assessed each site, 240 for 120 trees assessed each site#
for(i in 1:nrow(allscoresdf)){if(allsitesymptom[i]>0 & allvol2[i]>0 & allvol1[i]>0){se11[i] <- 1}else{
  se11[i] <- 0
}}
se21 <- c() #Symptom 2 is TP#
se21 <- 1:160 #Adjust value here if sampling more than 80 trees i.e. 200 for 100 trees assessed each site, 240 for 120 trees assessed each site#
for(i in 1:nrow(allscoresdf)){if(allsitesymptom[i]>0 & allvol2[i]>0){se21[i] <- 1}else{
  se21[i] <- 0
}}
se12 <- c() #Symptom 1 is TP#
se12 <- 1:160 #Adjust value here if sampling more than 80 trees i.e. 200 for 100 trees assessed each site, 240 for 120 trees assessed each site#
for(i in 1:nrow(allscoresdf)){if(allsitesymptom[i]>0 & allvol1[i]>0){se12[i] <- 1}else{
  se12[i] <- 0
}}

cseaa <- c()
csea <- c()
cseb <- c()


cseaa <- (sum(se11))/(sum(allsitesymptom)) #Probability of both symptoms giving a TP#
csea <- (sum(se12))/(sum(allsitesymptom)) #Probability of symptom one giving a TP#
cseb <- (sum(se21))/(sum(allsitesymptom)) #Probability of symptom two giving a TP#


#therefore to calculate the covariance between tests when the host is positive
actualcovdp[j] <- cseaa - (csea*cseb) #CovD+
maxcovdp[j] <- min(cseb,csea) - (csea*cseb) #the maximum possible CovD+ value
mincovdp[j] <- (csea-1)*(1-cseb) #the minimum possible CovD+ value


#Specificity#
sp11 <- c() #Both symptom observations are TN#
sp11 <- 1:160 #Adjust value here if sampling more than 80 trees i.e. 200 for 100 trees assessed each site, 240 for 120 trees assessed each site#
for(i in 1:nrow(allscoresdf)){if(allsitesymptom[i]<1 & allvol2[i]<1 & allvol1[i]<1){sp11[i] <- 1}else{
  sp11[i] <- 0
}}
sp21 <- c() #Symptom 2 is TN#
sp21 <- 1:160 #Adjust value here if sampling more than 80 trees i.e. 200 for 100 trees assessed each site, 240 for 120 trees assessed each site#
for(i in 1:nrow(allscoresdf)){if(allsitesymptom[i]<1 & allvol2[i]<1){sp21[i] <- 1}else{
  sp21[i] <- 0
}}
sp12 <- c() #Symptom 1 is TN#
sp12 <- 1:160 #Adjust value here if sampling more than 80 trees i.e. 200 for 100 trees assessed each site, 240 for 120 trees assessed each site#
for(i in 1:nrow(allscoresdf)){if(allsitesymptom[i]<1 & allvol1[i]<1){sp12[i] <- 1}else{
  sp12[i] <- 0
}}

cspaa <- (sum(sp11))/(160-(sum(allsitesymptom))) #Probability of both symptoms giving a TN: Adjust value here if sampling more than 80 trees i.e. 200 for 100 trees assessed each site, 240 for 120 trees assessed each site#
cspa <- (sum(sp12))/(160-(sum(allsitesymptom))) #Probability of symptom one giving a TN: Adjust value here if sampling more than 80 trees i.e. 200 for 100 trees assessed each site, 240 for 120 trees assessed each site#
cspb <- (sum(sp21))/(160-(sum(allsitesymptom))) #Probability of symptom two giving a TN: Adjust value here if sampling more than 80 trees i.e. 200 for 100 trees assessed each site, 240 for 120 trees assessed each site#


#therefore to calculate the covariance between tests when the host is negative

actualcovdn[j] <- cspaa - (cspa*cspb) #CovD-
maxcovdn[j] <- min(cspb,cspa) - (cspa*cspb) #the maximum possible CovD- value
mincovdn[j] <- (cspa-1)*(1-cspb) #the minimum possible CovD- value


#The sensitivity and specificity value used from the simulated distribution to simulate the survey dataset#
#acsimsevolunteer2[j] <- simvolse[2] # add in if we want to check for symptom 2 when no covariance is simulated
acsimsevolunteer1[j] <- simvolse[1]
#acsimspvolunteer2[j] <- simvolsp[2] # add in if we want to check for symptom 2 when no covariance is simulated
acsimspvolunteer1[j] <- simvolsp[1]

#Save simulated survey datasets here#
surveysymptom[((j*length(allsitesymptom)+1)-length(allsitesymptom)):(j*length(allsitesymptom))] <- allsitesymptom
surveyvol1[((j*length(allvol1)+1)-length(allvol1)):(j*length(allvol1))] <- allvol1
surveyvol2[((j*length(allvol2)+1)-length(allvol2)):(j*length(allvol2))] <- allvol2

#
#######Simulated observation dataset organised, now feed through the model#########


#########################Setting Priors#########################

##################Poor prior information both symptoms' sensitivity and specificity###################
spBetat1 <- epi.betabuster(median(bleedspec$BleedSpecificity)-0.2, 0.95, TRUE, 0.10, 0.95) #For specificity: Provides alpha and beta parameters from beta distribution with 'poor' prior information#


alphaspt1 <- c(rep(spBetat1$shape1,time=1)) #Saving the specificity alpha parameter prior information for symptom 1#
betaspt1 <- c(rep(spBetat1$shape2, time=1)) #Saving the specificity beta parameter prior information for symptom 1#

#Precision
precspt1 <- c(rep(1/spBetat1$variance,time=1)) #Saving the specificity precision prior value for symptom 1#

#Symptom two uses the same priors as symptom one
alphaspt2 <- c(rep(spBetat1$shape1,time=1)) #Saving the specificity alpha parameter prior information for symptom 2#
betaspt2 <- c(rep(spBetat1$shape2, time=1)) #Saving the specificity beta parameter prior information for symptom 2#
precspt2 <- c(rep(1/spBetat1$variance,time=1)) #Saving the specificity precision prior information for symptom 2#


seBetat1 <- epi.betabuster(median(bleedsens$BleedSensitivity)-0.2, 0.95, TRUE, 0.10, 0.95) #For sensitivity: Provides alpha and beta parameters from beta distribution with 'poor' prior information#


alphaset1 <- c(rep(seBetat1$shape1, time=1)) #Saving the sensitivity alpha parameter prior information for symptom 1#
betaset1 <- c(rep(seBetat1$shape2, time=1)) #Saving the sensitivity beta parameter prior information for symptom 1#

#Precision
precset1 <- c(rep(1/seBetat1$variance, time=1)) #Saving the sensitivity precision parameter prior information for symptom 1#

#Symptom two uses the same priors as symptom one
alphaset2 <- c(rep(seBetat1$shape1, time=1)) #Saving the sensitivity alpha parameter prior information for symptom 2#
betaset2 <- c(rep(seBetat1$shape2, time=1)) #Saving the sensitivity beta parameter prior information for symptom 2#

#Next, precision
precset2 <- c(rep(1/seBetat1$variance, time=1)) #Saving the sensitivity precision parameter prior information for symptom 2#

##################Use code below for poor prior information symptom one sensitivity and specificity, very good information symptom two sensitivity and specificity###################

#spBetat1 <- epi.betabuster(median(bleedspec$BleedSpecificity)-0.2, 0.95, TRUE, 0.10, 0.95) #For specificity: Provides alpha and beta parameters from beta distribution with 'poor' prior information#

#alphaspt1 <- c(rep(spBetat1$shape1,time=1)) #Saving the specificity alpha parameter prior information for symptom 1#
#betaspt1 <- c(rep(spBetat1$shape2, time=1)) #Saving the specificity beta parameter prior information for symptom 1#

#Precision
#precspt1 <- c(rep(1/spBetat1$variance,time=1)) #Saving the specificity precision prior value for symptom 1#


#alphaspt2 <- trialsp$estimate[[1]] #For specificity: symptom 2, taking the alpha parameter the data of the distribution data were simulated from
#betaspt2 <- trialsp$estimate[[2]] #For specificity: symptom 2, taking the beta parameter the data of the distribution data were simulated from

#Variance parameter for specificity from this distribution for symptom 2
#spvar2 <- (alphaspt2 * betaspt2) / (((alphaspt2 + betaspt2)^2) * (alphaspt2 + betaspt2 + 1))

#Next, calculate the precision prior parameter for specificity symptom 2
#precspt2 <- c(rep(1/spvar2,time=1))

#seBetat1 <- epi.betabuster(median(bleedsens$BleedSensitivity)-0.2, 0.95, TRUE, 0.10, 0.95) #For sensitivity: Provides alpha and beta parameters from beta distribution with 'poor' prior information#


#alphaset1 <- c(rep(seBetat1$shape1, time=1)) #Saving the sensitivity alpha parameter prior information for symptom 1#
#betaset1 <- c(rep(seBetat1$shape2, time=1)) #Saving the sensitivity beta parameter prior information for symptom 1#

#Precision
#precset1 <- c(rep(1/seBetat1$variance, time=1)) #Saving the sensitivity precision prior value for symptom 1#

#Test two uses informed priors
#alphaset2 <- trialse$estimate[[1]] #For sensitivity: symptom 2, taking the alpha parameter the data of the distribution data were simulated from
#betaset2 <- trialse$estimate[[2]] #For sensitivity: symptom 2, taking the beta parameter the data of the distribution data were simulated from

#Variance parameter for sensitivity from this distribution for symptom 2
#sevar2 <- (alphaset2 * betaset2) / (((alphaset2 + betaset2)^2) * (alphaset2 + betaset2 + 1))

#Next, calculate the precision prior parameter for sensitivity symptom 2
#precset2 <- c(rep(1/sevar2,time=1))


##################Use code below for good information both symptoms' sensitivity and specificity###################

#alphaspt1 <- c(rep(4.383722,time=1)) #Saving the specificity alpha parameter prior information for symptom 1#
#betaspt1 <- c(rep(1.728984, time=1)) #Saving the specificity beta parameter prior information for symptom 1#

#Variance parameter for specificity from this distribution for symptom 1
#spvar1 <- (alphaspt1 * betaspt1) / (((alphaspt1 + betaspt1)^2) * (alphaspt1 + betaspt1 + 1))

#Calculate the precision prior parameter for symptom 1
#precspt1 <- c(rep(1/spvar1,time=1))

#Symptom two uses the same priors as symptom one
#alphaspt2 <- c(rep(4.383722,time=1)) #Saving the specificity alpha parameter prior information for symptom 2#
#betaspt2 <- c(rep(1.728984, time=1)) #Saving the specificity beta parameter prior information for symptom 2#

#Variance parameter for specificity from this distribution for symptom 2
#spvar2 <- (alphaspt2 * betaspt2) / (((alphaspt2 + betaspt2)^2) * (alphaspt2 + betaspt2 + 1))

#Calculate the precision prior parameter for specificity symptom 2
#precspt2 <- c(rep(1/spvar2,time=1))


#alphaset1 <- c(rep(3.154012, time=1)) #Saving the sensitivity alpha parameter prior information for symptom 1#
#betaset1 <- c(rep(2.146211, time=1)) #Saving the sensitivity beta parameter prior information for symptom 1#

#Variance parameter for sensitivity from this distribution for symptom 1
#sevar1 <- (alphaset1 * betaset1) / (((alphaset1 + betaset1)^2) * (alphaset1 + betaset1 + 1))

#Calculate the precision prior parameter for sensitivity symptom 1
#precset1 <- c(rep(1/sevar1,time=1))

#Symptom two uses the same priors as symptom one
#alphaset2 <- c(rep(3.154012, time=1)) #Saving the sensitivity alpha parameter prior information for symptom 2#
#betaset2 <- c(rep(2.146211, time=1)) #Saving the sensitivity beta parameter prior information for symptom 2#

#Variance parameter for sensitivity from this distribution for symptom 2
#sevar2 <- (alphaset2 * betaset2) / (((alphaset2 + betaset2)^2) * (alphaset2 + betaset2 + 1))

#Calculate the precision prior parameter for sensitivity symptom 2
#precset2 <- c(rep(1/sevar2,time=1))


#################################Site prevalence priors#####################################################

nsite <- c(1:2) #Defining the number of sites#

#Extracting the alpha and beta parameters for the priors
alphaprev <- c(3,6) #alpha prior parameter for site 1 and 2 respectively
betaprev <- c(8, 4) #beta prior parameter for site 1 and 2 respectively

#Calculate variance for site 1 and 2 respectively#
varprev <- c((alphaprev[1] * betaprev[1]) / (((alphaprev[1] + betaprev[1])^2) * (alphaprev[1] + betaprev[1] + 1)),(alphaprev[2] * betaprev[2]) / (((alphaprev[2] + betaprev[2])^2) * (alphaprev[2] + betaprev[2] + 1)))

precprev <- c(1/varprev[1], 1/varprev[2]) #Provide precision prior parameters for site 1 and 2 respectively


############Covariance prior information#######################
#Generalised beta distribution: define a uniform beta distribution, then change the bounds of the distribution in the model#
alphacovdn <- 1 #alpha parameter covariance between 
betacovdn <- 1 #beta parameter covariance between 
alphacovdp <- 1 #alpha parameter covariance between 
betacovdp <- 1 #beta parameter covariance between 

#Variance for covariance between disease negatives#
varcovdn <- (alphacovdn * betacovdn) / (((alphacovdn + betacovdn)^2) * (alphacovdn + betacovdn + 1))

#Precision for covariance between disease negatives#
preccovdn <- c(rep(1/varcovdn,time=1))

#Variance for covariance between disease negatives#
varcovdp <- (alphacovdp * betacovdp) / (((alphacovdp + betacovdp)^2) * (alphacovdp + betacovdp + 1))

#Precision for covariance between disease negatives#
preccovdp <- c(rep(1/varcovdp,time=1))

########################Ok priors are set#######################################

#Formatting the observation data so it is consistent with number of observations for each site:
#Symptom 1 positive, symptom 2 positive,
#Symptom 1 positive, symptom 2 negative,
#Symptom 1 negative, symptom 2 positive,
#Symptom 1 negative, symptom 2 negative.

#Also specify the number of hosts inspected (in this case 80)
rm(vol1,vol2,vol1b,vol2b)
s1allpcomb <- datasite1 %>% mutate(pcombaa=vol1==1&vol2==1,
                                   pcombab=vol1==1&vol2==0,
                                   pcombba=vol1==0&vol2==1,
                                   pcombbb=vol1==0&vol2==0)

s1allpcomb[,c(8:11)]<- as.integer(c(s1allpcomb$pcombaa,s1allpcomb$pcombab,
                                    s1allpcomb$pcombba,s1allpcomb$pcombbb))

s1allpcomb <- s1allpcomb %>% mutate(pcombaa=sum(pcombaa), pcombab=sum(pcombab),
                                    pcombba=sum(pcombba), pcombbb=sum(pcombbb))

nhosts <- length(s1allpcomb$site1tree)

s2allpcomb <- datasite2 %>% mutate(pcombaa=vol1==1&vol2==1,
                                   pcombab=vol1==1&vol2==0,
                                   pcombba=vol1==0&vol2==1,
                                   pcombbb=vol1==0&vol2==0)

s2allpcomb[,c(8:11)]<- as.integer(c(s2allpcomb$pcombaa,s2allpcomb$pcombab,
                                    s2allpcomb$pcombba,s2allpcomb$pcombbb))

s2allpcomb <- s2allpcomb %>% mutate(pcombaa=sum(pcombaa), pcombab=sum(pcombab),
                                    pcombba=sum(pcombba), pcombbb=sum(pcombbb))

nhosts[2] <- length(s2allpcomb$site2tree)

y = matrix(c(s1allpcomb$pcombaa[1], s2allpcomb$pcombaa[1],
             s1allpcomb$pcombab[1], s2allpcomb$pcombab[1],
             s1allpcomb$pcombba[1], s2allpcomb$pcombba[1],
             s1allpcomb$pcombbb[1], s2allpcomb$pcombbb[1]),nrow=2,ncol=4)

#####Ready to go running models#################

#Specificy data for model#

datalist <- list(alphaprev,alphaset1,alphaset2,alphaspt1,alphaspt2,betaprev,
                 betaset1,betaset2,betaspt1,betaspt2,nhosts,nsite,
                 precprev,precset1,precset2,precspt1,precspt2,y)
#,alphacovdp,betacovdp,alphacovdn,betacovdn,preccovdp,preccovdn) #Add these parameters in if running model with covariance#

#name data for model#

names(datalist) <- c("alphaprev","alphase1","alphase2","alphasp1","alphasp2","betaprev",
                     "betase1","betase2","betasp1","betasp2","nhosts","nsite",
                     "precprev","precse1","precse2","precsp1","precsp2","y")
#,"alphacovdp","betacovdp","alphacovdn","betacovdn","preccovdp","preccovdn") #Add these parameters in if running model with covariance#


runjags.options(jagspath=#Insert location of JAGS on computer,
                force.summary=TRUE)
# Set starting values with 10 chains so we can assess convergence
inits1 = list(".RNG.name" ="base::Super-Duper"
              ,".RNG.seed" = 100022)
inits2 = list(".RNG.name" ="base::Super-Duper"
              ,".RNG.seed" = 200022)
inits3 = list(".RNG.name" ="base::Super-Duper"
              ,".RNG.seed" = 300022)
inits4 = list(".RNG.name" ="base::Super-Duper"
              ,".RNG.seed" = 400022)
inits5 = list(".RNG.name" ="base::Super-Duper"
              ,".RNG.seed" = 500022)
inits6 = list(".RNG.name" ="base::Super-Duper"
              ,".RNG.seed" = 600022)
inits7 = list(".RNG.name" ="base::Super-Duper"
              ,".RNG.seed" = 700022)
inits8 = list(".RNG.name" ="base::Super-Duper"
              ,".RNG.seed" = 800022)
inits9 = list(".RNG.name" ="base::Super-Duper"
              ,".RNG.seed" = 900022)
inits10 = list(".RNG.name" ="base::Super-Duper"
               ,".RNG.seed" = 1000022)
inits <- list(inits1, inits2,inits3,inits4,inits5,inits6,inits7,inits8,inits9,inits10)

##################################Model with covariance################################
model_string <- "model{
#Likelihood function for the multinomial distribution of the probability of each test outcome#
for (j in 1:length(nsite)){
y[j,1:4] ~ dmulti(pcomb[j,1:4],nhosts[j])
pcomb[j,1] <- prev[j]*se1*se2+(1-prev[j])*(1-sp1)*(1-sp2)
pcomb[j,2] <- prev[j]*se1*(1-se2)+(1-prev[j])*(1-sp1)*sp2
pcomb[j,3] <- prev[j]*(1-se1)*se2+(1-prev[j])*sp1*(1-sp2)
pcomb[j,4] <- prev[j]*(1-se1)*(1-se2)+(1-prev[j])*sp1*sp2

#Use the below lines rather than above to specificy the model with covariance parameters#
#pcomb[j,1] <- prev[j]*(se1*se2+covdp)+(1-prev[j])*((1-sp1)*(1-sp2)+covdn)
#pcomb[j,2] <- prev[j]*(se1*(1-se2)-covdp)+(1-prev[j])*((1-sp1)*sp2-covdn)
#pcomb[j,3] <- prev[j]*((1-se1)*se2-covdp)+(1-prev[j])*(sp1*(1-sp2)-covdn)
#pcomb[j,4] <- prev[j]*((1-se1)*(1-se2)+covdp)+(1-prev[j])*(sp1*sp2+covdn)



#Priors and reparameterising

prev[j] ~ dbeta(alphaprevpost[j], betaprevpost[j])T(0.01,0.99)
alphaprevpost[j] <- muprev[j]*phiprev[j]
betaprevpost[j] <- (1-muprev[j])*phiprev[j] 
muprev[j] ~ dbeta(alphaprev[j],betaprev[j])T(0.01,0.99)
phiprev[j] ~dgamma(precprev[j],1)}

se1 ~dbeta(alphasep1,betasep1)T(0.01,0.99)
alphasep1 <- muse1*phise1
betasep1 <- (1-muse1)*phise1
muse1 ~ dbeta(alphase1,betase1)T(0.01,0.99)
phise1 ~dgamma(precse1,1)

se2 ~dbeta(alphasep2,betasep2)T(0.01,0.99)
alphasep2 <- muse2*phise2
betasep2 <- (1-muse2)*phise2
muse2 ~ dbeta(alphase2,betase2)T(0.01,0.99)
phise2 ~dgamma(precse2,1)

sp1 ~ dbeta(alphaspp1,betaspp1)T(0.01,0.99)
alphaspp1 <- musp1*phisp1
betaspp1 <- (1-musp1)*phisp1
musp1 ~ dbeta(alphasp1,betasp1)T(0.01,0.99)
phisp1 ~dgamma(precsp1,1)

sp2 ~ dbeta(alphaspp2,betaspp2)T(0.01,0.99)
alphaspp2 <- musp2*phisp2
betaspp2 <- (1-musp2)*phisp2
musp2 ~ dbeta(alphasp2,betasp2)T(0.01,0.99)
phisp2 ~dgamma(precsp2,1)

#Add in below code if running covariance model#
#Specified as in Dendukuri and Joseph, 2001 (see Branscum et al., 2005 too), where we have a beta distribution stretched, or shrunk to correspond to a range other than the standrard 0,1#
#covdp <- simcovdp*((min(se1,se2)-se1*se2)-((se1-1)*(1-se2)))+((se1-1)*(1-se2))
#simcovdp~dbeta(alphacovdpp,betacovdpp)T(0.01,0.99)
#alphacovdpp <- mucovdp*phicovdp
#betacovdpp <- (1-mucovdp)*phicovdp
#mucovdp ~ dbeta(alphacovdp,betacovdp)T(0.01,0.99)
#phicovdp ~dgamma(preccovdp,1)

#Add in below code if running covariance model#
Specified as in Dendukuri and Joseph, 2001 (see Branscum et al., 2005 too), where we have a beta distribution stretched, or shrunk to correspond to a range other than the standrard 0,1#
#covdn <- simcovdn*((min(sp1,sp2)-sp1*sp2)-((sp1-1)*(1-sp2)))+((sp1-1)*(1-sp2))
#simcovdn~dbeta(alphacovdnn,betacovdnn)T(0.01,0.99)
#alphacovdnn <- mucovdn*phicovdn
#betacovdnn <- (1-mucovdn)*phicovdn
#mucovdn ~ dbeta(alphacovdn,betacovdn)T(0.01,0.99)
#phicovdn ~dgamma(preccovdn,1)

}"

#Monitor parameters of interest, and specificy the run details#
result <- run.jags(model_string, monitor=c('alphasep1','betasep1','alphaspp1', 'betaspp1',
                                           'alphasep2','betasep2','alphaspp2', 'betaspp2','prev'
                                           #,'alphacovdnn','betacovdnn','alphacovdpp','betacovdpp' Add in if want to check covariance parameters
                                           ),data=c(datalist), inits=inits,
                   n.chains=length(inits),
                   burnin = 10000,
                   sample = 500000,
                   adapt = 2000,confidence=c(0.95), method="parallel")
#Save the output#
s <- summary(result, autocorr.lags=50,confidence=c(0.95))

#Specify that if a chain has a psrf value of less than 1.1 it should be discarded#
#NOTE: This code only discards the first run if the psrf value was less than 1.1 - will need to check psrf values in the
#output excel file and remove runs as appropriate
if (result$psrf[[1]][1]>1.1|result$psrf[[1]][2]>1.1|result$psrf[[1]][3]>1.1|result$psrf[[1]][4]>1.1
    |result$psrf[[1]][4]>1.1|result$psrf[[1]][5]>1.1|result$psrf[[1]][6]>1.1|result$psrf[[1]][7]>1.1
    |result$psrf[[1]][8]>1.1|result$psrf[[1]][9]>1.1|result$psrf[[1]][10]>1.1){
  rm(result)}else{trialnew <- as.mcmc(result)}

#Keep the first chain, thin at every 10th value after the initial 12,000 iterations
#This means each runs provides 48,800 alpha and beta value estimates for the beta distributions of sensitivity and specificity
if (length(trialnew)>1){trialneww <- window(trialnew, start = 1, end = 500000, thin=10)}

if(length(trialfinal)>1){
  trialfinal <- mcmc(rbind(trialfinal, trialneww))}else{trialfinal <- trialneww}

trialnew <- c()


#Saving important values from each run#
actualsevol1[j] <- (datasite1$scoresevol1[1]+datasite2$scoresevol1[1])/2 #Sensitivity of symptom 1 across both sites
actualsevol2[j] <- (datasite1$scoresevol2[1]+datasite2$scoresevol2[1])/2 #Sensitivity of symptom 2 across both sites
actualspvol1[j] <- (datasite1$scorespvol1[1]+datasite2$scorespvol1[1])/2 #Specificity of symptom 1 across both sites
actualspvol2[j] <- (datasite1$scorespvol2[1]+datasite2$scorespvol2[1])/2 #Specificity of symptom 2 across both sites

actualprev1[j] <- sum(site1symptom)/length(site1symptom) #The true disease prevalence of sampled trees on site 1
actualprev2[j] <- sum(site2symptom)/length(site2symptom) #The true disease prevalence of sampled trees on site 2

estsealpha1[j] <- s[1,2] #The most likely alpha parameter for sensitivity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)
estsealpha1psrf[j] <- s[1,11] #The psrf value for alpha parameter sensitivity of symptom 1
estsealpha1MCofSD[j] <- s[1,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for alpha parameter sensitivity of symptom 1
estsealpha1AC50[j] <- s[1,10] #The autocorrelation after 50 iterations for alpha parameter sensitivity of symptom 1
estsealpha1u95[j] <- s[1,3] #The upper 95 estimate of alpha parameter for sensitivity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)
estsealpha1l95[j] <- s[1,1] #The lower 95 estimate of alpha parameter for sensitivity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)

estsealpha2[j] <- s[5,2] #The most likely alpha parameter for sensitivity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)
estsealpha2psrf[j] <- s[5,11] #The psrf value for alpha parameter sensitivity of symptom 2
estsealpha2MCofSD[j] <- s[5,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for alpha parameter sensitivity of symptom 2
estsealpha2AC50[j] <- s[5,10] #The autocorrelation after 50 iterations for alpha parameter sensitivity of symptom 2
estsealpha2u95[j] <- s[5,3] #The upper 95 estimate of alpha parameter for sensitivity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)
estsealpha2l95[j] <- s[5,1] #The lower 95 estimate of alpha parameter for sensitivity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)

estsebeta1[j] <- s[2,2] #The most likely beta parameter for sensitivity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)
estsebeta1psrf[j] <- s[2,11] #The psrf value for beta parameter sensitivity of symptom 1
estsebeta1MCofSD[j] <- s[2,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for beta parameter sensitivity of symptom 1
estsebeta1AC50[j] <- s[2,10] #The autocorrelation after 50 iterations for beta parameter sensitivity of symptom 1
estsebeta1u95[j] <- s[2,3] #The upper 95 estimate of alpha parameter for sensitivity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)
estsebeta1l95[j] <- s[2,1] #The lower 95 estimate of alpha parameter for sensitivity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)

estsebeta2[j] <- s[6,2] #The most likely beta parameter for sensitivity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)
estsebeta2psrf[j] <- s[6,11] #The psrf value for beta parameter sensitivity of symptom 2
estsebeta2MCofSD[j] <- s[6,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for beta parameter sensitivity of symptom 2
estsebeta2AC50[j] <- s[6,10] #The autocorrelation after 50 iterations for beta parameter sensitivity of symptom 2
estsebeta2u95[j] <- s[6,3] #The upper 95 estimate of alpha parameter for sensitivity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)
estsebeta2l95[j] <- s[6,1] #The lower 95 estimate of alpha parameter for sensitivity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)

estspalpha1[j] <- s[3,2] #The most likely alpha parameter for specificity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)
estspalpha1psrf[j] <- s[3,11] #The psrf value for alpha parameter specificity of symptom 1
estspalpha1MCofSD[j] <- s[3,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for alpha parameter specificity of symptom 1
estspalpha1AC50[j] <- s[3,10] #The autocorrelation after 50 iterations for alpha parameter specificity of symptom 1
estspalpha1u95[j] <- s[3,3] #The upper 95 estimate of alpha parameter for specificity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)
estspalpha1l95[j] <- s[3,1] #The lower 95 estimate of alpha parameter for specificity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)

estspalpha2[j] <- s[7,2] #The lower 95 estimate of alpha parameter for sensitivity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)
estspalpha2psrf[j] <- s[7,11] #The psrf value for alpha parameter specificity of symptom 2
estspalpha2MCofSD[j] <- s[7,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for alpha parameter specificity of symptom 2
estspalpha2AC50[j] <- s[7,10] #The autocorrelation after 50 iterations for alpha parameter specificity of symptom 2
estspalpha2u95[j] <- s[7,3] #The upper 95 estimate of alpha parameter for specificity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)
estspalpha2l95[j] <- s[7,1] #The lower 95 estimate of alpha parameter for specificity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)

estspbeta1[j] <- s[4,2] #The most likely beta parameter for specificity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)
estspbeta1psrf[j] <- s[4,11] #The psrf value for beta parameter specificity of symptom 1
estspbeta1MCofSD[j] <- s[4,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for beta parameter specificity of symptom 1
estspbeta1AC50[j] <- s[4,10] #The autocorrelation after 50 iterations for beta parameter specificity of symptom 1
estspbeta1u95[j] <- s[4,3] #The upper 95 estimate of beta parameter for specificity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)
estspbeta1l95[j] <- s[4,1] #The lower 95 estimate of beta parameter for specificity of symptom 1 according to runjags summary statistics (I do not use this value in analyses)

estspbeta2[j] <- s[8,2] #The most likely beta parameter for specificity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)
estspbeta2psrf[j] <- s[8,11] #The psrf value for beta parameter specificity of symptom 2
estspbeta2MCofSD[j] <- s[8,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for beta parameter specificity of symptom 2
estspbeta2AC50[j] <- s[8,10] #The autocorrelation after 50 iterations for beta parameter specificity of symptom 2
estspbeta2u95[j] <- s[8,3] #The upper 95 estimate of beta parameter for specificity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)
estspbeta2l95[j] <- s[8,1] #The lower 95 estimate of beta parameter for specificity of symptom 2 according to runjags summary statistics (I do not use this value in analyses)

estprev1[j] <- s[9,2] #The most likely site one true disease prevalence parameter according to runjags summary statistics (I do not use this value in analyses)
estprev1psrf[j] <- s[9,11] #The psrf value for site one true disease prevalence parameter
estprev1MCofSD[j] <- s[9,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for site one true disease prevalence parameter
estprev1AC50[j] <- s[9,10] #The autocorrelation after 50 iterations for site one true disease prevalence parameter
estprev1u95[j] <- s[9,3] #The upper 95 estimate of site one true disease prevalence parameter according to runjags summary statistics (I do not use this value in analyses)
estprev1l95[j] <- s[9,1] #The lower 95 estimate of site one true disease prevalence parameter according to runjags summary statistics (I do not use this value in analyses)

estprev2[j] <- s[10,2] #The most likely site one true disease prevalence parameter according to runjags summary statistics (I do not use this value in analyses)
estprev2psrf[j] <- s[10,11] #The psrf value for true disease prevalence parameter
estprev2MCofSD[j] <- s[10,8] #The Monte Carlo standard error expressed as a percentage of the standard deviation for site two true disease prevalence parameter
estprev2AC50[j] <- s[10,10] #The autocorrelation after 50 iterations for site two true disease prevalence parameter
estprev2u95[j] <- s[10,3] #The upper 95 estimate of site one true disease prevalence parameter according to runjags summary statistics (I do not use this value in analyses)
estprev2l95[j] <- s[10,1] #The lower 95 estimate of site one true disease prevalence parameter according to runjags summary statistics (I do not use this value in analyses)

}


final <- data.frame(cbind(estsealpha1, estsealpha1u95, estsealpha1l95, estsebeta1, estsebeta1u95, estsebeta1l95,
                          estspalpha1, estspalpha1u95, estspalpha1l95, estspbeta1, estspbeta1u95, estspbeta1l95, 
                          estsealpha2, estsealpha2u95, estsealpha2l95, estsebeta2, estsebeta2u95, estsebeta2l95,
                          estspalpha2, estspalpha2u95, estspalpha2l95, estspbeta2, estspbeta2u95, estspbeta2l95,
                          actualsevol1,
                          actualsevol2, actualspvol1, actualspvol2,
                          estprev1, estprev1u95, estprev1l95,
                          estprev2, estprev2u95, estprev2l95,
                          actualprev1, actualprev2,
                          estsealpha1psrf, estsealpha1MCofSD, estsealpha1AC50,
                          estsebeta1psrf, estsebeta1MCofSD, estsebeta1AC50,
                          estspalpha1psrf, estspalpha1MCofSD, estspalpha1AC50,
                          estspbeta1psrf, estspbeta1MCofSD, estspbeta1AC50,
                          estsealpha2psrf, estsealpha2MCofSD, estsealpha2AC50,
                          estsebeta2psrf, estsebeta2MCofSD, estsebeta2AC50,
                          estspalpha2psrf, estspalpha2MCofSD, estspalpha2AC50,
                          estspbeta2psrf, estspbeta2MCofSD, estspbeta2AC50,
                          #actualcovdn,maxcovdn,mincovdn,actualcovdp,maxcovdp,mincovdp, #Add in line if want to check actual covaraince from simulated data - and minimum and maximum possible covariance values#
                          acsimsevolunteer1,acsimspvolunteer1))

final$nsamples <- 80 #change if sample more trees

write.csv(final, file ="")#Insert location to write csv file of results summary table

saveRDS(trialfinal,) #Insert location to write rds file of the mcmc chain results

