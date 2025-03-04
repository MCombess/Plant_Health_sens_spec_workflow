
rm(list=ls())
unloadNamespace("car")
library(tidyverse)
library(dplyr)
library(DescTools)
library(ggplot2)
library(grid)
library(gridExtra)

#Read in the AOD data set and get exit hole sens and spec#

#########Part 1: the sensitivity and specificity plots regardless of symptom number#######

#Read in the volunteer AOD volunteer data
vdata <- read.csv()

#Correct the types of data
vdata$Site <- as.factor(vdata$Site)
vdata$Date <- as.Date(vdata$Date, "%d/%m/%Y",  tz="GMT")
vdata$SurveyorID <- as.factor(vdata$SurveyorID)
vdata$TreeNumber <- as.factor(vdata$TreeNumber)

#Change 20+ to 20 for StemBleeds and BarkCracks. Then convert to integer data.
vdata <- vdata %>% mutate(StemBleeds=recode(StemBleeds, '20+'='20')) %>% mutate(BarkCracks=recode(BarkCracks, '20+'='20'))

vdata$StemBleeds <- as.integer(vdata$StemBleeds)
vdata$BarkCracks <- as.integer(vdata$BarkCracks)

#Present or absent for volunteer data
vdata <- vdata %>% mutate(Stembleeds=StemBleeds, Barkcracks=BarkCracks, Exitholes=ExitHoles)

vdata$Stembleeds <- vdata$Stembleeds>0
vdata$Barkcracks <- vdata$Barkcracks>0
vdata$Exitholes <- vdata$Exitholes!=0

#Read in the reference data collected by expert
rdata <- read.csv()
#Correct the types of data
rdata$Site <- as.factor(rdata$Site)
rdata$TreeNumber <- as.factor(rdata$TreeNumber)

#Again, change 20+ to 20 for BarkCracks in expert dataset. Then convert to integer data.
rdata <- rdata %>%  mutate(BarkCracks=recode(BarkCracks, '20+'='30')) #We will use 30 as a replacement for the category 20+ for statistical analyses purposes
rdata <- rdata %>%  mutate(ExitHoles=recode(ExitHoles, '20+'='30')) #We will use 30 as a replacement for the category 20+ for statistical analyses purposes
rdata <- rdata %>%  mutate(ExitHoles=recode(ExitHoles, '50+'='70')) #We will use 70 as a replacement for the category 50+ for statistical analyses purposes
rdata$StemBleeds <- as.integer(rdata$StemBleeds)
rdata$BarkCracks <- as.integer(rdata$BarkCracks)
rdata$ExitHoles <- as.integer(rdata$ExitHoles)

#Present or absent data for expert data
rdata <- rdata %>% mutate(Stembleeds=StemBleeds, Barkcracks=BarkCracks, Exitholes=ExitHoles)

rdata$Stembleeds <- rdata$Stembleeds>0
rdata$Barkcracks <- rdata$Barkcracks>0
rdata$Exitholes <- rdata$Exitholes>0

colnames(vdata)
colnames(rdata)


#remove irrelevant columns and merge

#From vdata remove column 2 (Date)
vdata <- vdata[,-2]



#Convert the rdata (expert) assessments into new columns
rdata <- rdata %>% rename(rStemBleeds=StemBleeds, rStembleeds=Stembleeds, rBarkCracks=BarkCracks,
                          rBarkcracks=Barkcracks, rExitHoles=ExitHoles, rExitholes=Exitholes)


#merge expert and volunteer data frames
all <- merge(vdata,rdata)

#1) Calculate the individual's sensitivity
#True positives/(True positives+False Negatives)

#True positives = volunteer TRUE and expert TRUE
#False Negatives = volunteer FALSE and expert TRUE

#Create a column TP and FN for each volunteer
#with TRUE (1) or FALSE (0) for TP and FN columns.
#Then sum for each volunteer and calculate sensivitity

# 2) Calculate the specificity (Ability to correctly identify absence of disease) for individuals
# Specificity = True negatives/(True negatives+False Positives)
#calculate the number of true negatives (TNs) and false positives (FPs)

#Calculate true positives
all <- all %>% mutate(TPbleed=Stembleeds=="TRUE"&rStembleeds=="TRUE",TPholes=rExitholes=="TRUE"&Exitholes=="TRUE",TPcracks=rBarkcracks=="TRUE"&Barkcracks=="TRUE")

#Convert to integer data
all$TPbleed <- as.integer(all$TPbleed)
all$TPholes <- as.integer(all$TPholes)
all$TPcracks <- as.integer(all$TPcracks)


#Calculate  false negatives for symptoms
all <- all %>% mutate(FNbleed=Stembleeds=="FALSE"&rStembleeds=="TRUE",FNholes=Exitholes=="FALSE"&rExitholes=="TRUE",FNcracks=Barkcracks=="FALSE"&rBarkcracks=="TRUE")

#convert to integer data
all$FNbleed <- as.integer(all$FNbleed)
all$FNholes <- as.integer(all$FNholes)
all$FNcracks <- as.integer(all$FNcracks)


#Calculate true negatives
all <- all %>% mutate(TNbleed=Stembleeds=="FALSE"&rStembleeds=="FALSE",TNholes=Exitholes=="FALSE"&rExitholes=="FALSE",TNcracks=Barkcracks=="FALSE"&rBarkcracks=="FALSE")

#convert to integer data
all$TNbleed <- as.integer(all$TNbleed)
all$TNholes <- as.integer(all$TNholes)
all$TNcracks <- as.integer(all$TNcracks)

#Calculate false positives
all <- all %>% mutate(FPbleed=Stembleeds=="TRUE"&rStembleeds=="FALSE",FPholes=Exitholes=="TRUE"&rExitholes=="FALSE",FPcracks=Barkcracks=="TRUE"&rBarkcracks=="FALSE")

#convert to integer data
all$FPbleed <- as.integer(all$FPbleed)
all$FPholes <- as.integer(all$FPholes)
all$FPcracks <- as.integer(all$FPcracks)


#Total sample number for sensitivity (true positive trees) and specificity (true negative trees) to calculate confidence intervals

#True positive trees
all <- all %>% group_by(SurveyorID) %>% mutate(nposbleed=sum(rStembleeds),nposcrack=sum(rBarkcracks),nposholes=sum(rExitholes))


#True negative trees
all <- all %>% group_by(SurveyorID) %>% mutate(nnegbleed=length(rStembleeds)-sum(rStembleeds),nnegcrack=length(rBarkcracks)-sum(rBarkcracks),nnegholes=length(rExitholes)-sum(rExitholes))


#Summing the true positives, true negatives, false positives and false negatives for each volunteer
allbleedsite <- all %>% group_by(SurveyorID) %>% mutate(vTPbleed=sum(TPbleed),vTPholes=sum(TPholes),vTPcracks=sum(TPcracks),vFNbleed=sum(FNbleed),vFNholes=sum(FNholes),vFNcracks=sum(FNcracks),
                                                        vTNbleed=sum(TNbleed),vTNholes=sum(TNholes),vTNcracks=sum(TNcracks),
                                                        vFPbleed=sum(FPbleed),vFPholes=sum(FPholes),vFPcracks=sum(FPcracks))

# 1) Calculate the sensitivity (Ability to correctly identify presence of disease) for individuals
# Sensitivity = True positives/(True positives+False Negatives)
allbleedsite <- allbleedsite %>% mutate(BleedSensitivity=vTPbleed/(vTPbleed+vFNbleed),
                                        CrackSensitivity=vTPcracks/(vTPcracks+vFNcracks),
                                        HoleSensitivity=vTPholes/(vTPholes+vFNholes))

# 2) Calculate the specificity (Ability to correctly identify absence of disease) for individuals
# Specificity = True negatives/(True negatives+False Positives)
#calculate the number of true negatives (TNs) and false positives (FPs)
allbleedsite <- allbleedsite %>% mutate(BleedSpecificity=vTNbleed/(vTNbleed+vFPbleed),
                                        CrackSpecificity=vTNcracks/(vTNcracks+vFPcracks),
                                        HoleSpecificity=vTNholes/(vTNholes+vFPholes))
colnames(allbleedsite)

#Now just get the sensitvity and specificity data we need for each volunteer to plot the sensitivity and specificity of each
#symptom for each individual, with 95% CIs
allsite <- allbleedsite[,c(3,28:51)]

##########Calculating the confidence intervals of the proportions

dataerror <- allsite

dataerror$Hseuci <- NA
dataerror$Hselci <- NA
dataerror$Hspuci <- NA
dataerror$Hsplci <- NA

for(i in 1:nrow(dataerror)){
  dataerror$Hseuci[i] <- BinomCI(dataerror$vTPholes[i],dataerror$nposholes[i],
                                 conf.level=0.95, method="clopper-pearson")[1,3]
  
  dataerror$Hselci[i] <- BinomCI(dataerror$vTPholes[i],dataerror$nposholes[i],
                                 conf.level=0.95, method="clopper-pearson")[1,2]
  
  dataerror$Hspuci[i] <- BinomCI(dataerror$vTNholes[i],dataerror$nnegholes[i],
                                 conf.level=0.95, method="clopper-pearson")[1,3]
  
  dataerror$Hsplci[i] <- BinomCI(dataerror$vTNholes[i],dataerror$nnegholes[i],
                                 conf.level=0.95, method="clopper-pearson")[1,2]}


dataerror$Bseuci <- NA
dataerror$Bselci <- NA
dataerror$Bspuci <- NA
dataerror$Bsplci <- NA

for(i in 1:nrow(dataerror)){
  dataerror$Bseuci[i] <- BinomCI(dataerror$vTPbleed[i],dataerror$nposbleed[i],
                                 conf.level=0.95, method="clopper-pearson")[1,3]
  
  dataerror$Bselci[i] <- BinomCI(dataerror$vTPbleed[i],dataerror$nposbleed[i],
                                 conf.level=0.95, method="clopper-pearson")[1,2]
  
  dataerror$Bspuci[i] <- BinomCI(dataerror$vTNbleed[i],dataerror$nnegbleed[i],
                                 conf.level=0.95, method="clopper-pearson")[1,3]
  
  dataerror$Bsplci[i] <- BinomCI(dataerror$vTNbleed[i],dataerror$nnegbleed[i],
                                 conf.level=0.95, method="clopper-pearson")[1,2]}


dataerror$Cseuci <- NA
dataerror$Cselci <- NA
dataerror$Cspuci <- NA
dataerror$Csplci <- NA

for(i in 1:nrow(dataerror)){
  dataerror$Cseuci[i] <- BinomCI(dataerror$vTPcracks[i],dataerror$nposcrack[i],
                                 conf.level=0.95, method="clopper-pearson")[1,3]
  
  dataerror$Cselci[i] <- BinomCI(dataerror$vTPcracks[i],dataerror$nposcrack[i],
                                 conf.level=0.95, method="clopper-pearson")[1,2]
  
  dataerror$Cspuci[i] <- BinomCI(dataerror$vTNcracks[i],dataerror$nnegcrack[i],
                                 conf.level=0.95, method="clopper-pearson")[1,3]
  
  dataerror$Csplci[i] <- BinomCI(dataerror$vTNcracks[i],dataerror$nnegcrack[i],
                                 conf.level=0.95, method="clopper-pearson")[1,2]}



###################Now plot data with the CIs


hseplot <- ggplot(dataerror, aes(x = SurveyorID, y = HoleSensitivity)) +
  geom_point() +                         # Points
  geom_errorbar(aes(ymin = Hselci, ymax = Hseuci), width = 0.2) +
  labs(title = "Exit holes",x = NULL, y = NULL)+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme(plot.title = element_text(hjust = 0.45), text = element_text(size = 10),
        axis.text.x = element_blank())

hseplot


hspplot <- ggplot(dataerror, aes(x = SurveyorID, y = HoleSpecificity)) +
  geom_point() +
  geom_errorbar(aes(ymin = Hsplci, ymax = Hspuci), width = 0.2) +
  labs(title = NULL,x = NULL, y = NULL)+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme(plot.title = element_text(hjust = 0.45), text = element_text(size = 10),
        axis.text.x = element_blank())


bseplot <- ggplot(dataerror, aes(x = SurveyorID, y = BleedSensitivity)) +
  geom_point() +                        
  geom_errorbar(aes(ymin = Bselci, ymax = Bseuci), width = 0.2) +  
  labs(title = "Stem bleeds",x = NULL, y = "Sensitivity")+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme(plot.title = element_text(hjust = 0.45), text = element_text(size = 10),
        axis.text.x = element_blank())





bspplot <- ggplot(dataerror, aes(x = SurveyorID, y = BleedSpecificity)) +
  geom_point() +
  geom_errorbar(aes(ymin = Bsplci, ymax = Bspuci), width = 0.2) +
  labs(title = NULL,x = NULL, y = "Specificity")+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme(plot.title = element_text(hjust = 0.45), text = element_text(size = 10),
        axis.text.x = element_blank())


cseplot <- ggplot(dataerror, aes(x = SurveyorID, y = CrackSensitivity)) +
  geom_point() +                        
  geom_errorbar(aes(ymin = Cselci, ymax = Cseuci), width = 0.2) +  
  labs(title = "Bark cracks",x = NULL, y = NULL)+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme(plot.title = element_text(hjust = 0.45), text = element_text(size = 10),
        axis.text.x = element_blank())


cspplot <- ggplot(dataerror, aes(x = SurveyorID, y = CrackSpecificity)) +
  geom_point() +
  geom_errorbar(aes(ymin = Csplci, ymax = Cspuci), width = 0.2) + 
  labs(title = NULL,x = NULL, y = NULL)+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme(plot.title = element_text(hjust = 0.45), text = element_text(size = 10),
        axis.text.x = element_blank())

colnames(dataerror)


library(cowplot)

combined_plot <- plot_grid(hseplot, hspplot, bseplot, bspplot, cseplot, cspplot, ncol = 3, nrow=2)
combined_plot <- plot_grid(bseplot, hseplot, cseplot, bspplot, hspplot, cspplot, ncol = 3, nrow=2, align = 'v')
combined_plot

x.grob <- textGrob("Surveyors", 
                   gp=gpar(fontsize=12))


combined_plot <- grid.arrange(arrangeGrob(combined_plot, bottom = x.grob))
combined_plot



#ROC Plot
dataerror$HFP <- 1-dataerror$HoleSpecificity
dataerror$HFPlci <- 1-dataerror$Hsplci
dataerror$HFPuci <- 1-dataerror$Hspuci

dataerror$BFP <- 1-dataerror$BleedSpecificity
dataerror$BFPlci <- 1-dataerror$Bsplci
dataerror$BFPuci <- 1-dataerror$Bspuci

dataerror$CFP <- 1-dataerror$CrackSpecificity
dataerror$CFPlci <- 1-dataerror$Csplci
dataerror$CFPuci <- 1-dataerror$Cspuci


#Plots

hroc <- ggplot(dataerror, aes(x = HFP, y = HoleSensitivity)) +
  geom_point(size = 0.5) +                         
  geom_errorbar(aes(ymin = Hselci, ymax = Hseuci), linewidth = 0.1, width = 0.2, alpha=0.3) +
  geom_errorbar(aes(xmin = HFPlci, xmax = HFPuci), linewidth = 0.1, width = 0.2, alpha=0.3) +
  geom_abline(intercept=0, slope=1, linetype=2)+
  labs(title = "Exit holes", x = "", y = "")+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), minor_breaks = NULL)+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), minor_breaks = NULL)+ theme(plot.title = element_text(hjust = 0.5),
                                                                                            text = element_text(size = 8))



hroc

broc <- ggplot(dataerror, aes(x = BFP, y = BleedSensitivity)) +
  geom_point(size = 0.5) +                         
  geom_errorbar(aes(ymin = Bselci, ymax = Bseuci), linewidth = 0.1, width = 0.2, alpha=0.3) +
  geom_errorbar(aes(xmin = BFPlci, xmax = BFPuci), linewidth = 0.1, width = 0.2, alpha=0.3) +
  geom_abline(intercept=0, slope=1, linetype=2)+
  labs(title = "Stem bleeds", x = "", y = "")+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), minor_breaks = NULL)+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), minor_breaks = NULL)+ theme(plot.title = element_text(hjust = 0.5),
                                                                                            text = element_text(size = 8))



broc

croc <- ggplot(dataerror, aes(x = CFP, y = CrackSensitivity)) +
  geom_point(size = 0.5) +                         
  geom_errorbar(aes(ymin = Cselci, ymax = Cseuci), linewidth = 0.1, width = 0.2, alpha=0.3) +
  geom_errorbar(aes(xmin = CFPlci, xmax = CFPuci), linewidth = 0.1, width = 0.2, alpha=0.3) +
  geom_abline(intercept=0, slope=1, linetype=2)+
  labs(title = "Bark cracks", x = "", y = "")+
  theme_minimal()+
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), minor_breaks = NULL)+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), minor_breaks = NULL)+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 8))



croc

#True positive rate on the y axes
#False positive rate on the x axes

combined_plot <- plot_grid(broc, hroc, croc, ncol = 3, nrow=1)
combined_plot

x.grob <- textGrob("False positive rate", 
                   gp=gpar(fontsize=12), vjust = -0.5)

y.grob <- textGrob("True positive rate", 
                   gp=gpar(fontsize=12), rot=90)

combined_plot <- grid.arrange(arrangeGrob(combined_plot, bottom = x.grob, left = y.grob))
combined_plot

#################################Part 2: Statistical analyses################################################################################

#Remove duplicate values#
colnames(allbleedsite)
#Dataframe with Treenumber, surveyor ID and counts of true positives, true negatives, false positives, false negatives
allsitef <- allbleedsite[,c(2,3,16:27)]
colnames(allsitef)
#We will then need to remove the duplicate values
allsite <-allsitef[!duplicated(allsitef), ]
library(car)
library(DHARMa)

#Fitting generalised linear models with a binomial distribution and logit function
#to then perform type II likelihood ratio chi-square tests to assess variable significance given data visualisation


#FOR SPECIFICITY#


colnames(allsite)
Spdfbleed <- allsite[,c(2,9,12)]
Spdfbleed$type <- "bleed"
colnames(Spdfbleed)[2] <- "TN"
colnames(Spdfbleed)[3] <- "FP"

Spdfhole <- allsite[,c(2,10,13)]
Spdfhole$type <- "hole"
colnames(Spdfhole)[2] <- "TN"
colnames(Spdfhole)[3] <- "FP"

Spdfcrack <- allsite[,c(2,11,14)]
Spdfcrack$type <- "crack"
colnames(Spdfcrack)[2] <- "TN"
colnames(Spdfcrack)[3] <- "FP"

Spdf <- rbind(Spdfcrack, Spdfhole, Spdfbleed)
Spdf$type <- as.factor(Spdf$type)

#We need to remove data where there is not a TN or an FP (i.e. no observation on a given tree for a given symptom)!!!!
Spdf <- Spdf[!(Spdf$TN == 0 & Spdf$FP == 0),]

#Specificity model#
m1sp <- glm(cbind(TN,FP)~type*SurveyorID,family=binomial(link="logit"), data=Spdf)
m1sp

Anova(m1sp, test.statistic = "LR")

#> Anova(m1sp, test.statistic = "LR")
#Analysis of Deviance Table (Type II tests)

#Response: cbind(TN, FP)
#LR Chisq Df Pr(>Chisq)    
#type              297.87  2  < 2.2e-16 ***
#  SurveyorID        311.74 22  < 2.2e-16 ***
#  type:SurveyorID   192.19 44  < 2.2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Check model fit#
simulationOutput <- simulateResiduals(fittedModel = m1sp, plot = T)
testDispersion(simulationOutput)
#One, possible outlier, but looks fine


#######################Repeat for sensitivity including symptom frequency######################


#Remove duplicate values#
colnames(allbleedsite)
#Dataframe with Treenumber, surveyor ID and counts of true positives, true negatives, false positives, false negatives, and expert counts
allsitef <- allbleedsite[,c(2,3,10:12,16:27)]
colnames(allsitef)
#We will then need to remove any duplicate values
allsite <-allsitef[!duplicated(allsitef), ]
library(car)
library(DHARMa)

#Fitting generalised linear models with a binomial distribution and logit function
#to then perform type II likelihood ratio chi-square tests to assess variable significance based on visualisation of data

#SENSITIVITY#
colnames(allsite)
Sedfbleed <- allsite[,c(2,3,6,9)]
Sedfbleed$type <- "bleed"
colnames(Sedfbleed)[3] <- "TP"
colnames(Sedfbleed)[4] <- "FN"
colnames(Sedfbleed)[2] <- "Expert_Count"

Sedfhole <- allsite[,c(2,5,7,10)]
Sedfhole$type <- "hole"
colnames(Sedfhole)[3] <- "TP"
colnames(Sedfhole)[4] <- "FN"
colnames(Sedfhole)[2] <- "Expert_Count"

Sedfcrack <- allsite[,c(2,4,8,11)]
Sedfcrack$type <- "crack"
colnames(Sedfcrack)[3] <- "TP"
colnames(Sedfcrack)[4] <- "FN"
colnames(Sedfcrack)[2] <- "Expert_Count"

Sedf <- rbind(Sedfcrack, Sedfhole, Sedfbleed)
Sedf$type <- as.factor(Sedf$type)

#We need to remove rows where there is not a TP or an FN (i.e. no observation for a given symptom)!!!!
Sedf <- Sedf[!(Sedf$TP == 0 & Sedf$FN == 0),]

#Sensitivity model#
#We cannot examine any interaction between the frequency of symptoms and the individual surevyor due to a comparably limited dataset#
sens_model <- glm(cbind(TP,FN)~type*Expert_Count+SurveyorID*type,family=binomial(link="logit"), data=Sedf, na.action = "na.fail")

summary(sens_model)
Anova(sens_model)

#> Anova(sens_model, test.statistic = "LR")
#Analysis of Deviance Table (Type II tests)

#Response: cbind(TP, FN)
#LR Chisq Df Pr(>Chisq)    
#type                      132.889  2  < 2.2e-16 ***
#  scale(Expert_Count)       145.003  1  < 2.2e-16 ***
 # SurveyorID                104.705 22  9.617e-13 ***
 # type:scale(Expert_Count)   49.803  2  1.533e-11 ***
 # type:SurveyorID           106.582 44  4.140e-07 ***
#  ---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Check model assumptions
simulationOutput <- simulateResiduals(fittedModel = m1se, plot = T)
testDispersion(simulationOutput)
#No issues#




############################################Part 3: The sensitivity by symptom number plots##############################################

rm(list=ls())
unloadNamespace("car")

#Read in the AOD data set and get exit hole sens and spec

#Read in the data collected by volunteers
vdata <- read.csv()

#Correct the types of data
vdata$Site <- as.factor(vdata$Site)
vdata$Date <- as.Date(vdata$Date, "%d/%m/%Y",  tz="GMT")
vdata$SurveyorID <- as.factor(vdata$SurveyorID)
vdata$TreeNumber <- as.factor(vdata$TreeNumber)

#Change 20+ to 20 for StemBleeds and BarkCracks#

vdata <- vdata %>% mutate(StemBleeds=recode(StemBleeds, '20+'='20')) %>% mutate(BarkCracks=recode(BarkCracks, '20+'='20'))

#For now, change the 10+, 20+, 50+ categorical volunteer observations to 10, 20, 51 (there is already a 50 observation)
vdata <- vdata %>%  mutate(ExitHoles=recode(ExitHoles, '10+'='10'))
vdata <- vdata %>%  mutate(ExitHoles=recode(ExitHoles, '50+'='51'))
vdata <- vdata %>%  mutate(ExitHoles=recode(ExitHoles, '20+'='20'))
vdata <- vdata %>%  mutate(ExitHoles=recode(ExitHoles, '5+'='10'))

vdata$StemBleeds <- as.integer(vdata$StemBleeds)
vdata$BarkCracks <- as.integer(vdata$BarkCracks)
vdata$ExitHoles <- as.integer(vdata$ExitHoles)

#present or absent data
vdata <- vdata %>% mutate(Stembleeds=StemBleeds, Barkcracks=BarkCracks, Exitholes=ExitHoles)

vdata$Stembleeds <- vdata$Stembleeds>0
vdata$Barkcracks <- vdata$Barkcracks>0
vdata$Exitholes <- vdata$Exitholes>0

#Read in the reference data collected by expert
rdata <- read.csv()
#Correct the types of data
rdata$Site <- as.factor(rdata$Site)
rdata$TreeNumber <- as.factor(rdata$TreeNumber)

#Creating a dataframe with just trees and symptom numbers to decide if#

#For now, change the 20+ to 22
rdata <- rdata %>%  mutate(BarkCracks=recode(BarkCracks, '20+'='22'))

#For exit holes change 20+ to 20, 50+ to 50; later put in separate categories#
rdata <- rdata %>%  mutate(ExitHoles=recode(ExitHoles, '20+'='22'))
rdata <- rdata %>%  mutate(ExitHoles=recode(ExitHoles, '50+'='24'))

rdata$StemBleeds <- as.integer(rdata$StemBleeds)
rdata$BarkCracks <- as.integer(rdata$BarkCracks)
rdata$ExitHoles <- as.integer(rdata$ExitHoles)

#present or absent data
rdata <- rdata %>% mutate(Stembleeds=StemBleeds, Barkcracks=BarkCracks, Exitholes=ExitHoles)

rdata$Stembleeds <- rdata$Stembleeds>0
rdata$Barkcracks <- rdata$Barkcracks>0
rdata$Exitholes <- rdata$Exitholes>0

colnames(vdata)
colnames(rdata)


#remove irrelevant columns and merge
#Only need, Site, SurveyorID, TreeNumber, StemBleeds and Stembleeds, BarkCracks and Barkcracks,
#ExitHoles and Exitholes,

#From vdata remove column 2 (Date)

vdata <- vdata[,-c(2)]

#Convert the rdata (expert) assessments into new columns

rdata <- rdata %>% rename(rStemBleeds=StemBleeds, rStembleeds=Stembleeds, rBarkCracks=BarkCracks,
                          rBarkcracks=Barkcracks, rExitHoles=ExitHoles, rExitholes=Exitholes)


#merge data frames
all <- merge(vdata,rdata)

#1) Calculate the individual's sensitivity
#True positives/(True positives+False Negatives)

#True positives = volunteer TRUE and r(expert) TRUE
#False Negatives = volunteer FALSE and expert TRUE
#Therefore create a column TP and FN for each volunteer
#with TRUE = 1, FALSE 0 for TP and FN. Then sum for each volunteer and calculate sensitivity

# 2) Calculate the specificity (Ability to correctly identify absence of disease) for individuals
# Specificity = True negatives/(True negatives+False Positives)
#calculate the number of true negatives (TNs) and false positives (FPs)

#Calculate true positives for symptoms
all <- all %>% mutate(TPbleed=Stembleeds=="TRUE"&rStembleeds=="TRUE",TPholes=rExitholes=="TRUE"&Exitholes=="TRUE",TPcracks=rBarkcracks=="TRUE"&Barkcracks=="TRUE")

#Convert to integer data
all$TPbleed <- as.integer(all$TPbleed)
all$TPholes <- as.integer(all$TPholes)
all$TPcracks <- as.integer(all$TPcracks)


#Calculate  false negatives for symptoms
all <- all %>% mutate(FNbleed=Stembleeds=="FALSE"&rStembleeds=="TRUE",FNholes=Exitholes=="FALSE"&rExitholes=="TRUE",FNcracks=Barkcracks=="FALSE"&rBarkcracks=="TRUE")

#convert to integer data
all$FNbleed <- as.integer(all$FNbleed)
all$FNholes <- as.integer(all$FNholes)
all$FNcracks <- as.integer(all$FNcracks)


#Calculate true negatives
all <- all %>% mutate(TNbleed=Stembleeds=="FALSE"&rStembleeds=="FALSE",TNholes=Exitholes=="FALSE"&rExitholes=="FALSE",TNcracks=Barkcracks=="FALSE"&rBarkcracks=="FALSE")

#convert to integer data
all$TNbleed <- as.integer(all$TNbleed)
all$TNholes <- as.integer(all$TNholes)
all$TNcracks <- as.integer(all$TNcracks)

#Calculate false positives
all <- all %>% mutate(FPbleed=Stembleeds=="TRUE"&rStembleeds=="FALSE",FPholes=Exitholes=="TRUE"&rExitholes=="FALSE",FPcracks=Barkcracks=="TRUE"&rBarkcracks=="FALSE")

#convert to integer data
all$FPbleed <- as.integer(all$FPbleed)
all$FPholes <- as.integer(all$FPholes)
all$FPcracks <- as.integer(all$FPcracks)


#All symptoms all sites#

#True positive trees assessed for each individual (for CIs)
all <- all %>% group_by(SurveyorID) %>% mutate(nposbleed=sum(rStembleeds),nposcrack=sum(rBarkcracks),nposholes=sum(rExitholes))


#True negative trees assessed for each individual (for CIs)
all <- all %>% group_by(SurveyorID) %>% mutate(nnegbleed=length(rStembleeds)-sum(rStembleeds),nnegcrack=length(rBarkcracks)-sum(rBarkcracks),nnegholes=length(rExitholes)-sum(rExitholes))

#Summing the total nimber of true positives, true negatives, false positives, false negatives for each individual
allbleedsite <- all %>% group_by(SurveyorID) %>% mutate(vTPbleed=sum(TPbleed),vTPholes=sum(TPholes),vTPcracks=sum(TPcracks),vFNbleed=sum(FNbleed),vFNholes=sum(FNholes),vFNcracks=sum(FNcracks),
                                                        vTNbleed=sum(TNbleed),vTNholes=sum(TNholes),vTNcracks=sum(TNcracks),
                                                        vFPbleed=sum(FPbleed),vFPholes=sum(FPholes),vFPcracks=sum(FPcracks))

# 1) Calculate the sensitivity for individuals

# Sensitivity = True positives/(True positives+False Negatives)

allbleedsite <- allbleedsite %>% mutate(BleedSensitivity=vTPbleed/(vTPbleed+vFNbleed),
                                        CrackSensitivity=vTPcracks/(vTPcracks+vFNcracks),
                                        HoleSensitivity=vTPholes/(vTPholes+vFNholes))

# 2) Calculate the specificity for individuals
# Specificity = True negatives/(True negatives+False Positives)

allbleedsite <- allbleedsite %>% mutate(BleedSpecificity=vTNbleed/(vTNbleed+vFPbleed),
                                        CrackSpecificity=vTNcracks/(vTNcracks+vFPcracks),
                                        HoleSpecificity=vTNholes/(vTNholes+vFPholes))
colnames(allbleedsite)

allsitef <- allbleedsite[,c(2,3,10:12,16:27)]
allsite <-allsitef[!duplicated(allsitef), ]

#Making the sensitivity dataframe#

#We need to remove data where there is not a TP or an FN!!!!

names(allsite)
Sedfbleed <- allsite[,c(2,3,6,9)]
Sedfbleed$type <- "bleed"
colnames(Sedfbleed)[3] <- "TP"
colnames(Sedfbleed)[4] <- "FN"

Sedfhole <- allsite[,c(2,5,7,10)]
Sedfhole$type <- "hole"
colnames(Sedfhole)[3] <- "TP"
colnames(Sedfhole)[4] <- "FN"

Sedfcrack <- allsite[,c(2,4,8,11)]
Sedfcrack$type <- "crack"
colnames(Sedfcrack)[3] <- "TP"
colnames(Sedfcrack)[4] <- "FN"

#Remove those rows that were not scored for a given symptom (i.e. no true positives or false negatives)
Sedfbleed <- Sedfbleed[!(Sedfbleed$TP == 0 & Sedfbleed$FN == 0),]
Sedfhole <- Sedfhole[!(Sedfhole$TP == 0 & Sedfhole$FN == 0),]
Sedfcrack <- Sedfcrack[!(Sedfcrack$TP == 0 & Sedfcrack$FN == 0),]

#Plot to see if there is any visual trend of sensitivity by symptom number, pooling individual data together#
Sedfbleed <- Sedfbleed %>% group_by(rStemBleeds) %>% mutate(STP = sum(TP))
Sedfbleed <- Sedfbleed %>% group_by(rStemBleeds) %>% mutate(SFN = sum(FN))
Sedfbleed <- Sedfbleed %>% group_by(rStemBleeds) %>% mutate(sens = STP/(STP+SFN))
Sedfbleed <- Sedfbleed %>% group_by(rStemBleeds) %>% mutate(npos = length(rStemBleeds))


#Making a nice plot, because of the lack of replicate trees we pool the data together for all individuals#

#First make confidence intervals#
Sedfbleed$uci <- NA
Sedfbleed$lci <- NA
library(DescTools)
for(i in 1:nrow(Sedfbleed)){
  Sedfbleed$uci[i] <- BinomCI(Sedfbleed$STP[i],Sedfbleed$npos[i],
                              conf.level=0.95, method="clopper-pearson")[1,3]
  
  Sedfbleed$lci[i] <- BinomCI(Sedfbleed$STP[i],Sedfbleed$npos[i],
                              conf.level=0.95, method="clopper-pearson")[1,2]}


bplot <- ggplot(Sedfbleed, aes(x = rStemBleeds, y = sens)) +
  geom_point() +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.2) + 
  labs(title = NULL,x = "Number of bleeds", y = "")+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  scale_x_continuous(breaks = c(1:14),minor_breaks = seq(1, 14, by = 1))+
  theme(plot.title = element_text(hjust = 0.45), text = element_text(size = 10),axis.text.x = element_text(angle =45, hjust = 1))
bplot



#Repeat for exit holes
Sedfhole <- Sedfhole %>% group_by(rExitHoles) %>% mutate(STP = sum(TP))
Sedfhole <- Sedfhole %>% group_by(rExitHoles) %>% mutate(SFN = sum(FN))
Sedfhole <- Sedfhole %>% group_by(rExitHoles) %>% mutate(sens = STP/(STP+SFN))
Sedfhole <- Sedfhole %>% group_by(rExitHoles) %>% mutate(npos = length(rExitHoles))

Sedfhole$uci <- NA
Sedfhole$lci <- NA
library(DescTools)
for(i in 1:nrow(Sedfhole)){
  Sedfhole$uci[i] <- BinomCI(Sedfhole$STP[i],Sedfhole$npos[i],
                             conf.level=0.95, method="clopper-pearson")[1,3]
  
  Sedfhole$lci[i] <- BinomCI(Sedfhole$STP[i],Sedfhole$npos[i],
                             conf.level=0.95, method="clopper-pearson")[1,2]}


hplot <- ggplot(Sedfhole, aes(x = rExitHoles, y = sens)) +
  geom_point() +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.2) +
  labs(title = NULL,x = "Number of exit holes", y = "")+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  scale_x_continuous(breaks = c(1,5,10,15,20,22,24),minor_breaks = seq(1,20, by = 1), labels = c(1,5,10,15,20, "20+", "50+"))+
  theme(plot.title = element_text(hjust = 0.45), text = element_text(size = 10)
        ,panel.grid.minor.x = element_blank(),axis.text.x = element_text(angle =75, hjust = 1))
hplot


########Now bark cracks

Sedfcrack <- Sedfcrack %>% group_by(rBarkCracks) %>% mutate(STP = sum(TP))
Sedfcrack <- Sedfcrack %>% group_by(rBarkCracks) %>% mutate(SFN = sum(FN))
Sedfcrack <- Sedfcrack %>% group_by(rBarkCracks) %>% mutate(sens = STP/(STP+SFN))
Sedfcrack <- Sedfcrack %>% group_by(rBarkCracks) %>% mutate(npos = length(rBarkCracks))

Sedfcrack$uci <- NA
Sedfcrack$lci <- NA
library(DescTools)
for(i in 1:nrow(Sedfcrack)){
  Sedfcrack$uci[i] <- BinomCI(Sedfcrack$STP[i],Sedfcrack$npos[i],
                              conf.level=0.95, method="clopper-pearson")[1,3]
  
  Sedfcrack$lci[i] <- BinomCI(Sedfcrack$STP[i],Sedfcrack$npos[i],
                              conf.level=0.95, method="clopper-pearson")[1,2]}


cplot <- ggplot(Sedfcrack, aes(x = rBarkCracks, y = sens)) +
  geom_point() +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.2) +
  labs(title = NULL,x = "Number of cracks", y = "")+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  scale_x_continuous(breaks = c(1,5,10,15,20,22),minor_breaks = seq(1, 23, by = 1), labels = c(1,5,10,15,20, "21+"))+
  theme(plot.title = element_text(hjust = 0.45), text = element_text(size = 10),axis.text.x = element_text(angle =45, hjust = 1))
cplot



library(cowplot)
combined_plot <- plot_grid(bplot,
                           hplot, cplot,ncol = 3, nrow=1)
combined_plot

y.grob <- textGrob("Sensitivity", 
                   gp=gpar(fontsize=12), rot=90)


combined_plot <- grid.arrange(arrangeGrob(combined_plot, left = y.grob))

combined_plot


