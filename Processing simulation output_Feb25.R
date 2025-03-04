
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

######################Looking at the combined distributions############################################

#Load in rds workflow output file from the workflow of interest
blah <- readRDS()
blahdf <- as.data.frame(blah)

#Sample randomly from all the possible 1,220,000 distributions from the combined 25 mcmc chains from the workflow
#Note: if convergence did not happen on the first run, this chain will be discarded, chains will only be saved
#from the first chain that converged. Each chain has 48,800 value of alpha and beta parameters of the beta distributions
#of sensitivity and specificity of symptom one. Will need to check psrf values from the output summary csv file and 
#remove chains which did not converge

randomse <- c()
for(i in 1:nrow(blahdf)){
  randomse[i] <- rbeta(1,blahdf$alphasep1[i], blahdf$betasep1[i])
}

randomsp <- c()
for(i in 1:nrow(blahdf)){
  randomsp[i] <- rbeta(1,blahdf$alphaspp1[i], blahdf$betaspp1[i])
}


#Get the actual sensitivity and specificity data from the summary output csv file from the workflow#
actblah <- read.csv()#
actse <- data.frame(actblah$actualsevol1)

n <- 0
actse <- cbind(actse,n)

actsp <- data.frame(actblah$actualspvol1)
actsp <- cbind(actsp,n)

names(actse)[1] <- "actualse"
names(actsp)[1] <- "actualsp"

#Kernel density estimates from the workflow estimated sensitivity values
epdfse <- density(randomse, from = 0, to = 1)
epdfse_df <- data.frame(x = epdfse$x, y = epdfse$y)

#Make the density plot for the sensitivity output#
sensitivity_density_output_plot <- ggplot() +
  geom_ribbon(data=epdfse_df,aes(x = x, ymin = 0, ymax = y), fill = "blue", alpha = 0.3)+
  geom_point(data = actse, aes(actualse, n), size=8,alpha=0.8, shape=4, color = "red")+
  geom_vline(xintercept=0.677,colour="red",alpha=1)+ #0.677 is the median value from the distribution data was simulated from
  geom_vline(xintercept=0.47,colour="red", linetype = "dashed",alpha=1)+ #0.47 is the lower 5 value from the distribution data was simulated from
  geom_vline(xintercept=0.845,colour="red", linetype = "dashed",alpha=1)+ #0.845 is the upper 95 value from the distribution data was simulated from
  labs(title = "(a)",
       x = "", y = "") +
  theme_minimal()+
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.1))+
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5,1))+
  theme(plot.title = element_text(hjust = 0), text = element_text(size = 20))

sensitivity_density_output_plot

#Kernel density estimates from the workflow estimated sensitivity values
epdfsp <- density(randomsp, from = 0, to = 1)
epdfsp_df <- data.frame(x = epdfsp$x, y = epdfsp$y)

specificity_density_output_plot <- ggplot(epdfsp_df, aes(x, y)) +
  geom_ribbon(data=epdfsp_df,aes(x = x, ymin = 0, ymax = y), fill = "blue", alpha = 0.3)+
  geom_point(data = actsp, aes(actualsp, n), size=8,alpha=0.8, shape=4, color = "red")+
  geom_vline(xintercept=0.884,colour="red",alpha=1)+
  geom_vline(xintercept=0.695,colour="red", linetype = "dashed",alpha=1)+
  geom_vline(xintercept=0.977,colour="red", linetype = "dashed",alpha=1)+
  labs(title = "(b)",
    x = "", y = "") +
  theme_minimal()+
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.1))+
  scale_y_continuous(limits = c(0, 5.5), breaks = seq(0, 5.5, 1))+
  theme(text = element_text(size = 20))

specificity_density_output_plot

saveRDS(randomse,)#Save the output of the workflow estimated sensitivity and specificity values
saveRDS(randomsp,)#Save the output of the workflow estimated sensitivity and specificity values

#######Statistics#########

#Sensitivity#
#Median from simulated surveyor data = 
median(nocovactse$actualse)
#Median from workflow output = 
median(nocovrandomse)

#min from simulated surveyor data = 
min(nocovactse$actualse)
#0.05 from workflow output = 
quantile(nocovrandomse, probs = c(0.05))[[1]]

#max from simulated surveyor data = 
max(nocovactse$actualse)
#0.95 from workflow output = 
quantile(nocovrandomse, probs = c(0.95))[[1]]

#Specificity#
#Median from simulated surveyor data = 
median(nocovactsp$actualsp)
#Median from workflow output = 
median(nocovrandomsp)

#min from simulated surveyor data = 
min(nocovactsp$actualsp)
#0.05 from workflow output = 
quantile(nocovrandomsp, probs = c(0.05))[[1]]

#max from simulated surveyor data = 
max(nocovactsp$actualsp)
#0.95 from workflow output = 
quantile(nocovrandomsp, probs = c(0.95))[[1]]


#############################################Individuals' estimates###############################################################


rm(list=ls())
#Load in workflow estimated sensitivity and specificity values
randomse <- readRDS()
randomsp <- readRDS()


#Divide this up for each of the 25 surveyors: 48800 values per surveyor#
#Adjust as appropriate if have less surveyors due to model not converging#
randomsedf <- data.frame(randomse)
randomspdf <- data.frame(randomsp)

randomsedf <- randomsedf %>%
  mutate(label = case_when(
    row_number() <= 48800 ~ 1,
    row_number() <= (48800*2) ~ 2,
    row_number() <= (48800*3) ~ 3,
    row_number() <= (48800*4) ~ 4,
    row_number() <= (48800*5) ~ 5,
    row_number() <= (48800*6) ~ 6,
    row_number() <= (48800*7) ~ 7,
    row_number() <= (48800*8) ~ 8,
    row_number() <= (48800*9) ~ 9,
    row_number() <= (48800*10) ~ 10,
    row_number() <= (48800*11) ~ 11,
    row_number() <= (48800*12) ~ 12,
    row_number() <= (48800*13) ~ 13,
    row_number() <= (48800*14) ~ 14,
    row_number() <= (48800*15) ~ 15,
    row_number() <= (48800*16) ~ 16,
    row_number() <= (48800*17) ~ 17,
    row_number() <= (48800*18) ~ 18,
    row_number() <= (48800*19) ~ 19,
    row_number() <= (48800*20) ~ 20,
    row_number() <= (48800*21) ~ 21,
    row_number() <= (48800*22) ~ 22,
    row_number() <= (48800*23) ~ 23,
    row_number() <= (48800*24) ~ 24,
    row_number() <= (48800*25) ~ 25,
    TRUE ~ NA_real_
  ))

#Obtain median, upper 95, lower 5 estimates for each surveyor
randomsedf <- randomsedf %>% group_by(label) %>% mutate(vmed=quantile(randomse, probs = c(0.5))[[1]],
                                                        lmed=quantile(randomse, probs = c(0.05))[[1]],
                                                        umed=quantile(randomse, probs = c(0.95))[[1]])

#Only keep median, upper 95, lower 5 estimate data for each surveyor
randomsedf <- randomsedf[,c(2,3,4,5)]
randomsedf <-randomsedf[!duplicated(randomsedf), ]

#Get the actual sensitivity data from the summary output csv file from the workflow#
actblah <- read.csv()
actse <- data.frame(actblah$actualsevol1)
colnames(actse)[1] <- "actsens"

#Combine the sensitivity estimates from the workflow with the actual data
randomsedf <- cbind(randomsedf,actse)

write.csv(randomsedf,) #save data as csv file

#Plotting actual surveyor sensitivity values on the x axis, and the workflow estimated
#value on the y axis, with error bars representing upper 95 and lower 5 estimate from the workflow
sensitivity_individual_output_plot <- ggplot(randomsedf, aes(x = actsens, y = vmed)) +
  geom_point(colour="blue") +
  geom_errorbar(aes(ymin = lmed, ymax = umed),colour="blue", width = 0.05, alpha = 0.2) +
  geom_abline(intercept=0, slope=1, linetype=2, colour="red")+
  labs(title = "(a)",x = "Sensitivity", y = "")+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2))+
  theme(text = element_text(size = 20))
sensitivity_individual_output_plot
mse(randomsedf$actsens,randomsedf$vmed) #Obtain mean squared error value#


#Repeat process for specificity
#Divide workflow estimates up for each of the 25 surveyors: 48800 values per surveyor#
#Adjust as appropriate if have less surveyors due to model not converging#
randomspdf <- randomspdf %>%
  mutate(label = case_when(
    row_number() <= 48800 ~ 1,
    row_number() <= (48800*2) ~ 2,
    row_number() <= (48800*3) ~ 3,
    row_number() <= (48800*4) ~ 4,
    row_number() <= (48800*5) ~ 5,
    row_number() <= (48800*6) ~ 6,
    row_number() <= (48800*7) ~ 7,
    row_number() <= (48800*8) ~ 8,
    row_number() <= (48800*9) ~ 9,
    row_number() <= (48800*10) ~ 10,
    row_number() <= (48800*11) ~ 11,
    row_number() <= (48800*12) ~ 12,
    row_number() <= (48800*13) ~ 13,
    row_number() <= (48800*14) ~ 14,
    row_number() <= (48800*15) ~ 15,
    row_number() <= (48800*16) ~ 16,
    row_number() <= (48800*17) ~ 17,
    row_number() <= (48800*18) ~ 18,
    row_number() <= (48800*19) ~ 19,
    row_number() <= (48800*20) ~ 20,
    row_number() <= (48800*21) ~ 21,
    row_number() <= (48800*22) ~ 22,
    row_number() <= (48800*23) ~ 23,
    row_number() <= (48800*24) ~ 24,
    row_number() <= (48800*25) ~ 25,
    TRUE ~ NA_real_
  ))

#Obtain median, upper 95, lower 5 estimates for each surveyor
randomspdf <- randomspdf %>% group_by(label) %>% mutate(vmed=quantile(randomsp, probs = c(0.5))[[1]],
                                                        lmed=quantile(randomsp, probs = c(0.05))[[1]],
                                                        umed=quantile(randomsp, probs = c(0.95))[[1]])

#Only keep median, upper 95, lower 5 estimate data for each surveyor
randomspdf <- randomspdf[,c(2,3,4,5)]
randomspdf <-randomspdf[!duplicated(randomspdf), ]

#Get the actual specificity data from the summary output csv file from the workflow#
actsp <- data.frame(actblah$actualspvol1)
colnames(actsp)[1] <- "actspec"

#Combine the specificity estimates from the workflow with the actual data
randomspdf <- cbind(randomspdf,actsp)

write.csv(randomspdf,) ##save data as csv file


specificity_individual_output_plot <- ggplot(randomspdf, aes(x = actspec, y = vmed)) +
  geom_point(colour="blue") +
  #geom_line()+
  #geom_hline(yintercept=0.8876754, linetype="dashed",colour="red")+
  #geom_hline(yintercept=0.9914714, linetype="dashed",colour="red")+
  geom_errorbar(aes(ymin = lmed, ymax = umed),colour="blue", width = 0.05, alpha = 0.2) +
  geom_abline(intercept=0, slope=1, linetype=2, colour="red")+
  labs(title = "(b)",x = "Specificity", y = "")+
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  #scale_x_continuous(limits = c(0, 13.9999), breaks = seq(0, 13.9999, 1))+
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2))+
  theme(text = element_text(size = 20))
specificity_individual_output_plot
mse(randomspdf$actspec,randomspdf$vmed) #Obtain mean squared error value#

