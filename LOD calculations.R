########################################################
#Programmer: Wenqiao He
#Last update: May 2024
#Purpose: P. ovale real-time PCR assay LOD calculations
#        (adapted from Jeff Laux and Jonathan Parr)
########################################################

#Loading packages for use
library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)
# Set working directory
setwd("FILE_PATH")

#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------

#loading data 
LOD_data <- read_excel("FILE_NAME") ##Input data file name

##Singleplex P.ovale curtisi assay LOD calculation#############################################################################
Singleplex_Poc_LOD_data <- subset(LOD_data,LOD_data$Assay=="Singleplex" & LOD_data$Species=="P. ovale curtisi")
Singleplex_Poc_LOD_data <- Singleplex_Poc_LOD_data |> dplyr::mutate(Diagnosis = dplyr::case_when(`Cq Value` < 40 ~ "Positive",
                                                           .default = "Negative"))
Singleplex_Poc_LOD_data <- Singleplex_Poc_LOD_data |> dplyr::mutate(Detection = dplyr::case_when(`Cq Value` < 40 ~ 1,
                                                                                                 .default = 0))

mod.1_Poc = glm(Detection~Concentration, data=Singleplex_Poc_LOD_data, family=binomial(link="probit"))

summary(mod.1_Poc)

x.seq    = seq(from=0, to=10000, by=0.1)
preds_c    = predict(mod.1_Poc, newdata=data.frame(Concentration=x.seq), 
                     type="link", se.fit=TRUE)
preds_c.df = data.frame(Concentration=x.seq, eta=preds_c[[1]], se=preds_c[[2]])

preds_c.df$pred.prob = pnorm(preds_c.df$eta)
preds_c.df$lower.lim = pnorm(preds_c.df$eta - 1.96*preds_c.df$se)
preds_c.df$upper.lim = pnorm(preds_c.df$eta + 1.96*preds_c.df$se)
View(preds_c.df)

preds_c.df$Concentration[min(which(preds_c.df$upper.lim>.95))]
preds_c.df$Concentration[max(which(preds_c.df$lower.lim<.95))]
preds_c.df$Concentration[min(which(preds_c.df$pred.prob>.95))]

### Plotting Singleplex P.ovale curtisi assay LOD figure#######################################################################
windows() 
plot(Detection~Concentration, data=Singleplex_Poc_LOD_data,cex.axis=1.25,xlim=c(0,10),cex.lab=1.5,
     xlab = 'Concentration (Parasites/ul)', 
     ylab = 'Probability of detection',
     family = 'serif',font.lab=2)  
lines(x.seq, preds_c.df$pred.prob, col="blue",lwd=5) 
abline(h=.95, col="gray",lwd=2)
lines(x.seq, preds_c.df$lower.lim, col="lightblue",lwd=5) 
lines(x.seq, preds_c.df$upper.lim, col="lightblue",lwd=5)    

##Singleplex P.ovale wallikeri assay LOD calculation###########################################################################
Singleplex_Pow_LOD_data <- subset(LOD_data,LOD_data$Assay=="Singleplex"& LOD_data$Species=="P. ovale wallikeri")
Singleplex_Pow_LOD_data <- Singleplex_Pow_LOD_data |> dplyr::mutate(Diagnosis = dplyr::case_when(`Cq Value` < 45 ~ "Positive",
                                                                                                 .default = "Negative"))
Singleplex_Pow_LOD_data <- Singleplex_Pow_LOD_data |> dplyr::mutate(Detection = dplyr::case_when(`Cq Value` < 45 ~ 1,
                                                                                                 .default = 0))

mod.1_Pow = glm(Detection~Concentration, data=Singleplex_Pow_LOD_data, family=binomial(link="probit"))

summary(mod.1_Pow)

x.seq    = seq(from=0, to=10000, by=0.1)
preds_w    = predict(mod.1_Pow, newdata=data.frame(Concentration=x.seq), 
                        type="link", se.fit=TRUE)
preds_w.df = data.frame(concentration=x.seq, eta=preds_w[[1]], se=preds_w[[2]])

preds_w.df$pred.prob = pnorm(preds_w.df$eta)
preds_w.df$lower.lim = pnorm(preds_w.df$eta - 1.96*preds_w.df$se)
preds_w.df$upper.lim = pnorm(preds_w.df$eta + 1.96*preds_w.df$se)
View(preds_w.df)

preds_w.df$concentration[min(which(preds_w.df$upper.lim>.95))]
preds_w.df$concentration[max(which(preds_w.df$lower.lim<.95))]
preds_w.df$concentration[max(which(preds_w.df$pred.prob<.95))]

### Plotting Singleplex P.ovale wallikeri assay LOD figure#####################################################################
windows()
plot(Detection~Concentration, data=Singleplex_Pow_LOD_data,cex.axis=1.25,xlim=c(0,40),cex.lab=1.5,
     xlab = 'Concentration (Parasites/ul)', ylab = 'Probability of detection',family = 'serif',font.lab=2) 
lines(x.seq, preds_w.df$pred.prob, col="#fc8d59",lwd=5) 
abline(h=.95, col="gray",lwd=2)
lines(x.seq, preds_w.df$lower.lim, col="navajowhite",lwd=5) 
lines(x.seq, preds_w.df$upper.lim, col="navajowhite",lwd=5) 

##Duplex P.ovale curtisi assay LOD calculation#################################################################################
Duplex_Poc_LOD_data <- subset(LOD_data,LOD_data$Assay=="Duplex" & LOD_data$Species=="P. ovale curtisi")
Duplex_Poc_LOD_data <- Duplex_Poc_LOD_data |> dplyr::mutate(Diagnosis = dplyr::case_when(`Cq Value` < 45 ~ "Positive",
                                                                                                 .default = "Negative"))
Duplex_Poc_LOD_data <- Duplex_Poc_LOD_data |> dplyr::mutate(Detection = dplyr::case_when(`Cq Value` < 45 ~ 1,
                                                                                                 .default = 0))

mod.2_Poc = glm(Detection~Concentration, data=Duplex_Poc_LOD_data, family=binomial(link="probit"))

summary(mod.2_Poc)

x.seq    = seq(from=0, to=10000, by=0.1)
preds_dc    = predict(mod.2_Poc, newdata=data.frame(Concentration=x.seq), 
                     type="link", se.fit=TRUE)
preds_dc.df = data.frame(Concentration=x.seq, eta=preds_dc[[1]], se=preds_dc[[2]])

preds_dc.df$pred.prob = pnorm(preds_dc.df$eta)
preds_dc.df$lower.lim = pnorm(preds_dc.df$eta - 1.96*preds_dc.df$se)
preds_dc.df$upper.lim = pnorm(preds_dc.df$eta + 1.96*preds_dc.df$se)
View(preds_dc.df)

preds_dc.df$Concentration[min(which(preds_dc.df$upper.lim>.95))]
preds_dc.df$Concentration[max(which(preds_dc.df$lower.lim<.95))]
preds_dc.df$Concentration[min(which(preds_dc.df$pred.prob>.95))]

### Plotting Duplex P.ovale curtisi assay LOD figure###########################################################################
windows() 
plot(Detection~Concentration, data=Duplex_Poc_LOD_data,cex.axis=1.25,xlim=c(0,10),cex.lab=1.5,
     xlab = 'Concentration (Parasites/ul)', 
     ylab = 'Probability of detection',
     family = 'serif',font.lab=2)  
lines(x.seq, preds_dc.df$pred.prob, col="blue",lwd=5) 
abline(h=.95, col="gray",lwd=2)
lines(x.seq, preds_dc.df$lower.lim, col="lightblue",lwd=5) 
lines(x.seq, preds_dc.df$upper.lim, col="lightblue",lwd=5)  

##Duplex P.ovale wallikeri assay LOD calculation###############################################################################
Duplex_Pow_LOD_data <- subset(LOD_data,LOD_data$Assay=="Duplex"& LOD_data$Species=="P. ovale wallikeri")
Duplex_Pow_LOD_data <- Duplex_Pow_LOD_data |> dplyr::mutate(Diagnosis = dplyr::case_when(`Cq Value` < 45 ~ "Positive",
                                                                                                 .default = "Negative"))
Duplex_Pow_LOD_data <- Duplex_Pow_LOD_data |> dplyr::mutate(Detection = dplyr::case_when(`Cq Value` <45 ~ 1,
                                                                                                 .default = 0))

mod.2_Pow = glm(Detection~Concentration, data=Duplex_Pow_LOD_data, family=binomial(link="probit"))

summary(mod.2_Pow)

x.seq    = seq(from=0, to=10000, by=0.1)
preds_dw    = predict(mod.2_Pow, newdata=data.frame(Concentration=x.seq), 
                     type="link", se.fit=TRUE)
preds_dw.df = data.frame(concentration=x.seq, eta=preds_dw[[1]], se=preds_dw[[2]])

preds_dw.df$pred.prob = pnorm(preds_dw.df$eta)
preds_dw.df$lower.lim = pnorm(preds_dw.df$eta - 1.96*preds_dw.df$se)
preds_dw.df$upper.lim = pnorm(preds_dw.df$eta + 1.96*preds_dw.df$se)
View(preds_w.df)

preds_dw.df$concentration[min(which(preds_dw.df$upper.lim>.95))]
preds_dw.df$concentration[max(which(preds_dw.df$lower.lim<.95))]
preds_dw.df$concentration[max(which(preds_dw.df$pred.prob<.95))]

### Plotting Duplex P.ovale wallikeri assay LOD figure#########################################################################
windows()
plot(Detection~Concentration, data=Duplex_Pow_LOD_data,cex.axis=1.25,xlim=c(0,100),cex.lab=1.5,
     xlab = 'Concentration (Parasites/ul)', ylab = 'Probability of detection',family = 'serif',font.lab=2) 
lines(x.seq, preds_dw.df$pred.prob, col="#fc8d59",lwd=5) 
abline(h=.95, col="gray",lwd=2)
lines(x.seq, preds_dw.df$lower.lim, col="navajowhite",lwd=5) 
lines(x.seq, preds_dw.df$upper.lim, col="navajowhite",lwd=5) 
