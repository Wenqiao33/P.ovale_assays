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
     xlab = 'Concentration (Parasites/µl)', 
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
     xlab = 'Concentration (Parasites/µl)', ylab = 'Probability of detection',family = 'serif',font.lab=2) 
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
     xlab = 'Concentration (Parasites/µl)', 
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
     xlab = 'Concentration (Parasites/µl)', ylab = 'Probability of detection',family = 'serif',font.lab=2) 
lines(x.seq, preds_dw.df$pred.prob, col="#fc8d59",lwd=5) 
abline(h=.95, col="gray",lwd=2)
lines(x.seq, preds_dw.df$lower.lim, col="navajowhite",lwd=5) 
lines(x.seq, preds_dw.df$upper.lim, col="navajowhite",lwd=5) 

##Plotting Duplex P.ovale curtisi and P.ovale wallikeri assay LOD#############################################################
Duplex_Poc_Pow_LOD_data <- subset(LOD_data,LOD_data$Assay=="Duplex")
Duplex_Poc_Pow_LOD_data <- Duplex_Poc_Pow_LOD_data |> dplyr::mutate(Diagnosis = dplyr::case_when(`Cq Value` < 45 ~ "Positive",
                                                                                                 .default = "Negative"))
LOD <-ggplot(Duplex_Poc_Pow_LOD_data, aes(x=Species, y=Concentration, colour=Diagnosis))+
  geom_point(position=position_jitter(width = 0.23,height = 0),size=3.5, alpha=0.5)+
  guides(colour=guide_legend(title=" "))+
  scale_fill_brewer()+
  labs(title=" ")+
  ylab("Parasite density (parasites/µl)")+
  xlab("Species")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(text=element_text(family="serif"),axis.title=element_text(size=17,face="bold"))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5,size = 17, face="italic"),legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 17),axis.title.y = element_text(size = 17))+
  geom_segment(aes(y=4.2,x=0.7,yend=4.2,xend=1.3),linetype = "dashed",size=0.8)+
  geom_segment(aes(y=41.2,x=1.7,yend=41.2,xend=2.3),linetype = "dashed",size=0.8)+
  scale_y_log10()+
  scale_x_discrete(labels=c("P. ovale curtisi","P. ovale wallikeri"))+
  geom_text(aes(y=2.60,x=1.46, label='P. ovale curtisi \n LOD: 4.2 \n parasites/µl'),size=4.5,family="serif",color="#2c7bb6")+
  geom_text(aes(y=72.623,x=1.5, label='P. ovale wallikeri \n LOD: 41.2 \n parasites/µl'),size=4.5,family="serif",color="#2c7bb6")+
  scale_color_manual(breaks=c("Positive","Negative"),values=c("#2c7bb6","#d7191c"))

print(LOD)
