########################################################################################################################################
#Programmer: Wenqiao He
#Last update: May 2024
#Purpose: P. ovale real-time PCR assays Analysis
########################################################################################################################################

#Loading packages for use
library(readxl)
library(ggplot2)

# Set working directory
setwd("FILE_PATH")

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------

# Field sample analysis and plotting##################################################################################################### 

##Loading data
DRC_Po <- read_excel("FILE_NAME")

##Calculating average Cq values
DRC_Po$Cq_mean_18S <- rowMeans(DRC_Po[, c("Cq1_18S", "Cq2_18S")], na.rm = TRUE)

DRC_Po$Cq_mean_Poc <- rowMeans(DRC_Po[, c("Cq1_Poc", "Cq2_Poc")], na.rm = TRUE)

DRC_Po$Cq_mean_Pow <- rowMeans(DRC_Po[, c("Cq1_Pow", "Cq2_Pow")], na.rm = TRUE)

##Determine P. ovale species 
DRC_Po <- DRC_Po |> dplyr::mutate(Diagnosis = dplyr::case_when(Cq_mean_Poc < 45 & Cq_mean_Pow < 45 ~ "Mixed",
                                                     Cq_mean_Poc < 45 ~ "P. ovale curtisi",
                                                     Cq_mean_Pow < 45 ~ "P. ovale wallikeri",
                                                     .default = "Not determined"))

##Making figure
Po_diagnosis <-ggplot(DRC_Po, aes(y=`Parasite_density_18S`, x=Species, colour=Diagnosis,fill=Diagnosis))+
  geom_point(position=position_jitter(w=0.15,height = 0),size=4,shape=21,alpha=0.85)+
  guides(fill=guide_legend(title=" "))+
  labs(title=" ",shape=" ",colour=" ")+
  xlab("Species")+
  ylab("Parasite density (parasites/µl)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(text=element_text(family="serif"),title=element_text(size=16,face="bold"))+
  geom_segment(aes(y=4.2,x=0.7,yend=4.2,xend=4.2),linetype = "dashed",size=0.65, colour ="#2c7bb6")+
  geom_segment(aes(y=41.2,x=0.7,yend=41.2,xend=4.2),linetype = "dashed",size=0.65, colour ="#fc8d62")+
  scale_x_discrete(labels=c("P.ovale","P.ovale and \nP.falciparum","P.ovale and \nP.malariae","P.ovale, P.falciparum,\n and P.malariae"))+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.7,size=12),axis.text.y = element_text(size=12),legend.text = element_text(size = 16))+
  geom_text(aes(x=1.2,y=3.2, label='P. ovale curtisi LOD'),size=3,family="serif",color="#2c7bb6")+
  geom_text(aes(x=3.2,y=53.6, label='P. ovale wallikeri LOD'),size=3,family="serif",color="#fc8d62")+
  scale_color_manual(breaks=c("P. ovale curtisi","P. ovale wallikeri","Mixed","Not determined"),values=c("#2c7bb6", "#fc8d62", "#66c2a5", "black"))+
  scale_fill_manual(breaks=c("P. ovale curtisi","P. ovale wallikeri","Mixed","Not determined"), values=c("#2c7bb6", "#fc8d62", "#66c2a5", "white"))

print(Po_diagnosis)
