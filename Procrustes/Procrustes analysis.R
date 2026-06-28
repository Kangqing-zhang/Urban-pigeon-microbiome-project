---
  title: "Seasonal dynamics of the gut microbiome in urban feral pigeons are associated with environmental conditions, not with diet shifts"
author: "Kangqing Zhang"
date: "October 09, 2025"
---

  #####################################################################################

#Procrustes analysis on genus level
#R 4.2.3
rm(list = ls())
#set working environment
setwd("D:/pigeon/Microbiom_diet/Procrustes analysis")
getwd()

library("vegan");packageVersion("vegan")
library("ggplot2")
library("cowplot")
library(dplyr)
library(tidyr)

#jaccard similarity of diet and microbiome, Mantel test
metadata = read.csv(file = "metadata_FP.csv", header = TRUE,sep="", check.names = FALSE, row.names = 1)

Phylum_abundance = read.csv(file = "Phylum_abundance.csv", header= TRUE,sep=";",check.names = FALSE, row.names = 1)
genus_abundance = read.csv(file = "genus_abundance.csv", header= TRUE,sep=";",check.names = FALSE, row.names = 1)

Phylum_pa <- (data.matrix(Phylum_abundance) > 0) + 0L
Genus_pa  <- (data.matrix(genus_abundance)  > 0) + 0L

write.csv(Phylum_pa,"Phylum_pa.csv", row.names = T)
write.csv(Genus_pa,"Genus_pa.csv", row.names = T)

#diet
diet_abundance = read.csv(file = "diet_abundance.csv", header= TRUE,sep=";",check.names = FALSE, row.names = 1)
SampleID= read.csv(file = "SampleID.csv", header= TRUE,sep=";",check.names = FALSE)

merged_diet_abundance<-merge(diet_abundance,SampleID, by.x=c("Sample"), by.y=c("Sample"))
merged_diet_abundance<-merged_diet_abundance[,-1]
merged_diet_abundance <- merged_diet_abundance[, c("SampleID", setdiff(names(merged_diet_abundance), "SampleID"))]
write.csv(merged_diet_abundance,"merged_diet_abundance.csv", row.names = T)

barcode_abundance <- merged_diet_abundance %>%
  group_by(SampleID, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups="drop") %>%
  pivot_wider(names_from = Genus,
              values_from = Abundance,
              values_fill = 0)

barcode_abundance <- as.data.frame(barcode_abundance)
rownames(barcode_abundance) <- barcode_abundance$SampleID
barcode_abundance <- as.matrix(barcode_abundance[,-1])
write.table(barcode_abundance, "barcode_abundance.csv", sep=";", quote=FALSE)

barcode_pa <- (data.matrix(barcode_abundance)  > 0) + 0L
write.csv(barcode_pa,"barcode_pa.csv", row.names = T)


#Procrustes analysis on genus level
#data used: Genus_pa (microbe) and barcode_pa (diet)

df1 <- barcode_pa
df2 <- Genus_pa

#Presence/Absence data- jaccard dissimilarity
dist.abund <- vegdist(df1, method = "jaccard")
dist.abund <- as.dist(dist.abund)
mdist.abund = vegdist(df2, method = "jaccard")
mdist.abund <- as.dist(mdist.abund)
#make pcoas
dpcoa <- as.data.frame(cmdscale(dist.abund))

mpcoa <- as.data.frame(cmdscale(mdist.abund))

#procrustes analysis
pro <- procrustes(X = dpcoa, Y = mpcoa, scale = TRUE,symmetric = TRUE)

pro_test <- protest(dpcoa,mpcoa,perm=9999)

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
#beta_pro$type <- "Diet(Bray_curtis)"
beta_pro$type <- "Diet(jaccard)"
seasons=metadata[,20]
beta_pro=cbind(beta_pro,seasons)

trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Microbiome"
seasons=metadata[,20]
trans_pro=cbind(trans_pro,seasons)

colnames(trans_pro) <- colnames(beta_pro)
pval <- pro_test$signif
plot <- rbind(beta_pro, trans_pro)


colslocal <- c("darkseagreen", "darkgreen","saddlebrown")
colsSeason <- c("#DA5724","#508578")

#season
diet_microbiome_season <- ggplot(plot) +
  geom_point(size = 4, alpha=0.75, aes(x = V1, y = V2, color = seasons,shape=type))+
  scale_color_manual(values = colsSeason) +
  theme_classic() +
  scale_x_continuous()+
  scale_y_continuous()+
  geom_line(aes(x= V1, y=V2, group=UserName), col = "darkgrey", alpha = 0.6,linewidth=0.2) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10,colour="black"),
        legend.position = 'bottom',
        axis.text = element_text(size=10,colour="black"),
        axis.title = element_text(size=13,colour="black"),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  #annotate("text", x = 0.07, y = -0.13, label = paste0("p-value=",pval), size = 4) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]"))
diet_microbiome_season
diet_microbiome_season_leg <- get_legend(diet_microbiome_season)

diet_microbiome_season + theme(legend.position = "right")

###location
Location=metadata[,6]
beta_pro=cbind(beta_pro,Location)
trans_pro=cbind(trans_pro,Location)

colnames(trans_pro) <- colnames(beta_pro)
pval <- pro_test$signif
plot <- rbind(beta_pro, trans_pro)

#locations
diet_microbiome_location <- ggplot(plot) +
  geom_point(size = 4, alpha=0.75, aes(x = V1, y = V2, color = Location,shape=type))+
  scale_color_manual(values = colslocal) +
  theme_classic() +
  scale_x_continuous()+
  scale_y_continuous()+
  geom_line(aes(x= V1, y=V2, group=UserName), col = "darkgrey", alpha = 0.6,size=0.2) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10,colour="black"),
        legend.position = 'bottom',
        axis.text = element_text(size=10,colour="black"),
        axis.title = element_text(size=13,colour="black"),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  #annotate("text", x = 0.07, y = -0.13, label = paste0("p-value=",pval), size = 4) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]"))

diet_microbiome_location
diet_microbiome_location_leg <- get_legend(diet_microbiome_location)

diet_microbiome_location + theme(legend.position = "right")



















