---
  title: "Seasonal dynamics of the gut microbiome in urban feral pigeons are associated with environmental conditions, not with diet shifts"
author: "Kangqing Zhang"
date: "March 09, 2026"
---

  #####################################################################################

#R 4.4.0
rm(list=ls())
setwd("D:/pigeon/dbRDA")
library("phyloseq");packageVersion('phyloseq')
library(permute)
library(lattice)
library(carData)
library(ggplot2)
library(openxlsx)
library(vegan)
library(car)
library(dplyr)
library(ape)

# load phyloseq object
noncontam_physeq_rarefied2 <- readRDS("noncontam_physeq_rarefied2.rds")
taxa_noncontam_physeq_rarefied2 = as(tax_table(noncontam_physeq_rarefied2), "matrix")
write.table(
  taxa_noncontam_physeq_rarefied2,
  "taxa_noncontam_physeq_rarefied2.csv",
  sep = ";",
  quote = TRUE,
  row.names = TRUE,
  fileEncoding = "UTF-8"
)
ASV_noncontam_physeq_rarefied2= as(otu_table(noncontam_physeq_rarefied2), "matrix")
write.table(
  ASV_noncontam_physeq_rarefied2,
  "ASV_noncontam_physeq_rarefied2.csv",
  sep = ";",
  quote = TRUE,
  row.names = TRUE,
  fileEncoding = "UTF-8"
)
metadata_noncontam_physeq_rarefied2 <- data.frame(sample_data(noncontam_physeq_rarefied2))
write.table(
  metadata_noncontam_physeq_rarefied2,
  "metadata_noncontam_physeq_rarefied2.csv",
  sep = ";",
  quote = TRUE,
  row.names = TRUE,
  fileEncoding = "UTF-8"
)

# JACCARD DISTANCES (presence-absence)
FP_rarefied1_jacc <- phyloseq::distance(noncontam_physeq_rarefied2, "jaccard")

metadata <- as(sample_data(noncontam_physeq_rarefied2), "data.frame")
head(metadata)

metadata$BodyMass <- as.numeric(gsub(",", ".", metadata$BodyMass))
metadata$Wing <- as.numeric(gsub(",", ".", metadata$Wing))
metadata$HeadBill <- as.numeric(gsub(",", ".", metadata$HeadBill))
metadata$Daylength <- as.numeric(gsub(",", ".", metadata$Daylength))

### 1. BCI
cor(metadata$Wing, metadata$HeadBill) #0.49
library(car)
vif(lm(BodyMass ~ Wing + HeadBill, data = metadata)) #Wing (1.322011) HeadBill ( 1.322011 )
metadata$BCI <- resid(lm(BodyMass ~ Wing + HeadBill, data = metadata))

### 2. environment factor(PCA)
metadata <- metadata %>%
  mutate(
    TG = TG / 10,
    TN = TN / 10,
    TX = TX / 10)

climate_vars <- metadata[, c("TG", "TN", "TX", "DR", "RH")]
cor_matrix <- cor(climate_vars, use = "complete.obs") #"TG", "TN", "TX" have high correlation (r > 0.8),"DR", "RH" (r = 0.5)

##temperature
temp_mat <- metadata[, c("TG", "TN", "TX")]
temp_pca <- prcomp(temp_mat, scale. = TRUE)

summary(temp_pca)
temp_pca$rotation
plot(temp_pca, type = "l", main = "Temperature PCA scree")
metadata$Temp_PC1 <- temp_pca$x[, 1]

##rainfall
prec_mat <- metadata[, c("DR", "RH")]
prec_pca <- prcomp(prec_mat, scale. = TRUE)

summary(prec_pca)   #PC1(75.21%)
prec_pca$rotation

metadata$Prec_PC1 <- prec_pca$x[, 1]

write.csv(metadata,"metadata.csv", row.names = T)
write.xlsx(metadata, file = file.path("D:/pigeon/dbRDA", "metadata.xlsx"),rowNames= T)

### 3. diet
diet_genus = read.csv(file = "diet_genus.csv", header = TRUE,sep=";", check.names = FALSE, row.names = 1)
diet_genus<- data.matrix(diet_genus)
Observed<- specnumber(diet_genus)

metadata$Diet_Observed<- Observed[rownames(metadata)]

# Jaccard（presence/absence）
dist_jac <- vegdist(diet_genus, method = "jaccard", binary = TRUE)

pcoa <- cmdscale(dist_jac, k = 2, eig = TRUE)

eigvals <- pcoa$eig
variance_explained <- eigvals / sum(eigvals) * 100
variance_explained[1:2]
dist_pcoa2 <- dist(pcoa$points[, 1:2])
mantel(vegdist(diet_genus, method = "jaccard", binary = TRUE), dist_pcoa2)

# Cailliez
diet_pcoa <- pcoa(dist_jac, correction = "cailliez")

names(diet_pcoa$values)
diet_pcoa$values$Rel_corr_eig[1:2] * 100

metadata$Diet_PC1 <- diet_pcoa$vectors[,1]
metadata$Diet_PC2 <- diet_pcoa$vectors[,2]

write.csv(metadata,"metadata.csv", row.names = T)
write.xlsx(metadata, file = file.path("D:/pigeon/dbRDA", "metadata.xlsx"),rowNames= T)

#
dist_micro <- as.dist(FP_rarefied1_jacc)
#dist_micro <- as.dist(FP_rarefied1_bray)
#dist_micro <- as.dist(unweighted_unif_rar1)
dist_micro <- as.dist(weighted_unif_rar1)

m_full <- capscale(dist_micro ~ Temp_PC1  + Prec_PC1 +Location+
                     Diet_PC1 + Diet_Observed + BCI +Sex,
                   data = metadata)

anova(m_full, by="margin", permutations=999)
RsquareAdj(m_full)
vif.cca(m_full)

m_full_1 <- capscale(dist_micro ~ Temp_PC1  + Prec_PC1 +
                     Diet_PC1 + Diet_Observed + BCI +Sex+
                     Condition(Location),
                   data = metadata)
anova(m_full_1, by="margin", permutations=999)
RsquareAdj(m_full_1)

###variation partitioning
otu <- as(otu_table(noncontam_physeq_rarefied2), "matrix")
if (taxa_are_rows(noncontam_physeq_rarefied2)) otu <- t(otu)
otu_hel <- decostand(otu, "hellinger")

X_env  <- metadata[, c("Temp_PC1", "Prec_PC1")]
X_diet <- metadata[, c("Diet_PC1", "Diet_Observed")]
X_host <- metadata[, c("BCI", "Sex")]
X_host$Sex <- factor(X_host$Sex)

X_space <- data.frame(Location = factor(metadata$Location))
rownames(X_space) <- rownames(metadata)
stopifnot(all(rownames(X_space) == rownames(otu_hel)))

vp <- varpart(otu_hel, X_env, X_diet, X_host, X_space)
vp
plot(vp)

rda_env  <- rda(otu_hel ~ ., data = X_env,  Condition = cbind(X_diet, X_host, X_space))
rda_diet <- rda(otu_hel ~ ., data = X_diet, Condition = cbind(X_env,  X_host, X_space))
rda_host <- rda(otu_hel ~ ., data = X_host, Condition = cbind(X_env,  X_diet, X_space))
rda_spa  <- rda(otu_hel ~ ., data = X_space,Condition = cbind(X_env,  X_diet, X_host))

anova(rda_env,  permutations=999)
anova(rda_diet, permutations=999)
anova(rda_host, permutations=999)
anova(rda_spa,  permutations=999)


# PCoA
pcoa <- cmdscale(dist_micro, k=2, eig=TRUE)
scores <- as.data.frame(pcoa$points)
colnames(scores) <- c("PCoA1","PCoA2")
scores$Season <- metadata$Season
scores$Temp_PC1 <- metadata$Temp_PC1

# axis variance explained
var_exp <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)

# envfit
fit <- envfit(pcoa$points ~ Temp_PC1 + Prec_PC1 + Diet_PC1 + Diet_Observed + BCI, data=metadata, permutations=999)

# extract arrow
arrow <- as.data.frame(scores(fit, display="vectors"))
arrow$var <- rownames(arrow)

ggplot(scores, aes(PCoA1, PCoA2, color=Season)) +
  geom_point(size=2, alpha=0.8) +
  stat_ellipse(level=0.95) +
  geom_segment(data=arrow,
               aes(x=0, y=0, xend=Dim1, yend=Dim2),
               inherit.aes = FALSE,
               arrow=arrow(length=unit(0.25,"cm"))) +
  geom_text(data=arrow,
            aes(x=Dim1, y=Dim2, label=var),
            inherit.aes = FALSE,
            vjust=-0.5) +
  labs(x=paste0("PCoA1 (", var_exp[1], "%)"),
       y=paste0("PCoA2 (", var_exp[2], "%)")) +
  theme_classic()

m <- m_full_1
sites <- as.data.frame(scores(m, display = "sites", choices = 1:2))
colnames(sites) <- c("CAP1", "CAP2")

sites$Season   <- factor(metadata[rownames(sites), "Season"])
sites$Location <- factor(metadata[rownames(sites), "Location"])

eig <- m$CCA$eig
cap1_exp <- round(100 * eig[1] / sum(eig), 1)
cap2_exp <- round(100 * eig[2] / sum(eig), 1)

bp <- as.data.frame(scores(m, display = "bp", choices = 1:2))
bp$var <- rownames(bp)

bp <- bp[bp$var %in% c("Temp_PC1"), , drop = FALSE]

arrow_mul <- min(
  (max(abs(sites$CAP1)) / max(abs(bp$CAP1))),
  (max(abs(sites$CAP2)) / max(abs(bp$CAP2)))
) * 0.8
bp$CAP1 <- bp$CAP1 * arrow_mul
bp$CAP2 <- bp$CAP2 * arrow_mul

p <- ggplot(sites, aes(CAP1, CAP2, color = Season, shape = Location)) +
  geom_point(size = 2.2, alpha = 0.85) +
  stat_ellipse(aes(group = Season, color = Season), level = 0.95, linewidth = 0.6) +
  stat_ellipse(aes(group = Season, color = Season), level = 0.95, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4) +
  geom_segment(
    data = bp,
    aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.6,
    color = "black"
  ) +
  geom_text(
    data = bp,
    aes(x = CAP1, y = CAP2, label = var),
    inherit.aes = FALSE,
    vjust = -0.6,
    size = 3,
    color = "black"
  ) +
  labs(
    x = paste0("CAP1 (", cap1_exp, "%)"),
    y = paste0("CAP2 (", cap2_exp, "%)"),
    color = "Season",
    shape = "Location"
  ) +
  theme_classic()


ggsave("Figure 5.tif", p, width = 7.5, height = 5.5)



