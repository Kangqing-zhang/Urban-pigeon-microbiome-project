---
  title: "Seasonal dynamics of the gut microbiome in urban feral pigeons are associated with environmental conditions, not with diet shifts"
author: "Kangqing Zhang"
date: "October 09, 2025"
---

  #####################################################################################

#R 4.4.0
rm(list=ls())
setwd("D:/pigeon/Microbiom_diet")
#install.packages('devtools')
# install.packages("igraph")
#devtools::install_github('gastonstat/plspm')
library(plspm);packageVersion("plspm")
library(ggplot2)
library(openxlsx)
library(vegan)
library(ape)
library(dplyr)
library(igraph)

metadata = read.csv(file = "metadata_FP.csv", header = TRUE,sep="", check.names = FALSE, row.names = 1)

#host traits: BCI index (bodymass ~ wing + headbill), age
#environment factors: Pca (daily mean temp, max temp, min temp, precipitation duration, daily precipitation amount), day length
#diet: Observed and Jaccard (Presence/Absence data)
#microbiome: Observed and Weighted UniFrac

metadata$BodyMass <- as.numeric(gsub(",", ".", metadata$BodyMass))
metadata$Wing <- as.numeric(gsub(",", ".", metadata$Wing))
metadata$HeadBill <- as.numeric(gsub(",", ".", metadata$HeadBill))
metadata$Daylength <- as.numeric(gsub(",", ".", metadata$Daylength))

### 1. BCI (body condication index)
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

summary(temp_pca)            #Proportion of Variance / Cumulative PC1(95.96%)
temp_pca$rotation            #loading
plot(temp_pca, type = "l", main = "Temperature PCA scree")
metadata$Temp_PC1 <- temp_pca$x[, 1]

##rainfall
prec_mat <- metadata[, c("DR", "RH")]
prec_pca <- prcomp(prec_mat, scale. = TRUE)

summary(prec_pca)   #PC1(75.21%)
prec_pca$rotation

metadata$Prec_PC1 <- prec_pca$x[, 1]

write.csv(metadata,"metadata.csv", row.names = T)
write.xlsx(metadata, file = file.path("D:/pigeon/Microbiom_diet", "metadata.xlsx"),rowNames= T)

### 3. diet
plspm_diet_genus = read.csv(file = "plspm_diet_genus.csv", header = TRUE,sep=";", check.names = FALSE, row.names = 1)
plspm_diet_genus<- data.matrix(plspm_diet_genus)
Observed<- specnumber(plspm_diet_genus)

metadata$Diet_Observed<- Observed[rownames(metadata)]

# Jaccard（presence/absence）
dist_jac <- vegdist(plspm_diet_genus, method = "jaccard", binary = TRUE)

pcoa <- cmdscale(dist_jac, k = 2, eig = TRUE)

eigvals <- pcoa$eig
variance_explained <- eigvals / sum(eigvals) * 100
variance_explained[1:2]
dist_pcoa2 <- dist(pcoa$points[, 1:2])
mantel(vegdist(plspm_diet_genus, method = "jaccard", binary = TRUE), dist_pcoa2)

# Cailliez
diet_pcoa <- pcoa(dist_jac, correction = "cailliez")

names(diet_pcoa$values)
diet_pcoa$values$Rel_corr_eig[1:2] * 100

metadata$Diet_PC1 <- diet_pcoa$vectors[,1]
metadata$Diet_PC2 <- diet_pcoa$vectors[,2]

write.csv(metadata,"metadata.csv", row.names = T)
write.xlsx(metadata, file = file.path("D:/pigeon/Microbiom_diet", "metadata.xlsx"),rowNames= T)


##########################################pls-pm###############################
dat<-read.delim('plspm_9.txt',sep='\t')
dat[] <- lapply(dat, function(x) as.numeric(gsub(",", ".", x)))
str(dat)

#Specifies latent variables and relationships, storing the relationship between variables and latent variables in R as a list
#Latent variable ordering is consistent with relation matrix order *dat_blocks

dat_blocks<-list(Environment = c('Temp_PC1'),Diet=c('Diet_PCoA1','Diet_Observed'),Host =c('BCI','Sex'), Microbiome=c('Microbiome_PCoA1','Microbiome_Observed'))

#Correlations between latent variables are described by a 0-1 matrix, where 0 represents no correlation between variables and 1 represents correlation
Environment<- c(0,0,0,0)
Diet<-c(1,0,0,0)
Host <- c(1,1,0,0)
Microbiome<-c(1,1,1,0)

#Merge data
dat_path <- rbind(Environment,Diet,Host,Microbiome)
colnames(dat_path) <- rownames(dat_path)
dat_path

#Specify causation, either A or B
dat_modes <- rep('A', 4)
dat_modes
dat_modes <- rep('B', 4)
dat_modes <- c('B', 'B','A','A')

#PLS-PM
dat_pls<-plspm(dat, dat_path, dat_blocks, modes = dat_modes)
dat_pls

#View parameter estimates for path coefficients and related statistics
dat_pls$path_coefs
dat_pls$inner_model
dat_pls$unidim

#View the cause-and-effect path map, see ?innerplot
innerplot(dat_pls, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray', box.lwd = 0)

#To view the relationship between observed variables and latent variables, the outerplot() can be used to map a structure similar to the path diagram, see?outerplot
dat_pls$outer_model
outerplot(dat_pls, what = 'loadings', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')
outerplot(dat_pls, what = 'weights', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')

#The goodness of the model can be evaluated by the goodness-of-fit
dat_pls$gof #0.26
dat_pls$effects

#View the state of direct or indirect influence between variables
dat_pls$effects
summary(dat_pls)
plot(dat_pls)

boot_res <- plspm::plspm(dat, dat_path, dat_blocks, modes = dat_modes, boot.val = TRUE, br = 5000)
summary(boot_res)

#View the latent variable score, which can be understood as the value of the standardized latent variable
dat_pls$scores
dat_pls$data
# Create a data frame containing the path name, path coefficient, and standard error
out.long <- data.frame(
  paths = c("Environment -> Diet", "Environment -> Host", "Diet -> Host",
            "Environment -> Microbiome", "Diet -> Microbiome", "Host -> Microbiome"),
  coef = c(0.12, -0.27, -0.10, -0.49,0.01, 0.18),
  SE = c(0.13, 0.12, 0.12, 0.12, 0.11, 0.12)
)


# Generate a path coefficient graph
#Path Coefficients with Standard Errors
ggplot(data = out.long, aes(x = paths, y = coef)) +
  geom_point(position = position_dodge(0.25), size = 3) +
  geom_errorbar(aes(ymin = coef - SE, ymax = coef + SE), width = 0, position = position_dodge(0.25)) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c("#FEB24C")) +
  labs(title = "", x = "Paths", y = "Coefficient")

#
edges <- data.frame(
  from   = c("Environment","Environment","Environment","Diet","Diet","Host"),
  to     = c("Diet","Host","Microbiome","Host","Microbiome","Microbiome"),
  beta   = c( 0.12, -0.26, -0.48, -0.08, 0.08, 0.06),
  p      = c(0.35, 0.04, 0.0002, 0.55, 0.53, 0.60)
)

#
b <- abs(edges$beta); rng <- diff(range(b))
w <- if (rng == 0) rep(5, length(b)) else 1 + 12*(b - min(b))/rng

g <- graph_from_data_frame(edges, directed = TRUE)

#
coords <- rbind(
  Environment = c(0, 1),
  Diet        = c(-1, 0),
  Host        = c(0, -1),
  Microbiome  = c(1, 0)
)


stars <- ifelse(E(g)$p < 0.0001, "****",
                ifelse(E(g)$p < 0.001, "***",
                       ifelse(E(g)$p < 0.01,  "**",
                              ifelse(E(g)$p < 0.05,  "*", ""))))

plot(g,
     layout      = coords[V(g)$name, ],
     vertex.size = 34,
     vertex.color= "grey90",
     vertex.label.cex = 1.1,
     vertex.label.color = "black",

     edge.width  = w,
     edge.color  = ifelse(E(g)$beta > 0, "red", "blue"),
     edge.arrow.size = 1.5,
     edge.label  = paste0("β = ", sprintf("%.2f", E(g)$beta),
                          "\n", "P = ", signif(E(g)$p, 2), " ", stars),
     edge.label.cex = 0.9,
     edge.label.color = "black"
)






