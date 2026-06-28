---
  title: "Seasonal dynamics of the gut microbiome in urban feral pigeons are associated with environmental conditions, not with diet shifts"
author: "Kangqing Zhang"
date: "October 09, 2025"
---

  ####################################################################################
rm(list=ls())
setwd("D:/pigeon/seq/01052024/phyloseq")

packageVersion('phyloseq')
library("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library(ape)
theme_set(theme_bw())
library(phyloseq)
library(microbiome)
library(decontam)
library(vsn)
library(ggplot2)
library(vegan);packageVersion('vegan')
library(plyr)
library(dplyr)
library(ape)
library(venn)
library(directlabels)
library(doBy)
library(Rmisc)
library(lme4)
library(nlme);packageVersion("nlme")
library(nortest)
library(data.table)
library(cluster)
library(rockchalk)
library(reshape2)
library(blmeco)
library(glmmTMB)
library(bbmle)
library(MuMIn)
library(rcompanion)
library(qqplotr)
library(car)
library(btools);packageVersion('btools')
library(multcomp)
library(multcompView)
library(cowplot)
library(devtools)
library(pairwiseAdonis) # for posthoc ADONIS analysis
library(emmeans)
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc)
library(effects)
library(sjstats) #use for r2 functions
library(ggeffects)
library(ade4) # v 1.7-16
library(ICC)
library(ggpubr)
library(openxlsx)
library(tidyr)
library(pacman)
colslocal <- c("darkseagreen", "darkgreen","saddlebrown")
colsSeason <- c("#DA5724","#508578")
colsSex <- c("#CC79A7","#56B4E9")
Sex_labels <- c("F" = "females",
                "M" = "males")
Local_labels<-c("Molukkenpad","Noorderplantsoen","Vismarkt")


## load in otu table
otu_able= read.csv("otu_matrix.csv", sep=";", row.names = 1)
otu_table=as.matrix(otu_able)
##load in taxonomy
taxonomy = read.csv("taxonomy.csv", sep=";", row.names=1)
taxonomy= as.matrix(taxonomy)

metadata = read.csv("pigeons_metadata.csv", sep=";", row.names=1)
phy_tree = read_tree("tree.nwk")

class(otu_table)
#####------------------------------------------------------------------------------
## create phyloseq object
OTU = otu_table(otu_table, taxa_are_rows=TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)

taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

sample_names(OTU)
sample_names(META)

physeq = phyloseq(OTU, TAX, META, phy_tree)
physeq
rank_names(physeq)

physeq #1119 taxa and 74 samples
sum(sample_sums(physeq))
sample_sums(physeq)
range(taxa_sums(physeq))
range(sample_sums(physeq))
summary(sample_sums(physeq))
sample_variables(physeq)
table(tax_table(physeq)[,"Phylum"], exclude=NULL)
table(tax_table(physeq)[,"Kingdom"], exclude=NULL)
table(tax_table(physeq)[,"Genus"], exclude=NULL)
## remove Unassigned
## remove ASVs with 0 reads
physeq_clean <- subset_taxa(physeq, Kingdom != "Unassigned")
physeq_clean <- prune_taxa(taxa_sums(physeq_clean)>0, physeq_clean)
# check removal
table(tax_table(physeq_clean)[,"Kingdom"], exclude=NULL)

# Check data again
physeq_clean #1116 taxa and 74 samples
sum(sample_sums(physeq_clean))
sample_sums(physeq_clean)
range(taxa_sums(physeq_clean))
range(sample_sums(physeq_clean))
summary(sample_sums(physeq_clean))
table(tax_table(physeq_clean)[,"Phylum"], exclude=NULL)

# Make taxa and ASV table en write to file
taxa_physeq_clean = as(tax_table(physeq_clean), "matrix")
write.csv(taxa_physeq_clean, file='taxa_physeq_clean.csv')
ASV_physeq_clean = as(otu_table(physeq_clean), "matrix")
write.csv(ASV_physeq_clean, file='ASV_physeq_clean.csv')
metadata_physeq_clean <- data.frame(sample_data(physeq_clean))
write.csv(metadata_physeq_clean, "metadata_physeq_clean.csv")
tree_physeq_clean <- write.tree(phy_tree(physeq_clean))
write(tree_physeq_clean, file = "tree_physeq_clean.newick")

#save phyloseq object
saveRDS(physeq_clean, "physeq_clean.rds")
# load phyloseq object
#physeq_clean <- readRDS("physeq_clean.rds")

# Remove contamination from the data by using the package decontam.
# All options were compared: 1) the relationship between frequency & DNA concentration,
# 2) prevalence in NC samples.
# 13 of the 77 samples are negative controls, of which 5 are negative controls fro kit 1 to kit 5.
# 2 are PCR negative controls and the others (6) are for swabs.

# First inspect library sizes (i.e. the number of reads) in each sample
# as a function of whether that sample was a true positive sample or a negative control:
df_physeq_clean <- as.data.frame(sample_data(physeq_clean))
df_physeq_clean$LibrarySize <- sample_sums(physeq_clean)
df_physeq_clean <- df_physeq_clean[order(df_physeq_clean$LibrarySize),]
df_physeq_clean$Index <- seq(nrow(df_physeq_clean))

colsGroup3 <- c("#CC79A7","#56B4E9","#636363")
ggplot(data=df_physeq_clean, aes(x=Index, y=LibrarySize, color=Group)) +
  geom_point()+
  scale_colour_manual(values=colsGroup3)
#NCs have lowest read counts

## Identify contaminants - frequency based
# In this method, the distribution of the frequency of each sequence feature as a
# function of the input DNA concentration is used to identify contaminants.
seqtab <- as(otu_table(physeq_clean), "matrix")
if (taxa_are_rows(physeq_clean)) seqtab <- t(seqtab)
seqtab <- data.matrix(seqtab)
meta <- data.frame(sample_data(physeq_clean))
stopifnot("DNAconcentration" %in% names(meta))
conc_raw <- meta$DNAconcentration
conc <- suppressWarnings(readr::parse_number(as.character(conc_raw)))
names(conc) <- rownames(meta)
conc <- conc[rownames(seqtab)]
keep <- !is.na(conc) & conc > 0
seqtab <- seqtab[keep, , drop = FALSE]
conc   <- conc[keep]

contamdf.freq_FP_clean <- isContaminant(
  seqtab,
  method = "frequency",
  conc   = conc
)

# Inspect the results
table(contamdf.freq_FP_clean$contaminant)
head(which(contamdf.freq_FP_clean$contaminant))
write.csv(contamdf.freq_FP_clean, file='contamdf.freq_FP_clean_10092025.csv')
write.xlsx(contamdf.freq_FP_clean, file = file.path("D:/pigeon/seq/01052024/phyloseq", "contamdf.freq_FP_clean.xlsx"), rowNames = T)
# 5 out of the 1115 ASVs were classified as contaminant.
#find the contaminants ASVs
seqtab <- as(otu_table(physeq_clean), "matrix")
if (taxa_are_rows(physeq_clean)) seqtab <- t(seqtab)
meta <- data.frame(sample_data(physeq_clean))
meta$Group <- trimws(as.character(meta$Group))
meta$is_neg <- meta$Group == "NC"
meta <- meta[rownames(seqtab), , drop = FALSE]
neg  <- meta$is_neg
contam_res <- data.frame(
  ASV         = rownames(contamdf.freq_FP_clean),
  p           = contamdf.freq_FP_clean$p,
  contaminant = contamdf.freq_FP_clean$contaminant,
  row.names   = NULL
)
pollutant_ASVs <- contam_res %>% filter(contaminant) %>% pull(ASV)
prev_neg <- colSums(seqtab[neg, pollutant_ASVs, drop = FALSE] > 0, na.rm = TRUE)
prev_pos <- colSums(seqtab[!neg, pollutant_ASVs, drop = FALSE] > 0, na.rm = TRUE)

summary_tbl <- tibble(
  ASV = pollutant_ASVs,
  prev_in_neg = as.integer(prev_neg[ASV]),
  prev_in_pos = as.integer(prev_pos[ASV])
) %>%
  mutate(
    total_prev    = prev_in_neg + prev_in_pos,
    neg_share     = ifelse(total_prev > 0, prev_in_neg / total_prev, NA_real_),
    any_in_neg    = prev_in_neg > 0,
    mostly_in_neg = neg_share > 0.5
  ) %>%
  left_join(contam_res, by = "ASV") %>%
  arrange(desc(any_in_neg), desc(neg_share))
head(summary_tbl)
# However, both of them were not present in the NCs!

## Identify contaminants  - prevalence
# In this method, the prevalence (presence/absence across samples) of each sequence feature
# in true positive samples is compared to the prevalence in negative controls to identify contaminants.
sample_data(physeq_clean)$is.neg <- sample_data(physeq_clean)$Group == "NC"
contamdf.prev_physeq_clean <- isContaminant(physeq_clean, method="prevalence", neg="is.neg")
table(contamdf.prev_physeq_clean$contaminant)
head(contamdf.prev_physeq_clean)
write.csv(contamdf.prev_physeq_clean, file='contamdf.prev_physeq_clean_v5_10092025.csv')

contamdf.prev_physeq_clean_05 <- isContaminant(physeq_clean, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev_physeq_clean_05$contaminant)
write.csv(contamdf.prev_physeq_clean_05, file="contamdf_prev_physeq_clean_05_V5_10092025.csv")

contamdf.prev_physeq_clean_06<- isContaminant(physeq_clean, method="prevalence", neg="is.neg", threshold=0.6)
table(contamdf.prev_physeq_clean_06$contaminant)
write.csv(contamdf.prev_physeq_clean_06, file="contamdf_prev_physeq_clean_06_V5_10092025.csv")

# Make phyloseq object of presence-absence in negative controls
physeq_clean.neg <- prune_samples(sample_data(physeq_clean)$Group == "NC", physeq_clean)
physeq_clean.neg.presence <- transform_sample_counts(physeq_clean.neg, function(abund) 1*(abund>0))

# Make phyloseq object of presence-absence in true positive samples
physeq_clean.pos <- prune_samples(sample_data(physeq_clean)$Group == "True sample", physeq_clean)
physeq_clean.pos.presence <- transform_sample_counts(physeq_clean.pos, function(abund) 1*(abund>0))
# Make data.frame of prevalence in positive and negative samples
df.pres_physeq_clean <- data.frame(prevalence.pos=taxa_sums(physeq_clean.pos.presence), prevalence.neg=taxa_sums(physeq_clean.neg.presence),
                                   contam.prev=contamdf.prev_physeq_clean_06$contaminant)
colsGroup2 <- c("#CC79A7","#56B4E9")
ggplot(data=df.pres_physeq_clean, aes(x=prevalence.neg, y=prevalence.pos, color=contam.prev)) +
  geom_point() +
  scale_colour_manual(values=colsGroup2)+
  xlab("Prevalence (Negative Controls)")+
  ylab("Prevalence (True Samples)")

# CHECK ASVs OF NCs AND SAMPLES
physeq_clean.neg
physeq_clean.pos

# Write taxa and count information to csv file
# Extract abundance matrix from the phyloseq objects
physeq_clean.neg_otu <- as(otu_table(physeq_clean.neg), "matrix")
head(physeq_clean.neg_otu)
physeq_clean.pos_otu <- as(otu_table(physeq_clean.pos), "matrix")
head(physeq_clean.pos_otu)

# Coerce to data.frame
physeq_clean.neg_otu_df <- as.data.frame(physeq_clean.neg_otu)
head(physeq_clean.neg_otu_df)
physeq_clean.pos_otu_df <- as.data.frame(physeq_clean.pos_otu)
head(physeq_clean.pos_otu_df)

#Save files as csv
write.table(physeq_clean.neg_otu_df, file = "physeq_clean.neg_otu_df.csv")
write.table(physeq_clean.pos_otu_df, file = "physeq_clean.pos_otu_df.csv")
# Note that all ASVs occur in the NCs (and that all features are in the table),

#
contam_res <- data.frame(
  ASV         = rownames(contamdf.prev_physeq_clean_06),
  contaminant = contamdf.prev_physeq_clean_06$contaminant
)
pollutant_ASVs <- contam_res %>% filter(contaminant) %>% pull(ASV)

seqtab <- as(otu_table(physeq_clean), "matrix")
if (taxa_are_rows(physeq_clean)) seqtab <- t(seqtab)

meta <- data.frame(sample_data(physeq_clean))
meta$Group <- trimws(as.character(meta$Group))
meta <- meta[rownames(seqtab), , drop = FALSE]
meta$SampleID <-rownames(meta)
poll_counts <- seqtab[, pollutant_ASVs, drop = FALSE]
merged_poll_counts<-cbind(poll_counts, meta)

write.csv(merged_poll_counts, file='merged_poll_counts_10092025.csv')
write.xlsx(merged_poll_counts, file = file.path("D:/pigeon/seq/01052024/phyloseq", "merged_poll_counts_10092025.xlsx"), rowNames = T)
#most of them have high number of read counts in NCs and a low number of read counts in the samples.

noncontam_physeq_clean_06<-prune_taxa(!contamdf.prev_physeq_clean_06$contaminant, physeq_clean)
noncontam_physeq_clean_06   #1085 taxa and 74 samples
sum(sample_sums(noncontam_physeq_clean_06))
sample_sums(noncontam_physeq_clean_06)
range(taxa_sums(noncontam_physeq_clean_06))
range(sample_sums(noncontam_physeq_clean_06))
summary(sample_sums(noncontam_physeq_clean_06))
table(tax_table(noncontam_physeq_clean_06)[,"Phylum"], exclude=NULL)

# The NC samples were removed from the data set, as they should not be included in the downstream analyses.
noncontam_physeq_clean_noNC <- subset_samples(noncontam_physeq_clean_06, Group != "NC")
noncontam_physeq_clean_noNC #1085 taxa and 64 samples
sum(sample_sums(noncontam_physeq_clean_noNC))
sample_sums(noncontam_physeq_clean_noNC)
range(taxa_sums(noncontam_physeq_clean_noNC))
range(sample_sums(noncontam_physeq_clean_noNC))
summary(sample_sums(noncontam_physeq_clean_noNC))
table(tax_table(noncontam_physeq_clean_noNC)[,"Phylum"], exclude=NULL)

noncontam_physeq_clean_noNC <- prune_taxa(taxa_sums(noncontam_physeq_clean_noNC)>0, noncontam_physeq_clean_noNC)
noncontam_physeq_clean_noNC <- prune_samples(sample_sums(noncontam_physeq_clean_noNC)>0, noncontam_physeq_clean_noNC)
noncontam_physeq_clean_noNC #1052 taxa and 64 samples

# Make taxa and ASV table en write to file
taxa_noncontam_physeq_clean_noNC = as(tax_table(noncontam_physeq_clean_noNC), "matrix")
write.csv(taxa_noncontam_physeq_clean_noNC, file='taxa_noncontam_physeq_clean_noNC.csv')
ASV_noncontam_physeq_clean_noNC = as(otu_table(noncontam_physeq_clean_noNC), "matrix")
write.csv(ASV_noncontam_physeq_clean_noNC, file='ASV_noncontam_physeq_clean_noNC.csv')
metadata_noncontam_physeq_clean_noNC <- data.frame(sample_data(noncontam_physeq_clean_noNC))
write.csv(metadata_noncontam_physeq_clean_noNC, "metadata_noncontam_physeq_clean_noNC.csv")
tree_noncontam_physeq_clean_noNC <- write.tree(phy_tree(noncontam_physeq_clean_noNC))
write(tree_noncontam_physeq_clean_noNC, file = "tree_noncontam_physeq_clean_noNC.newick")

#save phyloseq object
saveRDS(noncontam_physeq_clean_noNC, "noncontam_physeq_clean_noNC.rds")
# load phyloseq object
#noncontam_physeq_clean_noNC <- readRDS("noncontam_physeq_clean_noNC.rds")
# save phyloseq to biom file


# Extract abundance matrix from the phyloseq object
noncontam_physeq_clean_noNC_abund <- as(otu_table(noncontam_physeq_clean_noNC), "matrix")
head(noncontam_physeq_clean_noNC_abund)
# Coerce to data.frame
noncontam_physeq_clean_noNC_abund_df <- as.data.frame(noncontam_physeq_clean_noNC_abund)
head(noncontam_physeq_clean_noNC_abund_df)
#Save files as csv
write.table(noncontam_physeq_clean_noNC_abund_df, file = "noncontam_physeq_clean_noNC_abund_df.csv")

# Remove zeros that have occured due to removal of the NCs from the data set
noncontam_physeq_clean_noNC1 <-prune_taxa(taxa_sums(noncontam_physeq_clean_noNC)>0, noncontam_physeq_clean_noNC)
noncontam_physeq_clean_noNC1 #1052 taxa and 64 samples
sum(sample_sums(noncontam_physeq_clean_noNC1))
sample_sums(noncontam_physeq_clean_noNC1)
range(taxa_sums(noncontam_physeq_clean_noNC1))
range(sample_sums(noncontam_physeq_clean_noNC1))
summary(sample_sums(noncontam_physeq_clean_noNC1))
table(tax_table(noncontam_physeq_clean_noNC1)[,"Phylum"], exclude=NULL)
# Extract abundance matrix and write to file
noncontam_physeq_clean_noNC_min1 <- as(otu_table(noncontam_physeq_clean_noNC), "matrix")
head(noncontam_physeq_clean_noNC_min1)
write.table(noncontam_physeq_clean_noNC_min1, file = "noncontam_physeq_clean_noNC_min1.csv")


##### INSPECT DATA FOR RARE ASVs AND LOW PREVALENCE ASVs

# One way to inspect the data for rare ASVs is to check how many ASVs have low read counts.
# and check which ASvs have the highest abundance in the unfiltered data set.
# or you remove low read count ASVs from the data. E.g. remove the ones with read counts <2,
# and check the effects. Etc.

# remove singletons
noncontam_physeq_clean_noNC_min2 <-prune_taxa(taxa_sums(noncontam_physeq_clean_noNC1)>1, noncontam_physeq_clean_noNC1)
noncontam_physeq_clean_noNC_min2 #1052 taxa and 64 samples
sum(sample_sums(noncontam_physeq_clean_noNC_min2))
sample_sums(noncontam_physeq_clean_noNC_min2)
range(taxa_sums(noncontam_physeq_clean_noNC_min2))
range(sample_sums(noncontam_physeq_clean_noNC_min2))
summary(sample_sums(noncontam_physeq_clean_noNC_min2))

noncontam_physeq_clean_noNC2 <- prune_samples(
  sample_names(noncontam_physeq_clean_noNC2) != "FPS018",
  noncontam_physeq_clean_noNC2
)
noncontam_physeq_clean_noNC2 #1052 taxa and 63 samples

# Another way to inspect rare ASVs is to look at the prevalence of the ASVs.
# Generate a table with prevelance (the number of samples each taxa occurs in)
# and total abundance (total number of read counts) for each taxa.
noncontam_physeq_prevelance_df = apply(X = otu_table(noncontam_physeq_clean_noNC1),
                                       MARGIN = 1,
                                       FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
noncontam_physeq_prevelance_df = data.frame(Prevalence = noncontam_physeq_prevelance_df,
                                            TotalAbundance = taxa_sums(noncontam_physeq_clean_noNC1),
                                            tax_table(noncontam_physeq_clean_noNC1))
noncontam_physeq_prevelance_df[1:10,]
write.table(noncontam_physeq_prevelance_df, file = "noncontam_physeq_prevelance_df_09052024.csv")

# Next explore prevalence (and abundance) per Phylum to gain more insight into the data
plyr::ddply(noncontam_physeq_prevelance_df, "Phylum", function(df1)(cbind(mean(df1$Prevalence), sum(df1$Prevalence), mean(df1$TotalAbundance), sum(df1$TotalAbundance))))
# Make a prevalence plot per phylum with indicator line at 5%.
ggplot(noncontam_physeq_prevelance_df, aes(TotalAbundance, Prevalence / nsamples(noncontam_physeq_clean_noNC1),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum, nrow = 5) + theme(legend.position="none")

## MAKE RAREFACTION CURVES
# basic script from  https://github.com/joey711/phyloseq/issues/143 adapted for the data
# Rarefaction curves are required to determine to what level to rarefy the data before
# calculation alpha-diversities, etc. Note that the data must be rarefied as sample sequence depth
# has a strong effect on alpha diversity, but is also a useful method to normalize the data (Weis et al 2015 & 2017, Callahan et al 2016 and McKnight et al 2018)

# data used: noncontam_physeq_clean_noNC1
FP_data <- noncontam_physeq_clean_noNC1
range(taxa_sums(FP_data))

# calculate alpha-diversity at various rarefaction levels
set.seed(42)

calculate_rarefaction_curves <- function(physeq, measures, depths) {

  estimate_rarified_richness <- function(physeq, measures, depth) {
    if(max(sample_sums(physeq)) < depth) return()
    physeq <- prune_samples(sample_sums(physeq) >= depth, physeq)

    physeq_rar <- rarefy_even_depth(physeq, depth, verbose = FALSE)

    alpha_diversity <- estimate_richness(physeq_rar, measures = measures)

    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('X.SampleID', 'Measure'), value.name = 'Alpha_diversity')

    molten_alpha_diversity
  }
  names(depths) <- depths
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, physeq = physeq, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))

  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]

  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(FP_data, c('Observed', 'Shannon'),
                                                       rep(c(1,10,100,500,1000, 1:100 * 1000, 1:1000 * 10000), each = 10))
# This calculation yields error messages due to the fact that the data have no singletons.
# This is caused by using DADA2 in QIIME, which removes singletons from the data.
# Hence the error message can be ignored.
summary(rarefaction_curve_data)
# summarize alpha diversity
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'X.SampleID', 'Measure'),
                                        summarise, Alpha_diversity_mean = mean(Alpha_diversity),
                                        Alpha_diversity_sd = sd(Alpha_diversity))
# add sample data
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary,
                                                data.frame(sample_data(FP_data)),
                                                by.x = 'X.SampleID', by.y = 'row.names')

# make plots per season & sex
theme_set(theme_bw())
pseason <- ggplot(data = rarefaction_curve_data_summary_verbose,
                  mapping = aes(x = Depth, y = Alpha_diversity_mean,
                                ymin = Alpha_diversity_mean - 2*Alpha_diversity_sd,
                                ymax = Alpha_diversity_mean + 2*Alpha_diversity_sd,
                                colour = Season,
                                group = X.SampleID)) +
  geom_point(size=0.25) +
  geom_line() +
  geom_pointrange() +
  geom_dl(aes(label = X.SampleID), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.5)) +
  scale_x_continuous(limits=c(0,15540),
                     breaks=c(0,2000,4000,6000,8000,10000,12000,14000,16000))+
  scale_color_manual(values = colsSeason,
                     labels=c("summer", "winter"))+
  facet_wrap(facets = ~Measure, scales = 'free_y') +
  theme(legend.position="none",
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=12),
        strip.text = element_text(size=10, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_rect(fill = "lightgrey",colour="lightgrey"))
pseason

summary(rarefaction_curve_data_summary_verbose)
rarefaction_curve_data_summary_verbose_1 <- as.data.frame(rarefaction_curve_data_summary_verbose)
write.csv(rarefaction_curve_data_summary_verbose_1, file='rarefaction_curve_data_summary_verbose_1.csv')

###########_____________________________________________________________________
# Apart from FPS62 and FPS73 , Richness levels off around 3000 reads.
# The Shannon index levels off around 2000 reads.
# When inspecting the number of reads per sample, the minimum number of reads per sample is 411 (FPS28).and then following 1110 (FPS24), 1306 (FPS63), 2140 (FPS33)
# Four samples (FPS28, FPS24, FPS63, FPS33) that will be lost when rarefying to the number of reads of
# the fifth lowest reads sample (FPS04, 3295 read counts), which matches the leveling off of Richness.
# Hence the data were rarefied to 3295 reads, as at this point the Richness rarefaction curves level off,
# and four samples were lost.

## RAREFY
# Data are rarefied at 3295 reads
# data used: Phyloseq object noncontam_physeq_clean_noNC2
set.seed(48)
noncontam_physeq_rarefied <- rarefy_even_depth(noncontam_physeq_clean_noNC2, 3295, replace=FALSE, trimOTUs = FALSE)
noncontam_physeq_rarefied #1052 taxa and 59 samples
sum(sample_sums(noncontam_physeq_rarefied))
sample_sums(noncontam_physeq_rarefied)
range(taxa_sums(noncontam_physeq_rarefied))
range(sample_sums(noncontam_physeq_rarefied))
summary(sample_sums(noncontam_physeq_rarefied))
table(tax_table(noncontam_physeq_rarefied)[,"Phylum"], exclude=NULL)
table(tax_table(noncontam_physeq_rarefied)[,"Genus"], exclude=NULL)

# remove taxa with zero counts
noncontam_physeq_rarefied1 <-prune_taxa(taxa_sums(noncontam_physeq_rarefied)>0, noncontam_physeq_rarefied)
noncontam_physeq_rarefied1 #981 taxa and 59 samples
sum(sample_sums(noncontam_physeq_rarefied1))
sample_sums(noncontam_physeq_rarefied1)
range(taxa_sums(noncontam_physeq_rarefied1))
range(sample_sums(noncontam_physeq_rarefied1))
summary(sample_sums(noncontam_physeq_rarefied1))
table(tax_table(noncontam_physeq_rarefied1)[,"Phylum"], exclude=NULL)
table(tax_table(noncontam_physeq_rarefied1)[,"Genus"], exclude=NULL)
# Note that after removal of the taxa with zero counts, there are 1001 taxa left,
# so 51 AVSs were excluded from the data after rarefying.

# Make taxa and ASV table en write to file
taxa_noncontam_physeq_rarefied1 = as(tax_table(noncontam_physeq_rarefied1), "matrix")
write.csv(taxa_noncontam_physeq_rarefied1, file='taxa_noncontam_physeq_rarefied1.csv')
ASV_noncontam_physeq_rarefied1 = as(otu_table(noncontam_physeq_rarefied1), "matrix")
write.csv(ASV_noncontam_physeq_rarefied1, file='ASV_noncontam_physeq_rarefied1.csv')
metadata_noncontam_physeq_rarefied1 <- data.frame(sample_data(noncontam_physeq_rarefied1))
write.csv(metadata_noncontam_physeq_rarefied1, "metadata_noncontam_physeq_rarefied1.csv")
tree_noncontam_physeq_rarefied1 <- write.tree(phy_tree(noncontam_physeq_rarefied1))
write(tree_noncontam_physeq_rarefied1, file = "tree_noncontam_physeq_rarefied1.newick")

#save phyloseq object
saveRDS(noncontam_physeq_rarefied1, "noncontam_physeq_rarefied1.rds")
# load phyloseq object
#noncontam_physeq_rarefied1 <- readRDS("noncontam_physeq_rarefied1.rds")


## INSPECT AGAIN THE DATA FOR RARE ASVs AND LOW PREVALENCE ASVs
# data used phyloseq object: noncontam_physeq_rarefied1

# check effect of removing singletons
noncontam_physeq_rarefied2 <-prune_taxa(taxa_sums(noncontam_physeq_rarefied1)>1, noncontam_physeq_rarefied1)
noncontam_physeq_rarefied2 #926 taxa and 59 samples
sum(sample_sums(noncontam_physeq_rarefied2))
sample_sums(noncontam_physeq_rarefied2)
range(taxa_sums(noncontam_physeq_rarefied2))
range(sample_sums(noncontam_physeq_rarefied2))
summary(sample_sums(noncontam_physeq_rarefied2))
sample_names(noncontam_physeq_rarefied2)

# Make taxa and ASV table en write to file
taxa_noncontam_physeq_rarefied2 = as(tax_table(noncontam_physeq_rarefied2), "matrix")
write.csv(taxa_noncontam_physeq_rarefied2, file='taxa_noncontam_physeq_rarefied2.csv')
ASV_noncontam_physeq_rarefied2 = as(otu_table(noncontam_physeq_rarefied2), "matrix")
write.csv(ASV_noncontam_physeq_rarefied2, file='ASV_noncontam_physeq_rarefied2.csv')
metadata_noncontam_physeq_rarefied2 <- data.frame(sample_data(noncontam_physeq_rarefied2))
write.csv(metadata_noncontam_physeq_rarefied2, "metadata_noncontam_physeq_rarefied2.csv")
tree_noncontam_physeq_rarefied2 <- write.tree(phy_tree(noncontam_physeq_rarefied2))
write(tree_noncontam_physeq_rarefied2, file = "tree_noncontam_physeq_rarefied2.newick")

#save phyloseq object
saveRDS(noncontam_physeq_rarefied2, "noncontam_physeq_rarefied2.rds")
# load phyloseq object
#noncontam_physeq_rarefied2 <- readRDS("noncontam_physeq_rarefied2.rds")


FP_prevelance_rar1_df = apply(X = otu_table(noncontam_physeq_rarefied1),
                              MARGIN = 1,
                              FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
FP_prevelance_rar1_df = data.frame(Prevalence = FP_prevelance_rar1_df,
                                   TotalAbundance = taxa_sums(noncontam_physeq_rarefied1),
                                   tax_table(noncontam_physeq_rarefied1))
FP_prevelance_rar1_df[1:10,]
write.table(FP_prevelance_rar1_df, file = "FP_prevelance_rar1_df.csv")

# Next plot prevalence (and abundance) per Phylum to gain more insight into the data
plyr::ddply(FP_prevelance_rar1_df, "Phylum", function(df1)(cbind(mean(df1$Prevalence), sum(df1$Prevalence), mean(df1$TotalAbundance), sum(df1$TotalAbundance))))
# Make a prevalence plot per phylum with indicator line at 5%.
ggplot(FP_prevelance_rar1_df, aes(TotalAbundance, Prevalence / nsamples(noncontam_physeq_rarefied1),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum, nrow = 5) + theme(legend.position="none")


###########________________________________________________________________________
## ALPHA DIVERSITY - CALCULATING RICHNESS, SHANNON & FAITH's PHYLOGENETIC DIVERSITY
# data used: noncontam_physeq_rarefied1

# Calculate alpha-diversities
# Estimate richness on all sample types
alpha_FP_rarefied1 <- estimate_richness(noncontam_physeq_rarefied1)
# link to mapping file
alpha_FP_rarefied1$SampleID <- rownames(alpha_FP_rarefied1)
metadata_FP = read.csv("metadata_FP.csv", sep="", row.names=1)

rownames(metadata_FP)
rownames(alpha_FP_rarefied1)

data_alpha_FP_rarefied1 <- merge(metadata_FP, alpha_FP_rarefied1, by = "row.names", all.x = TRUE)
data_alpha_FP_rarefied1 <- data_alpha_FP_rarefied1[, -1]

## FPS018 NO ID and body mass should be remove from data_alpha_HP_rarefied1
data_alpha_FP_rarefied1 <- subset(data_alpha_FP_rarefied1, Group == "True sample")
data_alpha_FP_rarefied2 <- subset(data_alpha_FP_rarefied1, SampleID != "FPS018")
data_alpha_FP_rarefied2 <- data_alpha_FP_rarefied2[, c("SampleID", setdiff(names(data_alpha_FP_rarefied2), "SampleID"))]
data_alpha_FP_rarefied2
write.table(data_alpha_FP_rarefied2, file = "data_alpha_FP_rarefied2.csv")
write.xlsx(data_alpha_FP_rarefied2, file = file.path("D:/pigeon/seq/01052024/phyloseq", "data_alpha_FP_rarefied2.xlsx"))

## TEST FOR SEASONAL VARIATION IN THE HOST PARAMETERS
# data used: data_alpha_FP_rarefied2
# Convert BodyMass to numeric

data_alpha_FP_rarefied2 <- data_alpha_FP_rarefied2[, c("SampleID", setdiff(names(data_alpha_FP_rarefied2), "SampleID"))]
data_alpha_FP_rarefied2$BodyMass <- as.numeric(gsub(",", ".", data_alpha_FP_rarefied2$BodyMass))
data_alpha_FP_rarefied2$Wing <- as.numeric(gsub(",", ".", data_alpha_FP_rarefied2$Wing))
data_alpha_FP_rarefied2$HeadBill <- as.numeric(gsub(",", ".", data_alpha_FP_rarefied2$HeadBill))
# Convert Daylength to numeric (The format is in minutes)
data_alpha_FP_rarefied2$Daylength <- as.numeric(gsub(",", ".", data_alpha_FP_rarefied2$Daylength))
# Check the structure of the data again
str(data_alpha_FP_rarefied2)

# Convert Season and Sex to factors
data_alpha_FP_rarefied2$Season <- as.factor(data_alpha_FP_rarefied2$Season)
data_alpha_FP_rarefied2$Sex <- as.factor(data_alpha_FP_rarefied2$Sex)
data_alpha_FP_rarefied2$Age <- as.factor(data_alpha_FP_rarefied2$Age)
data_alpha_FP_rarefied2$Location <- as.factor(data_alpha_FP_rarefied2$Location)
# Check the levels again
levels(data_alpha_FP_rarefied2$Season)
levels(data_alpha_FP_rarefied2$Sex)
levels(data_alpha_FP_rarefied2$Age)
levels(data_alpha_FP_rarefied2$Location)

# BODY MASS
BM.1 =lme(BodyMass ~ Season + Sex + Season:Sex, random= ~1 | BirdID, na.action=na.omit,
          data=data_alpha_FP_rarefied2)
# Inspect the results
anova(BM.1, type = "marginal")
summary(BM.1)
# inspect the residuals: Check for linearity, homoscedaticity and normality by plotting residuals, a histogram and a Q-plot:
plot(fitted(BM.1),residuals(BM.1))
qqnorm(residuals(BM.1))
res=residuals(BM.1)
plotNormalHistogram(res)
shapiro.test(data_alpha_FP_rarefied2$BodyMass)
shapiro.test(resid(BM.1))
# Body mass varied between season and sex.
# Males were heavier in winter, in females this seasonal variation is less clear.

#Body mass with no Season:Sex
BM.2 =lme(BodyMass ~ Season + Sex , random= ~1 | BirdID, na.action=na.omit,
           data=data_alpha_FP_rarefied2)
# Inspect the results
anova(BM.2, type = "marginal")
summary(BM.2)
anova(BM.2, BM.2)

library(emmeans)
em1 <- emmeans(BM.2, ~Season)
em1
em2 <- emmeans(BM.2, ~Sex)
em2

#Wing length
BM.3 =lme(Wing ~ Season + Sex + Season:Sex, random= ~1 |BirdID, na.action=na.omit,
          data=data_alpha_FP_rarefied2)
# Inspect the results
anova(BM.3, type = "marginal")
summary(BM.3)

BM.3a =lme(Wing ~  Sex , random= ~1 | BirdID, na.action=na.omit,
           data=data_alpha_FP_rarefied2)
anova(BM.3a, type = "marginal")
em3 <- emmeans(BM.3a, ~Sex)
em3
#HeadBill
BM.4 =lme(HeadBill ~ Season + Sex + Season:Sex , random= ~1 | BirdID, na.action=na.omit,
          data=data_alpha_FP_rarefied2)
anova(BM.4, type = "marginal")
BM.4a =lme(HeadBill ~ Sex, random= ~1 | BirdID, na.action=na.omit,
          data=data_alpha_FP_rarefied2)
anova(BM.4a, type = "marginal")
em4 <- emmeans(BM.4a, ~Sex)
em4

#####____________________________________________________________________
## ANALYSIS OF RICHNESS (Observed)
# data used: data_alpha_FP_rarefied2
fit<-glm(Observed ~ BirdID,data=data_alpha_FP_rarefied2)
summary(fit)
anova(fit, test = "F")

fit<-glm(Observed ~ Season+Location+Season*Location,data=data_alpha_FP_rarefied2)
anova(fit, test = "F")
fit<-glm(Observed ~ Season,data=data_alpha_FP_rarefied2)
summary(fit)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)

emm <- emmeans(fit, ~ Season*Location)
emm_df <- as.data.frame(emm)
pairs(emm)

#violin plot
Observed <- ggplot() +
  geom_violin(
    data = data_alpha_FP_rarefied2,
    aes(x = Location, y = Observed, fill = Season),
    alpha = 0.3, trim = FALSE,
    position = position_dodge(width = 0.9)
  ) +
  geom_jitter(
    data = data_alpha_FP_rarefied2,
    aes(x = Location, y = Observed, color = Season),
    size = 1.8, alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9)
  ) +
  geom_point(
    data = emm_df,
    aes(x = Location, y = emmean, color = Season),
    size = 3,
    position = position_dodge(width = 0.9)
  ) +
  geom_errorbar(
    data = emm_df,
    aes(x = Location, y = emmean, color = Season,
        ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2, linewidth = 0.6,
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = colsSeason, labels = c("Summer", "Winter")) +
  scale_color_manual(values = colsSeason, labels = c("Summer", "Winter")) +
  scale_x_discrete(labels = c("Molukkenpad","Noorderplantsoen","Vismarkt"))+
  theme_classic(base_size = 18) +
  labs(y = "Richness", x = "Location", fill = "Season", color = "Season") +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black")
  ) +
  annotate("text", x = 2.8, y = 200,
           label = "Season: P = 0.01\nSeason * Location: P = 0.04", size = 5)

Observed


ymax_v <- max(data_alpha_FP_rarefied2$Observed[data_alpha_FP_rarefied2$Location == "VISMARKT"], na.rm = TRUE)
y_bracket <- ymax_v + 30
h <- 5

#
x_vis <- 3
dx <- 0.9/4

Observed <- Observed +
  annotate("segment",
           x = x_vis - dx, xend = x_vis + dx,
           y = y_bracket, yend = y_bracket, linewidth = 0.7) +
  annotate("segment",
           x = x_vis - dx, xend = x_vis - dx,
           y = y_bracket - h, yend = y_bracket, linewidth = 0.7) +
  annotate("segment",
           x = x_vis + dx, xend = x_vis + dx,
           y = y_bracket - h, yend = y_bracket, linewidth = 0.7) +
  annotate("text",
           x = x_vis, y = y_bracket + 6,
           label = "*", size = 5)

Observed



## SHANNON INDEX
# data used: data_alpha_FP_rarefied2
fit<-glm(Shannon ~ BirdID,data=data_alpha_FP_rarefied2)
summary(fit)
anova(fit, test = "F")

fit<-glm(Shannon ~ Season+Location+Season*Location,data=data_alpha_FP_rarefied2)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)

## FAITH'S PHYLOGENETIC DIVERSITY
# data used: noncontam_physeq_rarefied1
# https://rdrr.io/github/twbattaglia/btools/src/R/estimate_pd.R
# https://rdrr.io/github/twbattaglia/btools/f/README.md

library(btools)
# calculate Faiths phylogenetic differences (PD) and Species richness (SR) and put them in a table
PDtable_FP_S <- estimate_pd(noncontam_physeq_rarefied1)
PDtable_FP_S
# write OTU table to file for picrust analysis
#write.table(PDtable_FP_S, file='PDtable_FP_S.csv')

mapping_FP_noNC <- subset(metadata_FP, Group == "True sample")
mapping_FP_noNC

# link to mapping file
PDtable_FP_S1 <- merge(mapping_FP_noNC, PDtable_FP_S, by = "row.names", all.x = TRUE)
PDtable_FP_S1 <- PDtable_FP_S1 %>% rename_at( 1 , ~'SampleID')
PDtable_FP_S2<-copy(PDtable_FP_S1)

col1 <- PDtable_FP_S1[, 1]
rownames(PDtable_FP_S1) <- col1

PDtable_FP_S1
rownames(mapping_FP_noNC)
rownames(PDtable_FP_S1)

write.table(PDtable_FP_S2, file = "PDtable_FP_S2.csv")

PDtable_FP_S1$Sex <- mapping_FP_noNC$Sex
PDtable_FP_S1$Season <- mapping_FP_noNC$Season
PDtable_FP_S1$Year <- mapping_FP_noNC$Year
PDtable_FP_S1$Location <- mapping_FP_noNC$Location
PDtable_FP_S1$BirdID <- mapping_FP_noNC$BirdID
PDtable_FP_S1

PDtable_FP_S1 <- subset(PDtable_FP_S1, Group == "True sample")
PDtable_FP_S1 <- subset(PDtable_FP_S1, SampleID != "FPS018")
PDtable_FP_S1

# GLMs: PD, same full models as for Observed and Shannon
# data used: PDtable_FP_S1
fit<-glm(PD ~ BirdID,data=PDtable_FP_S1)
summary(fit)
anova(fit, test = "F")

fit<-glm(PD ~ Season+Location+Season*Location,data=PDtable_FP_S1)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)

######### BETA -DIVERSITY
# The beta-diversity analysis is using rarefied data
# Weiss et al 2015 & 2017, Callahan et al 2016 and McKnight et al 2018
# consider rarefied or proportional data the best data to use for beta-diversity analysis.

# JACCARD DISTANCES (presence-absence)
FP_rarefied1_jacc <- phyloseq::distance(noncontam_physeq_rarefied2, "jaccard")
df_jacc_rar1 <- as(sample_data(noncontam_physeq_rarefied2), "data.frame")
head(df_jacc_rar1)
# test for differences between seasons and Locations.
# Likewise BirdID will give even more problems
# Note that since ADONIS is a sequential analysis,
# the model should be run with all potential orders of the vars.
# # see also https://rdrr.io/github/vegandevs/vegan/man/adonis.html
# Following Baniel et al 2020 Microbiome: First test for individual effect, then include this effect in ADONIS via Strata function
# Effect of individuals
set.seed(22)
adonis2(FP_rarefied1_jacc ~ BirdID, df_jacc_rar1, permutations = 999)
# BirdID is not significant

set.seed(22)
adonis2(FP_rarefied1_jacc ~ Season + Location + Season*Location, df_jacc_rar1,by = "terms",permutations = 999)
set.seed(22)
adonis2(FP_rarefied1_jacc ~ Season + Location , df_jacc_rar1,by = "terms",permutations = 999)

# use Season as factor for betadisper
beta_jacc3 <- betadisper(FP_rarefied1_jacc, df_jacc_rar1$Season)
permutest(beta_jacc3)
# betadisper was not significant, thus differences are not due to differences in group dispersions

# use Location as factor for betadisper
beta_jacc3 <- betadisper(FP_rarefied1_jacc, df_jacc_rar1$Location)
permutest(beta_jacc3)
# betadisper was not significant, thus differences are not due to differences in group dispersions

# Posthoc analysis
# cite package as Martinez Arbizu, P. (2020). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.4
# note that the PairWiseAdonais does not allow for interaction terms
pairwise.adonis(FP_rarefied1_jacc,factors=df_jacc_rar1$Season)
pairwise.adonis(FP_rarefied1_jacc,factors=df_jacc_rar1$Location)


# Make plot with mean +/- inter-quantile range added to plot
# make a seperate panel per sex
otu_rar1 <- t(otu_table(noncontam_physeq_rarefied2))
set.seed(22)
jacc_rar1 <- ordinate(otu_rar1, method="PCoA", distance ="jaccard")
# plot the PCoA scores to obtain axe-labels for final plot
Jacc_pcoa2b <- plot_ordination(noncontam_physeq_rarefied2,jacc_rar1, type="samples", color="Season")+
  theme_bw() +
  geom_point(size=5) +
  scale_color_manual(values=colsSeason,labels = c("Summer", "Winter"))+
  facet_wrap(~Location) +
  theme(legend.position.inside = c(.1, .1),
        legend.background = element_rect(),
        legend.text = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Jacc_pcoa2b
# Yields: Y-axe: [6.5%], X-axe [12.3%]; use this for the multi-panel plot below

# create data from of PCoA scores along the first two axes
pjacc_rar1_df = plot_ordination(justDF=TRUE, noncontam_physeq_rarefied2, jacc_rar1,
                                axes= c(1,2), type="samples", color="Season")
pjacc_rar1_df2 = plot_ordination(justDF=TRUE, noncontam_physeq_rarefied2, jacc_rar1,
                                 axes= c(1,2), type="samples", color="Location")
# write jaccard axis data to file
write.csv(pjacc_rar1_df, file="pjacc_rar1_df.csv")
# create a summary of PCoA scores to calculate median, quantiles and interquantile range
# summary of PCoA df for plotting median and IQR across experimental groups (Season)
summary_jacc_pcoa_rar1 <- summaryBy(Axis.1 + Axis.2 ~ Season,
                                    data = pjacc_rar1_df,
                                    FUN = function(x) { c(m = median(x), quant_25 = quantile(x, 0.25), quant_75 = quantile(x, 0.75),
                                                          s = sd(x), n=length(x), sem= sd(x)/sqrt(length(x)), n=length(x))} )
colnames(summary_jacc_pcoa_rar1)[3] <- "Axis.1.quant_25"
colnames(summary_jacc_pcoa_rar1)[4] <- "Axis.1.quant_75"
colnames(summary_jacc_pcoa_rar1)[10] <- "Axis.2.quant_25"
colnames(summary_jacc_pcoa_rar1)[11] <- "Axis.2.quant_75"
colnames(summary_jacc_pcoa_rar1)[2] <- "Axis.1"
colnames(summary_jacc_pcoa_rar1)[9] <- "Axis.2"

# calculate IQR
summary_jacc_pcoa_rar1[16] <- with(summary_jacc_pcoa_rar1, Axis.1 - Axis.1.quant_25)
summary_jacc_pcoa_rar1[17] <- with(summary_jacc_pcoa_rar1, Axis.1.quant_75 - Axis.1 )
summary_jacc_pcoa_rar1[18] <- with(summary_jacc_pcoa_rar1, Axis.2 - Axis.2.quant_25)
summary_jacc_pcoa_rar1[19] <- with(summary_jacc_pcoa_rar1, Axis.2.quant_75 - Axis.2 )
colnames(summary_jacc_pcoa_rar1)[16] <- "ax1_25"
colnames(summary_jacc_pcoa_rar1)[17] <- "ax1_75"
colnames(summary_jacc_pcoa_rar1)[18] <- "ax2_25"
colnames(summary_jacc_pcoa_rar1)[19] <- "ax2_75"
summary_jacc_pcoa_rar1
# make plot
p <- ggplot(data = pjacc_rar1_df, aes(x = Axis.1, y = Axis.2))
jacc <- p + theme_bw()+
  geom_point(aes(colour = Season), size = 3, alpha=0.5) +
  geom_point(data = summary_jacc_pcoa_rar1,
             aes(x = Axis.1, y = Axis.2, colour = Season), size = 4) +
  geom_errorbarh(data=summary_jacc_pcoa_rar1,
                 aes(xmin = Axis.1 - ax1_25, xmax = Axis.1 + ax1_75, colour = Season), height=0) +
  geom_errorbar(data=summary_jacc_pcoa_rar1,
                aes(ymin = Axis.2 - ax2_25, ymax = Axis.2 + ax2_75, colour = Season), width=0) +
  scale_color_manual(values=colsSeason,"",
                     labels=c(" Summer", " Winter"))+
  facet_wrap(~Location, labeller=labeller(Location= Local_labels)) +
  xlab("Axis 1 [12.3%]") +
  ylab("Axis 2 [6.5%]") +
  theme(legend.position.inside = c(0.9, 0.92),
        legend.text = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour = "black"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  coord_fixed(ratio = 1)
jacc

ggsave(
  filename = "jacc.pdf",
  plot = jacc,
  width = 8,
  height = 2.5,
  units = "in")

# BRAY-CURTIS (relative abundance/presence-absence)
# data used: noncontam_physeq_rarefied2

# Stats calculate Bray Curtis
FP_rarefied1_bray <- distance(noncontam_physeq_rarefied2, "bray")  # vegdist bray
df_bray_rar1 <- as(sample_data(noncontam_physeq_rarefied2), "data.frame")
head(df_bray_rar1)

# Following Baniel et al 2020 Microbiome: First test for individual effect, then include this effect in ADONIS via Strata function
#Effect of individuals
set.seed(22)
adonis2(FP_rarefied1_bray ~ BirdID, df_bray_rar1, permutations = 999)
# BirdID is not significant

set.seed(22)
adonis2(FP_rarefied1_bray ~ Season + Location + Season*Location, data=df_bray_rar1, by = "terms",permutations = 999)
# Location, Season both significant
set.seed(22)
adonis2(FP_rarefied1_bray ~ Season, data=df_bray_rar1, by = "terms",permutations = 999)


#betadisper
beta_bray3 <- betadisper(FP_rarefied1_bray, df_bray_rar1$Location)
permutest(beta_bray3)
# Perform test
anova(beta_bray3)
# betadisper was not significant, thus differences are not due to differences in group dispersions
# Posthoc analysis on Season
pairwise.adonis(FP_rarefied1_bray,factors=df_bray_rar1$Season)

# Make plot with mean +/- inter-quantile range added to plot
# make a seperate panel per sex
otu_rar1 <- t(otu_table(noncontam_physeq_rarefied2))
set.seed(22)
bray_rar1 <- ordinate(otu_rar1, method="PCoA", distance ="bray")

theme_set(theme_bw())
# plot the PCoA scores to obtain axe-labels for final plot
BC_pcoa2b <- plot_ordination(noncontam_physeq_rarefied2,bray_rar1, type="samples", color="Season")+
  theme_bw() +
  geom_point(size=5) +
  scale_color_manual(values=colsSeason)+
  facet_wrap(~Location) +
  theme(legend.position.inside = c(.1, .1),
        legend.background = element_rect(),
        legend.text = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
BC_pcoa2b
# Y-axe: [8.8%], X-axe [19.8%]

# create data from of PCoA scores along the first two axes
pbray_rar1_df = plot_ordination(justDF=TRUE, noncontam_physeq_rarefied2, bray_rar1,
                                axes= c(1,2), type="samples", color="Season")
write.csv(pbray_rar1_df, file="pbray_rar1_df.csv")

# create a summary of PCoA scores to calculate median, quantiles and interquantile range
# summary of PCoA df for plotting median and IQR across experimental groups (Season)
summary_bray_pcoa_rar1 <- summaryBy(Axis.1 + Axis.2 ~ Season,
                                    data = pbray_rar1_df,
                                    FUN = function(x) { c(m = median(x), quant_25 = quantile(x, 0.25), quant_75 = quantile(x, 0.75),
                                                          s = sd(x), n=length(x), sem= sd(x)/sqrt(length(x)), n=length(x))} )
colnames(summary_bray_pcoa_rar1)[3] <- "Axis.1.quant_25"
colnames(summary_bray_pcoa_rar1)[4] <- "Axis.1.quant_75"
colnames(summary_bray_pcoa_rar1)[10] <- "Axis.2.quant_25"
colnames(summary_bray_pcoa_rar1)[11] <- "Axis.2.quant_75"
colnames(summary_bray_pcoa_rar1)[2] <- "Axis.1"
colnames(summary_bray_pcoa_rar1)[9] <- "Axis.2"

# calculate IQR
summary_bray_pcoa_rar1[16] <- with(summary_bray_pcoa_rar1, Axis.1 - Axis.1.quant_25)
summary_bray_pcoa_rar1[17] <- with(summary_bray_pcoa_rar1, Axis.1.quant_75 - Axis.1 )
summary_bray_pcoa_rar1[18] <- with(summary_bray_pcoa_rar1, Axis.2 - Axis.2.quant_25)
summary_bray_pcoa_rar1[19] <- with(summary_bray_pcoa_rar1, Axis.2.quant_75 - Axis.2 )
colnames(summary_bray_pcoa_rar1)[16] <- "ax1_25"
colnames(summary_bray_pcoa_rar1)[17] <- "ax1_75"
colnames(summary_bray_pcoa_rar1)[18] <- "ax2_25"
colnames(summary_bray_pcoa_rar1)[19] <- "ax2_75"
summary_bray_pcoa_rar1

# make plot
p <- ggplot(data = pbray_rar1_df, aes(x = Axis.1, y = Axis.2))
bc <- p + theme_bw()+
  geom_point(aes(colour = Season), size = 3, alpha=0.5) +
  geom_point(data = summary_bray_pcoa_rar1,
             aes(x = Axis.1, y = Axis.2, colour = Season), size = 4) +
  geom_errorbarh(data=summary_bray_pcoa_rar1,
                 aes(xmin = Axis.1 - ax1_25, xmax = Axis.1 + ax1_75, colour = Season), height=0) +
  geom_errorbar(data=summary_bray_pcoa_rar1,
                aes(ymin = Axis.2 - ax2_25, ymax = Axis.2 + ax2_75, colour = Season), width=0) +
  xlab("Axis 1 [19.8%]") +
  ylab("Axis 2 [8.8%]") +
  scale_color_manual(values=colsSeason,"",
                     labels=c(" Summer", " Winter"))+
  #facet_wrap(~Location, labeller=labeller(Location= Local_labels)) +
  theme(legend.position="right",
        legend.text = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour = "black"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
bc


# UNWEIGHTED UNIFRAC (presence-absence with phylogenetic relationships taken into account)
# data used: noncontam_physeq_rarefied2

# calculate unweighted unifracs - PCoA
unweighted_unif_rar1 <- UniFrac(noncontam_physeq_rarefied2, weighted = FALSE)
uwunif_rar1 <- ordinate(noncontam_physeq_rarefied2, method="PCoA", distance ="unifrac")
uwunif_df.rar1 <- as(sample_data(noncontam_physeq_rarefied2), "data.frame")

# Following Baniel et al 2020 Microbiome: First test for individual effect, then include this effect in ADONIS via Strata function
#Effect of individuals
set.seed(22)
adonis2(unweighted_unif_rar1 ~ BirdID, uwunif_df.rar1, permutations = 999)
# BirdID is significant

# use blocks of individuals in the ADONIS:
perm <- how(nperm = 999)
setBlocks(perm) <- with(uwunif_df.rar1, BirdID)
set.seed(22)
adonis2(unweighted_unif_rar1 ~ Season + Location + Season*Location, data=uwunif_df.rar1,by= "term", permutations = perm)

# use Location as factor for betadisperser
beta_uwunif2s <- betadisper(unweighted_unif_rar1, uwunif_df.rar1$Location)
permutest(beta_uwunif2s)
# Perform test
anova(beta_uwunif2s)
# Permutation test for F
permutest(beta_uwunif2s, pairwise = TRUE)
# Tukey's Honest Significant Differences
(beta_uwunif2s.HSD <- TukeyHSD(beta_uwunif2s))
plot(beta_uwunif2s.HSD)
# plot betadispersions
plot(beta_uwunif2s, cex.lab=1, label= FALSE , cex = 1, col=Local_cols)
# betadisper was not significant

# use Season as grouping factor for betadisperser
beta_uwunif3 <- betadisper(unweighted_unif_rar1, uwunif_df.rar1$Season)
permutest(beta_uwunif3)
# Perform test
anova(beta_uwunif3)
# Permutation test for F
permutest(beta_uwunif3, pairwise = TRUE)
# Tukey's Honest Significant Differences
(beta_uwunif3.HSD <- TukeyHSD(beta_uwunif3))
plot(beta_uwunif3.HSD)
# plot betadispersions
plot(beta_uwunif3, cex.lab=1, label= FALSE , cex = 1, col=colsSeason)
# betadisper was not significant

# Make plot with mean +/- inter-quantile range added to plot
# Create new distance matrix with normalized = TRUE
uwunif_rar1 <- UniFrac(noncontam_physeq_rarefied2, weighted = FALSE, normalized = TRUE)
set.seed(22)
pcoa_uwunif_rar1 <- ordinate(noncontam_physeq_rarefied2, method="PCoA", distance ="unifrac")

# plot PCoA scores to obtain axe-labels for final plot
uwunif_pcoa2b <- plot_ordination(noncontam_physeq_rarefied2, pcoa_uwunif_rar1, type="samples", color="Season")+
  theme_bw() +
  geom_point(size=5) +
  #labs(title = "Bray-Curtis Deseq2 vsd2 normalized data")+
  #annotate(geom="text", x = 0.5, y = 0.5, label = "A")+
  scale_color_manual(values=colsSeason)+
  facet_wrap(~Location) +
  theme(legend.position="right",
        legend.background = element_rect(),
        legend.text = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
uwunif_pcoa2b
# Y-axe: [16.1%], X-axe [41.3%]

# create data from of PCoA scores along the first two axes
pcoa_scores_uwunif_rar1_df <- plot_ordination(justDF = TRUE, noncontam_physeq_rarefied2,
                                              pcoa_uwunif_rar1, axes = c(1,2), type="samples",
                                              color = "Season")
write.csv(pcoa_scores_uwunif_rar1_df, file="pcoa_scores_uwunif_rar1_df.csv")

# create a summary of PCoA scores to calculate median, quantiles and interquantile range
# summary of PCoA df for plotting median and IQR across experimental groups (Season)
summary_pcoa_uwunif_rar1 <- summaryBy(Axis.1 + Axis.2 ~ Season, data = pcoa_scores_uwunif_rar1_df,
                                      FUN = function(x) { c(m = median(x), quant_25 = quantile(x, 0.25), quant_75 = quantile(x, 0.75),
                                                            s = sd(x), n=length(x), sem= sd(x)/sqrt(length(x)), n=length(x))} )
colnames(summary_pcoa_uwunif_rar1)[3] <- "Axis.1.quant_25"
colnames(summary_pcoa_uwunif_rar1)[4] <- "Axis.1.quant_75"
colnames(summary_pcoa_uwunif_rar1)[10] <- "Axis.2.quant_25"
colnames(summary_pcoa_uwunif_rar1)[11] <- "Axis.2.quant_75"
colnames(summary_pcoa_uwunif_rar1)[2] <- "Axis.1"
colnames(summary_pcoa_uwunif_rar1)[9] <- "Axis.2"

# calculate IQR
str(summary_pcoa_uwunif_rar1)
summary_pcoa_uwunif_rar1[16] <- with(summary_pcoa_uwunif_rar1, Axis.1 - Axis.1.quant_25)
summary_pcoa_uwunif_rar1[17] <- with(summary_pcoa_uwunif_rar1, Axis.1.quant_75 - Axis.1 )
summary_pcoa_uwunif_rar1[18] <- with(summary_pcoa_uwunif_rar1, Axis.2 - Axis.2.quant_25)
summary_pcoa_uwunif_rar1[19] <- with(summary_pcoa_uwunif_rar1, Axis.2.quant_75 - Axis.2 )
colnames(summary_pcoa_uwunif_rar1)[16] <- "ax1_25"
colnames(summary_pcoa_uwunif_rar1)[17] <- "ax1_75"
colnames(summary_pcoa_uwunif_rar1)[18] <- "ax2_25"
colnames(summary_pcoa_uwunif_rar1)[19] <- "ax2_75"
summary_pcoa_uwunif_rar1

p <- ggplot(data = pcoa_scores_uwunif_rar1_df, aes(x = Axis.1, y = Axis.2))
uwunif <- p + theme_bw()+
  geom_point(aes(colour =Season), size = 3, alpha=0.5) +
  geom_point(data = summary_pcoa_uwunif_rar1,
             aes(x = Axis.1, y = Axis.2, colour = Season), size = 4) +
  geom_errorbarh(data=summary_pcoa_uwunif_rar1,
                 aes(xmin = Axis.1 - ax1_25, xmax = Axis.1 + ax1_75, colour = Season), height=0) +
  geom_errorbar(data=summary_pcoa_uwunif_rar1,
                aes(ymin = Axis.2 - ax2_25, ymax = Axis.2 + ax2_75, colour = Season), width=0) +
  xlab("Axis 1 [41.3%]") +
  ylab("Axis 2 [16.1%]") +
  scale_color_manual(values=colsSeason,"",
                     labels=c(" Summer", " Winter"))+
  facet_wrap(~Location, labeller=labeller(Location= Local_labels)) +
  theme(legend.position="right",
        legend.text = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour = "black"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
uwunif


# WEIGHTED UNIFRAC
# data used: noncontam_physeq_rarefied2
# calculate unweighted UNifracs - PCoA
weighted_unif_rar1 <- UniFrac(noncontam_physeq_rarefied2, weighted = TRUE)
wunif_rar1 <- ordinate(noncontam_physeq_rarefied1, method="PCoA", distance ="wunifrac")
wunif_df.rar1 <- as(sample_data(noncontam_physeq_rarefied2), "data.frame")

# Following Baniel et al 2020 Microbiome: First test for individual effect, then include this effect in ADONIS via Strata function
#Effect of individuals
set.seed(22)
adonis2(weighted_unif_rar1 ~ BirdID, wunif_df.rar1, permutations = 999)
# BirdID is significant

# use blocks of individuals in the ADONIS:
perm_WU <- how(nperm = 999)
setBlocks(perm_WU) <- with(uwunif_df.rar1, BirdID)
set.seed(22)
adonis2(weighted_unif_rar1 ~ Season + Location + Season*Location, data=wunif_df.rar1,by = "term", permutations = perm_WU)

set.seed(22)
adonis2(weighted_unif_rar1 ~ Season + Location, data=wunif_df.rar1,by = "term", permutations = perm_WU)


# use Season as factor for betadisperser
beta_wunif2s <- betadisper(weighted_unif_rar1, wunif_df.rar1$Season)
permutest(beta_wunif2s)
# Perform test
anova(beta_wunif2s)
# Permutation test for F
permutest(beta_wunif2s, pairwise = TRUE)
# Tukey's Honest Significant Differences
(beta_wunif2s.HSD <- TukeyHSD(beta_wunif2s))
plot(beta_wunif2s.HSD)
# plot betadispersions
plot(beta_wunif2s, cex.lab=1, label= FALSE , cex = 1, col=colsSeason)
# betadisper was not significant

# use Location as grouping factor for betadisperser
beta_wunifL <- betadisper(weighted_unif_rar1, wunif_df.rar1$Location)
permutest(beta_wunifL)
# Perform test
anova(beta_wunifL)
# Permutation test for F
permutest(beta_wunifL, pairwise = TRUE)
# Tukey's Honest Significant Differences
(beta_wunifL.HSD <- TukeyHSD(beta_wunifL))
plot(beta_wunifL.HSD)
# plot betadispersions
plot(beta_wunifL, cex.lab=1, label= FALSE , cex = 1, col=Local_cols)
# betadisper both were not significant, thus differences are not due to differences in group dispersions

# Make plot with mean +/- inter-quantile range added to plot
# Create new distance matrix with normalized = TRUE
wunif_rar1 <- UniFrac(noncontam_physeq_rarefied2, weighted = TRUE, normalized = TRUE)
set.seed(22)
pcoa_wunif_rar1 <- ordinate(noncontam_physeq_rarefied2, method="PCoA", distance ="wunifrac")

# plot PCoA scores to obtain axe-labels for final plot
wunif_pcoa2b <- plot_ordination(noncontam_physeq_rarefied2,pcoa_wunif_rar1, type="samples", color="Season")+
  theme_bw() +
  geom_point(size=5) +
  #labs(title = "Bray-Curtis Deseq2 vsd2 normalized data")+
  #annotate(geom="text", x = 0.5, y = 0.5, label = "A")+
  scale_color_manual(values=colsSeason)+
  facet_wrap(~Location) +
  theme(legend.position="right",
        legend.background = element_rect(),
        legend.text = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
wunif_pcoa2b
# Y-axe: [12.6%], X-axe [35.7%]

# create data from of PCoA scores along the first two axes
pcoa_scores_wunif_rar1_df <- plot_ordination(justDF = TRUE, noncontam_physeq_rarefied2,
                                             pcoa_wunif_rar1, axes = c(1,2), type="samples",
                                             color = "Season")
write.csv(pcoa_scores_wunif_rar1_df, file="pcao_scores_wunif_rar1_df.csv")
write.xlsx(pcoa_scores_wunif_rar1_df, file = file.path("D:/pigeon/seq/01052024/phyloseq", "pcoa_scores_wunif_rar1_df.xlsx"),rowNames= T)



# create a summary of PCoA scores to calculate median, quantiles and interquantile range
# summary of PCoA df for plotting median and IQR across experimental groups (Season)
summary_pcoa_wunif_rar1 <- summaryBy(Axis.1 + Axis.2 ~ Season, data = pcoa_scores_wunif_rar1_df,
                                     FUN = function(x) { c(m = median(x), quant_25 = quantile(x, 0.25), quant_75 = quantile(x, 0.75),
                                                           s = sd(x), n=length(x), sem= sd(x)/sqrt(length(x)), n=length(x))} )
colnames(summary_pcoa_wunif_rar1)[3] <- "Axis.1.quant_25"
colnames(summary_pcoa_wunif_rar1)[4] <- "Axis.1.quant_75"
colnames(summary_pcoa_wunif_rar1)[10] <- "Axis.2.quant_25"
colnames(summary_pcoa_wunif_rar1)[11] <- "Axis.2.quant_75"
colnames(summary_pcoa_wunif_rar1)[2] <- "Axis.1"
colnames(summary_pcoa_wunif_rar1)[9] <- "Axis.2"

# calculate IQR
summary_pcoa_wunif_rar1[16] <- with(summary_pcoa_wunif_rar1, Axis.1 - Axis.1.quant_25)
summary_pcoa_wunif_rar1[17] <- with(summary_pcoa_wunif_rar1, Axis.1.quant_75 - Axis.1 )
summary_pcoa_wunif_rar1[18] <- with(summary_pcoa_wunif_rar1, Axis.2 - Axis.2.quant_25)
summary_pcoa_wunif_rar1[19] <- with(summary_pcoa_wunif_rar1, Axis.2.quant_75 - Axis.2 )
colnames(summary_pcoa_wunif_rar1)[16] <- "ax1_25"
colnames(summary_pcoa_wunif_rar1)[17] <- "ax1_75"
colnames(summary_pcoa_wunif_rar1)[18] <- "ax2_25"
colnames(summary_pcoa_wunif_rar1)[19] <- "ax2_75"
summary_pcoa_wunif_rar1


# make plot
p <- ggplot(data = pcoa_scores_wunif_rar1_df, aes(x = Axis.1, y = Axis.2))
wunif <- p + theme_bw()+
  geom_point(aes(colour = Season), size = 3, alpha=0.5) +
  geom_point(data = summary_pcoa_wunif_rar1,
             aes(x = Axis.1, y = Axis.2, colour = Season), size = 4) +
  geom_errorbarh(data=summary_pcoa_wunif_rar1,
                 aes(xmin = Axis.1 - ax1_25, xmax = Axis.1 + ax1_75, colour = Season), height=0) +
  geom_errorbar(data=summary_pcoa_wunif_rar1,
                aes(ymin = Axis.2 - ax2_25, ymax = Axis.2 + ax2_75, colour = Season), width=0) +
  xlab("Axis 1 [35.7%]") +
  ylab("Axis 2 [12.6%]") +
  scale_color_manual(values=colsSeason,"",
                     labels=c(" Summer", " Winter"))+
  facet_wrap(~Location, labeller=labeller(Location= Local_labels)) +
  theme(legend.position="right",
        legend.text = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=10, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour = "black"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
wunif



###################################
setwd("D:/pigeon/seq/01052024/phyloseq/phylum")
library(reshape2)
library(tidyr)
library(dplyr)

colslocal <- c("darkseagreen", "darkgreen","saddlebrown")
colsSeason <- c("#DA5724","#508578")
Sex_labels <- c("F" = "females",
                "M" = "males")
Local_labels<-c("Molukkenpad","Noorderplantsoen","Vismarkt")
#phylum
phylum <- read.csv("Phylum_abundance.csv", sep=";", row.names = 1)
head(phylum)

phylum.ave <- colMeans(phylum)

phylum.0 <- phylum[, order(-phylum.ave)]
phylum.1 <- phylum[, order(-phylum.ave)[1:10]]

phylum_per <- as.data.frame(lapply(phylum, function(x) x / sum(x)))
rownames(phylum_per)<- row.names(phylum)

phylum.2 <- phylum_per[, order(-phylum.ave)[1:10]]
phylum.2 <- cbind(phylum.2, others = 1 - rowSums(phylum.2))
phylum.2 <- cbind(SampleID = rownames(phylum.2), phylum.2)


write.csv(phylum.0,"phylum.0.csv", row.names = T)
write.csv(phylum.1,"phylum.1.csv", row.names = T)
write.csv(phylum.2,"phylum.2.csv", row.names = T)


phylum.2$SampleID <- rownames(phylum.2)
phylum.2 <- phylum.2[, c("SampleID", setdiff(names(phylum.2), "SampleID"))]


phylum_long <- phylum.2 %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "PhylumID",
    values_to = "Abundance"
  )

head(phylum_long)

write.csv(phylum_long,file = "phylum_long.csv",row.names = F)
metadata_FP = read.csv("metadata_FP.csv", sep="", row.names=1)

metadata_FP$SampleID <- rownames(metadata_FP)
metadata_FP <- metadata_FP[, c("SampleID", setdiff(names(metadata_FP), "SampleID"))]

phylum.gg<-merge(phylum_long,metadata_FP, by.x=c("SampleID"), by.y=c("SampleID"))
write.table(
  phylum.gg,
  "phylum.gg.csv",
  sep = ";",
  quote = TRUE,
  row.names = TRUE,
  fileEncoding = "UTF-8"
)

##----------------------------------------------------------------------------------
#phylum.gg <- read.table("phylum.gg.csv", sep=";", header=T, row.names = 1)
head(phylum.gg)

library(wesanderson)
library(RColorBrewer)
library(ggpubr)
library("dplyr")

phylum_reabundance<-ggbarplot(phylum.gg, x = "SampleID", y="Abundance", color="black", fill="PhylumID",
                              legend="middle",
                              legend.title="Phylum", main="Relative abundance per Phylum",
                              font.main = c(14,"bold", "black"), font.x = c(14, "bold"),
                              font.y=c(14,"bold")) +
  theme_bw() +
  rotate_x_text() +
  scale_fill_brewer(palette="Set3")+
  facet_grid(~ Season, scales = "free_x", space='free') +
  labs(x = "Samples", y = "Relative abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
phylum_reabundance
phylum_reabundance + theme(
  plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5))


##
# select for Firmicutes
phylum_Firm_df_L = subset.data.frame(phylum.gg, PhylumID=='Firmicutes')

# logit transform of proportions according to Warton & Hui Ecology 2011:
# They have proposed log(y/1-y) as a transformation to consider in the analysis of non-binomial proportions.
# They suggested dealing with sample proportions equal to exactly zero or one using an approach analogous to the
# empirical logistic transform (Collett, 2002) - by adding some small value e to the numerator and denominator
# of the logit transform, i.e., log(y+e /1-y+e) (note log = natural log). They proposed using as e the smallest non-zero proportion (or if
# counts are near one, the minimum non-zero value of 1-y).
# lowest proportion for Firmicutes is 0.02, hence e was considered 0.02.
phylum_Firm_df_L<-phylum_Firm_df_L %>%
  mutate(
    LogitProp = log((Abundance+0.02)/(1-Abundance+0.02))
  )


## FIRMICUTES
# Plot proportion Firmicutes vs Season
# boxplot season vs sex with season as colour
phylum_Firm_df_L$Season <- factor(phylum_Firm_df_L$Season)
theme_set(theme_bw())
phylum_Firm_Season_Location<-ggplot(data = phylum_Firm_df_L,
                                    aes(x=Season, y=LogitProp, fill= Season))+
  geom_boxplot(alpha=0.9) +
  scale_fill_manual(values = colsSeason) +
  facet_wrap(Location~., labeller=labeller(Location= Local_labels)) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=2, position=position_jitterdodge(0.3), color="black") +
  ylab("Proportion Firmicutes\n")+
  xlab("Season") +
  ylim(-3.2,1.8)+
  scale_y_continuous(breaks=seq(-3.2,1.8,0.2))+
  #scale_x_discrete(labels=c("summer 2013", "winter 2014", "summer 2014", "winter 2015"))+
  theme(axis.text.y = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black", angle = 90),
        strip.text = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_rect(fill = "lightgrey",colour="lightgrey"),
        legend.position = "none")
phylum_Firm_Season_Location

## GLMs TO TEST IF FIRMICUTES LOGIT TRANSFORMED PROPORTIONS VARY WITH SAME FULL MODEL AS FOR ALPHA-DOVERSITY
# data used: phylum_Firm_df_L
# see also David I. Warton and Francis K. C. Hui 2011 about using logit transformed proportions

# Model a)
str(phylum_Firm_df_L)
phylum_Firm_df_L$Season<-as.factor(phylum_Firm_df_L$Season)
phylum_Firm_df_L$Location<-as.factor(phylum_Firm_df_L$Location)
phylum_Firm_df_L$Sex<-as.factor(phylum_Firm_df_L$Sex)
phylum_Firm_df_L$Age<-as.factor(phylum_Firm_df_L$Age)
# Check the levels again
levels(phylum_Firm_df_L$Season)
levels(phylum_Firm_df_L$Sex)
levels(phylum_Firm_df_L$Age)
levels(phylum_Firm_df_L$Location)

# Check the structure of the data again
str(phylum_Firm_df_L)
summary(phylum_Firm_df_L)
class(phylum_Firm_df_L)

phylum_Firm_df_L<-as.data.frame(phylum_Firm_df_L)

fit<-glm(LogitProp ~Season*Location,data=phylum_Firm_df_L)
summary(fit)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)

stepAIC(fit,direction = "backward")


#install.packages("glmulti")
library(glmulti)
fit.model <- glmulti(fit, level = 1, crit = "aicc")
summary(fit.model)
summary(fit.model)$icvalue

##------------------------------
# select for Actinobacteria
phylum_Actin_df_L = subset.data.frame(phylum.gg, PhylumID=='Actinobacteria')

phylum_Actin_df_L<-phylum_Actin_df_L %>%
  mutate(
    LogitProp = log((Abundance+0.06)/(1-Abundance+0.06))
  )

## Actinobacteria
# Plot proportion Actinobacteria vs Season
# boxplot season vs sex with season as colour
phylum_Actin_df_L$Season <- factor(phylum_Actin_df_L$Season)
theme_set(theme_bw())
phylum_Actin_Season_Location<-ggplot(data = phylum_Actin_df_L,
                                     aes(x=Season, y=Abundance, fill= Season))+
  geom_boxplot(alpha=0.95) +
  scale_fill_manual(values = colsSeason) +
  #facet_wrap(Location~., labeller=labeller(Location= Local_labels)) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=2, position=position_jitterdodge(0.3), color="black") +
  ylab("Proportion Actinobacteria\n")+
  xlab("Season") +
  ylim(0,1)+
  scale_y_continuous()+
  #scale_x_discrete(labels=c("summer 2013", "winter 2014", "summer 2014", "winter 2015"))+
  theme(axis.text.y = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black", angle = 90),
        strip.text = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_rect(fill = "lightgrey",colour="lightgrey"),
        legend.position = "none")
phylum_Actin_Season_Location

## GLMs TO TEST IF FIRMICUTES LOGIT TRANSFORMED PROPORTIONS VARY WITH SAME FULL MODEL AS FOR ALPHA-DOVERSITY
# data used: phylum_Firm_df_L
# see also David I. Warton and Francis K. C. Hui 2011 about using logit transformed proportions

# Model a)
str(phylum_Actin_df_L)
phylum_Actin_df_L$Season<-as.factor(phylum_Actin_df_L$Season)
phylum_Actin_df_L$Location<-as.factor(phylum_Actin_df_L$Location)
phylum_Actin_df_L$Sex<-as.factor(phylum_Actin_df_L$Sex)
phylum_Actin_df_L$Age<-as.factor(phylum_Actin_df_L$Age)
# Check the levels again
levels(phylum_Actin_df_L$Season)
levels(phylum_Actin_df_L$Sex)
levels(phylum_Actin_df_L$Age)
levels(phylum_Actin_df_L$Location)

# Check the structure of the data again
str(phylum_Actin_df_L)
summary(phylum_Actin_df_L)
class(phylum_Actin_df_L)

phylum_Actin_df_L<-as.data.frame(phylum_Actin_df_L)

std <- function(x) sd(x)/sqrt(length(x))
aggregate(phylum_Actin_df_L$Abundance,by=list(phylum_Actin_df_L$Season),FUN=mean)
aggregate(phylum_Actin_df_L$Abundance,by=list(phylum_Actin_df_L$Season),FUN=std)
shapiro.test(phylum_Actin_df_L$Abundance)
qqnorm(phylum_Actin_df_L$Abundance)
qqline(phylum_Actin_df_L$Abundance)
bartlett.test(Abundance~Season,phylum_Actin_df_L)


## Generalize linear Model
ggplot(data = phylum_Actin_df_L) +
  geom_histogram(aes(x = LogitProp), bins = 30, color = 'gray30', fill = '#FFE8A2') +
  #geom_histogram(aes(x = length), binwidth = 2500, color = 'gray30', fill = '#FFE8A2') +
  theme_bw() +
  labs(x = 'LogitProp', y = 'Number of Sequences')


fit<-glm(LogitProp ~Season*Location,data=phylum_Actin_df_L)
summary(fit)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)
stepAIC(fit,direction = "backward")


#install.packages("glmulti")
library(glmulti)
fit.model <- glmulti(fit, level = 1, crit = "aicc")
summary(fit.model)
summary(fit.model)$icvalue
fit<-glm(LogitProp ~ Season,data=phylum_Actin_df_L)
summary(fit)

##------------------------------
# select for Proteobacteria
phylum_Prote_df_L = subset.data.frame(phylum.gg, PhylumID=='Proteobacteria')
# lowest proportion for Proteobacteria is 0.009, hence e was considered 0.009.
phylum_Prote_df_L<-phylum_Prote_df_L %>%
  mutate(
    LogitProp = log((Abundance+0.009)/(1-Abundance+0.009))
  )

## Proteobacteria
# Plot proportion Proteobacteria vs Season
# boxplot season vs sex with season as colour
phylum_Prote_df_L$Season <- factor(phylum_Prote_df_L$Season)
theme_set(theme_bw())

## Season
phylum_Prote_Season_Location<-ggplot(data = phylum_Prote_df_L,
                                     aes(x=Season, y=LogitProp, fill= Season))+
  geom_boxplot(alpha=0.9, outliers.size = 2) +
  scale_fill_manual(values = colsSeason) +
  facet_wrap(Location~., labeller=labeller(Location= Local_labels)) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=2, position=position_jitterdodge(0.3), color="black") +
  ylab("Proportion Proteobacteria\n")+
  xlab("Season") +
  ylim(-4.2,1.4)+
  scale_y_continuous()+
  #scale_x_discrete(labels=c("summer 2013", "winter 2014", "summer 2014", "winter 2015"))+
  theme(axis.text.y = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black", angle = 90),
        strip.text = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_rect(fill = "lightgrey",colour="lightgrey"),
        legend.position = "none")
phylum_Prote_Season_Location

## GLMs TO TEST IF FIRMICUTES LOGIT TRANSFORMED PROPORTIONS VARY WITH SAME FULL MODEL AS FOR ALPHA-DOVERSITY
# data used: phylum_Prote_df_L
# see also David I. Warton and Francis K. C. Hui 2011 about using logit transformed proportions

# Model a)
str(phylum_Prote_df_L)
phylum_Prote_df_L$Season<-as.factor(phylum_Prote_df_L$Season)
phylum_Prote_df_L$Location<-as.factor(phylum_Prote_df_L$Location)
phylum_Prote_df_L$Sex<-as.factor(phylum_Prote_df_L$Sex)
phylum_Prote_df_L$Age<-as.factor(phylum_Prote_df_L$Age)
# Check the levels again
levels(phylum_Prote_df_L$Season)
levels(phylum_Prote_df_L$Sex)
levels(phylum_Prote_df_L$Age)
levels(phylum_Prote_df_L$Location)

# Check the structure of the data again
str(phylum_Prote_df_L)
summary(phylum_Prote_df_L)
class(phylum_Prote_df_L)

phylum_Prote_df_L<-as.data.frame(phylum_Prote_df_L)


## Generalize linear Model
ggplot(data = phylum_Prote_df_L) +
  geom_histogram(aes(x = LogitProp), bins = 30, color = 'gray30', fill = '#FFE8A2') +
  #geom_histogram(aes(x = length), binwidth = 2500, color = 'gray30', fill = '#FFE8A2') +
  theme_bw() +
  labs(x = 'Abundance', y = 'Number of Sequences')

fit<-glm(LogitProp ~Season+Location+Season*Location,data=phylum_Prote_df_L)
summary(fit)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)
stepAIC(fit,direction = "backward")

##------------------------------
# select for Tenericutes
phylum_Tener_df_L = subset.data.frame(phylum.gg, PhylumID=='Tenericutes')
# lowest proportion for Tenericutes is 0.0006, hence e was considered 0.0006.
phylum_Tener_df_L<-phylum_Tener_df_L %>%
  mutate(
    LogitProp = log((Abundance+0.0006)/(1-Abundance+0.0006))
  )

phylum_Tener_df_L$Season <- factor(phylum_Tener_df_L$Season)
theme_set(theme_bw())
phylum_Tener_Season_Location<-ggplot(data = phylum_Tener_df_L,
                                     aes(x=Season, y=Abundance, fill= Season))+
  geom_boxplot(alpha=0.9, outliers.size = 2) +
  scale_fill_manual(values = colsSeason) +
  facet_wrap(Location~., labeller=labeller(Location= Local_labels)) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=2, position=position_jitterdodge(0.3), color="black") +
  ylab("Proportion Tenericutes\n")+
  xlab("Season") +
  ylim(-4.2,1.4)+
  scale_y_continuous()+
  #scale_x_discrete(labels=c("summer 2013", "winter 2014", "summer 2014", "winter 2015"))+
  theme(axis.text.y = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black", angle = 90),
        strip.text = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_rect(fill = "lightgrey",colour="lightgrey"),
        legend.position = "none")
phylum_Tener_Season_Location

# Model a)
str(phylum_Tener_df_L)
phylum_Tener_df_L$Season<-as.factor(phylum_Tener_df_L$Season)
phylum_Tener_df_L$Location<-as.factor(phylum_Tener_df_L$Location)
phylum_Tener_df_L$Sex<-as.factor(phylum_Tener_df_L$Sex)
phylum_Tener_df_L$Age<-as.factor(phylum_Tener_df_L$Age)
# Check the levels again
levels(phylum_Tener_df_L$Season)
levels(phylum_Tener_df_L$Sex)
levels(phylum_Tener_df_L$Age)
levels(phylum_Tener_df_L$Location)

# Check the structure of the data again
str(phylum_Tener_df_L)
summary(phylum_Tener_df_L)
class(phylum_Tener_df_L)

phylum_Tener_df_L<-as.data.frame(phylum_Tener_df_L)


## Generalize linear Model
ggplot(data = phylum_Tener_df_L) +
  geom_histogram(aes(x = LogitProp), bins = 30, color = 'gray30', fill = '#FFE8A2') +
  #geom_histogram(aes(x = length), binwidth = 2500, color = 'gray30', fill = '#FFE8A2') +
  theme_bw() +
  labs(x = 'LogitProp', y = 'Number of Sequences')

fit<-glm(LogitProp ~Season*Location,data=phylum_Tener_df_L)
summary(fit)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)
stepAIC(fit,direction = "backward")
##--------------------------------------------------------------------------------------
##genus
setwd("D:/pigeon/seq/01052024/phyloseq/genus")
#genus
genus<- read.csv("genus_abundance.csv", sep=";", row.names = 1)
head(genus)

genus.ave <- colMeans(genus)

genus.0 <- genus[, order(-genus.ave)]
genus.1 <- genus[, order(-genus.ave)[1:10]]

genus_per <- as.data.frame(lapply(genus, function(x) x / sum(x)))
rownames(genus_per)<- row.names(genus)

genus.2 <- genus_per[, order(-genus.ave)[1:10]]
genus.2 <- cbind(genus.2, others = 1 - rowSums(genus.2))
genus.2 <- cbind(SampleID = rownames(genus.2), genus.2)


write.csv(genus.0,"genus.0.csv", row.names = T)
write.csv(genus.1,"genus.1.csv", row.names = T)
write.csv(genus.2,"genus.2.csv", row.names = T)

library(reshape2)
library(tidyr)
library(dplyr)
genus.2$SampleID <- rownames(genus.2)
genus.2 <- genus.2[, c("SampleID", setdiff(names(genus.2), "SampleID"))]


genus_long <- genus.2 %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "GenusID",
    values_to = "Abundance"
  )

head(genus_long)

write.csv(genus_long,file = "genus_long.csv",row.names = F)
metadata_FP = read.csv("metadata_FP.csv", sep="", row.names=1)
metadata_FP$SampleID <- rownames(metadata_FP)
metadata_FP <- metadata_FP[, c("SampleID", setdiff(names(metadata_FP), "SampleID"))]

genus.gg<-merge(genus_long,metadata_FP, by.x=c("SampleID"), by.y=c("SampleID"))
write.csv(genus.gg,"genus.gg.csv", row.names = T)
##----------------------------------------------------------------------------------
#genus.gg <- read.table("genus.gg.csv", sep=";", header=T, row.names = 1)
head(genus.gg)

library(wesanderson)
library(RColorBrewer)
library(ggpubr)
genus_reabundance<-ggbarplot(genus.gg, x = "SampleID", y="Abundance", color="black", fill="GenusID",
                             legend="middle",
                             legend.title="Genus", main="Relative abundance per Genus",
                             font.main = c(14,"bold", "black"), font.x = c(14, "bold"),
                             font.y=c(14,"bold")) +
  theme_bw() +
  rotate_x_text() +
  scale_fill_brewer(palette="Set3")+
  facet_grid(~ Season, scales = "free_x", space='free') +
  labs(x = "Samples", y = "Relative abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
genus_reabundance
genus_reabundance + theme(
  plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5))
ggsave(filename = "genus_reabundance_Season.pdf", device="pdf", width=8, height=6)
## Location_season
genus_reabundance<-ggbarplot(genus.gg, x = "SampleID", y="Abundance", color="black", fill="GenusID",
                             legend="middle",
                             legend.title="Genus", main="Relative abundance per Genus",
                             font.main = c(14,"bold", "black"), font.x = c(14, "bold"),
                             font.y=c(14,"bold")) +
  theme_bw() +
  rotate_x_text() +
  scale_fill_brewer(palette="Set3")+
  facet_grid(Season~Location , scales = "free_x", space='free') +
  labs(x = "Samples", y = "Relative abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
genus_reabundance
genus_reabundance + theme(
  plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5))
ggsave(filename = "genus_reabundance_season_Location.pdf", device="pdf", width=8, height=6)

##-----------------------------
# select for Lactobacillus
genus_Lacto_df_L = subset.data.frame(genus.gg, GenusID=='Lactobacillus')

summary(genus_Lacto_df_L)

# lowest proportion for Firmicutes is 0.02, hence e was considered 0.02.
genus_Lacto_df_L<-genus_Lacto_df_L %>%
  mutate(
    LogitProp = log((Abundance+0.02)/(1-Abundance+0.02))
  )


## Lactobacillus
# Plot proportion Lactobacillus vs Season
# boxplot season vs sex with season as colour
genus_Lacto_df_L$Season <- factor(genus_Lacto_df_L$Season)
theme_set(theme_bw())
genus_Lacto_Season_Location<-ggplot(data = genus_Lacto_df_L,
                                    aes(x=Season, y=Abundance, fill= Season))+
  geom_boxplot(alpha=0.9) +
  scale_fill_manual(values = colsSeason) +
  #facet_wrap(Location~., labeller=labeller(Location= Local_labels)) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=2, position=position_jitterdodge(0.3), color="black") +
  ylab("Proportion Lactobacillus\n")+
  xlab("Season") +
  ylim(0,10)+
  scale_y_continuous()+
  theme(axis.text.y = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black", angle = 90),
        strip.text = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_rect(fill = "lightgrey",colour="lightgrey"),
        legend.position = "none")
genus_Lacto_Season_Location
ggsave(filename = "genus_Lacto_Season_Location.pdf", device="pdf", width=8, height=6)


## GLMs TO TEST IF FIRMICUTES LOGIT TRANSFORMED PROPORTIONS VARY WITH SAME FULL MODEL AS FOR ALPHA-DOVERSITY
# see also David I. Warton and Francis K. C. Hui 2011 about using logit transformed proportions

# Model a)
str(genus_Lacto_df_L)
genus_Lacto_df_L$Season<-as.factor(genus_Lacto_df_L$Season)
genus_Lacto_df_L$Location<-as.factor(genus_Lacto_df_L$Location)
genus_Lacto_df_L$Sex<-as.factor(genus_Lacto_df_L$Sex)
genus_Lacto_df_L$Age<-as.factor(genus_Lacto_df_L$Age)
# Check the levels again
levels(genus_Lacto_df_L$Season)
levels(genus_Lacto_df_L$Sex)
levels(genus_Lacto_df_L$Age)
levels(genus_Lacto_df_L$Location)

# Check the structure of the data again
str(genus_Lacto_df_L)
summary(genus_Lacto_df_L)
class(genus_Lacto_df_L)

genus_Lacto_df_L<-as.data.frame(genus_Lacto_df_L)
fit<-glm(LogitProp ~Season*Location,data=genus_Lacto_df_L)
anova(fit, test = "F")
fit<-glm(LogitProp ~Season,data=genus_Lacto_df_L)
anova(fit, test = "F")

logit.step<-step(fit,direction = c("both"))
summary(logit.step)
stepAIC(fit,direction = "backward")

ggplot(data = genus_Lacto_df_L) +
  geom_histogram(aes(x = LogitProp), bins = 30, color = 'gray30', fill = '#FFE8A2') +
  #geom_histogram(aes(x = length), binwidth = 2500, color = 'gray30', fill = '#FFE8A2') +
  theme_bw() +
  labs(x = 'LogitProp', y = 'Number of Sequences')
##------------------------------
# select for Enterococcus
Genus_Enter_df_L = subset.data.frame(genus.gg, GenusID=='Enterococcus')

Genus_Enter_df_L <-Genus_Enter_df_L  %>%
  mutate(
    LogitProp = log((Abundance+0.002)/(1-Abundance+0.002))
  )

Genus_Enter_df_L $Season <- factor(Genus_Enter_df_L $Season)
theme_set(theme_bw())
Genus_Enter_Season_Location<-ggplot(data = Genus_Enter_df_L,
                                    aes(x=Season, y=Abundance, fill= Season))+
  geom_boxplot(alpha=0.9, outliers.size = 2) +
  scale_fill_manual(values = colsSeason) +
  facet_wrap(Location~., labeller=labeller(Location= Local_labels)) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=2, position=position_jitterdodge(0.3), color="black") +
  ylab("Proportion Enterococcus\n")+
  xlab("Season") +
  ylim(0,1)+
  scale_y_continuous()+
  scale_x_discrete(labels=c("summer 2013", "winter 2014", "summer 2014", "winter 2015"))+
  theme(axis.text.y = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black", angle = 90),
        strip.text = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_rect(fill = "lightgrey",colour="lightgrey"),
        legend.position = "none")
Genus_Enter_Season_Location

## GLMs TO TEST IF FIRMICUTES LOGIT TRANSFORMED PROPORTIONS VARY WITH SAME FULL MODEL AS FOR ALPHA-DOVERSITY
# data used: phylum_Firm_df_L
# see also David I. Warton and Francis K. C. Hui 2011 about using logit transformed proportions

# Model a)
str(Genus_Enter_df_L)
Genus_Enter_df_L$Season<-as.factor(Genus_Enter_df_L$Season)
Genus_Enter_df_L$Location<-as.factor(Genus_Enter_df_L$Location)
Genus_Enter_df_L$Sex<-as.factor(Genus_Enter_df_L$Sex)
Genus_Enter_df_L$Age<-as.factor(Genus_Enter_df_L$Age)
# Check the levels again
levels(Genus_Enter_df_L$Season)
levels(Genus_Enter_df_L$Sex)
levels(Genus_Enter_df_L$Age)
levels(Genus_Enter_df_L$Location)


# Check the structure of the data again
str(Genus_Enter_df_L)
summary(Genus_Enter_df_L)
class(phylum_Actin_df_L)

Genus_Enter_df_L<-as.data.frame(Genus_Enter_df_L)

fit<-glm(LogitProp ~Season+Location+Season*Location,data=Genus_Enter_df_L)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)

##------------------------------
# select for Corynebacterium
genus_Cory_df_L = subset.data.frame(genus.gg, GenusID=='Corynebacterium')
# lowest proportion for Proteobacteria is 0.009, hence e was considered 0.009.
genus_Cory_df_L<-genus_Cory_df_L %>%
  mutate(
    LogitProp = log((Abundance+0.009)/(1-Abundance+0.009))
  )


genus_Cory_df_L$Season <- factor(genus_Cory_df_L$Season)
theme_set(theme_bw())

## Season
genus_Cory_Location<-ggplot(data = genus_Cory_df_L,
                            aes(x=Location, y=LogitProp, fill= Location))+
  geom_boxplot(alpha=0.9, outliers.size = 2) +
  scale_fill_manual(values = colslocal) +
  #facet_wrap(Location~., labeller=labeller(Location= Local_labels)) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=2, position=position_jitterdodge(0.3), color="black") +
  ylab("Proportion Corynebacterium\n")+
  xlab("Location") +
  ylim(-4.2,1.4)+
  scale_y_continuous()+
  #scale_x_discrete(labels=c("summer 2013", "winter 2014", "summer 2014", "winter 2015"))+
  theme(axis.text.y = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black", angle = 90),
        strip.text = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_rect(fill = "lightgrey",colour="lightgrey"),
        legend.position = "none")
genus_Cory_Location

## GLMs TO TEST IF FIRMICUTES LOGIT TRANSFORMED PROPORTIONS VARY WITH SAME FULL MODEL AS FOR ALPHA-DOVERSITY
# data used: phylum_Firm_df_L
# see also David I. Warton and Francis K. C. Hui 2011 about using logit transformed proportions

# Model a)
str(genus_Cory_df_L)
genus_Cory_df_L$Season<-as.factor(genus_Cory_df_L$Season)
genus_Cory_df_L$Location<-as.factor(genus_Cory_df_L$Location)
genus_Cory_df_L$Sex<-as.factor(genus_Cory_df_L$Sex)
genus_Cory_df_L$Age<-as.factor(genus_Cory_df_L$Age)
# Check the levels again
levels(genus_Cory_df_L$Season)
levels(genus_Cory_df_L$Sex)
levels(genus_Cory_df_L$Age)
levels(genus_Cory_df_L$Location)

# Check the structure of the data again
str(genus_Cory_df_L)
summary(genus_Cory_df_L)
class(genus_Cory_df_L)

genus_Cory_df_L<-as.data.frame(genus_Cory_df_L)

fit<-glm(LogitProp ~Season*Location,data=genus_Cory_df_L)
anova(fit, test = "F")
fit<-glm(LogitProp ~Location,data=genus_Cory_df_L)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)
stepAIC(fit,direction = "backward")

emm  <- emmeans(fit, ~ Location)
pairs(emm, adjust = "tukey")


##------------------------------
# select for Shigella
genus_Shige_df_L = subset.data.frame(genus.gg, GenusID=='Shigella')
# lowest proportion for Tenericutes is 0.0006, hence e was considered 0.0006.
genus_Shige_df_L<-genus_Shige_df_L %>%
  mutate(
    LogitProp = log((Abundance+0.002)/(1-Abundance+0.002))
  )

genus_Shige_df_L$Season <- factor(genus_Shige_df_L$Season)
theme_set(theme_bw())
genus_Shige_Season_Location<-ggplot(data = genus_Shige_df_L,
                                    aes(x=Season, y=Abundance, fill= Season))+
  geom_boxplot(alpha=0.9, outliers.size = 2) +
  scale_fill_manual(values = colsSeason) +
  facet_wrap(Location~., labeller=labeller(Location= Local_labels)) +
  stat_summary(fun.y=mean, geom="point", shape=15, size=2, position=position_jitterdodge(0.3), color="black") +
  ylab("Proportion Shigella\n")+
  xlab("Season") +
  ylim(-4.2,1.4)+
  scale_y_continuous()+
  #scale_x_discrete(labels=c("summer 2013", "winter 2014", "summer 2014", "winter 2015"))+
  theme(axis.text.y = element_text(size=10, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(size=10, colour="black", angle = 90),
        strip.text = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_rect(fill = "lightgrey",colour="lightgrey"),
        legend.position = "none")
genus_Shige_Season_Location


## GLMs TO TEST IF FIRMICUTES LOGIT TRANSFORMED PROPORTIONS VARY WITH SAME FULL MODEL AS FOR ALPHA-DOVERSITY
# data used: phylum_Firm_df_L
# see also David I. Warton and Francis K. C. Hui 2011 about using logit transformed proportions

# Model a)
str(genus_Shige_df_L)
genus_Shige_df_L$Season<-as.factor(genus_Shige_df_L$Season)
genus_Shige_df_L$Location<-as.factor(genus_Shige_df_L$Location)
genus_Shige_df_L$Sex<-as.factor(genus_Shige_df_L$Sex)
genus_Shige_df_L$Age<-as.factor(genus_Shige_df_L$Age)
# Check the levels again
levels(genus_Shige_df_L$Season)
levels(genus_Shige_df_L$Sex)
levels(genus_Shige_df_L$Age)
levels(genus_Shige_df_L$Location)

# Check the structure of the data again
str(genus_Shige_df_L)
summary(genus_Shige_df_L)
class(genus_Shige_df_L)

genus_Shige_df_L<-as.data.frame(genus_Shige_df_L)

fit<-glm(LogitProp ~Season+Location+Season*Location,data=genus_Shige_df_L)
anova(fit, test = "F")
logit.step<-step(fit,direction = c("both"))
summary(logit.step)
stepAIC(fit,direction = "backward")

## Generalize linear Model
ggplot(data = genus_Shige_df_L) +
  geom_histogram(aes(x = LogitProp), bins = 30, color = 'gray30', fill = '#FFE8A2') +
  #geom_histogram(aes(x = length), binwidth = 2500, color = 'gray30', fill = '#FFE8A2') +
  theme_bw() +
  labs(x = 'LogitProp', y = 'Number of Sequences')

