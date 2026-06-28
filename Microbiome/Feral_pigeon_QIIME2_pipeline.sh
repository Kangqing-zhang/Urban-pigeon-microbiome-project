---
#title: "Environmental variation, not diet, shapes seasonal gut microbiome dynamics in urban feral pigeons"
#author: "Kangqing Zhang"
#date: "October 09, 2025"
---


module load QIIME2/2022.11
module list

cd /mnt/d/pigeon/seq/Rawsequence
awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} \
      NR>1{print $1"\t$PWD/"$1"_R1_001.fastq.gz\t$PWD/"$1"_R2_001.fastq.gz"}' \
      metadata.txt > manifest
head -n3 manifest

## import demultiplexed data
## manually make "manifest" file (.CSV); https://docs.qiime2.org/2018.6/tutorials/importing/?highlight=import%20data

	qiime tools import \
	  --type 'SampleData[PairedEndSequencesWithQuality]' \
	  --input-path manifest \
	  --output-path demux.qza \
	  --input-format PairedEndFastqManifestPhred33V2
cp demux.qza /mnt/d/pigeon/seq/01052024
cp manifest /mnt/d/pigeon/seq/01052024
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization pigeon-demux.qzv

## DADA2 denoising, including primer trimming, truncating based on quality plots, chimera removal, PhiX removal, and MERGING of R1 and R2 reads
### one could use --verbose flag to get stats on the different steps how many reads were filtered, denoised, merged, and non-chimeric
time qiime dada2 denoise-paired \
	  --i-demultiplexed-seqs demux.qza \
	  --p-n-threads 12 \
	  --p-trim-left-f 19 --p-trim-left-r 20 \
	  --p-trunc-len-f 0 --p-trunc-len-r 260 \
	  --o-table dada2-table.qza \
	  --o-representative-sequences dada2-rep-seqs.qza \
	  --o-denoising-stats denoising-stats.qza

## Create summary of FeatureTable and FeatureData
qiime feature-table summarize \
  --i-table dada2-table.qza \
  --o-visualization dada2-table.qzv \
  --m-sample-metadata-file metadata.txt

qiime feature-table tabulate-seqs \
  --i-data dada2-rep-seqs.qza \
  --o-visualization dada2-rep-seqs.qzv

cp dada2-table.qza table.qza
cp dada2-rep-seqs.qza rep-seqs.qza

##train dataset
## Silva 138 99% OTUs full-length sequences 
wget -c https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza

##Greengenes2 2022.10 full length sequences
wget -c ftp://download.nmdc.cn/tools/amplicon/silva/silva-138-99-nb-classifier.qza

wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ./gg_13_8_otus/rep_set/99_otus.fasta \
  --output-path 99_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ./gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza

##Train the classifier
time qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads 99_otus.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier_gg_13_8_99.qza

time qiime feature-classifier extract-reads \
  --i-sequences 99_otus.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --o-reads GG_99_ref-seqs.qza

# train classifier
  qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads GG_99_ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier GG_99_classifier.qza

# Filter non-target DNA from dataset; removes everything that differes more than 50% from GG/Silva db = ref-seqs.qza
qiime quality-control exclude-seqs \
  --i-query-sequences dada2-rep-seqs.qza \
  --i-reference-sequences GG_99_ref-seqs.qza \
  --p-method blast \
  --p-perc-identity 0.50 \
  --p-perc-query-aligned 0.50 \
  --o-sequence-hits pigeon-rep-seqs_min12S.qza \
  --o-sequence-misses pigeon12S.qza

 qiime feature-table filter-features \
  --i-table dada2-table.qza \
  --m-metadata-file pigeon12S.qza \
  --o-filtered-table pigeon-table_min12S.qza \
  --p-exclude-ids

## Summarize non-target-filtered FeatureTable and FeatureData and non-target sequences
qiime feature-table summarize \
  --i-table pigeon-table_min12S.qza \
  --o-visualization pigeon-table_min12S.qzv \
  --m-sample-metadata-file metadata.txt
	
qiime feature-table tabulate-seqs \
  --i-data pigeon-rep-seqs_min12S.qza \
  --o-visualization pigeon-rep-seqs_min12S.qzv
	
qiime feature-table tabulate-seqs \
  --i-data pigeon12S.qza \
  --o-visualization pigeon12S.qzv
	
# test classifier (use 150G mem)
qiime feature-classifier classify-sklearn \
  --i-classifier GG_99_classifier.qza \
  --i-reads pigeon-rep-seqs_min12S.qza \
  --o-classification pigeon_taxonomy_min12S.qza

### summarize classifier
qiime metadata tabulate \
  --m-input-file pigeon_taxonomy_min12S.qza \
  --o-visualization pigeon_taxonomy_min12S.qzv

## FILTER mtDNA/cpDNA/Archaea from dataset (also possible in Phyloseq)
## https://docs.qiime2.org/2017.10/tutorials/filtering/#taxonomy-based-filtering-of-tables-and-sequences

qiime taxa filter-table \
  --i-table pigeon-table_min12S.qza \
  --i-taxonomy pigeon_taxonomy_min12S.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-table pigeon-table_min12S_no-chloro-mito-arch.qza
qiime taxa filter-seqs \
  --i-sequences pigeon-rep-seqs_min12S.qza \
  --i-taxonomy pigeon_taxonomy_min12S.qza \
  --p-exclude mitochondria,chloroplast,archaea \
  --o-filtered-sequences pigeon-rep-seqs_min12S_no-chloro-mito-arch.qza

# Summarize non-mito-chloro-arch-filtered FeatureTable and FeatureData
qiime feature-table summarize \
  --i-table pigeon-table_min12S_no-chloro-mito-arch.qza \
  --o-visualization pigeon-table_min12S_no-chloro-mito-arch.qzv \
  --m-sample-metadata-file pigeons_metadata.txt

qiime feature-table tabulate-seqs \
  --i-data pigeon-rep-seqs_min12S_no-chloro-mito-arch.qza \
  --o-visualization pigeon-rep-seqs_min12S_no-chloro-mito-arch.qzv

## Generate a tree for phylogenetic diversity analyses

## QIIME supports several phylogenetic diversity metrics, including Faiths Phylogenetic Diversity and 
## weighted and unweighted UniFrac. In addition to counts of features per sample (i.e., the data in the 
## FeatureTable[Frequency] QIIME 2 artifact), these metrics require a rooted phylogenetic tree relating 
## the features to one another. This information will be stored in a Phylogeny[Rooted] QIIME 2 artifact. 
## The following steps will generate this QIIME 2 artifact.

## First, we perform a multiple sequence alignment of the sequences in our FeatureData[Sequence] to create 
## a FeatureData[AlignedSequence] QIIME 2 artifact. Here we do this with the mafft program.
## follow: https://docs.qiime2.org/2017.9/tutorials/moving-pictures/

qiime alignment mafft \
  --i-sequences pigeon-rep-seqs_min12S_no-chloro-mito-arch.qza  \
  --o-alignment aligned-pigeon-rep-seqs_min12S_no-chloro-mito-arch.qza

## Next, we mask (or filter) the alignment to remove positions that are highly variable. These positions 
## are generally considered to add noise to a resulting phylogenetic tree.

qiime alignment mask \
  --i-alignment aligned-pigeon-rep-seqs_min12S_no-chloro-mito-arch.qza \
  --o-masked-alignment masked-aligned-pigeon-rep-seqs_min12S_no-chloro-mito-arch.qza

## Next, we will apply FastTree to generate a phylogenetic tree from the masked alignment.

qiime phylogeny fasttree \
  --i-alignment masked-aligned-pigeon-rep-seqs_min12S_no-chloro-mito-arch.qza\
  --o-tree unrooted-tree_pigeon.qza

## The FastTree program creates an unrooted tree, so in the final step in this section we apply midpoint 
## rooting to place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree. 

qiime phylogeny midpoint-root \
  --i-tree unrooted-tree_pigeon.qza \
  --o-rooted-tree rooted-tree_pigeon.qza

# Export the feature table to biom format
    qiime tools export \
      --input-path pigeon-table_min12S_no-chloro-mito-arch.qza \
      --output-path feature-table
# Convert to tsv format
    biom convert -i feature-table/feature-table.biom \
      -o feature-table/feature-table.txt \
      --to-tsv
# Removing comment lines
    sed -i '/# Const/d' feature-table/feature-table.txt
    
   # Export representative sequences
    qiime tools export \
      --input-path pigeon-rep-seqs_min12S_no-chloro-mito-arch.qza \
      --output-path rep-seqs
    
    # Export species annotations
    qiime tools export \
      --input-path pigeon_taxonomy_min12S.qza \
      --output-path taxonomy

##Export each horizontal OTU table
qiime taxa collapse\
  --i-table pigeon-table_min12S_no-chloro-mito-arch.qza \
  --i-taxonomy pigeon_taxonomy_min12S.qza \
  --p-level 6\
  --o-collapsed-table table-l6.qza
  
qiime tools export\
  --input-path table-l6.qza\
  --output-path exported-table

biom convert -i exported-table/feature-table.biom\
  -o exported-table/silva_l6.txt --to-tsv

#Export ASV table
qiime tools export\
  --input-path pigeon-table_min12S_no-chloro-mito-arch.qza\
  --output-path exported-table
  
biom convert -i exported-table/feature-table.biom\
  -o exported-table/asv_table.txt --to-tsv

####Export ASV table
qiime tools export\
  --input-path table.qza\
  --output-path exported-table
  
biom convert -i exported-table/feature-table.biom\
  -o exported-table/asv_table.txt --to-tsv

##Export tree
qiime tools export\
  --input-path rooted-tree_pigeon.qza\
  --output-path exported-tree

qiime tools export\
  --input-path unrooted-tree_pigeon.qza\
  --output-path exported-unrooted-tree
