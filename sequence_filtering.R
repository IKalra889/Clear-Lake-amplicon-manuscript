#Set working directory
setwd("~/Desktop/Caron_lab_research/Clear_lake/Amplicon_analysis/final_analysis_2024/")

#load the necessary libraries
library("phyloseq")
library("ggplot2")      
library("dplyr")
library("tibble") 
library("tidyverse")
library(decontam)

#####16S data filtering####

#make matrix tables
otu_mat<- read.csv("ASV_counts_16S_exclude_chloro_mito.csv")
tax_mat<- read.csv("taxonomy_16S_exclude_chloro_mito.csv")
sample_info <- read.csv("sample_info_sorted.csv")

#define otu_id and sample as rownames
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu_id") 

tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu_id")

sample_info <- sample_info %>% 
  tibble::column_to_rownames("Sample") 

#Make OTU and tax table as matrix
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#sample info remains as dataframe

#Make phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(sample_info)

clear_bac <- phyloseq(OTU, TAX, samples) #phyloseq object
clear_bac

#check the phyloseq object
sample_names(clear_bac)
rank_names(clear_bac)
sample_variables(clear_bac)

#Identifying library sizes of the samples
df <- as.data.frame(sample_data(clear_bac)) #need dataframe for ggplot
df$LibrarySize <- sample_sums(clear_bac)
df <- df[order(df$LibrarySize),] #order according to librarysize
df$Index <- seq(nrow(df))

#plot library sizes
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Type)) + geom_point()

#using prevelance method to identify contaminants using decontam package

sample_data(clear_bac)$is.neg <- sample_data(clear_bac)$Type == "Control"
contamdf <- isContaminant(clear_bac, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf$contaminant) # 13 contaminants were identified
which(contamdf$contaminant) # decontaminating OTU indexes
#73  100  482 1639 1967 4633 4968 5407 5646 5738 7518 7569 8293

#make a vector with contaminant indices
i <- c(73,100,482,1639,1967,4633,4968,5407,5646,5738,7518,7569,8293)
#remove contaminating ASVs
bad_taxa <- row.names(contamdf[i,])
good_taxa <- setdiff(taxa_names(clear_bac), bad_taxa)
clear_bac2 <- prune_taxa(good_taxa, clear_bac)
#remove negative contol as no longer needed
clear_bac3 <- prune_samples(sample_names(clear_bac2) != c("C61", "C62"), clear_bac2)
clear_bac3 #60 samples

write.csv(as.data.frame(otu_table(clear_bac3)), "asv_16s_decontam.csv")
write.csv(as.data.frame(tax_table(clear_bac3)), "taxonomy_16s_decontam.csv")
write.csv(as.data.frame(sample_data(clear_bac3)), "sample_info_16s_decontam.csv")


#####18S data filtering#####

#make matrix tables
otu_mat<- read.csv("asv_18s_pr2.csv")
tax_mat<- read.csv("taxonomy_18s_pr2.csv")
sample_info_18s <- read.csv("sample_info_18s.csv")

#define otu_id and sample as rownames
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu_id") 

tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu_id")

sample_info_18s <- sample_info_18s %>% 
  tibble::column_to_rownames("Sample") 

#Make OTU and tax table as matrix
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#sample info remains as dataframe

#Make phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(sample_info_18s)

clear_euk <- phyloseq(OTU, TAX, samples) #phyloseq object
clear_euk

#check the phyloseq object
sample_names(clear_euk)
rank_names(clear_euk)
sample_variables(clear_euk)

#Identifying library sizes of the samples
df <- as.data.frame(sample_data(clear_euk)) #need dataframe for ggplot
df$LibrarySize <- sample_sums(clear_euk) 

df <- df[order(df$LibrarySize),] #order according to librarysize
df$Index <- seq(nrow(df)) 
#total 347,748 sequences

#plot library sizes
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Type)) + geom_point()

#using prevelance method to identify contaminants using decontam package
sample_data(clear_euk)$is.neg <- sample_data(clear_euk)$Type == "Control"
contamdf.euk <- isContaminant(clear_euk, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.euk$contaminant) # no contaminants were identified 

#remove negative contol as no longer needed
clear_euk2 <- prune_samples(sample_names(clear_euk) !="C62", clear_euk)
clear_euk2 #60 samples

#save the decontaminated data
write.csv(as.data.frame(otu_table(clear_euk2)), "asv_18s_decontam.csv")
write.csv(as.data.frame(tax_table(clear_euk2)), "taxonomy_18s_decontam.csv")
write.csv(as.data.frame(sample_data(clear_euk2)), "sample_info_18s_decontam.csv")

#total number of sequences before decontamination
sum(sample_sums(clear_bac)) #16S - 10908910
sum(sample_sums(clear_euk)) #18S - 347748

#total number of ASVs before decontamination
nrow(otu_table(clear_bac)) #16S - 8728 ASVs
nrow(otu_table(clear_euk)) #18S - 1265 ASVs

#total number of sequences after decontamination and removal of neg control
sum(sample_sums(clear_bac3)) #16S - 10,853,875
sum(sample_sums(clear_euk2)) #18S - 346,868

#total number of ASVs after decontamination and removal of neg control
nrow(otu_table(clear_bac3)) #16S - 8715 ASVs
nrow(otu_table(clear_euk2)) #18S - 1265 ASVs
