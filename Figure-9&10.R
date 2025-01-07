library(tidyverse)
library(microeco)
library(reshape2)
library(ggpubr)
library(mecodev)
library(edgeR)
#set working directory
setwd("~/Desktop/Caron_lab_research/Clear_lake/Amplicon_analysis/final_analysis_2024/")

#load normalized data
asv_16s <- read.csv("asv_16s.tss.csv")
tax_16s <- read.csv("taxonomy_16s_no10.csv") 
asv_18s <- read.csv("asv_18s.tss.csv")
tax_18s <- read.csv("taxonomy_18s_no10.csv")
tax1 <- colsplit(tax_18s$taxonomy, ";", c("domain", "supergroup", "phylum", "class", "order",
                                               "family", "genus", "species"))
tax_18s <- cbind(tax_18s, tax1) %>% select(-taxonomy)

sample_info <- read.csv("sample_metadata2.csv")

#setting the sample group and site levels
sample_info$Sample <- factor(sample_info$Sample, 
                             levels = c("Aug_2019", "Aug_2020", "July_2021", "Aug_2021", "Sep_2021", "Oct_2021"),
                             labels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"))

sample_info$Site <- factor(sample_info$Site,
                           levels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"), 
                           labels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"))

#making all column1 as rownames
row.names(asv_16s) <- NULL
asv_16s <- asv_16s %>%
  tibble::column_to_rownames("otu_id")

tax_16s <- tax_16s %>% 
  tibble::column_to_rownames("otu_id")

row.names(asv_18s) <- NULL
asv_18s <- asv_18s %>%
  tibble::column_to_rownames("otu_id")

tax_18s <- tax_18s %>% 
  tibble::column_to_rownames("otu_id")

sample_info <- sample_info %>% 
  tibble::column_to_rownames("Sample_ID")

#create a microtable object
proks <- microtable$new(sample_table = sample_info, otu_table = asv_16s, 
                            tax_table = tax_16s)
euks <- microtable$new(sample_table = sample_info, otu_table = asv_18s, 
                           tax_table = tax_18s)
#remove rows with no genera
euks$tax_table <- euks$tax_table %>% filter(genus!= "")
euks$tidy_dataset()

####Figure 9 - RDA plots####
##proks
#create new env
proks_env <- trans_env$new(dataset = proks, env_cols = 10:18, complete_na = TRUE)
#calculate RDA
proks_env$cal_ordination(method = "RDA", taxa_level = "genus")
#adjust ordination
proks_env$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
#plot
fig9a <- proks_env$plot_ordination(plot_color = "Sample", plot_shape = "Site")

##euks
#create new env
euks_env <- trans_env$new(dataset = euks, env_cols = 10:18, complete_na = TRUE)
#calculate RDA
euks_env$cal_ordination(method = "RDA", taxa_level = "genus")
#adjust ordination
euks_env$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
#plot
fig9b <- euks_env$plot_ordination(plot_color = "Sample", plot_shape = "Site")
fig9b
#save plots
ggsave("fig.9a.pdf", plot=fig9a, width = 9, height = 6)
ggsave("fig.9b.pdf", plot=fig9b, width = 9, height = 6)


####Correlation Matrix####

##proks
#cal correlation with the env object
proks_env$cal_cor(use_data = "genus", cor_method = "spearman", p_adjust_method = "fdr", p_adjust_type = "Env", use_taxa_num = 10) #
# default ggplot2 method with clustering
fig10a <- proks_env$plot_cor(filter_feature = c(""))

##euks
#cal correlation with the env object
euks_env$cal_cor(use_data = "genus", cor_method = "spearman", p_adjust_method = "fdr", p_adjust_type = "Env", use_taxa_num = 10) #
# default ggplot2 method with clustering
fig10b <- euks_env$plot_cor(filter_feature = c(""))

#save plots
ggsave("fig.10a.pdf", plot=fig10a, width = 7.5, height = 6)
ggsave("fig.10b.pdf", plot=fig10b, width = 7.5, height = 6)
