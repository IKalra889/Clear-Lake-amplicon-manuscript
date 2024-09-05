##Figure 3##
##Eukaryotic community diversity and composition at Clear Lake##

library(tidyverse)
library(ggplot2)
library(microeco)

#set working directory
setwd("~/Desktop/Caron_lab_research/Clear_lake/Amplicon_analysis/final_analysis_2024/")

# use pipe operator in magrittr package
library(magrittr)
# fix the random number generation to make the results repeatable
set.seed(123)
# set the plotting background
theme_set(theme_bw())

#load files
asv_18s_2021 <- read.csv("asv_18s_2021.csv") %>% select(-total)
tax_18s_2021 <- read.csv("tax_18s_2021.csv") 
asv_18s_aug <- read.csv("asv_18s_aug.csv") %>% select(-total)
tax_18s_aug <- read.csv("tax_18s_aug.csv")
sample_info <- read.csv("sample_info_18s_decontam.csv")

tax1 <- colsplit(tax_18s_2021$taxonomy, ";", c("Domain", "Supergroup", "Phylum", "Class", "Order",
                                              "Family", "Genus", "Species"))
tax_18s_2021 <- cbind(tax_18s_2021, tax1) %>% select(-taxonomy)

tax2 <- colsplit(tax_18s_aug$taxonomy, ";", c("Domain", "Supergroup", "Phylum", "Class", "Order",
                                               "Family", "Genus", "Species"))
tax_18s_aug <- cbind(tax_18s_aug, tax2) %>% select(-taxonomy)

#remove ASVs that don't have phyla level information
x <- tax_18s_aug %>% filter(Phylum != "") #lost 13 ASVs
y <- tax_18s_2021 %>% filter(Phylum != "") #lost 23 ASVs

#setting the sample group and site levels
sample_info$Group <- factor(sample_info$Group, 
                            levels = c("August_2019", "August_2020", "July_2021", "August_2021", "September_2021", "October_2021"),
                            labels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"))

sample_info$Site <- factor(sample_info$Site,
                           levels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"), 
                           labels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"))

#save
write.csv(sample_info, "sample_metadata.csv", row.names = FALSE)

#making all column1 as rownames
row.names(asv_18s_2021) <- NULL
asv_18s_2021 <- asv_18s_2021 %>%
  tibble::column_to_rownames("otu_id")

tax_18s_2021 <- tax_18s_2021 %>% 
  tibble::column_to_rownames("otu_id")

sample_info_2021 <- sample_info %>% filter(Year == 2021) %>% 
  tibble::column_to_rownames("Sample")

row.names(asv_18s_aug) <- NULL
asv_18s_aug <- asv_18s_aug %>%
  tibble::column_to_rownames("otu_id") 

tax_18s_aug <- tax_18s_aug %>% 
  tibble::column_to_rownames("otu_id")

sample_info_aug <- sample_info %>% filter(Month == "Aug") %>%
  tibble::column_to_rownames("Sample")

#create a microtable object
euks_2021 <- microtable$new(sample_table = sample_info_2021, otu_table = asv_18s_2021, 
                             tax_table = tax_18s_2021)
euks_aug <- microtable$new(sample_table = sample_info_aug, otu_table = asv_18s_aug, 
                            tax_table = tax_18s_aug)

###Venn diagrams for ASV shared and unique between different months and years###
##ASVs are total for all 5 sites within the sampling date##
# merge samples as one community for each group
df1 <- euks_2021$merge_samples(use_group = "Group")
df2 <- euks_aug$merge_samples(use_group = "Group")
# create trans_venn object
t1 <- trans_venn$new(df1, ratio = NULL)
t2 <- trans_venn$new(df2, ratio = NULL)
##venn diagrams
venn.2021 <- t1$plot_venn(color_circle = RColorBrewer::brewer.pal(8, "Set1"))
venn.aug <- t2$plot_venn()

###Phyla level box plots averaged for all sites###

# remove ASVs that don't have phyla level information
x <- tax_18s_2021 %>% filter(Phylum != "") #lost 13 ASVs
y <- tax_18s_aug %>% filter(Phylum != "") #lost 23 ASVs
#make new microeco objects
euks_2021.x <- microtable$new(sample_table = sample_info_2021, otu_table = asv_18s_2021, 
                            tax_table = x)
euks_aug.y <- microtable$new(sample_table = sample_info_aug, otu_table = asv_18s_aug, 
                           tax_table = y)

## normalisation of the dataset

#2021 samples
#make a temp object
df <- trans_norm$new(dataset = euks_2021.x)
#total sum scaling method for normalisation
euks_2021.tss <- df$norm(method = "TSS")

#august samples
df <- trans_norm$new(dataset = euks_aug.y)
euks_aug.tss <- df$norm(method = "TSS")

#boxplots
t3 <- trans_abund$new(dataset = euks_2021.tss, taxrank = "Phylum", ntaxa = 10)
box.2021 <- t3$plot_box(group = "Group", xtext_angle = 30, color_values = RColorBrewer::brewer.pal(8, "Set1"))

t4 <- trans_abund$new(dataset = euks_aug.tss, taxrank = "Phylum", ntaxa = 10)
box.aug <- t4$plot_box(group = "Group", xtext_angle = 30)

##plot figures
venn.2021 #Fig 3.A
box.2021 <- box.2021 + guides(fill=guide_legend(title="")) + 
  theme(legend.position = "top", legend.text = element_text(size=12)) #Fig.3B
venn.aug #Fig. 3C
box.aug <- box.aug + guides(fill=guide_legend(title="")) +
  theme(legend.position = "top", legend.text = element_text(size=12)) #Fig.3D


#save files
ggsave("fig.3a.pdf",plot=venn.2021, width=8, height=6)
ggsave("fig.3b.pdf",plot=venn.aug, width=6, height=6)
ggsave("fig.3c.pdf", plot=box.2021, width = 8, height=6)
ggsave("fig.3d.pdf", plot=box.aug, width = 8, height=6)

#make a combined figure
library(ggpubr)
figure3 <- ggarrange(venn.2021, box.2021, venn.aug,box.aug,
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2)
#save figure
ggsave("figure3.pdf", figure3, width = 12.5, height = 10)