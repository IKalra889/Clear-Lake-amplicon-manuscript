##Figure 2##
##Bacterial community diversity and composition at Clear Lake##

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
asv_16s_2021 <- read.csv("asv_16s_2021.csv") %>% select(-total)
tax_16s_2021 <- read.csv("tax_16s_2021.csv") 
asv_16s_aug <- read.csv("asv_16s_aug.csv") %>% select(-total)
tax_16s_aug <- read.csv("tax_16s_aug.csv")
sample_info <- read.csv("sample_info_16s_decontam.csv")
colnames(sample_info)[1] <- "Sample"
#setting the sample group and site levels
sample_info$Group <- factor(sample_info$Group, 
                                levels = c("August_2019", "August_2020", "July_2021", "August_2021", "September_2021", "October_2021"),
                                labels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"))

sample_info$Site <- factor(sample_info$Site,
                               levels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"), 
                               labels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"))

#making all column1 as rownames
row.names(asv_16s_2021) <- NULL
asv_16s_2021 <- asv_16s_2021 %>%
  tibble::column_to_rownames("otu_id")

tax_16s_2021 <- tax_16s_2021 %>% 
  tibble::column_to_rownames("otu_id")

sample_info_2021 <- sample_info %>% filter(Year == 2021) %>% 
  tibble::column_to_rownames("Sample")
  
row.names(asv_16s_aug) <- NULL
asv_16s_aug <- asv_16s_aug %>%
  tibble::column_to_rownames("otu_id") 

tax_16s_aug <- tax_16s_aug %>% 
  tibble::column_to_rownames("otu_id")

sample_info_aug <- sample_info %>% filter(Month == "Aug") %>%
  tibble::column_to_rownames("Sample")

#create a microtable object
proks_2021 <- microtable$new(sample_table = sample_info_2021, otu_table = asv_16s_2021, 
                             tax_table = tax_16s_2021)
proks_aug <- microtable$new(sample_table = sample_info_aug, otu_table = asv_16s_aug, 
                            tax_table = tax_16s_aug)

###Venn diagrams for ASV shared and unique between different months and years###
##ASVs are total for all 5 sites within the sampling date##
# merge samples as one community for each group
df1 <- proks_2021$merge_samples(use_group = "Group")
df2 <- proks_aug$merge_samples(use_group = "Group")
# create trans_venn object
t1 <- trans_venn$new(df1, ratio = NULL)
t2 <- trans_venn$new(df2, ratio = NULL)
##venn diagrams
venn.2021 <- t1$plot_venn(color_circle = RColorBrewer::brewer.pal(8, "Set1"))
venn.aug <- t2$plot_venn()

###Phyla level box plots averaged for all sites###
#normalisation of the dataset
#2021 samples
#make a temp object
df <- trans_norm$new(dataset = proks_2021)
#total sum scaling method for normalisation
proks_2021.tss <- df$norm(method = "TSS")

#august samples
df <- trans_norm$new(dataset = proks_aug)
proks_aug.tss <- df$norm(method = "TSS")

#boxplots
t3 <- trans_abund$new(dataset = proks_2021.tss, taxrank = "phylum", ntaxa = 10)
box.2021 <- t3$plot_box(group = "Group", xtext_angle = 30, color_values = RColorBrewer::brewer.pal(8, "Set1"))

t4 <- trans_abund$new(dataset = proks_aug.tss, taxrank = "phylum", ntaxa = 10)
box.aug <- t4$plot_box(group = "Group", xtext_angle = 30)

##plot figures
venn.2021 #Fig 2.A
box.2021 <- box.2021 + guides(fill=guide_legend(title="")) + theme(legend.position = "top", legend.text = element_text(size=12)) #Fig.2B
venn.aug #Fig. 2C
box.aug <- box.aug + guides(fill=guide_legend(title="")) +theme(legend.position = "top", legend.text = element_text(size=12)) #Fig.2D
box.2021

#save files
ggsave("fig.2a.pdf",plot=venn.2021, width=8, height=6)
ggsave("fig.2b.pdf",plot=venn.aug, width=6, height=6)
ggsave("fig.2c.pdf", plot=box.2021, width = 8, height=6)
ggsave("fig.2d.pdf", plot=box.aug, width = 8, height=6)

#make a combined figure
library(ggpubr)
figure2 <- ggarrange(venn.2021, box.2021, venn.aug,box.aug,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
#save figure
ggsave("figure2.pdf", figure2, width = 12.5, height = 10)
