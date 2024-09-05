library(tidyverse)
library(reshape2)

#set working directory
setwd("~/Desktop/Caron_lab_research/Clear_lake/Amplicon_analysis/final_analysis_2024/")

#####16S rRNA#####

#16S rRNA gene amplicon sequencing dataset; metadata table; data.frame
#all are previously processed by decontam package to remove neg control and contaminants
sample_info_16s <- read.csv("sample_info_16s_decontam.csv")
colnames(sample_info_16s)[1] <- "Sample"
# feature table; data.frame
asv_table_16s <- read.csv("asv_16s_decontam.csv")

# taxonomic assignment table; data.frame [silva +ncbi database used]
taxonomy_table_16s <- read.csv("taxonomy_16s_decontam_ncbi.csv")

#setting the sample group and site levels
sample_info_16s$Group <- factor(sample_info_16s$Group, 
                                levels = c("August_2019", "August_2020", "July_2021", "August_2021", "September_2021", "October_2021"),
                                labels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"))

sample_info_16s$Site <- factor(sample_info_16s$Site,
                               levels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"), 
                               labels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"))

#make a table with total reads
asv_total <- apply(asv_table_16s[2:61], 1, sum) #sum function applied on rows to get number of reads for an OTU/ASV

#remove rare ASVs (less than 10 sequences total)
asv_16s_no10 <- asv_table_16s[asv_total > 10, ]
colnames(asv_16s_no10)[1] <- "otu_id"

#subset taxonomy to remove taxa corresponding to rare ASVs
taxonomy_16s_no10 <- taxonomy_table_16s %>% filter(otu_id %in% asv_16s_no10$otu_id)

#save files for future analyses
write.csv(asv_16s_no10, "asv_16s_no10.csv")
write.csv(taxonomy_16s_no10, "taxonomy_16s_no10.csv")

#subset samples
asv_16s_long <- melt(asv_16s_no10)
colnames(asv_16s_long)[2] <- "Sample"
asv_16s_long <- inner_join(asv_16s_long, sample_info_16s, by = "Sample")

#subset august
#filter samples collected in August
df <- asv_16s_long %>% filter(Month == "Aug")
#select columns needed
asv_16s_aug <- df %>% 
  select(otu_id, Sample, value) %>% 
  pivot_wider(names_from = "Sample", values_from = "value") #3365 ASVs total
#calculate sum counts of ASVs for all samples
asv_16s_aug$total <- apply(asv_16s_aug[2:31], 1, sum)
#remove ASVs with 0 total counts
#August 2019,2020 and 2021 ASV count table
asv_16s_aug <- asv_16s_aug %>% filter(total != 0) #3365 ASVs total
#subset taxonomy for august
tax_16s_aug <- taxonomy_16s_no10 %>% filter(otu_id %in% asv_16s_aug$otu_id)

write.csv(asv_16s_aug, "asv_16s_aug.csv")
write.csv(tax_16s_aug, "tax_16s_aug.csv")

#subset 2021
#filter samples collected in 2021
df <- asv_16s_long %>% filter(Year == "2021")
#select columns needed
asv_16s_2021 <- df %>% 
  select(otu_id, Sample, value) %>% 
  pivot_wider(names_from = "Sample", values_from = "value")
#calculate sum counts of ASVs for all samples
asv_16s_2021$total <- apply(asv_16s_2021[2:41], 1, sum)
#remove ASVs with 0 total counts
#July, Aug, Sep and Oct 2021 ASV count table
asv_16s_2021 <- asv_16s_2021 %>% filter(total != 0) #4235 ASVs total
#subset taxonomy for 2021
tax_16s_2021 <- taxonomy_16s_no10 %>% filter(otu_id %in% asv_16s_2021$otu_id)
#save files
write.csv(asv_16s_2021, "asv_16s_2021.csv")
write.csv(tax_16s_2021, "tax_16s_2021.csv")

#####18S ASV filtering and subsetting####
#18S rRNA gene amplicon sequencing dataset; metadata table; data.frame
#all are previously processed by decontam package to remove neg control and contaminants
sample_info_18s <- read.csv("sample_info_18s.csv")
colnames(sample_info_18s)[1] <- "Sample"
# feature table; data.frame
asv_table_18s <- read.csv("asv_18s_decontam.csv")

# taxonomic assignment table; data.frame [silva +ncbi database used]
taxonomy_table_18s <- read.csv("taxonomy_18s_decontam.csv")
colnames(taxonomy_table_18s)[1] <- "otu_id"

#setting the sample group and site levels
sample_info_18s$Group <- factor(sample_info_18s$Group, 
                                levels = c("August_2019", "August_2020", "July_2021", "August_2021", "September_2021", "October_2021"),
                                labels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"))

sample_info_18s$Site <- factor(sample_info_18s$Site,
                               levels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"), 
                               labels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"))

#make a table with total reads
asv_total <- apply(asv_table_18s[2:61], 1, sum) #sum function applied on rows to get number of reads for an OTU/ASV

#remove rare ASVs (less than 10 sequences total)
asv_18s_no10 <- asv_table_18s[asv_total > 10, ]
colnames(asv_18s_no10)[1] <- "otu_id"

#subset taxonomy to remove taxa corresponding to rare ASVs
taxonomy_18s_no10 <- taxonomy_table_18s %>% filter(otu_id %in% asv_18s_no10$otu_id)

#save files for future analyses
write.csv(asv_18s_no10, "asv_18s_no10.csv")
write.csv(taxonomy_18s_no10, "taxonomy_18s_no10.csv")

#subset samples
asv_18s_long <- melt(asv_18s_no10)
colnames(asv_18s_long)[2] <- "Sample"
asv_18s_long <- inner_join(asv_18s_long, sample_info_18s, by = "Sample")

#subset august
#filter samples collected in August
df <- asv_18s_long %>% filter(Month == "Aug")
#select columns needed
asv_18s_aug <- df %>% 
  select(otu_id, Sample, value) %>% 
  pivot_wider(names_from = "Sample", values_from = "value")
#calculate sum counts of ASVs for all samples
asv_18s_aug$total <- apply(asv_18s_aug[2:31], 1, sum)
#remove ASVs with 0 total counts
#August 2019,2020 and 2021 ASV count table
asv_18s_aug <- asv_18s_aug %>% filter(total != 0) 
#subset taxonomy for august
tax_18s_aug <- taxonomy_18s_no10 %>% filter(otu_id %in% asv_18s_aug$otu_id)
#save files
write.csv(asv_18s_aug, "asv_18s_aug.csv")
write.csv(tax_18s_aug, "tax_18s_aug.csv")

#subset 2021
#filter samples collected in 2021
df <- asv_18s_long %>% filter(Year == "2021")
#select columns needed
asv_18s_2021 <- df %>% 
  select(otu_id, Sample, value) %>% 
  pivot_wider(names_from = "Sample", values_from = "value")
#calculate sum counts of ASVs for all samples
asv_18s_2021$total <- apply(asv_18s_2021[2:41], 1, sum)
#remove ASVs with 0 total counts
#July, Aug, Sep and Oct 2021 ASV count table
asv_18s_2021 <- asv_18s_2021 %>% filter(total != 0)
#subset taxonomy for 2021
tax_18s_2021 <- taxonomy_18s_no10 %>% filter(otu_id %in% asv_18s_2021$otu_id)
#save files
write.csv(asv_18s_2021, "asv_18s_2021.csv")
write.csv(tax_18s_2021, "tax_18s_2021.csv")
