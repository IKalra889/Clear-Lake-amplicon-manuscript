####Figure 4 & Figure 8####
##Monthly and yearly variations in microbial community and diversity at Clear Lake##

library(tidyverse)
library(ggplot2)
library(microeco)
library(reshape2)
library(ggpubr)
#set working directory
setwd("~/Desktop/Caron_lab_research/Clear_lake/Amplicon_analysis/final_analysis_2024/")

# use pipe operator in magrittr package
library(magrittr)
# fix the random number generation to make the results repeatable
set.seed(123)
# set the plotting background
theme_set(theme_bw())

####PCoA Analysis####

#load files
#2021
asv_16s_2021 <- read_csv("asv_16s_2021.csv") %>% select(-total)
tax_16s_2021 <- read_csv("tax_16s_2021.csv")
asv_18s_2021 <- read_csv("asv_18s_2021.csv") %>% select(-total)
tax_18s_2021 <- read_csv("tax_18s_2021.csv") 
tax1 <- colsplit(tax_18s_2021$taxonomy, ";", c("Domain", "Supergroup", "Phylum", "Class", "Order",
                                               "Family", "Genus", "Species"))
tax_18s_2021 <- cbind(tax_18s_2021, tax1) %>% select(-taxonomy)

#august-2019, 2020 and 2021
asv_16s_aug <- read_csv("asv_16s_aug.csv") %>% select(-total)
tax_16s_aug <- read_csv("tax_16s_aug.csv")
asv_18s_aug <- read_csv("asv_18s_aug.csv") %>% select(-total)
tax_18s_aug <- read_csv("tax_18s_aug.csv")
tax1 <- colsplit(tax_18s_aug$taxonomy, ";", c("Domain", "Supergroup", "Phylum", "Class", "Order",
                                               "Family", "Genus", "Species"))
tax_18s_aug <- cbind(tax_18s_aug, tax1) %>% select(-taxonomy)

sample_info <- read_csv("sample_metadata.csv", locale=locale(encoding="latin1"))
#setting the sample group and site levels
sample_info$Sample <- factor(sample_info$Sample, 
                            levels = c("Aug_2019", "Aug_2020", "July_2021", "Aug_2021", "Sep_2021", "Oct_2021"),
                            labels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"))

sample_info$Site <- factor(sample_info$Site,
                           levels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"), 
                           labels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"))

##2021 analyses##
#making all column1 as rownames
row.names(asv_16s_2021) <- NULL
asv_16s_2021 <- asv_16s_2021 %>%
  tibble::column_to_rownames("otu_id")

tax_16s_2021 <- tax_16s_2021 %>% 
  tibble::column_to_rownames("otu_id")

row.names(asv_18s_2021) <- NULL
asv_18s_2021 <- asv_18s_2021 %>%
  tibble::column_to_rownames("otu_id")

tax_18s_2021 <- tax_18s_2021 %>% 
  tibble::column_to_rownames("otu_id")

sample_info_2021 <- sample_info %>% filter(Year == 2021) %>% 
  tibble::column_to_rownames("Sample_ID")

#create a microtable object
proks_2021 <- microtable$new(sample_table = sample_info_2021, otu_table = asv_16s_2021, 
                            tax_table = tax_16s_2021)
euks_2021 <- microtable$new(sample_table = sample_info_2021, otu_table = asv_18s_2021, 
                           tax_table = tax_18s_2021)

##normalize data using total sum scaling
#make a temp object
df <- trans_norm$new(dataset = proks_2021)
#total sum scaling method for normalisation
proks_2021.tss <- df$norm(method = "TSS")

#make a temp object
df <- trans_norm$new(dataset = euks_2021)
#total sum scaling method for normalisation
euks_2021.tss <- df$norm(method = "TSS")

##proks##
#calculate beta diversity#
proks_2021.tss$cal_betadiv()
t1 <- trans_beta$new(dataset = proks_2021.tss, group = "Sample", measure = "bray")
# PcoA calculation
t1$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
proks.2021.pcoa <- t1$plot_ordination(plot_color = "Sample", plot_shape = "Site",
                                      color_values = RColorBrewer::brewer.pal(8, "Set1"),
                                      plot_type = c("point", "ellipse"))

##euks##
#calculate beta diversity
euks_2021.tss$cal_betadiv()
t1 <- trans_beta$new(dataset = euks_2021.tss, group = "Sample", measure = "bray")
# PcoA calculation
t1$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
euks.2021.pcoa <- t1$plot_ordination(plot_color = "Sample", plot_shape = "Site",
                                      color_values = RColorBrewer::brewer.pal(8, "Set1"),
                                      plot_type = c("point", "ellipse"))

#save figures
ggsave("fig.4a.pdf",plot=proks.2021.pcoa, width=6, height=4)
ggsave("fig.4b.pdf",plot=euks.2021.pcoa, width=6, height=4)


###August - 2019, 2020, 2021 analysis###
#making all column1 as rownames
row.names(asv_16s_aug) <- NULL
asv_16s_aug <- asv_16s_aug %>%
  tibble::column_to_rownames("otu_id")

tax_16s_aug <- tax_16s_aug %>% 
  tibble::column_to_rownames("otu_id")

row.names(asv_18s_aug) <- NULL
asv_18s_aug <- asv_18s_aug %>%
  tibble::column_to_rownames("otu_id")

tax_18s_aug <- tax_18s_aug %>% 
  tibble::column_to_rownames("otu_id")

sample_info_aug <- sample_info %>% filter(Month == "Aug") %>% 
  tibble::column_to_rownames("Sample_ID")

#create a microtable object
proks_aug <- microtable$new(sample_table = sample_info_aug, otu_table = asv_16s_aug, 
                             tax_table = tax_16s_aug)
euks_aug <- microtable$new(sample_table = sample_info_aug, otu_table = asv_18s_aug, 
                            tax_table = tax_18s_aug)

##normalize data using total sum scaling
#make a temp object
df <- trans_norm$new(dataset = proks_aug)
#total sum scaling method for normalisation
proks_aug.tss <- df$norm(method = "TSS")

#make a temp object
df <- trans_norm$new(dataset = euks_aug)
#total sum scaling method for normalisation
euks_aug.tss <- df$norm(method = "TSS")

##proks##
#calculate beta diversity#
proks_aug.tss$cal_betadiv()
t1 <- trans_beta$new(dataset = proks_aug.tss, group = "Sample", measure = "bray")
# PcoA calculation
t1$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
proks.aug.pcoa <- t1$plot_ordination(plot_color = "Sample", plot_shape = "Site",
                                      color_values = RColorBrewer::brewer.pal(8, "Dark2"),
                                      plot_type = c("point", "ellipse"))

##euks##
#calculate beta diversity
euks_aug.tss$cal_betadiv()
t1 <- trans_beta$new(dataset = euks_aug.tss, group = "Sample", measure = "bray")
# PcoA calculation
t1$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
euks.aug.pcoa <- t1$plot_ordination(plot_color = "Sample", plot_shape = "Site",
                                     color_values = RColorBrewer::brewer.pal(8, "Dark2"),
                                     plot_type = c("point", "ellipse"))
euks.aug.pcoa +theme(legend.position = "top")
#save figures
ggsave("fig.8a.pdf",plot=proks.aug.pcoa, width=6, height=4)
ggsave("fig.8b.pdf",plot=euks.aug.pcoa, width=6, height=4)


####Composition Bar Plots by Site and Month####

#bacteria
asv_16s <- read.csv("asv_16s_no10.csv") %>% select(-X)
tax_16s <- read.csv("taxonomy_16s_no10.csv")
colnames(tax_16s)[1] <- "otu_id"

#eukaryotes
asv_18s <- read.csv("asv_18s_no10.csv") %>% select(-X)
tax_18s <- read.csv("taxonomy_18s_no10.csv")

#sample metadata
sample_info <- read.csv("sample_metadata.csv")

##normalization of data - total sum scaling

#bacteria
total <- colSums(asv_16s[2:61])
df <- asv_16s$otu_id %>% as.data.frame()
for (i in 2:61){
  x <- asv_16s[i]/total[i-1]
  df <- cbind(df, x)
} 
asv_16s.tss <- df
colnames(asv_16s.tss)[1] <- "otu_id"
#save
write.csv(asv_16s.tss, "asv_16s.tss.csv", row.names = FALSE)

#eukaryotes
total <- colSums(asv_18s[2:61])
df <- asv_18s$otu_id %>% as.data.frame()
for (i in 2:61){
  x <- asv_18s[i]/total[i-1]
  df <- cbind(df, x)
} 
asv_18s.tss <- df
colnames(asv_18s.tss)[1] <- "otu_id"
#save
write.csv(asv_18s.tss, "asv_18s.tss.csv", row.names = FALSE)

##averaging the replicates##

##bacteria
tmp <- melt(asv_16s.tss)
colnames(tmp)[2] <- "Sample"
tmp <- inner_join(tmp, sample_info, by = "Sample")
tmp2 <- tmp %>% group_by(otu_id, Name) %>%
  summarise(avg = sum(value)/2)
asv_16s.avg <- tmp2

tmp <- asv_16s.avg %>% pivot_wider(names_from = Name, values_from = avg)

#make bacterial taxonomy clean
df <- tax_16s$otu_id %>% as.data.frame()
df$domain <- str_remove_all(tax_16s$domain, "d__")
df$phylum <- str_remove_all(tax_16s$phylum, "p__")
df$class <- str_remove_all(tax_16s$class, "c__")
df$order <- str_remove_all(tax_16s$order, "o__")
df$family <- str_remove_all(tax_16s$family, "f__")
df$genus <- str_remove_all(tax_16s$genus, "g__")
df$species <- str_remove_all(tax_16s$species, "s__")
colnames(df)[1] <- "otu_id"
tax_16s <- df

#join taxonomy with normalised and averaged asv counts
asv_16s.avg_wtax <- inner_join(asv_16s.avg, tax_16s, by = "otu_id")
#save
write.csv(asv_16s.avg_wtax %>% pivot_wider(names_from = Name, values_from = avg), "asv_16s_avg.csv", row.names = FALSE)

##eukaryotes
tmp <- melt(asv_18s.tss)
colnames(tmp)[2] <- "Sample"
tmp <- inner_join(tmp, sample_info, by = "Sample")
tmp2 <- tmp %>% group_by(otu_id, Name) %>%
  summarise(avg = sum(value)/2)
asv_18s.avg <- tmp2

#splitting taxonomy
df <- colsplit(tax_18s$taxonomy, ";", c("domain", "supergroup", "phylum", "class", "order", "family", "genus", "species"))
df$otu_id <- tax_18s$otu_id

#join taxonomy with normalised and averaged asv counts
asv_18s.avg_wtax <- inner_join(asv_18s.avg, df, by = "otu_id")

#save
write.csv(asv_18s.avg_wtax %>% pivot_wider(names_from = Name, values_from = avg), "asv_18s_avg.csv", row.names = FALSE)


####TAXONOMY BAR PLOTS####

#sample info with required columns
sample_info <- read.csv("sample_metadata.csv")
df <- sample_info %>% select(Name, Site, Date, Month, Year, Group) %>% unique()

#read averaged and normalized asv tables with taxonomy
asv_16s.avg_wtax <- read_csv("asv_16s_avg.csv") %>% melt(variable.name="Name", value.name="avg")
asv_18s.avg_wtax <- read_csv("asv_18s_avg.csv") %>% melt(variable.name="Name", value.name="avg")

df_16s <- inner_join(df, asv_16s.avg_wtax, by="Name")
df_18s <- inner_join(df, asv_18s.avg_wtax, by="Name")

#set tax color
tax_color<-c('firebrick4','indianred1','tomato3','forestgreen',
                          'yellowgreen','darkgrey','lightgrey','gold1','moccasin',
                          'lightblue','blue','ivory1','magenta',
                          'mediumvioletred',
                          'lightpink','#DDAD4B','tan2','tan3','tan4','hotpink4','darkblue','lightseagreen','lightskyblue1',
                          'sienna1','sienna3','sienna4','thistle','thisle4')



##Cyanobacteria bar plot##
#selecting only cyanos from the 16S table

cyanos <- df_16s %>% filter(phylum == "Cyanobacteria")
##simplifying taxonomy
cyanos$taxa="Other"
cyanos$taxa[cyanos$genus == "Lyngbya"]="Lyngbya"
cyanos$taxa[cyanos$genus == "Planktothrix_NIVA-CYA_15"]="Planktothrix"
cyanos$taxa[cyanos$genus == "Aphanizomenon_NIES81"]="Aphanizomenon"
cyanos$taxa[cyanos$genus == "Cyanobium_PCC-6307"]="Cyanobium"
cyanos$taxa[cyanos$genus == "Dolichospermum_NIES41"]="Dolichospermum"
cyanos$taxa[cyanos$genus == "Limnothrix"]="Limnothrix"
cyanos$taxa[cyanos$genus == "Microcystis_PCC-7914"]="Microcystis"
cyanos$taxa[cyanos$genus == "Nodularia_PCC-9350"]="Nodularia"
cyanos$taxa[cyanos$species == "Anabaenopsis_sp."]="Anabaenopsis"
cyanos$taxa[cyanos$genus == "Pseudanabaena_PCC-7429"]="Pseudanabaena"
cyanos$taxa[cyanos$genus == "Synechocystis_SAG_90.79"]="Synechocystis"

#aggregate the new taxonomy
NTC <- aggregate(cyanos$avg, by=list(site=cyanos$Site, sample=cyanos$Group, taxa=cyanos$taxa), sum)
NTC$Group <- NTC$sample
NTC <- NTC %>% separate(col=Group, into=c("month", "year"), sep="-")
NTC$rel <- NTC$x * 100
#factoring the variables
NTC$sample <- factor(NTC$sample, 
                     levels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"),
                     labels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"))
NTC$site <- factor(NTC$site,
                   levels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"), 
                   labels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"))
##phyla bar plots##
#2021
cyano_bar.2021 <- NTC %>% filter(year == "2021") %>%
  ggplot(aes(y=rel,x=site,fill=taxa))+geom_bar(stat="identity", position="fill", color="black")+labs(title="", x="",y="Relative abundance (%)")+
  theme_classic() +
  theme(legend.title=element_blank(),legend.position="right", axis.text = element_text(color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.text.x = element_text(angle = 90, hjust = 1, size = 11) )+
  scale_fill_manual(values=c("firebrick","lightblue2","pink","sienna1","lightgreen","gold1","tan","grey","darkblue","moccasin","hotpink4","ivory"))+
  facet_grid(.~sample, scales = "free", space = "free") +
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(labels = scales::percent)
cyano_bar.2021

#august
cyano_bar.aug <- NTC %>% filter(month == "August") %>%
  ggplot(aes(y=rel,x=site,fill=taxa))+geom_bar(stat="identity", position="fill", color="black")+labs(title="", x="",y="Relative abundance (%)")+
  theme_classic() +
  theme(legend.title=element_blank(),legend.position="right", axis.text = element_text(color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.text.x = element_text(angle = 90, hjust = 1, size = 11) )+
  scale_fill_manual(values=c("firebrick","lightblue2","pink","sienna1","lightgreen","gold1","tan","grey","darkblue","moccasin","hotpink4","ivory"))+
  facet_grid(.~sample, scales = "free", space = "free") +
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(labels = scales::percent)
cyano_bar.aug

##Eukaryotic Bar plots##
euks <- df_18s
euks$taxa <- "Other"
euks$taxa[is.na(euks$class)] <- "NA"
euks$taxa[euks$supergroup == "Amoebozoa"] <- "Amoebozoa"
euks$taxa[euks$phylum == "Apicomplexa"] <- "Apicomplexa"
euks$taxa[euks$class == "Bicoecea"] <- "Bicoecea"
euks$taxa[euks$phylum == "Cercozoa"] <- "Cercozoa"
euks$taxa[euks$phylum == "Chlorophyta"] <- "Chlorophytes"
euks$taxa[euks$class == "Chrysophyceae"] <- "Chrysophytes"
euks$taxa[euks$phylum == "Ciliophora"] <- "Ciliates"
euks$taxa[euks$phylum == "Metazoa"] <- "Metazoa"
euks$taxa[euks$phylum == "Fungi"] <- "Fungi"
euks$taxa[euks$phylum == "Cryptophyta"] <- "Cryptophytes"
euks$taxa[euks$class == "Bacillariophyta"] <- "Diatoms"
euks$taxa[euks$phylum == "Dinoflagellata"] <- "Dinoflagellates"
euks$taxa[euks$phylum == "Perkinsea"] <- "Perkinsea"

#aggregate the new taxonomy
NTE <- aggregate(euks$avg, by=list(site=euks$Site, sample=euks$Group, taxa=euks$taxa), sum)
NTE$Group <- NTE$sample
NTE <- NTE %>% separate(col=Group, into=c("month", "year"), sep="-")
NTE$rel <- NTE$x * 100
#factoring the variables
NTE$sample <- factor(NTE$sample, 
                     levels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"),
                     labels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"))
NTE$site <- factor(NTE$site,
                   levels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"), 
                   labels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"))
##phyla bar plots##

euks_bar.2021 <- NTE %>% filter(year == "2021") %>%
  ggplot(aes(y=rel,x=site,fill=taxa))+geom_bar(stat="identity", position="fill", color="black")+labs(title="", x="",y="Relative abundance (%)")+
  theme_classic() +
  theme(legend.title=element_blank(),legend.position="right", axis.text = element_text(color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.text.x = element_text(angle = 90, hjust = 1, size = 11) )+
  scale_fill_manual(values=c("firebrick2","lightblue","moccasin","tomato3","yellow3","magenta","pink","sienna1","forestgreen","gold1","tan4","darkblue","darkgrey","lightgrey","ivory"))+
  facet_grid(.~sample, scales = "free", space = "free") +
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(labels = scales::percent)


euks_bar.aug <- NTE %>% filter(month == "August") %>%
  ggplot(aes(y=rel,x=site,fill=taxa))+geom_bar(stat="identity", position="fill", color="black")+labs(title="", x="",y="Relative abundance (%)")+
  theme_classic() +
  theme(legend.title=element_blank(),legend.position="right", axis.text = element_text(color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line(), axis.text.x = element_text(angle = 90, hjust = 1, size = 11) )+
  scale_fill_manual(values=c("firebrick2","lightblue","moccasin","tomato3","yellow3","magenta","pink","sienna1","forestgreen","gold1","tan4","darkblue","darkgrey","lightgrey","ivory"))+
  facet_grid(.~sample, scales = "free", space = "free") +
  theme(strip.text.x = element_text(size = 12))+
  scale_y_continuous(labels = scales::percent)

#figure 4c and 4d
cyano_bar.2021 #fig.4c
euks_bar.2021 #fig.4d
#save figures
ggsave("fig.4c.pdf", plot=cyano_bar.2021, width = 8, height=5)
ggsave("fig.4d.pdf", plot=euks_bar.2021, width = 8, height=5)

#figure 8c and 8d
cyano_bar.aug
euks_bar.aug
#save figures
ggsave("fig.8c.pdf", plot=cyano_bar.aug, width = 8, height=5)
ggsave("fig.8d.pdf", plot=euks_bar.aug, width = 8, height=5)

