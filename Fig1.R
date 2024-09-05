setwd("/Users/ishakalra/Desktop/Caron_lab_research/Clear_lake/manuscripts/amplicon_manuscript/R_figure_files/")
##Fig. 1###

##Fig. 1A - Clear Lake Map
library(leaflet)

# Set the latitude and longitude coordinates of Clear Lake
clear_lake_lat <- 39.0086
clear_lake_long <- -122.7610

# Create the map using leaflet
clear_lake_map <- leaflet() %>%
  addProviderTiles("Stamen.Watercolor") %>%
  addProviderTiles("Stamen.TonerHybrid") %>%
  setView(lng = clear_lake_long, lat = clear_lake_lat, zoom = 11)

# Display the map
clear_lake_map

##Fig. 1B - Microcystin levels
library(tidyverse)
library(ggplot2)
library(ggbreak)
theme_bw()
#read metadata sample file
sample_info <- read.csv("sample_info_decontam.csv")

#setting the sample group and site levels
sample_info$Group <- factor(sample_info$Group, 
                            levels = c("August_2019", "August_2020", "July_2021", "August_2021", "September_2021", "October_2021"),
                            labels = c("August-2019", "August-2020", "July-2021", "August-2021", "September-2021", "October-2021"))

sample_info$Site <- factor(sample_info$Site,
                           levels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"), 
                           labels = c("Upper Arm", "Soda Bay", "The Narrows", "Oaks Arm", "Lower Arm"))


#microcystin plots
mc <- sample_info %>%
  ggplot(aes(x= Site, y = Microcystin, fill = Group)) +
  geom_bar(stat="identity", position = "dodge", color = "black", show.legend = FALSE) +
  labs(y=expression("Microcystin ("~mu~"g/L)"), fill = "Sample", tag = "B") +
  scale_y_break(c(22, 50)) +
  geom_hline(yintercept = c(0.8, 6, 20), linetype = 2, color = "firebrick")+
  annotate("text", x= as.character("Upper Arm") ,y=1.50, label="Caution", color = "red")+
  annotate("text", x= as.character("Upper Arm") ,y=7.0, label="Warning", color = "red")+
  annotate("text", x= as.character("Upper Arm") ,y=21.0, label="Danger", color = "red")+
  theme_bw()
mc+chl

##Fig. 1C - chlorophyll graph

#chlorophyll plot
chl <- sample_info %>% group_by(Name) %>% mutate(chl = mean(Chlorophyll.a)) %>%
  ggplot(aes(x= Site, y = chl, fill = Group)) +
  geom_bar(stat="identity", position = "dodge", color = "black", show.legend = FALSE) +
  labs(y=expression("Chlorophyll-a ("~mu~"g/L)"), fill = "Sample", tag = "C") + 
  geom_hline(yintercept = 73, linetype = 2, color = "firebrick")+
  scale_y_continuous(breaks = c(0, 73, 100, 200, 300, 400))+
  theme_bw()

chl

#combine plots and display
fig1 <- mc+chl
ggsave("fig1-bc.pdf", plot=fig1)
fig1
