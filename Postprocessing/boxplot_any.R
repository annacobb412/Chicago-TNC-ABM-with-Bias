library(tidyverse)
library(dplyr)
library(colorBlindness)
# Date: 07/09/24
# Author: Anna Cobb
# Description: Generates a nice boxplot showing the distribution of wait times 
# for a single "coordinate file" (processed simulation results file). Can take
# any number of rider races.

# Inputs to modify for your singular file
date_of_interest = "01_13_22"
characteristic = "extra_dropoffs" # leave blank if no characteristic in file name
bias_ratio = "0.0" # fraction of biased drivers (written as decimal)\
dir_coord = "Extra South Side Dropoffs Coord Files/" # need slash

# Reading in coordinate file
c_coord_file = paste(dir_coord,date_of_interest,"_",characteristic,"_",bias_ratio,"_coord.csv",sep="")
customer_coord = read.csv(c_coord_file)

# Function to create boxplot for specific racial groups
# Inputs
# races: list of races of interest (should be capitalized--ex: c("Black","White))
# customer_coord: the read-in csv file in line 18
create_boxplot = function(races, customer_coord) {
  # Processing
  customers = customer_coord %>% 
    filter(race_word %in% races) %>% 
    mutate(WT_min = round(WT_sec/60,2)) %>% 
    rename(race = race_word)
  grouped_WT = customers %>% 
    group_by(race)
  WT_avgs = grouped_WT %>%  
    summarize(WT_avg = mean(WT_min))
  # Generic boxplot of customer wait times
  dodge = position_dodge(width=0.75)  
  title = paste0("Simulated Rider Wait Times: ", gsub("_","-",date_of_interest))
  ggplot() +
    geom_boxplot(data = grouped_WT, aes(x = race, y = WT_min, fill = race),outlier.shape = NA) +
    ylab("Wait Time [sec]") +
    xlab("Customer Race") +
    coord_cartesian(y = c(0,20)) +
    geom_point(data = WT_avgs, aes(x = race, y = WT_avg), position = dodge) +
    xlab("") +
    labs(fill = "Customer\nRace") +
    ggtitle(title) +
    theme_bw() +
    theme(axis.text=element_text(size=11), legend.text = element_text(size = 11)) +
    theme(text = element_text(family = "Times New Roman")) +
    theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"), 
          axis.title = element_text(size = 13)) 
  
}

create_boxplot(c("Black","White"),customer_coord)

