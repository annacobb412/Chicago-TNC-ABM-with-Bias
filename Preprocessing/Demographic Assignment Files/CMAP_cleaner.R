library(dplyr)
library(miscTools)
library(tidyr)
library(ggplot2)
library(ggmap)

data_file = 'Chicago_Community_Areas_2022.csv'
demo_data = read.csv(data_file)
demo_data = demo_data %>% select(GEOID,GEOG,WHITE,HISP,BLACK,ASIAN,OTHER)

# get % of each race for each neighborhood -- there has got to be a prettier way but today is not the day to find out
demo_data$total_pops = rowSums(demo_data[,3:7]) 
demo_data = demo_data %>% 
  mutate(WHITE = WHITE/total_pops) %>% 
  mutate(HISP = HISP/total_pops) %>% 
  mutate(BLACK = BLACK/total_pops) %>% 
  mutate(ASIAN = ASIAN/total_pops) %>% 
  mutate(OTHER = OTHER/total_pops)

output_file = 'Community_Area_Demographics.csv'
write_csv(demo_data,output_file)
