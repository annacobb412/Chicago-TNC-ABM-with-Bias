library(tidyverse)
library(dplyr)
library(choroplethrZip)
library(sf)
library(scales)
library(lubridate)
library(viridis)
library(plotrix)
library(extrafont)
library(sfheaders)
library(transformr)
theme_set(theme_bw())
# Date: 07/09/2024
# Author: Anna Cobb
# Description: creates map showing average simulated wait times by community
# area. Also has some code for adding circles centered on downtown to map. 

# inputs to modify
date = "09_12_22"
bias_level = 0.0
characteristic = "extra_dropoffs" # if one is specified in coordinate file
char_for_plot = "Wait Times"
file_dir = "Extra South Side Dropoffs Coord Files"

# reading in coordinate file
if (bias_level == 0) {
  input_file = paste0(file_dir,"/",date,"_",characteristic,"_0.0_coord.csv")
} else {
  input_file = paste0(file_dir,"/",date,"_",characteristic,"_B_",bias_level,"_coord.csv")
}

# reading in coordinate file
results = read.csv(input_file)

# processing
results = results %>% 
  rename(race = race_word) %>% 
  mutate(RP_strategy = "no_bias",WT_min = WT_sec/60) %>% 
  filter(race == "Black" | race == "White")
year = paste("20",substring(date,7,8),sep="")
formatted_date = paste(substring(date,1,6),year,sep="")
exceed_twenty_min = results %>% filter(WT_sec > 1200) %>% mutate(WT_length = "> 20 min")
exceed_one_hour = results %>% filter(WT_sec > 3600) %>% mutate(WT_length = "> 1 hour")
exceed_1.5_hour = results %>% filter(WT_sec > 5400) %>% mutate(WT_length = "> 1.5 hours")
exceed_2_hour = results %>% filter(WT_sec > 7200) %>% mutate(WT_length = "> 2 hours")
upper = rbind(exceed_twenty_min,exceed_one_hour,exceed_1.5_hour,exceed_2_hour)
lower = results %>% filter(WT_sec < 1200) 
percent_over_20 = nrow(exceed_twenty_min)/(nrow(lower) + nrow(exceed_twenty_min))
CA_sf = st_read("Boundaries - Community Areas/geo_export_ca2499f2-5635-434c-ab80-5247bfba606e.shp") %>% 
  mutate(area_num_1 = as.numeric(area_num_1))
CA_demo_data = read.csv("Community_Area_Demographics.csv")
dt_CA = CA_sf %>% filter(area_num_1 %in% c(8,32,33))

results_sf = st_as_sf(results,coords=c('start_lon','start_lat'),crs = st_crs(CA_sf))
fake_row = CA_sf[75,]
results_CA = results_sf %>% 
  mutate(intersection = as.integer(st_intersects(geometry,CA_sf))) %>% 
  drop_na(intersection) %>% 
  mutate(start_CA = CA_sf$area_num_1[intersection]) 

WT_CA = results_CA %>% group_by(start_CA) %>% summarize(WT_avg = mean(WT_min)) 
to_plot = CA_sf %>% st_join(WT_CA) %>% filter(!(area_num_1 %in% c(56,76)))

demo_CA = CA_sf %>% left_join(CA_demo_data,by = c("area_num_1" = "GEOID"))
black_CA = demo_CA %>% filter(BLACK > 0.6) 
white_CA = demo_CA %>% filter(WHITE > 0.5) %>% filter(!(area_num_1 %in% c(56,76)))
hisp_CA = demo_CA %>% filter(HISP > 0.5) %>% filter(!(area_num_1 %in% c(56,76)))

my_title = paste(gsub("_","/",date),": Average Wait Times")

CA_B_west = CA_sf %>% filter(area_num_1 %in% c(25,26,27,29)) %>% st_union()
ggplot(data = to_plot) +
  geom_sf(aes(fill = WT_avg)) +
  scale_fill_viridis_c(option = "plasma", limits = c(0,max(to_plot$WT_avg))) +
  geom_sf(data = CA_B_west,linewidth = .6, alpha = 0, color = "black") +
  ggtitle(my_title) +
  #geom_sf(data = dt_CA,linewidth = .6, alpha = 0, color = "red") +
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(axis.text=element_text(size=11), legend.text = element_text(size = 11)) +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"), axis.title = element_text(size = 13)) +
  scale_y_continuous(breaks = c(41.7, 41.8, 41.9, 42), limits = c(41.64, 42.03)) +
  scale_x_continuous(breaks = seq(-87.85,-87.55,by = .1)) +
  labs(fill = "Average\nWait Time\n[min]")

createCircle <- function(LonDec, LatDec, radius_m) {#Corrected function
  #LatDec = latitude in decimal degrees of the center of the circle
  #LonDec = longitude in decimal degrees
  #radius_m = radius of the circle in kilometers
  ER <- 3959 #Mean Earth radius in kilometers. Change this to 3959 and you will have your function working in miles.
  AngDeg <- seq(1,362, by = .25) #angles in degrees 
  Lat1Rad <- LatDec*(pi/180)#Latitude of the center of the circle in radians
  Lon1Rad <- LonDec*(pi/180)#Longitude of the center of the circle in radians
  AngRad <- AngDeg*(pi/180)#angles in radians
  Lat2Rad <-asin(sin(Lat1Rad)*cos(radius_m/ER)+cos(Lat1Rad)*sin(radius_m/ER)*cos(AngRad)) #Latitude of each point of the circle rearding to angle in radians
  Lon2Rad <- Lon1Rad+atan2(sin(AngRad)*sin(radius_m/ER)*cos(Lat1Rad),cos(radius_m/ER)-sin(Lat1Rad)*sin(Lat2Rad))#Longitude of each point of the circle rearding to angle in radians
  Lat2Deg <- Lat2Rad*(180/pi)#Latitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
  Lon2Deg <- Lon2Rad*(180/pi)#Longitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
  df = data.frame(Lat = Lat2Deg, Lon = Lon2Deg)
  return(df)
  # polygon(Lon2Deg,Lat2Deg,lty=2)
}
dt_lon = -87.627758
dt_lat = 41.875718
radius_12 = createCircle(dt_lon,dt_lat,12)
radius_8 = createCircle(dt_lon,dt_lat,8)
radius_4 = createCircle(dt_lon,dt_lat,4)

ggplot(data = to_plot) +
  geom_sf(aes(fill = WT_avg)) +
  scale_fill_viridis_c(option = "plasma", limits = c(0,max(to_plot$WT_avg))) +
  # geom_sf(data = black_CA,linewidth = .6, alpha = 0, color = "black") +
  ggtitle(my_title) +
  #geom_sf(data = dt_CA,linewidth = .6, alpha = 0, color = "red") +
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(axis.text=element_text(size=11), legend.text = element_text(size = 11)) +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.title = element_text(size = 13)) +
  scale_y_continuous(breaks = c(41.7, 41.8, 41.9, 42), limits = c(41.64, 42.06)) +
  scale_x_continuous(breaks = seq(-87.85,-87.45,by = .2)) +
  labs(fill = "Average\nWait Time\n[min]") +
  geom_point(aes(x = dt_lon, y = dt_lat)) +
  geom_point(data = radius_12, aes(x = Lon, y = Lat), size = .1, color = "black") +
  geom_point(data = radius_8, aes(x = Lon, y = Lat), size = .1, color = "black") +
  geom_point(data = radius_4, aes(x = Lon, y = Lat), size = .1, color = "black")
  