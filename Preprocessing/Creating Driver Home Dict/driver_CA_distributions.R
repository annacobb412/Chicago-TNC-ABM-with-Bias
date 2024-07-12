library(tidyverse)
library(dplyr)
library(choroplethrZip)
library(zipcodeR)
data(zip.regions)
# Inputs: Chicago TNC driver data 
# Outputs: Driver_CA_distribution file(s)
# Description: for specified month and year, creates file showing distribution
# of driver residences across Chicago zip codes. Note that the month and year 
# designate the time period in which these drivers were reported as being 
# registered to drive for a TNC in Chicago.
# Author: Anna Cobb
# Last modified/checked: 07/12/2024

# Inputs to modify
dir = "Driver CA Distributions" # where output files are stored
driver_data_file = "TNC_drivers_2023-2024.csv"

# Reading in & cleaning:
driver_data = read.csv(driver_data_file)
chicago_zips = zip.regions %>% 
  filter(county.name == "cook" & state.name == "illinois")

lats = c()
lons = c()
for (i in 1:nrow(chicago_zips)) {
  zip = chicago_zips$region[i]
  lats[i] = geocode_zip(zip)$lat
  lons[i] = geocode_zip(zip)$lng
}  
chicago_zips$lat = lats
chicago_zips$lon = lons

distrib_maker = function(date_reported,dir) {
  
  driver_data = driver_data %>% filter(CITY == "Chicago") %>% filter(MONTH_REPORTED == date_reported)
  by_zip = driver_data %>% group_by(ZIP)
  driver_count_zip = summarize(by_zip, driver_count = n())
  
  joined_zips = chicago_zips %>% 
    rename(ZIP = region) %>% 
    inner_join(driver_count_zip, by = "ZIP") %>% 
    select(ZIP, lat, lon, driver_count)
  total_drivers = sum(joined_zips$driver_count)
  joined_zips = joined_zips %>% mutate(percent = driver_count/total_drivers)
  
  file_name = paste0(dir,"/","Driver_CA_Distribution_",date_reported,".csv")
  write.csv(joined_zips,file_name)
  
  return(joined_zips)
}

# Months/Years that we need:
# Example: distrib_maker("2018-11",dir)
