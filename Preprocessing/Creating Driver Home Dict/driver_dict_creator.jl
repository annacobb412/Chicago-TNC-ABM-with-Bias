using CSV 
using DataFrames
using Distributions
using JLD2
# INPUTS: all Driver_CA_Distribution files created in R (name of directory where these 
# are kept should stored in variable input_dir below)
# OUTPUTS: dictionary called driver_home_dict.jld
# DESCRIPTION:
# This code creates a dictionary [to be used in simulation] with two types of keys. One key is 
# the month and year for the date of trips being simulated and the other is the same month and
# year with "_latlon" appended. The plain month-year combo key yields a value that is the 
# distribution of driver homes across Chicago zip codes. The month-year + "_latlon" key yields
# a value that is another dictionary for that specific month-year. This second dictionary takes 
# as keys the numbers 1 - the total number of zip codes in Chicago. When indexed at a particular
# zip code number, the value it yields is a list of four coordinates, which represent a 0.04 by 
# 0.04 degrees box from which drivers from that zip code will appear in the simulation.

input_dir = "Driver CA Distributions"

driver_home_dict = Dict()

list_of_dist = readdir(input_dir)
list_of_dist = list_of_dist[2:end] # getting rid of .DS_Store

for date in list_of_dist
    println(date)
    driver_latlon_dict = Dict() # this dictionary resets each time, since it's different for each date

    date_key = match(r"\d\d\d\d-\d\d",date).match # this is the dictionary key "YYYY-MM"
    date_dist = DataFrame(CSV.File(input_dir * "/" * date))
    probabilities = date_dist[:,:percent]
    driver_home_dict[date_key] = Categorical(probabilities)

    # creation of the second dictionary, which stores the coordinates for each # associated with a zipcode
    for zip in 1:nrow(date_dist)
        driver_latlon_dict[zip] = [date_dist[zip,:lat] - 0.02,date_dist[zip,:lat] + 0.02, date_dist[zip,:lon] - 0.02, date_dist[zip,:lon] + 0.02]
    end

    # adding the secondary dictionary to the main dictionary
    latlon_dict_key = date_key * "_latlon"
    driver_home_dict[latlon_dict_key] = driver_latlon_dict
    println("complete")
end

@save "driver_home_dict.jld" driver_home_dict