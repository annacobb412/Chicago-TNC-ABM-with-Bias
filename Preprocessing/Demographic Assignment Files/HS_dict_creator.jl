using JLD
using Distributions
using CSV
using DataFrames
# INPUTS: csv with distribution of home statuses ("home_status_weekday.csv")
# OUTPUTS: dictionary with key = hour of the day (>23 means weekdend) &
# value = "home status" of trip distribution
# DESCRIPTION: creates dictionary from CSV

HS_day_file = "home_status_weekday.csv"
HS_day = DataFrame(CSV.File(HS_day_file)) # reading in community area level demographic info

# creating dictionariy
HS_day_dict = Dict()
day_hrs = HS_day[:,"start_hr"]
for i in 1:nrow(HS_day)
    local probs = [HS_day[i,"P_from"],HS_day[i,"P_to"],HS_day[i,"P_n"]]
    HS_day_dict[i-1] = Categorical(probs) # 1 = from home, 2 = to home, 3 = neither
end

@save "HS_day_dict.jld" HS_day_dict
