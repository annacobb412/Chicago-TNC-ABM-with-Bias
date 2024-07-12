using CategoricalDistributions
using CategoricalArrays
using JLD
using Distributions
using CSV
using DataFrames
# Inputs: Community_Area_Demographics.csv
# Outputs: dictionary file ("CA_dictionary.jld")
# Description: creates dictionary with key = community area number, value = racial composition 
# of that community area (1 = White, 2 = Black, 3 = Hispanic, 4 = Asian, 5 = other)

CA_demo_file = "Demographic Assignment Files/Community_Area_Demographics.csv"
CA_demo = DataFrame(CSV.File(CA_demo_file)) # reading in community area level demographic info
sort!(CA_demo,:GEOID)

# creating dictionary with key = community area number, value = demographic distribution
CA_dict = Dict()
for i in 1:nrow(CA_demo) + 1
    if i == nrow(CA_demo) + 1
        local probs = [0.331,0.292,0.287,0.068,0.022] # city of Chicago demographic breakdown
        CA_dict[i] = Categorical(probs) # 1 = White, 2 = Black, 3 = Hispanic, 4 = Asian, 5 = other
    else
        local probs = [CA_demo[i,"WHITE"],CA_demo[i,"BLACK"],CA_demo[i,"HISP"],CA_demo[i,"ASIAN"],CA_demo[i,"OTHER"]]
        CA_dict[i] = Categorical(probs)
    end
end

@save "CA_dictionary.jld" CA_dict
