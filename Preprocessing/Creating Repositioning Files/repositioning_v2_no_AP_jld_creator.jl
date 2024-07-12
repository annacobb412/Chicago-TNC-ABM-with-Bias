using CSV
using DataFrames
using Dates
using Pipe
using Statistics
using JLD2
using OpenStreetMapX
using LightGraphs
using Agents
using OpenStreetMapX
using Geodesics
using DataFrames
import DrWatson: @dict
using Base.Threads
using Distributions
using Distributed
using DataFramesMeta
using Random

# Description: This script takes in a month's worth of trips that have been grouped by pickup community area 
# and start trip_start_timestamp (presumably using the "Endpoint_data_retrieval.py" script in this directory). It
# outputs a file that lists the OSM coordinates of the (centroids of) the 5 community areas with the highest
# trip counts for each hour of the day for that month. The created files are stored in the directory "Repositioning
# Files".
# Author: Anna Cobb
# Date (last checked/modified): 07/12/2024
# Notes: chicago OSM file can be downloaded from: https://download.bbbike.org/osm/bbbike/Chicago/

# INPUTS YOU MAY WANT/NEED TO MODIFY
dir_name = "TNC Trip Counts" # this is the directory from which all the CA files will be created
chicago_osm_path = " " # path to Chicago.osm

# Reading in files
list_of_dates = readdir(dir_name) # can make this just specific files if you don't want to work on all of them
centroid_latlon = DataFrame(CSV.File("CA_centroids.csv"))

# ------------------------MODEL CODE (gotta have model instantiated)-------------------------------
if !(@isdefined model)

    @agent Myagents OSMAgent begin
        status::Int # (0 means available/waiting, 1 means trip, 2 means pickup allotted)
        paired::Int # for customers and drivers: ID of agent you are paired with
        revenue:: Float64 # for drivers: earned revenue, for customer: waiting time until start of their trip  
        fuel:: Float64 # for drivers: fuel level of the vehicle
        target:: Float64  # for customers or drivers: decision on when to exit the game (if driver has earned enough money or customer has waited too long)
        queue :: Float64  # for customers: amount of time spent waiting to be assigned taxi, for drivers: drivers that need route computation
        deadhead :: Float64 # for drivers: number of miles spent deadheading (driving w/o customer in the car)
        vmt :: Float64 # for drivers: vehicle miles travelled 
        race :: Int # for drivers: bias (1 = true, 0 = false), for customers: race
        cancel :: Vector{Int} # for drivers: customers that they cancel on due to bias,
        type::Symbol # customer or taxi
    end

    function initialize_model(; map_path = chicago_osm_path,maxtime = 10) 
        m = get_map_data(
            map_path,
            use_cache = true,
            trim_to_connected_graph = true,
        )
        properties = @dict maxtime
        properties[:tick] = 0 # every second
        properties[:ratio] = 1 
        model = AgentBasedModel(Myagents, OpenStreetMapSpace(map_path); properties)
        return model
    end

    model = initialize_model()
    show("Model initialized")
end

# ---------------------------------MODEL CODE ENDS---------------------------------------

sort!(centroid_latlon,[:Community_Area])
if list_of_dates[1] == ".DS_Store"
    list_of_dates = list_of_dates[2:end]
end

for cur_date in list_of_dates

    regex_date = r"\d\d-\d+"
    the_date = match(regex_date,cur_date).match
    
    # generates random point within a box given a set of lat lon locations
    function random_point_custom(bounds_specified) 
        return (
            rand() * (bounds_specified.max_y - bounds_specified.min_y) + bounds_specified.min_y,
            rand() * (bounds_specified.max_x - bounds_specified.min_x) + bounds_specified.min_x
        )
    end

    # converts community area centroids (lat lon) to OSM coordinates
    function OSM_converter(coord)
        coord_bounds = [coord[1] - 0.0005,coord[1] + 0.0005, coord[2] - 0.0005, coord[2] + 0.0005]
        small_box_bounds_drop = OpenStreetMapX.Bounds{LLA}(coord_bounds[1],coord_bounds[2],coord_bounds[3],coord_bounds[4])
        random_pt_zip_bounds = random_point_custom(small_box_bounds_drop)
        random_road = OSM.road(random_pt_zip_bounds,model)
        return(random_road)
    end

    yearly_data = DataFrame(CSV.File(dir_name * "/" * cur_date))
    dropmissing!(yearly_data)
    rename!(yearly_data,"trip_start_timestamp" => :start_stamp,"pickup_community_area" => :pickup_CA, "COUNT_trip_ID" => :trip_count)
    transform!(yearly_data, :start_stamp => ByRow(x -> Dates.value(Hour(x))) => :start_hr)
    transform!(yearly_data, :start_stamp => ByRow(x -> Dates.dayofweek(x)) => :trip_day)
    yearly_data = yearly_data[findall(in([1,2,3,4,5]),yearly_data[!,:trip_day]),:] # filter for week days only
    yearly_data = filter(:pickup_CA => !=(76),yearly_data) # remove O'Hare airport pickups
    yearly_data = filter(:pickup_CA => !=(56),yearly_data) # remove Midway airport pickups

    # group by community area & hour
    grouped_no_days = groupby(yearly_data,[:start_hr,:pickup_CA])
    counted_no_days = combine(grouped_no_days, :trip_count => sum, keepkeys = true)
    # intermediate_file = "TNC_trips_grouped_" * the_date * ".csv"
    # CSV.write(intermediate_file, counted_no_days)
    top_5_no_days = combine(groupby(counted_no_days,[:start_hr])) do sdf
        first(sort(sdf,[:trip_count_sum],rev = true),5) end
    meep = @pipe top_5_no_days |> filter(:pickup_CA => ==(76),_)
    transform!(top_5_no_days, :pickup_CA => ByRow(x -> tuple(centroid_latlon[x,:Latitude],centroid_latlon[x,:Longitude])) => :coord)
    transform!(top_5_no_days, :coord => ByRow(x -> OSM_converter(x)) => :OSM_coord)

    # for visualization purposes (this is truly ugly code)
    long_no_days = top_5_no_days[!, Not(r"count_mean")]
    vec = ["CA 1", "CA 2", "CA 3", "CA 4", "CA 5"]
    long_no_days[:,:CA_rank] = repeat(vec,24)
    wide_no_days = unstack(long_no_days,[:start_hr],:CA_rank,:OSM_coord)
    
    # writing to jld file
    jld_file = the_date * "_top_5_CAs_by_hour_no_AP.jld2"
    save(jld_file,"Repositioning Files/reposition_df",wide_no_days)
end