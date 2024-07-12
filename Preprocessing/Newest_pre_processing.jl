using OpenStreetMapX
using Geodesics
using DataFrames
import DrWatson: @dict
using CSV
using Agents
using Base.Threads
using Dates
using Pipe
using JLD
using LightGraphs
using Distributions
using Distributed
using DataFramesMeta
using Random
# Author: Aniruddh Mohan + modifications by Anna
# Date (last modified/checked): 07/11/24
# Description: takes in raw Chicago TNC trip data (26 hours worth) and converts to 
# ABM input file.
# Notes: chicago OSM file can be downloaded from: https://download.bbbike.org/osm/bbbike/Chicago/

# Inputs to modify:
dir_name = "2019 TNC Data Inputs" # directory where raw trip data is stored (can be for multiple days)
dir_out_name = "2019 Input Files" # directory where the created model input files will go
chicago_osm_path = " " # OpenStreetMap file of Chicago

list_of_files = readdir(dir_name)
if list_of_files[1] == ".DS_Store"
    list_of_files = list_of_files[2:end]
end

# ---------------------------------MODEL CODE---------------------------------------
# Model must be instantiated to do the lat/lon coordinate conversions to OpenStreetMap
# space used in simulation
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

# ---------------------------------MODEL CODE---------------------------------------

function get_full_time(start_timestamp)
    time_regex = r"\d+:\d\d"
    the_time = (match(time_regex,start_timestamp)).match
    if length(the_time) < 5
        the_time = "0" * the_time * ":00"
    else
        the_time = the_time * ":00"
    end
    return(the_time)
end

function get_full_date(start_timestamp)
    d_regex = r"\d+/\d+/\d\d"
    the_date = (match(d_regex,start_timestamp)).match
    if length(the_date) < 8
        the_date = "0" * the_date 
    end
    the_date = SubString(the_date,1,6) * "20" * SubString(the_date,7,8)
    return(the_date)
end

for i in eachindex(list_of_files)
    input_file = dir_name * "/" * list_of_files[i]
    
    # Loading in input file & naming output file
    date_regex = r"\d\d_\d\d_\d\d"
    date = (match(date_regex,input_file)).match
    day = parse(Int,date[4:5])
    month = parse(Int,date[1:2])
    year = parse(Int,("20" * date[end-1:end]))
    actual_date = Date(year,month,day)
    output_file = "Input_ready_" * date * "_new.csv"
    cd = DataFrame(CSV.File(input_file)) # cd is for chicago data

    # Loading in appropriate home status dictionary based on weekday vs. weekend
    if Dates.dayofweek(actual_date) < 6 # weekday (monday-friday = 1-5)
        HS_dict = load("Demographic Assignment Files/HS_day_dict.jld","HS_day_dict")
    else # weekend
        HS_dict = load("Demographic Assignment Files/HS_end_dict.jld","HS_end_dict")
    end

    # Reading in demographics dictionary
    CA_dict = load("Demographic Assignment Files/CA_dictionary.jld","CA_dict")

    # Removing rows without pickup or dropoff centroids (or community areas)
    dropmissing!(cd,[:"Pickup Centroid Location",:"Dropoff Centroid Location"])
    filter!(:"Pickup Centroid Location" => !=(""),cd)
    pre = nrow(cd)
    dropmissing!(cd,[:"Pickup Community Area",:"Dropoff Community Area"])
    post = nrow(cd)
    num_NAs = pre - post

    # Filter out airport trips
    filter!(:"Pickup Community Area" => !=(76),cd) # O'Hare
    filter!(:"Dropoff Community Area" => !=(76),cd) # O'Hare
    filter!(:"Pickup Community Area" => !=(56),cd) # Midway
    filter!(:"Dropoff Community Area" => !=(56),cd) # Midway


    # Creating start hour column (most of this is converting to military time)
    cd.Is_it_PM = ifelse.(occursin.("PM",cd."Trip Start Timestamp"),1,0)
    cd.Is_it_PM = ifelse.(occursin.(" 12:",cd."Trip Start Timestamp"),0,cd.Is_it_PM) 
    if length(cd."Trip Start Timestamp"[1]) < 19
        transform!(cd, "Trip Start Timestamp" => ByRow(ts -> get_full_time(ts)) => :full_time)
        transform!(cd, "Trip Start Timestamp" => ByRow(ts -> get_full_date(ts)) => :full_date)
        cd.start_timestamp = cd.full_date .* " " .* cd.full_time
    else
        cd.start_timestamp = SubString.(cd."Trip Start Timestamp",1,19)
    end
    cd.start_datetime = DateTime.(cd.start_timestamp,dateformat"mm/dd/yyyy HH:MM:SS")
    hour_vec = ["12:00:00 AM","12:15:00 AM","12:30:00 AM","12:45:00 AM"]
    for i in 1:length(hour_vec)
        cd.start_datetime = ifelse.(occursin.(hour_vec[i],cd."Trip Start Timestamp"),cd.start_datetime - Dates.Hour(12),cd.start_datetime)
    end
    cd.start_datetime = cd.start_datetime + Dates.Hour(12)*cd.Is_it_PM # adding 12 hrs to PM times 
    cd.start_hr = Dates.hour.(cd.start_datetime)

    # Removing unnecessary columns
    select!(cd, Not(:"Trip ID"), Not(:"Pickup Census Tract"), Not(:"Dropoff Census Tract"), Not(:"Shared Trip Authorized"), Not(:"Trips Pooled"), Not(:"Is_it_PM"), Not(:"Trip Start Timestamp"), Not(:"Trip End Timestamp"))

    # Assigning dispatch sec
    start_time = @pipe cd |> unique(_,:start_datetime) |> select(_,:start_datetime) 
    start_time.for_dispatch = 0:nrow(start_time) - 1
    cd = leftjoin(cd,start_time,on = :start_datetime)
    cd.dispatch_sec = cd.for_dispatch .* 900 
    transform!(cd, :dispatch_sec => ByRow(x -> x + rand(1:900)) => :dispatch_sec) # this is a thing of beauty

    # Calculating average speed of trip (currently unused)
    cd.speed = round.((cd."Trip Miles" ./ cd."Trip Seconds") * 1609.34)

    # Create box around centroids for selecting precise location (this happens later)
    transform!(cd, :"Pickup Centroid Latitude" => ByRow( x -> x + 0.01) => :pickuplatmax)
    transform!(cd, :"Pickup Centroid Latitude" => ByRow( x -> x - 0.01) => :pickuplatmin)
    transform!(cd, :"Pickup Centroid Longitude" => ByRow( x -> x + 0.01) => :pickuplongmax)
    transform!(cd, :"Pickup Centroid Longitude" => ByRow( x -> x - 0.01) => :pickuplongmin)
    transform!(cd, :"Dropoff Centroid Latitude" => ByRow( x -> x + 0.01) => :droplatmax)
    transform!(cd, :"Dropoff Centroid Latitude" => ByRow( x -> x - 0.01) => :droplatmin)
    transform!(cd, :"Dropoff Centroid Longitude" => ByRow( x -> x + 0.01) => :droplongmax)
    transform!(cd, :"Dropoff Centroid Longitude" => ByRow( x -> x - 0.01) => :droplongmin)

    function random_point_custom(bounds_specified) # generate random point within a box given a set of lat lon locations
        return (
            rand() * (bounds_specified.max_y - bounds_specified.min_y) + bounds_specified.min_y,
            rand() * (bounds_specified.max_x - bounds_specified.min_x) + bounds_specified.min_x
        )
    end

    df3_latlons = DataFrame(a = Any[], b = Any[])
    df_positions_map = DataFrame(a = Any[], b = Any[], c = Any[])
    full_df = DataFrame(start = Any[], finish = Any[], test = Any[])
    race_vec = Vector{Any}(nothing,nrow(cd))
    
    for i in 1:nrow(cd)
        pickuplatmin = cd[i,"pickuplatmin"]
        pickuplatmax = cd[i,"pickuplatmax"]
        pickuplonmax = cd[i,"pickuplongmax"]
        pickuplonmin = cd[i,"pickuplongmin"]
        small_box_bounds_pickup = OpenStreetMapX.Bounds{LLA}(pickuplatmin, pickuplatmax, pickuplonmin, pickuplonmax)
        randomized_location = random_point_custom(small_box_bounds_pickup)
        
        randomized_intersection = OSM.intersection(randomized_location, model)

        droplatmin = cd[i,"droplatmin"]
        droplatmax = cd[i,"droplatmax"]
        droplonmax = cd[i,"droplongmax"]
        droplonmin = cd[i,"droplongmin"]
        small_box_bounds_drop = OpenStreetMapX.Bounds{LLA}(droplatmin, droplatmax, droplonmin, droplonmax)
        randomized_location_drop = random_point_custom(small_box_bounds_drop)
        
        randomized_intersection_drop = OSM.intersection(randomized_location_drop, model)

        test_var = 0
        if randomized_intersection == randomized_intersection_drop
            test_var = 1
        end

        push!(full_df, [randomized_intersection randomized_intersection_drop test_var])
        push!(df3_latlons, [randomized_location randomized_location_drop])

        # ANNA'S RACE ASSIGNMENT CODE
        start_hr = cd[i,"start_hr"]
        HS_dist = HS_dict[start_hr]
        HS_trip = rand(HS_dist) # assign home status of trip based on that start hour's distribution
        pickup_CA = cd[i,"Pickup Community Area"]
        dropoff_CA = cd[i,"Dropoff Community Area"]
        if pickup_CA == 76 # if picked up in O'Hare CA (airport)
            race_dist = CA_dict[78] # use general Chicago demographics
            race_vec[i] = rand(race_dist)
        elseif HS_trip == 1 || dropoff_CA == 76 # if dropoff was in O'Hare CA or designated as "from home" Trip
            race_dist = CA_dict[pickup_CA]
            race_vec[i] = rand(race_dist)
        elseif HS_trip == 2 # to home, so look at dropoff CA
            race_dist = CA_dict[dropoff_CA]
            race_vec[i] = rand(race_dist)
        else # neither, so look at entire city
            race_dist = CA_dict[78] # entry 78 is demographic breakdown of Chicago
            race_vec[i] = rand(race_dist)
        end
    end

    # Now generate a vector of all of the routes between start and end positions 
    V = [Vector{Int64}() for _ in 1:nrow(full_df)]
    Threads.@threads for i in 1:nrow(full_df)   # USE MULTI THREADING TO SPEED IT UP BUT THIS WILL TAKE A WHILE
        startc =  full_df[i,1]
        finishc = full_df[i,2]
        V[i] = OSM.plan_route(startc, finishc, model)
        println(i)
    end

    full_df.routes = V
    full_df.dispatch_sec = cd[:,"dispatch_sec"]
    full_df.pickup_CA = cd[:,"Pickup Community Area"]
    full_df.dropoff_CA = cd[:,"Dropoff Community Area"]
    full_df.race = race_vec
    full_df.trip_distance = cd[:,"Trip Miles"]
    sort!(full_df,:dispatch_sec)
    CSV.write(output_file, full_df)
end