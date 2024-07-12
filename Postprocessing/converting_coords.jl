using DataFrames
using CSV
using Pipe
using OpenStreetMapX
using Geodesics
using Dates
using Agents
using LightGraphs
import DrWatson: @dict

# Description: takes a simulation results file and converts it to a smaller file,
# preserving info that is valuable for the project (e.g., wait times for customers and
# start positions for drivers). The new files created are referred to as "coordinate" 
# or "coord" files (they contain lat/lon info rather than the OSM coord used during 
# simulation).
# Date: 6-15-23
# Author: Anna Cobb

# to modify:
results_dir_name = "Extra South Side Dropoffs" # don't need slash at the end
coord_file_dir_name = "Extra South Side Dropoffs Coord Files"
chicago_osm_path = "/usr/path/to/Chicago.osm"
characteristic_output = "extra_dropoffs" # if running sensitivity case, can specify here
create_taxi_coord = false # output file with refined results for drivers
create_customer_coord = true # output file with refined results for customers

# ---------------------------------MODEL CODE---------------------------------------
# Note: model has to be instantiated to do the conversion from the ABM's OpenStreetMapSpace
# coordinate system used in the simulation to lat/lon for the coordinate file output.
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
    function initialize_model(; map_path = chicago_osm_path,maxtime = 10) # max is a 24hours period
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
    println("Model initialized")
end
# ---------------------------------MODEL CODE---------------------------------------
function model_step_to_time(step,time_step,date)
    beginning_hr = 2 # simulation starts at 2:30 AM but we've already cut off 2:30 - 3:30 AM
    beginning_min = 30
    start_hour = floor(step*time_step/3600) + beginning_hr
    start_min = floor(mod(step*time_step,3600)/60) + beginning_min
    start_sec = mod(step*time_step,3600)
    if start_min >= 60
        start_hour = start_hour + 1
        start_min = start_min - 60 
    end
    day = parse(Int,date[4:5])
    if start_hour > 23
        start_hour = mod(start_hour,24)
        day = day + 1
    end
    month = parse(Int,date[1:2])
    year = 2000 + parse(Int,date[7:8])
    new_date = Date(year,month,day)
    new_time = Time(start_hour,start_min,start_sec)
    return new_date,new_time
end

function model_tick_to_time(start_time,date)
    beginning_hr = 2 # simulation starts at 2:30 AM
    beginning_min = 30
    start_hour = floor(start_time/3600) + beginning_hr
    start_min = floor(mod(start_time,3600)/60) + beginning_min
    start_sec = mod(start_time,60)
    if start_min >= 60
        start_hour = start_hour + 1
        start_min = start_min - 60
    end
    day = parse(Int,date[4:5])
    if start_hour > 23
        start_hour = mod(start_hour,24)
        day = day + 1
    end
    month = parse(Int,date[1:2])
    year = 2000 + parse(Int,date[7:8])
    new_date = Date(year,month,day)
    new_time = Time(start_hour,start_min,start_sec)
    return new_date,new_time
end

function toOsmCoord(my_string)
    split_str = split(my_string,",")
    int1 = parse(Int,split_str[1][2:end])
    int2 = parse(Int,split_str[2])
    float3 = parse(Float64,split_str[3][1:end-1])
    OSM_coord = tuple(int1,int2,float3)
    return OSM_coord
end

function complex_filter(status)::Bool
    out = status == 1 || status == 2 || status == 3
    out 
end

list_of_files = readdir(results_dir_name)
if list_of_files[1] == ".DS_Store"
    list_of_files = list_of_files[2:end]
end

for i = 1:length(list_of_files)
    file = list_of_files[i]
    results = DataFrame(CSV.File(results_dir_name * "/" * file))
    date_reg = r"\d\d_\d\d_\d\d"
    date = (match(date_reg,file)).match
    bias_reg = r"\d.\d+_[A-Z]+"
    bias_long = (match(bias_reg,file)).match
    bias = (match(r"[A-Z]+",bias_long)).match
    bias_num = (match(r"\d.\d+",bias_long)).match
    time_step = 8

    # DAY PROCESSING: filter so we're only looking at 24 hours
    first_hour = floor(3600/time_step) # number of time steps in one hour
    last_hour = 25*first_hour # number of time steps in 25 hours
    filter!(row -> row.step >= first_hour && row.step < last_hour,results) # 0-450 = 2:30-3:30 AM, 11250-11700 = 3:30-4:30 AM
    filter!(row -> row.start >= first_hour*time_step && row.step < last_hour*time_step,results)

    ## TAXIS
    if create_taxi_coord
        taxis = @pipe results |> filter(:type => ==("taxi"),_)
        transform!(taxis, :step => ByRow(step -> step*time_step) => :time_step)
        
        date_and_time = model_tick_to_time.(taxis[:,:time_step],date)
        start_date_and_time = model_tick_to_time.(taxis[:,:start],date)
        start_time = []
        start_date = []
        current_time = []
        current_date = []
        for row in eachindex(date_and_time)
            push!(start_time,start_date_and_time[row,:][1,1][2])
            push!(start_date,start_date_and_time[row,:][1,1][1])
            push!(current_time,date_and_time[row,:][1,1][2])
            push!(current_date,date_and_time[row,:][1,1][1])
        end
        taxis[!,:start_time] = start_time
        taxis[!,:start_date] = start_date
        taxis[!,:current_time] = current_time
        taxis[!,:current_date] = current_date

        transform!(taxis, :pos => ByRow(pos -> toOsmCoord(pos)) => :pos_OSM)
        transform!(taxis, :pos_OSM => ByRow(coord -> OSM.latlon(coord,model)) => :pos_latlon)    
        transform!(taxis, :pos_latlon => ByRow(pair -> pair[1]) => :pos_lat)
        transform!(taxis, :pos_latlon => ByRow(pair -> pair[2]) => :pos_lon)

        transform!(taxis, :destination => ByRow(pos -> toOsmCoord(pos)) => :dest_OSM)
        transform!(taxis, :dest_OSM => ByRow(coord -> OSM.latlon(coord,model)) => :dest_latlon)    
        transform!(taxis, :dest_latlon => ByRow(pair -> pair[1]) => :dest_lat)
        transform!(taxis, :dest_latlon => ByRow(pair -> pair[2]) => :dest_lon)
        
        if bias == "NA"
            bias_output = ""
        else
            bias_output = "_" * bias
        end
        taxis_output = @pipe taxis |> select(_,[:pos_lat,:pos_lon,:dest_lat,:dest_lon,:start_time,:start_date,:current_time,:current_date,:status,:id,:step,:vmt])
        output_file_name = coord_file_dir_name * "/" * date * "_" * characteristic_output * bias_output * "_taxi_coord.csv"
        CSV.write(output_file_name,taxis_output)
        println("Done creating taxis file for " * date * "!")
    end

    ## CUSTOMERS
    if create_customer_coord
        customers_last = @pipe results |> filter(:type => ==("customer"),_) |> sort(_,:step,rev = true) |> unique(_, [:id]) |> sort(_,:id)
        customers_first = @pipe results |> filter(:type => ==("customer"),_) |> sort(_,:step) |> unique(_, [:id]) |> sort(_,:id)
        start_pos = customers_first[:,:pos]
        start_date_and_time = model_tick_to_time.(customers_first[:,:start],date)
        start_time = []
        start_date = []
        for row in eachindex(start_date_and_time)
            push!(start_time,start_date_and_time[row,:][1,1][2])
            push!(start_date,start_date_and_time[row,:][1,1][1])
        end
        customers_last[!,:start_str_OSM] = start_pos
        customers_last[!,:start_time] = start_time
        customers_last[!,:start_date] = start_date

        transform!(customers_last, :start_str_OSM => ByRow(pos -> toOsmCoord(pos)) => :start_OSM)
        transform!(customers_last, :start_OSM => ByRow(coord -> OSM.latlon(coord,model)) => :start_latlon)
        transform!(customers_last, :destination => ByRow(pos -> toOsmCoord(pos)) => :dest_OSM)
        transform!(customers_last, :dest_OSM => ByRow(coord -> OSM.latlon(coord,model)) => :finish_latlon)

        transform!(customers_last, :start_latlon => ByRow(pair -> pair[1]) => :start_lat)
        transform!(customers_last, :start_latlon => ByRow(pair -> pair[2]) => :start_lon)
        transform!(customers_last, :finish_latlon => ByRow(pair -> pair[1]) => :finish_lat)
        transform!(customers_last, :finish_latlon => ByRow(pair -> pair[2]) => :finish_lon)
        transform!(customers_last,:race => ByRow(race -> race == 1 ? "White" :
                                                race == 2 ? "Black" :
                                                race == 3 ? "Hispanic" :
                                                race == 4 ? "Asian" :
                                                race == 5 ? "Other" :
                                                string(x)) => 
                                                :race_word)
        transform!(customers_last,:revenue => (x -> x.*time_step) => :WT_sec)
        
        println("Number of customers total: " * string(nrow(customers_last)))
        customers_last = @pipe customers_last |> filter(:status => complex_filter,_)
        println("Number of customers served: " * string(nrow(customers_last)))
        customers_output = @pipe customers_last |> select(_,[:start_lat,:start_lon,:finish_lat,:finish_lon,:WT_sec,:race_word,:start_time,:start_date,:cancel,:id,:paired])
        if bias == "NA"
            bias_output = ""
        else
            bias_output = "_" * bias
        end
        output_file_name = coord_file_dir_name * "/" * date * "_" * characteristic_output * bias_output * "_" * bias_num * "_coord.csv"
        CSV.write(output_file_name,customers_output)
        println("Done creating customers file for " * date * "!")
    end
end