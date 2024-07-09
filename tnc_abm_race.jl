using Agents
using OpenStreetMapX
using Distributions
using Distributed
import DrWatson: @dict
using DataFrames
using CSV
using DataFramesMeta
using Base.Threads
using Random
using Dates
using JLD
using Pipe
using LightGraphs
using LinearAlgebra: dot
# Description: final version of ABM used in paper "Ridehailing Technologies Mitigate the 
# Effects of Driver Racial Discrimination, but Effects of Residential Segregation Persist"
# Date: July 9, 2024

# Inputs to modify:
list_of_dates = ["02_17_23","12_08_22"] # days of trips to simulate
input_dir = "Input Files/" # Chicago.osm file should also be in here
bias_ratio_vec = [0.1] # proportion of biased drivers
bias_vec = [2] # race of customer driver is biased against [1 = White, 2 = Black, 3 = Hispanic, 4 = Asian, 5 = other]
race_str_vec = ["B"] # initial of race of customer driver is biased agains (for output file name)
maxtime_in = 93600 # 26 hours

# Reading in input files
map_file = input_dir * "Chicago.osm"
locations_file = input_dir * "central_locations_input.csv"
dict_file = input_dir * "driver_home_dict.jld"

for date in list_of_dates
    for i in eachindex(bias_vec)

        # loading/configuring in additional files
        input_file = input_dir * "Input_ready_" * date * "_new.csv"
        date_regex = r"\d\d_\d\d_\d\d"
        date = match(date_regex,input_file).match
        year = ("20" * split(date,"_")[3])
        month = split(date,"_")[1]
        reposition_df = load(input_dir * "Repositioning V2 files/" * month * "-" * year * "_top_5_CAs_by_hour_no_AP.jld2","reposition_df")
        key = year * "-" * month
        key_lat_lon = key * "_latlon"
        driver_home_dict = load(dict_file,"driver_home_dict")
        monthly_driver_dist = driver_home_dict[key]
        latlon_dict = driver_home_dict[key_lat_lon]

        race_str = race_str_vec[i]
        bias_ratio = bias_ratio_vec[i]
        bias_race = bias_vec[i]
        SM1 = 1
        SM2 = 5
        code_start_time = Dates.format(now(),"yyyy-mm-dd HH:MM:SS")
        
        @agent Myagents OSMAgent begin
            status:: Int # (0 means available/waiting, 1 means trip, 2 means pickup allotted)
            paired:: Int # for customers and drivers: ID of agent you are paired with
            revenue:: Float64 # for drivers: earned revenue, for customer: waiting time until start of their trip  
            fuel:: Float64 # for drivers: fuel level of the vehicle, for customers: taxis that canceled
            target:: Float64  # for customers or drivers: decision on when to exit the game (if driver has earned enough money or customer has waited too long)
            queue:: Float64  # for customers: amount of time spent waiting to be assigned taxi, for drivers: drivers that need route computation
            deadhead:: Float64 # for drivers: number of miles spent deadheading (driving w/o customer in the car)
            vmt:: Float64 # for drivers: vehicle miles travelled 
            race:: Int # for drivers: bias (1 = true, 0 = false), for customers: race
            cancel:: Vector{Int} # for drivers: customers that they cancel on due to bias,
            start:: Int # for drivers: step at which they enter the game
            type:: Symbol # customer or taxi
        end

        max_soc = 14 # max gallons of fuel in a vehicle
        soc_dist = Uniform(0.3*max_soc, max_soc) # Uniform distribution of fuel level between 40% and 100% 
        bias_dist = Bernoulli(bias_ratio)
        earnings_dist = [220.0, 250.0, 65.0, 72.0, 60.0, 55.0, 80.0, 120.0, 150.0, 100.0]
        mpg = 29 * 0.9 # miles per gallon or fuel efficiency of vehicle, 25% reduction in winter, 10% in summer

        datainput = DataFrame(CSV.File(input_file)) # reading in trip data
        # creating output file name
        date_regex = r"\d\d_\d\d_\d\d"
        date = (match(date_regex,input_file)).match
        if length(string(bias_ratio)) < 3 # this is for post processing purposes; need bias ratio to have a "."
            str_BR = string(round(bias_ratio;digits = 2))
        else
            str_BR = string(bias_ratio)
        end
        output_file = "/ocean/projects/eng230003p/acobb/sensitivity results/Chicago_" * date * "_" * str_BR * "_" * race_str * "_no_AP_results.csv"
        if isfile(output_file) # change results file name in case it already exists
            output_file = output_file[1:end-4] * "_new.csv"
        end

        tripdata = datainput[datainput.test .< 1, :] # remove trips that are from same start to end points
        c_id_count = nrow(tripdata)
        if c_id_count > 220000
            model_ratio_lim = 4
        elseif c_id_count > 180000
            model_ratio_lim = 3
        else
            model_ratio_lim = 2
        end
        routes_filename = "routes_tracker.csv" # file that holds all precalculated routes

        Customer_m(id, pos, route, destination, status, paired, revenue, fuel, target, queue, deadhead, vmt, race, cancel, start) = Myagents(id, pos, route, destination, status, paired, revenue, fuel, target, queue, deadhead, vmt, race, cancel, start, :customer)
        Taxi_m(id, pos, route, destination, status, paired, revenue, fuel, target, queue, deadhead, vmt, race, cancel, start) = Myagents(id, pos, route, destination, status, paired, revenue, fuel, target, queue, deadhead, vmt, race, cancel, start, :taxi)
        
        function random_point_custom(bounds_specified) # generate random point within a box given a set of lat lon locations
            return (
                rand() * (bounds_specified.max_y - bounds_specified.min_y) + bounds_specified.min_y,
                rand() * (bounds_specified.max_x - bounds_specified.min_x) + bounds_specified.min_x
            )
        end

        function map_data_road(ll, m) # this is Anna's modified OSM.road function (uses map data instead of model as input)
            ll_enu = ENU(LLA(ll...), m.bounds)
            P = (ll_enu.east, ll_enu.north, ll_enu.up)
        
            idx = getindex(m.v, point_to_nodes(ll,m))
            idx_enu = m.nodes[m.n[idx]]
        
            candidates = Tuple{Tuple{Int,Int,Float64},Float64}[]
            
            if abs(ll_enu.east - idx_enu.east) > abs(ll_enu.north - idx_enu.north)
                # idx is the first position
                nps = LightGraphs.outneighbors(m.g,idx) # tells us all the... neighbors/vertices?
                isempty(nps) && return (idx, idx, 0.0)
                np_enus = map(np -> m.nodes[m.n[np]], nps) # extracting the ENU coord
                A = (idx_enu.east, idx_enu.north, idx_enu.up)
                short = nps
                if length(np_enus) < length(nps)
                    short = np_enus
                end
                for i in eachindex(short)
                    np_enu = np_enus[i]
                    np = nps[i]
                    B = (np_enu.east, np_enu.north, np_enu.up)
                    closest = OSM.orthognonal_projection(A, B, P)
                    candidate = (idx, np, distance(np_enu, ENU(closest...)))
                    push!(candidates, (candidate, sum(abs.(OpenStreetMapX.latlon(m,candidate[1]) .- ll))))
                end
            else
                # idx is the second position
                nps = LightGraphs.inneighbors(m.g,idx)
                isempty(nps) && return (idx, idx, 0.0)
                np_enus = map(np -> m.nodes[m.n[np]], nps)
                B = (idx_enu.east, idx_enu.north, idx_enu.up)
                short = nps
                if length(np_enus) < length(nps)
                    short = np_enus
                end
                for i in eachindex(short)
                    np_enu = np_enus[i]
                    np = nps[i]
                    A = (np_enu.east, np_enu.north, np_enu.up)
                    closest = OSM.orthognonal_projection(A, B, P)
                    candidate = (np, idx, distance(np_enu, ENU(closest...)))
                    push!(candidates, (candidate, sum(abs.(OpenStreetMapX.latlon(m,candidate[2]) .- ll))))
                end
            end
            bestidx = findmin(last.(candidates))[2]
            return first(candidates[bestidx])
        end

        function draw(input_distribution,latlon_dict,seedy) 
            Random.seed!(model.rng,seedy)
            zip = rand(input_distribution)
            zip_bounds = latlon_dict[zip]
            small_box_bounds_drop = OpenStreetMapX.Bounds{LLA}(zip_bounds[1], zip_bounds[2], zip_bounds[3], zip_bounds[4])
            random_pt_zip_bounds = random_point_custom(small_box_bounds_drop)
            random_road = OSM.road(random_pt_zip_bounds,model)
            return random_road
        end

        function draw_model(input_distribution,latlon_dict,seedy) 
            Random.seed!(model.rng,seedy)
            zip = rand(input_distribution)
            zip_bounds = latlon_dict[zip]
            small_box_bounds_drop = OpenStreetMapX.Bounds{LLA}(zip_bounds[1], zip_bounds[2], zip_bounds[3], zip_bounds[4])
            random_pt_zip_bounds = random_point_custom(small_box_bounds_drop)
            random_road = OSM.road(random_pt_zip_bounds,model)
            return random_road
        end

        #chicago map obtained from https://download.bbbike.org/osm/bbbike/Chicago/
        function initialize_model(; map_path = map_file, n_customers = 0, n_taxis =500, maxtime=maxtime_in) # max is a 24hours period
            m = get_map_data(
                map_path,
                use_cache = true,
                trim_to_connected_graph = true,
            )
            properties = @dict maxtime
            properties[:tick] = 0 # every second
            properties[:ratio] = round(Int, n_taxis/(n_customers+0.001)) 
            model = AgentBasedModel(Myagents, OpenStreetMapSpace(map_path); properties)
            print("model loaded")
            
            id = c_id_count
            for _ in 1:n_taxis
                id += 1
                start = draw(monthly_driver_dist,latlon_dict,id*SM1) # from driver data
                Random.seed!(model.rng,id*SM2)
                finish = OSM.random_position(model) # at a random position
                route = Int64[]
                Random.seed!(id)
                fuel = rand(soc_dist, 1)[1]
                Random.seed!(id)
                targets = rand(earnings_dist, 1)[1]
                Random.seed!(id)
                bias = rand(bias_dist)
                deadhead =  0       
                taxi = Taxi_m(id, start, route, finish, 0, 201, 0.0, fuel, targets , 0, deadhead, 0, Int(bias), [], 0)        # for taxi target is earning target for the day
                add_agent_pos!(taxi, model)
            end

            return model
        end

        model = initialize_model()
        show("Model initialized")

        space(agent) = OSM.map_coordinates(agent, model)
        level(agent) = agent.pos[3]

        function check_AP_CA(taxi) # check if a taxi is in the airport community area
            taxi_ll = OSM.latlon(taxi.pos,model)
            lat_up = 42.008295
            lat_low = 41.954232
            lon_up = -87.880861
            lon_low = -87.939216
            lat = taxi_ll[1]
            lon = taxi_ll[2]
            if lat > lat_low && lat < lat_up && lon > lon_low && lon < lon_up
                return true
            else
                return false
            end
        end 

        function near_customers_osm( # create a vector of agent IDs in nearby space that are NOT taxis
            agent,
            model::ABM{<:OpenStreetMapSpace},
            r, 
        )
            t1,t2 = OSM.map_coordinates(agent, model)
            # you first find location of agent, then search r in every directon and create a box
            # then in that box find all agents 
            x1 = t1 + r
            x2 = t1 - r
            y1 = t2 + r
            y2 = t2 - r

            n = zeros(0)

            allagentz = collect(allids(model))
            customers_near = filter!(id -> model[id].type == :customer && model[id].status == 0, allagentz)

            for idx in customers_near   #search for customers in box
                if space(model[idx])[1] < x1 && space(model[idx])[1] > x2 && space(model[idx])[2] >y2 && space(model[idx])[2] < y1 && !(model[idx].id in agent.cancel)
                # if space(model[idx])[1] < x1 && space(model[idx])[1] > x2 && space(model[idx])[2] >y2 && space(model[idx])[2] < y1
                    append!(n, idx)
                end
            end
            return n #vector of agent ids
        end

        function near_taxis_osm(   #create a vector of agent IDs in nearby space that are taxis
            agent,
            model::ABM{<:OpenStreetMapSpace},
            r, 
        )
            t1,t2 = OSM.map_coordinates(agent, model)
            # you first find location of agent, then search r in every directon and create a box
            # then in that box find all agents 
            x1 = t1 + r
            x2 = t1 - r
            y1 = t2 + r
            y2 = t2 - r

            n = zeros(0)
            allagentz = collect(allids(model))
            taxis_near = filter!(id -> model[id].type == :taxi && model[id].status == 0, allagentz)

            for idx in taxis_near   #search for customers in box
                if space(model[idx])[1] < x1 && space(model[idx])[1] > x2 && space(model[idx])[2] >y2 && space(model[idx])[2] < y1 && !(model[idx].id in agent.cancel)
                # if space(model[idx])[1] < x1 && space(model[idx])[1] > x2 && space(model[idx])[2] >y2 && space(model[idx])[2] < y1 
                    append!(n, idx)
                end
            end
            return n #vector of agent ids
        end

        function inverse(lon1, lat1, lon2, lat2, a, f)::Tuple{Float64,Float64,Float64}
            for lat in (lat1, lat2)
                abs(lat <= π/2) || throw(ArgumentError("Latitude ($lat) must be in range [-π/2, π/2]"))
            end
            a > 0 || throw(ArgumentError("Semimajor axis ($a) must be positive"))
            abs(f) < 1 || throw(ArgumentError("Magnitude of flattening ($f) must be less than 1"))
            lambda1, phi1, lambda2, phi2 = Float64(lon1), Float64(lat1), Float64(lon2), Float64(lat2)
            a, f = Float64(a), Float64(f)
            tol = 1.0e-8
            if (abs(phi2 - phi1) < tol) && (abs(lambda2 - lambda1) < tol)
                return 0.0, 0.0, 0.0
            end

            b = a*(1 - f)

            TanU1 = (1 - f)*tan(phi1)
            TanU2 = (1 - f)*tan(phi2)

            U1 = atan(TanU1)
            U2 = atan(TanU2)

            lambda = lambda2 - lambda1
            last_lambda = -4000000.0 # an impossibe value
            omega = lambda

            # Iterate the following equations until there is no significant change in lambda
            alpha, sigma, Sin_sigma, Cos2sigma_m, Cos_sigma, sqr_sin_sigma =
                -999999., -999999., -999999., -999999., -999999., -999999.
            while ((last_lambda < -3000000.0) || (lambda != 0)) &&
                    (abs((last_lambda - lambda)/lambda) > 1.0e-9)
                sqr_sin_sigma = (cos(U2)*sin(lambda))^2 +
                                ((cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda)))^2
                Sin_sigma = sqrt(sqr_sin_sigma)
                Cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda)
                sigma = atan(Sin_sigma, Cos_sigma)

                Sin_alpha = cos(U1)*cos(U2)*sin(lambda)/sin(sigma)

                if Sin_alpha >= 1
                    Sin_alpha = 1.0
                elseif Sin_alpha <= -1
                    Sin_alpha = -1.0
                end

                alpha = asin(Sin_alpha)
                Cos2sigma_m = cos(sigma) - 2*sin(U1)*sin(U2)/cos(alpha)^2
                C = (f/16)*cos(alpha)^2*(4 + f*(4 - 3*cos(alpha)^2))
                last_lambda = lambda
                lambda = omega + (1 - C)*f*sin(alpha)*(sigma +
                    C*sin(sigma)*(Cos2sigma_m + C*cos(sigma)*(-1 + 2*Cos2sigma_m^2)))
            end

            u2 = cos(alpha)^2*(a*a - b*b)/(b*b)
            A = 1 + (u2/16384)*(4096 + u2*(-768 + u2*(320 - 175*u2)))
            B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
            delta_sigma = B*Sin_sigma*(Cos2sigma_m + (B/4)*(
                Cos_sigma*(-1 + 2*Cos2sigma_m^2) -
                (B/6)*Cos2sigma_m*(-3 + 4*sqr_sin_sigma)*(-3 + 4*Cos2sigma_m^2)))
            s = b*A*(sigma - delta_sigma)

            alpha12 = atan((cos(U2)*sin(lambda)), ( cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda)))
            alpha21 = atan((cos(U1)*sin(lambda)), (-sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(lambda)))

            alpha12 = mod(alpha12, 2π)
            alpha21 = mod(alpha21 + π, 2π)

            return s, alpha12, alpha21
        end

        function surface_distance(lon0, lat0, lon1, lat1, a, f, degrees::Bool=true)
            if degrees
                lon0, lat0, lon1, lat1 = deg2rad(lon0), deg2rad(lat0), deg2rad(lon1), deg2rad(lat1)
            end
            distance, az, baz = inverse(lon0, lat0, lon1, lat1, a, f)
            distance
        end

        function distance_agents(a1,b1, model::ABM{<:OpenStreetMapSpace})
            agent1lat, agent1lon = OSM.latlon(a1, model)
            agent2lat, agent2lon = OSM.latlon(b1, model)
            a, f = 6.378137e6, 0.0033528106718309896
            dist_calc = surface_distance(agent1lon, agent1lat, agent2lon, agent2lat, a, f)
            return dist_calc
        end

        function nearest_customer_osm(
            agent::A,
            model::ABM{<:OpenStreetMapSpace,A},
            r;
        ) where {A}
            n = near_customers_osm(agent, model, r)
            length(n) == 0 && return nothing
            d, j = Inf, 1
            for i in 1:length(n)
                agentid = trunc(Int, n[i])
                @inbounds dnew = distance_agents(agent, model[agentid], model)
                optimal = round(dnew - model[agentid].revenue*8, digits =1)
                if optimal < d  
                    d, j = optimal, i
                end
            end
            nearestcustomer = trunc(Int, n[j])
            return model[nearestcustomer]
        end

        function nearest_taxi_osm(
            agent::A,
            model::ABM{<:OpenStreetMapSpace,A},
            r;
        ) where {A}
            n = near_taxis_osm(agent, model, r)
            length(n) == 0 && return nothing
            d, j = Inf, 1
            for i in 1:length(n)
                agentid = trunc(Int, n[i])
                @inbounds dnew = distance_agents(agent, model[agentid], model)
                optimal = round(dnew + model[agentid].revenue, digits =1)   ### factor in taxi's revenue when calculating nearest
                if optimal < d  
                    d, j = optimal, i
                end
            end
            nearesttaxi = trunc(Int, n[j])
            return model[nearesttaxi]
        end

        function modelid_unique(model,type)
            idset = collect(allids(model))
            idset_c = idset[idset .< c_id_count]
            idset_t = idset[idset .> c_id_count]
            if type == "customer"
                if idset_c == Int64[] # no customer agents have been created
                    newid = 1
                else
                    newid = maximum(idset_c) + 1
                end
            else # type = taxi
                newid = maximum(idset_t) + 1
            end
            return newid
        end

        function model_step!(model)
            model.tick += 8 # next step
            start_time = Time("02:30:00")
            mins_passed = floor(mod(model.tick,3600)/60)
            hours_passed = floor(model.tick/3600) # there should be a way to add an amount of time to a current time (e.g., 02:30:00 + 1.75 hours)
            #show(model.tick)
            if model.tick % 120 == 0
                print(start_time + Minute(mins_passed) + Hour(hours_passed))
                print(" ")
            end
            
            subset = tripdata[model.tick - 8 .< tripdata.dispatch_sec .<= model.tick, :] # only load in those trips that must be dispatched in these 8 seconds
            if model.tick % 30 == 0
                ids1 = collect(allids(model))
                agents_done = filter!(id -> model[id].status == 3, ids1)   # delete all agents who are status 3, trip finished or cancelled
                for a in agents_done
                    delete!(model.agents, a)
                end
            end
            new_activity!(model, subset)
        end

        function new_activity!(model, data)   # dispatch in new customers
            for i in 1:nrow(data)
                id = modelid_unique(model,"customer")
                startx = eval(Meta.parse(data[i,1])) # At an intersection 
                finishx = eval(Meta.parse(data[i,2])) # At an intersection  
                routex = eval(Meta.parse(data[i,4])) # route planned
                racex = data[i,"race"]
                wait_target = Inf
                customerx = Customer_m(id, startx, routex, finishx, 0, 201, 0.0, 0.0, wait_target, 0, 0.0, 0.0, racex, [], model.tick)
                add_agent_pos!(customerx, model)
            end

            ids = collect(allids(model))
            taxis = filter!(id -> model[id].type == :taxi && model[id].status == 0, ids)
            n_t = length(taxis) # Number of unmatched taxis

            ids = collect(allids(model))
            customers = filter!(id -> model[id].type == :customer && model[id].status == 0, ids)   
            n_c = length(customers) # Number of unmatched customers
            
            if model.tick % 120 == 0   # do this every 15 steps, or 2 mins, so it's not too sensitive and gives chance for taxis to come back online
                model.ratio = round(Int, n_t/(n_c+0.001))
                while model.ratio < model_ratio_lim   # taxi to customer is 1.5
                    id = modelid_unique(model,"taxi")
                    starty = draw_model(monthly_driver_dist,latlon_dict,id*SM1)
                    Random.seed!(model.rng,id*SM2)
                    finishy = OSM.random_position(model) # Somewhere on a road
                    routey = Int64[]
                    Random.seed!(id)
                    fuely = rand(soc_dist, 1)[1]
                    deadhead = 0
                    Random.seed!(id)
                    targetsy = rand(earnings_dist, 1)[1]
                    Random.seed!(id)
                    bias = rand(bias_dist)
                    taxiy = Taxi_m(id, starty, routey, finishy, 0, 201, 0.0, fuely, targetsy, 0, deadhead, 0.0, Int(bias), [], model.tick)
                    add_agent_pos!(taxiy, model)
                    n_t += 1
                    model.ratio = round(Int, n_t/(n_c+0.001))
                end
            end
        end

        function complex_step!(model)
            current_hour = mod(round(Int,floor(2.5 + model.tick/3600)),23) + 1 
            CAs_to_check = reposition_df[current_hour,["CA 1","CA 2","CA 3","CA 4","CA 5"]] 
            for id in model.scheduler(model)
                if model[id].type == :customer
                    customer_step!(model[id], model)
                elseif model[id].type == :taxi
                    try
                        taxi_step!(model[id], model, CAs_to_check)
                    catch e
                        print(id)
                        rethrow(e)
                    end
                end
            end
            allx = collect(allids(model))
            taxis = filter!(id -> model[id].type == :taxi && model[id].queue > 9900, allx)
            V = [Vector{Int64}() for _ in 1:length(taxis)]
            Threads.@threads for i in 1:length(taxis) # for those taxis that need a route calculation, calculate route and reset their queue number
                idx = trunc(Int, taxis[i])
                V[i] = OSM.plan_route(model[idx].pos, model[idx].destination, model)
                model[idx].route = V[i]
                model[idx].queue = model[idx].queue - 10000
            end
            model_step!(model)
        end


        function customer_step!(customer, model)
            prev_match = customer.paired # taxi ID that customer matched with previous step
            customer.queue += 1
            if customer.status == 0 # waiting to get paired with taxi
                customer.revenue += 1
                if customer.revenue > customer.target # longer than waiting target
                    customer.status = 6 # status 6 means they decided to abort plan to use a taxi
                end
                ans_c = nearest_taxi_osm(customer, model, 5000)
                if isnothing(ans_c) # no nearest taxi
                    ans_c = nearest_taxi_osm(customer, model, 10000) # expand search radius
                    if !isnothing(ans_c)
                        customer.paired = ans_c.id
                    end
                    if customer.paired != prev_match # this means we matched with a different taxi than previous round
                        customer.fuel = 1 # 1 will signify that taxi "cancelled"
                    end
                else # a nearest taxi is located
                    customer.paired = ans_c.id
                end
            elseif customer.status == 2
                customer.status = 3
            elseif customer.status == 4 # taxi is on the way to get customer
                customer.revenue += 1
                if customer.revenue > customer.target # longer than waiting target
                    customer.status = 6  # status 6 means they decided to abort plan to use a taxi
                end
            elseif customer.status == 6 # aborting taxi plan
                customer.deadhead = customer.deadhead + 1
                if customer.deadhead > 1
                    customer.status = 3
                end
            end
        end

        function taxi_step!(taxi, model, CAs_to_check)
            taxi.queue += 1
            if taxi.status == 0 # if available 
                if taxi.fuel <= 0.1 * max_soc # if less than 10% SOC you need to refuel
                    taxi.status = 5
                    taxi.fuel = 0 # so that each car takes 10 min to refuel, set this to a baseline fuel level
                    taxi.destination = taxi.pos
                    taxi.route = Int64[]
                else # you're not in emergency state
                    if taxi.revenue > taxi.target # revenue target reached! go home and refuel somewhere near home
                        taxi.status = 9 # time to go home
                        taxi.destination = taxi.pos
                        taxi.route = Int64[]
                    elseif taxi.queue >= 900 && (taxi.revenue*60)/(taxi.queue*8/60) < 10  # rev/hour is really low after 2 hours of working so just quit
                        taxi.status = 9 # time to go home
                        taxi.destination = taxi.pos
                        taxi.route = Int64[]
                    else # still in the game    
                        ans = nearest_customer_osm(taxi, model, 5000)   
                        if isnothing(ans) # did not find anyone in 5km radius
                            ans = nearest_customer_osm(taxi, model, 10000)  
                            if isnothing(ans) # did not find anyone in 10km radius
                                if taxi.route == Int64[]
                                    AP_CA = (123177, 123170, 214.11765590392517) # OSM coordinates
                                    cur_hour = floor(2.5 + model.tick/3600)
                                    CA_next = CAs_to_check[1]
                                    opt_distance = distance_agents(taxi,CA_next,model)
                                    for CA in 1:5
                                        distance = distance_agents(taxi,CAs_to_check[CA],model)
                                        if distance < opt_distance
                                            opt_distance = distance
                                            CA_next = CAs_to_check[CA]
                                        end
                                    end
                                    taxi.destination = CA_next
                                    # taxi.queue = taxi.queue + 10000
                                    #taxi.route = OSM.plan_route(taxi.pos, taxi.destination, model) 
                                else
                                    move_along_route!(taxi, model, 10*8) # you already have a plan so keep going
                                    taxi.fuel = taxi.fuel - (1/mpg)*0.0062*8 
                                    taxi.vmt = taxi.vmt + 0.0062*8
                                    taxi.deadhead = taxi.deadhead + 0.0062*8
                                end 
                            elseif  ans.paired == taxi.id  # 10 km search: you matched with someone!
                                if taxi.race == 1 && ans.race in bias_race  # driver is biased against customer
                                    push!(taxi.cancel, ans.id) # record customer ID in driver's cancel list
                                    push!(ans.cancel, taxi.id)
                                    taxi.status = 0 
                                else
                                    taxi.paired = ans.id # pair them up!
                                    taxi.destination = ans.pos
                                    # taxi.queue = taxi.queue + 10000
                                    taxi.status = 1 # taxi is now heading for pickup 
                                    ans.status = 4 # customer status changed to taxi assigned
                                end
                            elseif taxi.route == Int64[]  # you matched but with wrong customer
                                AP_CA = (123177, 123170, 214.11765590392517) # OSM coordinates
                                cur_hour = floor(2.5 + model.tick/3600)
                                CA_next = CAs_to_check[1]
                                opt_distance = distance_agents(taxi,CA_next,model)
                                for CA in 1:5
                                    distance = distance_agents(taxi,CAs_to_check[CA],model)
                                    if distance < opt_distance
                                        opt_distance = distance
                                        CA_next = CAs_to_check[CA]
                                    end
                                end
                                taxi.destination = CA_next
                                # taxi.queue = taxi.queue + 10000 
                            else 
                                move_along_route!(taxi, model, 10*8) # you already have a plan so keep going
                                taxi.fuel = taxi.fuel - (1/mpg)*0.0062*8
                                taxi.vmt = taxi.vmt + 0.0062*8
                                taxi.deadhead = taxi.deadhead + 0.0062*8
                            end
                        elseif ans.paired == taxi.id # 5km search: you matched with someone!
                            if taxi.race == 1 && ans.race in bias_race # driver is biased against customer
                                push!(taxi.cancel, ans.id)
                                push!(ans.cancel, taxi.id)
                                taxi.status = 0 
                            else
                                taxi.paired = ans.id # pair them up!
                                taxi.destination = ans.pos
                                # taxi.queue = taxi.queue + 10000
                                taxi.status = 1 # taxi is now heading for pickup or charging 
                                ans.status = 4 # customer status changed to taxi assigned
                            end
                        elseif taxi.route == Int64[] # you matched with wrong customer but head towards the nearest 5km customer anyway 
                            taxi.destination = ans.pos # go there anyway 
                            # taxi.queue = taxi.queue + 10000
                        else # you already have a plan, keep going 
                            move_along_route!(taxi, model, 10*8)
                            taxi.fuel = taxi.fuel - (1/mpg)*0.0062*8
                            taxi.vmt = taxi.vmt + 0.0062*8
                            taxi.deadhead = taxi.deadhead + 0.0062*8
                        end
                    end
                end
            elseif taxi.status == 1 # going to customer
                if model[taxi.paired].paired == taxi.id && model[taxi.paired].status == 4 # if customer is waiting for you 
                    if distance_agents(taxi, model[taxi.paired], model)<60 && distance_agents(taxi, model[taxi.paired], model)>2
                        move_along_route!(taxi, model, 10)
                        taxi.fuel = taxi.fuel - (1/mpg)*0.0062*8 
                        taxi.vmt = taxi.vmt + 0.0062
                        taxi.deadhead = taxi.deadhead + 0.0062*8
                    elseif distance_agents(taxi, model[taxi.paired], model)<=2
                        taxi.destination = model[taxi.paired].destination # pick up and journey to your customer destination 
                        model[taxi.paired].status = 1 # customer is picked up 
                        taxi.status = 2 # start of trip 
                    else   
                        move_along_route!(taxi, model, 10*8)
                        taxi.fuel = taxi.fuel - (1/mpg)*0.0062*8 
                        taxi.vmt = taxi.vmt + 0.0062*8
                        taxi.deadhead = taxi.deadhead + 0.0062*8
                    end
                else    
                    taxi.status = 0 # customer got picked up by someone else so new search
                    taxi.paired = 201 # reset pairing  
                    # taxi.deadhead = taxi.deadhead + 1 # JUST TRACK HOW MANY TIMES TAXIS HAD TO ABORT A PICKUP CUZ CUSTOMER CANCELLED ON THEM
                end
            elseif taxi.status == 2 # doing a trip for customer 
                if model[taxi.paired].pos == model[taxi.paired].destination  # if your taxi has reached the drop-off spot
                    taxi.status = 0 # taxi is now heading for pickup
                    model[taxi.paired].status = 2 # customer has been serviced 
                elseif taxi.route != model[taxi.paired].route # you've picked up customer and set on the way to destination!
                    xx = model[taxi.paired].route
                    taxi.route = xx
                else     # you're on the way
                    move_along_route!(taxi, model, 10*8)
                    model[taxi.paired].pos = taxi.pos
                    taxi.revenue = taxi.revenue + 0.0062*8*2.5   # 2.5 is revenue of ride per mile
                    taxi.fuel = taxi.fuel - (1/mpg)*0.0062*8
                    taxi.vmt = taxi.vmt + 0.0062*8
                end 
            elseif taxi.status == 3
                kill_agent!(taxi, model)
            elseif taxi.status == 5
                if taxi.fuel < max_soc
                    # taxi.fuel = taxi.fuel + ((max_soc - taxi.fuel)/600)*8 # regardless of SOC, will take 10 min to refuel
                    taxi.fuel = taxi.fuel + max_soc/75 # 75 = number of (8 sec) time steps in 10 min
                else
                    taxi.status = 0
                end
            elseif taxi.status == 9   # need to go home 
                taxi.status = 3   # kill the agent and post processing later
            end
        end

        steps = Int(floor(maxtime_in/8))
        avector = [0:5:steps;] # collects data on agents every so steps because we don't need to get data every step
        @time begin
            adata = [:type, :status, :pos, :destination, :paired, :revenue, :target, :deadhead, :vmt, :race, :cancel, :fuel, :queue, :start]
            mdata = [:ratio]
            data, modelxyz_data = run!(model, dummystep, complex_step!, steps; adata, mdata, when=avector) # we need 10800 steps to run full 24h period because 1 step = 8 seconds 
        end

        CSV.write(output_file, data)
        code_end_time = Dates.format(now(),"yyyy-mm-dd HH:MM:SS")
        println("Code began running: " * code_start_time)
        println("Code finished running: " * code_end_time)
    end
end

# @save "routes_dict.jld" routes_dict