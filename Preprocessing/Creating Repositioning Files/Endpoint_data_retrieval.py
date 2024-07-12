#!/usr/bin/env python
# Author: Anna Cobb + ChatGPT
# Date (last modified/checked): 07/11/24
# Description: this code interfaces with the Chicago Data Portal API to query, transform,
# and filter TNC Trip data to be used to inform driver agent repositioning behavior. The
# files it creates (and which are used as inputs to "repositioning_v2_no_AP_jld_creator.jl")
# are stored in the TNC Trip Counts directory. 

import pandas as pd
from sodapy import Socrata

# MODIFY:
my_year = '2023'
my_months = ['11','12']
end_days = ['30','31'] # number of days in the month

for ind in range(len(my_months)):
    my_year = my_year
    my_month = my_months[ind]
    end_day = end_days[ind]

    if my_year == '2018':
        my_endpoint = 'rpfj-eb3a'
    elif my_year == '2019':
        my_endpoint = 'iu3g-qa69'
    elif my_year == '2021':
        my_endpoint = 'unf9-2zu4'
    elif my_year == '2022':
        my_endpoint = 'm6dm-c72p'
    elif my_year == '2020':
        my_endpoint = 'rmc8-eqv4'
    elif my_year == '2023':
        my_endpoint = 'n26f-ihde'
    elif my_year == '2024':
        my_endpoint = 'aesv-xzh6'

    # Unauthenticated client only works with public data sets. Note 'None'
    # in place of application token, and no username or password:
    app_token = "NRZAZqgR2lywoiW7H1fALVgDa"
    client = Socrata("data.cityofchicago.org", app_token)
    client.timeout = 400
    start_timestamp = my_year + '-' + my_month + '-01T00:00:00'
    end_timestamp = my_year + '-' + my_month + '-' + end_day + 'T23:59:00'
    where_clause = f"trip_start_timestamp between '{start_timestamp}' and '{end_timestamp}'" # only get trips for one month at a time
    group_clause = "trip_start_timestamp,pickup_community_area" # group by time & community area
    select_clause = "COUNT(trip_ID),trip_start_timestamp,pickup_community_area" # get trip counts & other info
    results = client.get(my_endpoint, 
                        where= where_clause,
                        group= group_clause,
                        select= select_clause,
                        limit=10000000000)

    # Convert to pandas DataFrame
    results_df = pd.DataFrame.from_records(results)
    output_file = "TNC Trip Counts/TNC_trips_" + my_month + "-" + my_year + ".csv"
    results_df.to_csv(output_file)