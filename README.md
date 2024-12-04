# Chicago-TNC-ABM-with-Bias
Materials used in paper "Ridehailing Technology Mitigates the Effects of Driver Racial Discrimination, but Effects of Residential Segregation Persist".
> Anna Cobb, Aniruddh Mohan, Corey Harper, Destenie Nock, and Jeremy Michalek
> Proceedings of the National Academy of Sciences 2024 121 (41) 
> DOI: 10.1073/pnas.2408936121 
AgentX is an agent-based model for simulating TNC operations originally reported in the study "Life Cycle Air Pollution, Greenhouse Gas, and Traffic Externality Benefits and Costs of Electrifying Uber and Lyft”.
> Aniruddh Mohan, Matthew Bruchon, Jeremy Michalek, and Parth Vaishnav
> Environmental Science & Technology 2023 57 (23), 8524-8535
> DOI: 10.1021/acs.est.2c07030

This version of AgentX builds on the original formulation with modifications to driver starting locations and repositioning behavior; it also includes code for simulating biased driver behavior towards customers and some changes to preprocessing procedures. Runs on Julia version 1.7.0, which can be downloaded from [here](https://julialang.org/downloads/oldreleases/).

## What's Included
### V4 Julia Environment:
The provided julia environment called V4 contains Project and Manifest .toml files; this environment must be activated and instantiated to run the preprocessing code, model, and postprocessing code. There are duplicates of the environment inside the Preprocessing and Postprocessing directories.
### Preprocessing Files [Directory]:
Due to file size restrictions on GitHub, only one of the trip input files used in the paper is included in the Input Files directory. This folder includes all scripts and instructions necessary to process raw trip data into a trip input file ready to be run. 
### Input Files [Directory]: 
One cleaned & pre-processed day of trips (07/26/21) ready to be run. Other files in this folder are for determining driver starting positions ("driver_home_dict.jld") and determining driver repositioning behavior (files in the "Repositioning Files" folder) during simulation. 
### Model: tnc_abm_race.jl
This is the version of the model used for generating the core results in the paper.
### Postprocessing Files [Directory]:
Includes a julia script for transforming results files into smaller, more manageable files that preserve important info and R scripts for visualizations of wait times.

## Usage
#### To run the model:
1. Download an OSM file for Chicago to the Input Files directory; the OSM file can be downloaded from [here](https://download.bbbike.org/osm/bbbike/Chicago/) (select the OSM XML gzip'd 142M option).
2. Download Julia 1.7.0 (linked above). 
3. Make sure you are working in Julia 1.7.0, then activate and instantiate the V4 environment.
4. Check that all paths are correctly specified in the tnc_abm_race.jl model file.
5. Run the model. On a 2021 MacBook Pro (M1 chip), estimated run time is 7 hours.
   
#### To run the preprocessing code to simulate a new date:
1. Follow instructions in doc called "Process to simulate a new date" in the Preprocessing directory.

#### To run the postprocessing code on a results file:
1. Use the converting_coords.jl script to transform raw simulation results files into what are internally referred to as "coordinate files" (smaller, refined versions of the raw simulation results that preserve important info). You will need to have the V4 environment activated for this.
2. Create visualizations of wait time results using R scripts.
