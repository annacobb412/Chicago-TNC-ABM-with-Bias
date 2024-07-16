# Chicago-TNC-ABM-with-Bias
Materials used in paper "Ridehailing Technology Mitigates the Effects of Driver Racial Discrimination, but Effects of Residential Segregation Persist".
AgentX is an agent-based model for simulating TNC operations originally reported in the study "Life Cycle Air Pollution, Greenhouse Gas, and Traffic Externality Benefits and Costs of Electrifying Uber and Lyft”.
> Aniruddh Mohan, Matthew Bruchon, Jeremy Michalek, and Parth Vaishnav
> Environmental Science & Technology 2023 57 (23), 8524-8535
> DOI: 10.1021/acs.est.2c07030

This version of AgentX builds on the original formulation with modifications to driver starting locations and repositioning behavior; it also includes code for simulating biased driver behavior towards customers and some changes to preprocessing procedures. Runs on Julia version 1.7.0. 

## What's Included
### V4 Julia Environment:
The provided julia environment called V4 contains Project and Manifest .toml files; this environment must be activated and instantiated to run the preprocessing code, model, and postprocessing code. There are duplicates of the environment inside the Preprocessing and Postprocessing directories.
### Link to Chicago OpenStreetMap File:
To run the preprocessing code, model, and postprocessing code, you **must** download an OSM file for Chicago, which is linked [here](https://download.bbbike.org/osm/bbbike/Chicago/) (select the OSM XML gzip'd 142M option).
### Preprocessing Files:
Due to file size restrictions on GitHub, only one of the trip input files used in the paper is included in the Input Files directory. To recreate the trip input files used in the paper, or new ones, follow the instructions in the "Process to simulate a new date" doc. This document is detailed and explains a lot: read it and follow carefully. Contact me (annacobb@andrew.cmu.edu) if questions arise.
### Input Files: 
One cleaned & pre-processed day of trips (07/26/21) ready to be run. Other files in this folder are for determining driver starting positions ("driver_home_dict.jld") and determining driver repositioning behavior (files in the "Repositioning Files" folder) during simulation. 
### Model: tnc_abm_race.jl
This is the version of the model used for generating the core results in the paper.
### Postprocessing Files:
To run the postprocessing code, **simulation results files need to be converted to what are internally referred to as "coordinate files" (smaller, refined versions of the raw simulation results that preserve important info) using "converting_coords.jl"**. After this, the provided R scripts can create a variety of boxplot and map visualizations of wait times.
