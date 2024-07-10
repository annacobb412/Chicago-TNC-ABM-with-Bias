# Chicago-TNC-ABM-with-Bias
Agent-based model (ABM) used for running simulations of transportation network company (TNC) trips in Chicago with and without biased drivers.
Based on agent-based model created by Aniruddh Mohan and used in Life Cycle Air Pollution, Greenhouse Gas, and Traffic Externality Benefits and Costs of Electrifying Uber and Lyft.

## What's Included
Replication materials for Ridehailing Technology Mitigates the Effects of Driver Racial Discrimination, but Effects of Residential Segregation Persist
### Julia Environment:
The provided "V4" environment contains Project and Manifest .toml files; this environment must be activated and instantiated to run the preprocessing code, model, and postprocessing code.
### Preprocessing Files:
Due to file size restrictions, not all input files used in the paper are included in the Input Files directory. These, or new, input files can be created by following the instructions in the "Process to simulate a new date" doc. This document is very detailed and explains a lot: please read it and follow carefully. Contact me (annacobb@andrew.cmu.edu) if questions arise.
### Input Files: 
One cleaned & pre-processed day of trips (07/26/21) ready to be run. Other files in this folder are for determining driver starting positions ("driver_home_dict.jld") and determining driver repositioning behavior (directory "Repositioning V2 files") during simulation.
### Model:
This is the version of the model used for generating results in the paper.
### Postprocessing Files:
To run the postprocessing code, **simulation results files need to be converted to what are internally referred to as "coordinate files" (smaller, refined versions of the raw simulation results that preserve important info) using "converting_coords.jl"**. After this, the provided R scripts can create a variety of boxplot and map visualizations of wait times.
