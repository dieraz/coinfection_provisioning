# code README for the manuscript: 
## Parasite responses to resource provisioning can be altered by within-host co-infection interactions 

Overview

This repository contains R code and Mathematica scripts for parameter estimation and analysis of two models: a population dynamics model and a coinfection model. The structure is as follows:
Directory Structure
1.	population_model/:
o	Contains all necessary scripts and data for estimating parameters of the population dynamics model.
o	Key files:
- main_population.R: The master script that initiates the population dynamics analysis.
- popdyn_det.R: Implements the population dynamics model.
- analysis_population.R: Analyzes the output from main_population.R.
- Data_population.RData: Contains population data.

2.	coinfection_model/:
o	Contains all necessary scripts and data for estimating parameters of the coinfection model.
o	Key files:
- main_coinfection.R: The master script that initiates the coinfection analysis.
- coinf_det_I9.R: Implements the I9 coinfection dynamics.
- coinf_det_L9.R: Implements the L9 coinfection dynamics.
- analysis_coinfection.R: Analyzes the output from main_coinfection.R.
- Data_coinfection.RData: Contains coinfection data.

3.	functions.R:
o	A collection of functions used across both models, promoting code reusability and modularity.
4.	plots.nb:
o	A Mathematica notebook used to generate Figures 4 and 5 from the manuscript, visualizing the results of the models.
Summary of Workflow
•	In each model directory, main.R serves as the main entry point and calls the detailed model scripts (popdyn_det.R for population dynamics, and coinf_det_I9.R / coinf_det_L9.R for coinfection).
•	After running the simulations, analysis.R processes the output data and performs further analysis.
•	The plots.nb notebook is then used to create Figures 4 and 5, displaying the model results graphically.
