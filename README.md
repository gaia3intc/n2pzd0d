# n2pzd0d
N2PZD marine plankton food web model in 0-dimensional, chemostat setting

README file for N2PZD0D model
June 23, 2025
Katsumi Matsumoto

The N2PZD0D model describes food web interaction between two functional types of phytoplankton and one functional type for zooplankton in a 0-dimensional, chemostat setting. Unique features are flexible phytoplankton C:N:P stoichiometry, homeostatic zooplankton stoichiometry, and grazing based on the stoichiometric imbalance between predator/zooplankton and prey/phytoplankton. N2PZD0D model consists of five MATLAB scripts:
•	N2PZD0D.m: main script
•	N2PZD0D_params.m: sets the parameter values
•	N2PZD0D_init.m: initializes the model
•	N2PZD0D_eqs.m: a function called y N2PZD0D.m that contains the ODEs
•	N2PZD0D_plots.m: makes some basic plots

These scripts have been modified after Sergio Vallina at the Spanish Institute of Oceanography. He developed these scripts as part of his educational/outreach activities.
