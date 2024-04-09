# SWAP_regional

Author: Malve Heinz (malve.heinz@unibe.ch, malvemaria.heinz@agroscope.admin.ch)
Date: 02.04.2024

This repository is meant to help to run SWAP:

- from R  for multiple sites at once & (run_SWAP_in_R_for_multiple_sites.R, write_swapfile.swp)
- perform a global sensitivity analysis (Sensitivity_Analysis_SWAP.R, write_cropfile_for_SA.R)
- perform parameter optimization (DEoptim_parameter_optimization.R, write_cropfile_for_DEoptim.R)
- run SWAP on a regional scale and (regional_SWAP_application.R, write_swapfile_regional_application.R)
- prepare soil input data and (extract_Soilhydro_input_for_spatial_application.R)
- prepare cliamte input data for that (write_cliamte_inputfiles.R, template_climatedata.csv)

The basic structure is to execute the runswap.bat file and alter the input files (swap.swp and <CROP>.crp) from R.
The input files are altered by simply changing certain parameters and then writing the input files to your directory.

The set-up for the regional application is meant for a grid of cells, in my case potato fields in a specific catchment.
For each cell, first climate and soil input data and files are generated that can be called when SWAP is run for that cell.

Since the SWAP output is relative to the input data (like mm/m2) also the output is on that scale, so you might want to upscale the results. 

Beware that if your input data differs from mine, your units will to. check if you need to convert the to SWAPs input data requirements.

The setup for regional model application is meant to focus on the results for one year (but is run for a warm up period of three years + the year we are interest in).

You might consider the proposed folder structure, also provided in the repository (proposed_folder_structure.png)


