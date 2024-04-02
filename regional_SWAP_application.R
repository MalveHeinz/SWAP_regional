# This script provides the setup to run SWAP  on a regional scale by looping over gridcells
# and adjusting the main input file (swap.swp) to the pedo-cliamtic doncitions at each site/cell!

# In order to run this, the cliamte input files for each cell and the soil textural and hydraulic information must be gathered 
# and processed using the scripts "write_climate_inputfiles.R" and "extract_Soilhydro_input_for_spatial_application.R"

# The script: "write_swapfile_regional_application.R" is needed

# IMPORTANT: make sure, that you always use the same ID for one cell, make sure you have the soil and cliamte input data from the right location

# Author: Malve Heinz
# Date: 02.04.2024


#------------------ 01:  packages, paths, required data-------------------------

library(raster)
library(dplyr)
library(gdata)
library(tictoc)
library(ncdf4)
library(ggplot2)
library(tidyr)
library(profvis)

# _____________________________  01: set parameters ____________________________

area<-"EZG" # the mask area I am using (irr_area, EZG,...)
crop<-"potato"
yr<-"22"
ye<- "022"
year<-2022
bans<-0  #irrigation bans yes=1 or no=0 # I included the possibilit y of irrigation restriction by restricting the period where irrigation is possible. 0=no bans
path.base<-"C:your_path/"
climDir<- paste0(path.base,"meteofiles/")
resultDIR<- paste0(path.base,"Results_plots_ban=",bans,"/")

# _____________________________  02: read data  ____________________________

# soil textural and soil hydraulic input data -> preprocessed for the cells of interest
hydro.tex.all<-read.csv(paste0(path.base,"soil_data/soil_hydro_tex_",crop,yr,".csv"), sep=",")
kobo<-hydro.tex.all[,2:52]
subset_coords<-kobo[,50:51] # these are the coordinates for each cell
crop_coords<-subset_coords

#_______________ 03: containers to fill with results later _____________________

# Of course choose whatever is interesting for you!
# yield
sim_y<-matrix(NA, nrow = nrow(subset_coords), ncol=3)
colnames(sim_y)<-c("DVS","CWSO", "CPWSO")

# irrigation amount
sim_irr<-matrix(NA, nrow = nrow(subset_coords), ncol=2) 
colnames(sim_irr)<-c("count","amount(mm)")

# for transpiration
sim_Transp <- data.frame(matrix(NA, nrow=365, ncol=nrow(subset_coords)))

#  for drought stress
sim_Treddry <- data.frame(matrix(NA, nrow = 365, ncol = nrow(subset_coords)))


#_______________ 04: pre-defined input for SWAP main file  _____________________

# Cropsection
# I kept this static, but of course adjust to your needs!
cropmatrix<-as.data.frame(matrix(NA, ncol=6, nrow=1))
cropmatrix[1,]<-c(2, paste0("15-Apr-",year-3), paste0("10-Nov-",year-3), "'Potato'", "'PotatoD'",2)
cropmatrix[2,]<-c(2, paste0("15-Apr-",year-2), paste0("10-Nov-",year-2), "'Potato'", "'PotatoD'",2)
cropmatrix[3,]<-c(2, paste0("15-Apr-",year-1), paste0("10-Nov-",year-1), "'Potato'", "'PotatoD'",2)
cropmatrix[4,]<-c(2, paste0("15-Apr-",year), paste0("10-Nov-",year), "'Potato'", "'PotatoD'",2)
colnames(cropmatrix)<-c("INITCRP","CROPSTART","CROPEND", "CROPNAME","CROPFIL","CROPTYPE")

crp<- as.matrix(cropmatrix) 

#_____________________ 05: Loop to run over all cells  _________________________

##### define area and associated coordinates -> then loop through each pixel, row wise left to right (make sure the order is correct!)

tic()
#profvis({
for (x in 1:nrow(crop_coords)){  # this should correspond to the names I have in my climdata directory! very important

  
  soilmatrix<-as.data.frame(matrix(data=NA, nrow=5, ncol=5))
  colnames(soilmatrix)<-c("ISUBLAY","ISOILLAY","HSUBLAY","HCOMP","NCOMP")
  soilmatrix$ISUBLAY<-c(1:5)
  soilmatrix$ISOILLAY<-c(1,1,1,2,3)
  G<-kobo[x,49]
  
  if (kobo[x,49]=="120"){
    soilmatrix$HSUBLAY<-format(c(5,10,15,30,60), nsmall=1)
    soilmatrix$NCOMP<-c(5,5,3,4,6)
    soilmatrix$HCOMP<-format(c(1,2,5,7.5,10), nsmall=1)
  } else {
    soilmatrix$HSUBLAY<-format(c(5,10,15,30,30), nsmall=1)
    soilmatrix$NCOMP<-c(5,5,3,4,4)
    soilmatrix$HCOMP<-format(c(1,2,5,7.5,7.5), nsmall=1)
  }

  soil<- as.matrix(soilmatrix)
  
  source(paste0(path.base,"your_path/write_swapfile_regional_application.R"), local=TRUE)
  
  #----------------------------- 06: run swap ------------------------------------
  
  setwd(paste0(path.base,"regional_swap_default/"))
  # Delete old results
  system('DeleteOldResults.cmd')
  # Run SWAP 
  system.time(system('runswap.bat'))
  
  print(paste0("cell ",x," of ",nrow(crop_coords)))
  
  #-------------------------- 07: save results -----------------------------------
  
  ### yield
  # Read the text file
  data<-read.table(paste0(path.base,"regional_swap_default/result.crp"),skip=7, header=TRUE, sep=",")
  last_actual_yield <- tail(data$CWSO, 1) # Get the last day's yield -> automatically get the last year
  last_potential_yield <- tail(data$CPWSO, 1)
  DVS<-tail(data$DVS,1)
  sim_y[x,1]<- DVS
  sim_y[x,2]<- last_actual_yield/100 # kg in dt!
  sim_y[x,3]<- last_potential_yield/100 # kg in dt!
  
  ### Irrigation
  irr_sim<-read.table(paste0(path.base,"regional_swap_default/result.irg"), skip=9, header=TRUE, sep=",")
  # Extract the first four characters from the 'date' column
  first_four_digits <- substr(irr_sim$Date, 1, 4)
  # Filter rows where the first four digits match the target year
  irr_sim_year <- irr_sim[first_four_digits == year, ]
  sim_irr[x,1]<-nrow(irr_sim_year)
  sim_irr[x,2]<-sum(irr_sim_year$Irrigation) # I assume that relates to 1m2 (as also precipitation is given in mm/day/m2)
  
  ### transpiration
  sim_trsp<-read.table(paste0(path.base,"regional_swap_default/result.wba"), skip=6, header=TRUE, sep=",")
  first_four_digits <- substr(sim_trsp$Date, 1, 4)
  sim_trsp_year <- sim_trsp[first_four_digits == year, ]
  sim_trsp_year<-as.numeric(sim_trsp_year$Tact)
  sim_Transp[,x]<-sim_trsp_year[1:365] 
  
  ### drought stress
  sim_stress<-read.table(paste0(path.base,"regional_swap_default/result.str"), skip=6, header=TRUE, sep=",")
  first_four_digits <- substr(sim_stress$Date, 1, 4)
  sim_stress_year <- sim_stress[first_four_digits == year, ]
  DOY <-sim_stress_year$Day
  stress<-as.numeric(sim_stress_year$Treddry)
  start<-DOY[1]
  end<-start+nrow(sim_stress_year)-1
  sim_Treddry[start:end, x] <- stress

#--------------------------- 08. close loop ------------------------------------
  
  }

#--------------------------- 09. save results ----------------------------------
  
  # save results for further processing/visualization
  write.csv( sim_y, paste0(resultDIR, crop, "_yield_",year, "irrigationbans=", bans, ".csv"))
  write.csv( sim_irr, paste0(resultDIR, crop, "_irr_",year, "irrigationbans=", bans, ".csv"))
  write.csv( sim_Transp, paste0(resultDIR, crop, "_transp_",year, "irrigationbans=", bans, ".csv"))
  write.csv( sim_Treddry, paste0(resultDIR, crop, "_treddry_",year, "irrigationbans=", bans, ".csv"))
  
#  })

toc()


#----------------- 10: calc. irr demand per fraction area ----------------------

# results are in dt/ha dry yield and cm/m2 total amount 
# upscale e.g. seasonal irirgation amount from cm/m2 to mm per your grid cell area
irr_m<-sim_irr[,2]*10 # cm in mm -> zb 35mm pro m2 -> 35l/m2
irr_m_totalcell<-irr_m*900 # m2 of cell, in my case = 30x30 = 900

#------------------- 11: visualize and save results ----------------------------

write.csv(irr_m, paste0(resultDIR, "mean_irrigation_l_per_m2_",year,"bans=", bans,".csv"))
sim_y_df<-as.data.frame(sim_y)
write.csv(sim_y_df, paste0(resultDIR, "yield_dt_per_ha_",year,"bans=", bans,".csv"))

sim_y_df<-read.csv(paste0(resultDIR,"yield_dt_per_ha_",year,"bans=",bans,".csv"))
sim_y_df<-sim_y_df[,2:4]
sim_y_df$Ncell<-c(1:nrow(sim_y_df))
sim_y_df$wetyield<-(sim_y_df$CWSO/23)*100 # from dry to wet yield 
sim_y_df$wetpotyield<-(sim_y_df$CPWSO/23)*100 # from dry to wet yield 

colMeans(sim_y_df)

plot_yield<-ggplot(sim_y_df, aes(x=CWSO, y=CPWSO)) + geom_bin2d(bins = 20) + 
  ggtitle(title) + scale_fill_continuous(type = "viridis") + 
  theme_bw() + geom_abline(intercept = 0, slope = 1)  + ylim(60, 155) + xlim(60, 155) 
plot_yield

plotfile = paste0(resultDIR, "pot_act_yield_",year,"_",crop,"_bans=",bans,".png")
png(plotfile, width = 1800, height = 1400, res = 300)
print(plot_yield)
dev.off()

plot_actual_yield<-ggplot(sim_y_df, aes(x=Ncell, y=wetyield))+
  geom_point(alpha=0.3)+theme_bw()+ labs(x = "Ncell", y= "(wet) yield in dt/ha")
plot_actual_yield

plotf = paste0(resultDIR, "act_yield_",year,"_",crop,"_bans=",bans,".png")
png(plotf, width = 1800, height = 1200, res = 300)
print(plot_actual_yield)
dev.off()

irr_m<-read.csv(paste0(resultDIR, "/mean_irrigation_l_per_m2_",year,"bans=", bans,".csv"))
colnames(irr_m)<-c("Ncell", "mean_irr_l_m2")

plot_mean_irr<-ggplot(irr_m, aes(x=Ncell, y=mean_irr_l_m2))+
  geom_point(alpha=0.3)+theme_bw()+ labs(x = "Ncell", y= "mean irrigation amount in l/m2")
plot_mean_irr

plotfl = paste0(resultDIR, "mean_irr_",year,"_",crop,"_bans=",bans,".png")
png(plotfl, width = 1800, height = 1200, res = 300)
print(plot_mean_irr)
dev.off()

#---------------------------- END OF SCRIPT ------------------------------------
#_______________________________________________________________________________