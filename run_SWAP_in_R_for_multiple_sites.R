####   Run SWAP for multiple location from R

# using this script, SWAP can also very simply be run for one location or with fixes soil input etc.
# i also need the script : write_swapfile.R

# Author: Malve Heinz
# Date: 02.04.2024

# required packages:

library(devtools)
library(glue)
library(Rdpack)
library(euptf2)
library(dplyr)
library(reticulate)
library(scales)
library(rlang)
library(hydroGOF)
library(ggplot2)
library(cowplot)
library(tidyr)
library(lubridate)
library(qpcR)
library(reshape2)
library(ggplot2)
library(TDPanalysis)
library(dplyr) 

#------------------------------  01 : data prep --------------------------------#

homedir<- "C:/Users/giub/Malve/Agrohydro/SWAP/SWAP_SwissIrrigation/"
# you could define a treatment that you represent with changes in the mainfile or cropfile
# treatment <- "default_potato_irr"
path_crop <- paste0("'..\\crop_files\\",treatment,"\\'")
 

# input data for several sites ( in my case i had reference data on sowing and harvesting dates, yilds, irrigation applications and soil texture of several sites)
bew<-read.csv("path_to_reference_data/bew.csv", sep = ";", fileEncoding="latin1")
bew_amount<-read.csv("path_to_reference_data_on_irrigation/bew_amount.csv", sep = ";", fileEncoding="latin1")
soil<-read.csv("path_to_reference_data_on_soil_hydraulic_parameters/soil_hydro_tex.csv")
observed_yield<-read.csv2("path_to_reference_data_on_yield/yield.csv", sep=";")

# If you alter the sowing and harvest dates, make sure to match SWAP requirements, e.g. "01-Oct-2018"

#------------------------------  01b : folder prep --------------------------------#

# If i want to run SWAP for several locations (each with a unique ID), I created one folder for each location within my "treatment" folder
for(r in unique(bew$ID)) {
  folder_path<-file.path(paste0(homdir,treatment,"/",r))
  if (dir.exists(folder_path)){
    unlink(folder_path,recursive=TRUE)
  }
  dir.create(paste0(homedir,treatment,"/",r), recursive=TRUE)
}

#_______________ 04: pre-defined input for SWAP main file  _____________________

station<-bew

for (x in 1:nrow(station)){ # loop over all stations
  
  #---------------------------  03 : prep input data -----------------------------#
  site<-station[x,]
  filename<-site[,1] # ID
  year<-site$Jahr
  
  # Cropsection
  cropmatrix<-as.data.frame(matrix(NA, ncol=6, nrow=1))
  cropmatrix[1,]<-c(2, paste0("15-Apr-",year), paste0("10-Nov-",year), "'Potato'", "'PotatoD'",2)
  colnames(cropmatrix)<-c("INITCRP","CROPSTART","CROPEND", "CROPNAME","CROPFIL","CROPTYPE")
  crp<- as.matrix(cropmatrix) 
  
  #_______________________________________________________________________________
  
  # soil textural data : here you alter the discretizaion of your soil profile (depending on the available soil input data)
  soilmatrix<-as.data.frame(matrix(data=NA, nrow=5, ncol=5))
  colnames(soilmatrix)<-c("ISUBLAY","ISOILLAY","HSUBLAY","HCOMP","NCOMP")
  soilmatrix$ISUBLAY<-c(1:5)
  soilmatrix$ISOILLAY<-c(1,1,1,2,3)
  G<-soil[x,49]
  
  if (soil[x,49]=="120"){
    soilmatrix$HSUBLAY<-format(c(5,10,15,30,60), nsmall=1)
    soilmatrix$NCOMP<-c(5,5,3,4,6)
    soilmatrix$HCOMP<-format(c(1,2,5,7.5,10), nsmall=1)
  } else {
    soilmatrix$HSUBLAY<-format(c(5,10,15,30,30), nsmall=1)
    soilmatrix$NCOMP<-c(5,5,3,4,4)
    soilmatrix$HCOMP<-format(c(1,2,5,7.5,7.5), nsmall=1)
  }
  
  soil<- as.matrix(soilmatrix)
  
  #-------------------------  04 : write SWAP text file --------------------------#
  
  # path where to store the txt file
  source(paste0(homedir,treatment, "/write_swapfile.R"), local=TRUE) # here i can then change the sowing/harvest dates or soil conditions for each site!
  pathy <-paste0(homedir,treatment,"/Swap.swp")
  
  #----------------------------- 06: run swap ------------------------------------
  
  setwd(paste0(homedir,treatment, "/"))
  # Delete old results
  system('DeleteOldResults.cmd')
  # Run SWAP 
  system.time(system('runswap.bat'))
  # set the path to your cliamte inout data, maybe also dependent on the site, if yes, that name with unique ID of that site!
  path_clim<-paste0("'..\\climate_HAFL_grid\\",x,"_PAY_clim'")  
  
  #--------------------  06 : store results in right folder ---------------------#
  
  # Take all swap results and store them in the right folder in order to continue the script: 
  my_files <- list.files(paste0(homedir,treatment),pattern="result.")
  my_files<-c(my_files,list.files(paste0(homedir,treatment), pattern="swap.ok"),list.files(paste0(homedir,treatment),pattern="heatparam"),list.files(paste0(homedir,treatment), 
                                                                                                                                                     pattern="soilphysparam"),list.files(paste0(homedir,treatment),pattern="swap_swap"))
  file.copy(from = paste0(homedir,treatment,"/", my_files), to = paste0(homedir,treatment,"/",filename,"/", my_files))
  file.remove(from = paste0(homedir,treatment,"/", my_files))
  
  
} 
  
  
#-------------------  12 : calculate overall performance ----------------------#


metrics_d<-matrix(NA, nrow = nrow, ncol=2) 
metrics_r<-matrix(NA, nrow = nrow, ncol=2) 
metrics_rmse<-matrix(NA, nrow = nrow, ncol=2) 
metrics_nse<-matrix(NA, nrow = nrow, ncol=2) 
metrics_d[,1]<-bew$ID
metrics_r[,1]<-bew$ID
metrics_rmse[,1]<-bew$ID
metrics_nse[,1]<-bew$ID
colnames(metrics_d)<-c("site","d")
colnames(metrics_r)<-c("site","r")
colnames(metrics_rmse)<-c("site","rmse")
colnames(metrics_nse)<-c("site","nse")


# Here are just some outputs I looked at, probably not all relevant

#----------------- 13 : observed vs simulated yield ----------------------#

# observed yield 
obs_yield<-observed_yield
# SWAP gives only the dry yield, but e.g. potato contains 73.8 – 81 % Water (Elbatawi et al. 2008)
# Default before: I take the mean of 77.4%. So the yield I have here makes up only 22,6% of the total yield
# NEW: I take 80%. So the yield I have here makes up only 20% of the total yield # lets try the old way
obs_yield$yield.dt.ha<-(obs_yield$yield.dt.ha/100)*23

# simulated yield 
sim_y<-matrix(NA, nrow = nrow, ncol=5)
sim_y[,1]<-obs_yield$site
colnames(sim_y)<-c("site","DVS","harvestdate","CWSO","CWPSO")


# This way I can read in my model outputs:
for (x in 1:nrow(station)){
  tab<-read.table(paste0(homedir,treatment,"/",station[x,1],"/result.crp"),skip=11, header=TRUE, sep=",")
  sim_y[x,2]<-tab[nrow(tab),4]
  sim_y[x,3]<-tab[nrow(tab),1]
  sim_y[x,4]<-tab[nrow(tab),21] # CWSO
  sim_y[x,5]<-tab[nrow(tab),20]} # CPWSO

sim_y<- as.data.frame(sim_y)
sim_yield<-data.frame(ID=sim_y[,1], yield=sim_y[,4]) 
sim_yield$yield<-as.numeric(sim_yield$yield)
sim_yield$yield<-sim_yield$yield/100 # from kg/ha to dt/ha

mean(as.numeric(sim_yield$yield))

potyield<-as.numeric(sim_y$CWPSO)/100 # from kg/ha to dt/ha
potyielddf<-data.frame(sim_pot_yield=potyield, observed=obs_yield$yield.dt.ha)
hydroGOF::d(sim=potyielddf$sim_pot_yield, obs=potyielddf$observed) 

obs_sim_yield<-cbind(sim_yield$yield, obs_yield$yield.dt.ha)
obs_sim_yield<-as.data.frame(obs_sim_yield)
colnames(obs_sim_yield)<- c("simulated","observed")

# now performance measure for yield
gof_yield<-gof(obs_sim_yield$simulated,obs_sim_yield$observed)
gof_yield[4] #rmse
gof_yield[10] #nse
gof_yield[12] #d  == 0.45
gof_yield[16] #r

gof_pot_yield<-gof(potyielddf$sim_pot_yield, potyielddf$observed)
gof_pot_yield[4] #rmse
gof_pot_yield[10] #nse
gof_pot_yield[12] #d  == 0.45
gof_pot_yield[16] #r

plot_yield<-ggplot(obs_sim_yield, aes(x=simulated, y=observed)) + geom_bin2d(bins = 20) + 
  ggtitle("simulated actual vs. observed yield in dt/ha") + scale_fill_continuous(type = "viridis") + 
  theme_bw() + geom_abline(intercept = 0, slope = 1) +  ylim(0,230) + xlim(0,230) +
  geom_text(x=10, y=220, label=paste0("rmse=",gof_yield[4]))+
  geom_text(x=10, y=210, label=paste0("nse=",gof_yield[10]))+
  geom_text(x=10, y=200, label=paste0("d=",gof_yield[12]))+
  geom_text(x=10, y=190, label=paste0("r=",gof_yield[16]))
plot_yield

plot_pot_yield<-ggplot(potyielddf, aes(x=sim_pot_yield, y=observed)) + geom_bin2d(bins = 20) + 
  ggtitle("simulated potential vs. observed yield in dt/ha") + scale_fill_continuous(type = "viridis") + 
  theme_bw() + geom_abline(intercept = 0, slope = 1) + ylim(50,180) + xlim(50,180)+
  geom_text(x=60, y=175, label=paste0("rmse=",gof_pot_yield[4]))+
  geom_text(x=60, y=165, label=paste0("nse=",gof_pot_yield[10]))+
  geom_text(x=60, y=155, label=paste0("d=",gof_pot_yield[12]))+
  geom_text(x=60, y=145, label=paste0("r=",gof_pot_yield[16])) 
plot_pot_yield

#----------------- 14 : observed vs simulated irrigation ----------------------#

obs_irr<-data.frame(bew$Zuordnungscode, bew$irr_count, bew$irr_amount)
colnames(obs_irr)<-c("site","count_o","amount_o (mm)")

sim_irr<-matrix(NA, nrow = nrow, ncol=3) # 36
sim_irr[,1]<-bew$Zuordnungscode
colnames(sim_irr)<-c("site","count_s","amount_s (mm)")

for (i in 1:nrow(station)){
  irr_sim<-read.csv(paste0("C:/Users/giub/Malve/Agrohydro/SWAP/SWAP_SwissIrrigation/",treatment,"/",
                           station[i,1],"/result.irg"), skip=9, header=TRUE, sep=",")
  sim_irr[i,2]<-nrow(irr_sim)
  sim_irr[i,3]<-sum(irr_sim$Irrigation)*10
}

sim_irr<-as.data.frame(sim_irr)

obs_sim_count<-cbind(as.numeric(sim_irr$count_s), as.numeric(obs_irr$count_o))
obs_sim_count<-as.data.frame(obs_sim_count)
colnames(obs_sim_count)<- c("simulated","observed")

gof_count<-gof(obs_sim_count$simulated,obs_sim_count$observed)
gof_count[4] #rmse
gof_count[10] #nse
gof_count[12] #d
gof_count[16] #r


obs_sim_amount<-cbind(as.numeric(sim_irr$`amount_s (mm)`), as.numeric(obs_irr$`amount_o (mm)`))
obs_sim_amount<-as.data.frame(obs_sim_amount)
colnames(obs_sim_amount)<- c("simulated","observed")

gof_amount<-gof(obs_sim_amount$simulated,obs_sim_amount$observed)
gof_amount[4] #rmse
gof_amount[10] #nse
gof_amount[12]  #d
gof_amount[16]  #r

plot_obs2<-ggplot(obs_sim_count, aes(x=simulated, y=observed)) + geom_bin2d(bins = 20) + 
  ggtitle("simulated vs. observed counts of irrigation applications") + scale_fill_continuous(type = "viridis") + 
  theme_bw() + geom_abline(intercept = 0, slope = 1) + ylim(0,10) + xlim(0,10) + 
  geom_text(x=1, y=10, label=paste0("rmse=",gof_count[4]))+
  geom_text(x=1, y=9, label=paste0("nse=",gof_count[10]))+
  geom_text(x=1, y=8, label=paste0("d=",gof_count[12]))+
  geom_text(x=1, y=7, label=paste0("r=",gof_count[16]))
plot_obs2

plot_obs3<-ggplot(obs_sim_amount, aes(x=simulated, y=observed)) + geom_bin2d(bins = 20) + 
  scale_fill_continuous(type = "viridis") + 
  ylim(0,225) + xlim(0,225)+ theme_bw() + geom_abline(intercept = 0, slope = 1)+ 
  geom_text(x=15, y=220, label=paste0("d=",gof_amount[12]))+
  geom_text(x=15, y=200, label=paste0("r=",gof_amount[16]))
plot_obs3

ggsave(paste0(filename="C:/Users/giub/Malve/Agrohydro/SWAP/SWAP_SwissIrrigation/",treatment,"/plots/irr_amount.png"), 
       plot=plot_obs3, width=4, height=3.3, dpi = "retina")


Z<-sim_y$V3[sim_y$V3==2]
X<-length(Z)/nrow
paste0("DVS 2 is reached in ", X ,"% of all (",nrow,") cases")

# difference between actual and potential yield
sim_y[,4]
sim_y[,3]

#--------------------------------    end     ----------------------------------#
