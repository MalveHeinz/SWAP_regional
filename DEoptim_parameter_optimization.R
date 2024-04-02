######### SWAP parameter optimization (here only for the crop file) using DEoptim

# Script written by Maria Eliza Turek and Malve Heinz
# Date: 27.03.2024

# define a representative soil input and crop settings (sowing/harvest dates)
# in your main input file you define how long etc the model is run and where the meteofiles and cropfile are located!
# read in observational data that can be used to maximize the fit (in this case, data on crop yield and irrigation amounts were used)
# to run this script, you need this as well: crop_files/write_cropfile_for_DEoptim.R

rm(list = ls(all.names = TRUE))

#_______________________________________________________________________________
#_______________________________________________________________________________

# Packages needed:
library(tidyverse)
library(tidyr)
library(dplyr)
library(rlang)
library(devtools)
library(glue)
library(Rdpack)
library(DEoptim)
library(parallel)
library(lubridate)
library(reshape2)
library(TDPanalysis)
library(tseries)
library(qpcR)
library(tictoc)
library(hydroGOF)
library(ggplot2)
library(scales)

# DEoptim:
# Mullen, K. M., Ardia, D., Gil, D. L., Windover, D., & Cline, J. (2011). DEoptim: An R Package for Global Optimization by Differential Evolution. Journal of Statistical Software, 40(6). 


# DEoptim settings:
homedir<-"C:/YOUR_DIRECTORY"
result_plots = paste0(homedir,"plots/")

### with 9 parameters and 100 iterations this takes me 4 days (16GB RAM)
iterations<-200 # maximum number of iterations in DEOptim -> default is 200
np=10 # number of parameters
NP<-10*np # number of populations/generation

#_______________________________________________________________________________
######    read in reference data  
#_______________________________________________________________________________

# reference data on yield, dont forget, SWAP gives you dry yield!
# transform wet yield (reference data) to dry yield, e.g. by assuming a water content of 73% [mean of 73.8 – 81 % Water (Elbatawi et al. 2008)]

# read in your reference data , e.g. on yield
# example:
obs_yield<-matrix(NA, ncol=2, nrow=10)
obs_yield[,1]<- c(2012:2021)
obs_yield[,2]<-rep(300,10)
colnames(obs_yield)<-c("year", "yield")

# create list to fill later with the model fit to this reference data for each run
d_yield_<- c()

#_______________________________________________________________________________
######    define DEoptim function  
#_______________________________________________________________________________

# Define our DEoptim function! This governs how the model is run and evaluated after each run, as well as which parameters are changed.
# You can add as many parameters as you want, but they should be passed as a vector and the order here matters, must correspond to the ranges below!

crop_opt <- function(st){     
  TSUMEA=       st[1] # the name should exactly correspond to how you name that paramter in the crop file
  FLTBb=        st[2] 
  FOTBc=        st[3]
  ALPHACRIT=    st[4]
  SLATB=        st[5] 
  AMAXTB=       st[6]
  FLTBc=        st[7]
  CF=           st[8]
  ADCRL=        st[9]
  RDC=          st[10]
  
  
  ## running SWAP
  # to run SWAP you need to have the runswap batchfile and the deleteoldresults.cmd in one directory
  # in your homedir, you will need to have the source folder (which is downloaded when you get SWAP)
  source(paste0(homedir,"crop_files/write_cropfile_for_DEoptim.R"), local=TRUE)
  setwd(paste0(homedir,"cali/"))
  system('DeleteOldResults.cmd')
  system(paste0(homedir,"cali/runswap.bat"))
  
  # After each run, here you read in the results of swap and extract the output you want to fit our model onto, in this case the made-up yield.
  # Estimated values from SWAP #
  y.s = read.csv(file = paste0("result.crp"), header = TRUE, as.is = TRUE, skip = 7)
  y.s$Date = as.Date(y.s$Date) # ,format="%d-%b-%Y")
  y.s$Year = format(y.s$Date, "%Y")
  
  # get the maximum dry actual yield for each year you simulate SWAP for
  max_values <- y.s %>% group_by(Year = lubridate::year(Date)) %>% 
    summarize(Max_CWSO = (max(CWSO, na.rm = TRUE))/100)# kg in dt 
  yield_sim_obs<-as.data.frame(max_values)
  sim_yield<-subset(FAT_sim_obs, (Year %in% years_obsFAT))
  sim_obs_yield = merge(sim_yield, obs_yield, by='Year')
  
  # you can do the same thin for Irrigation
  
  # (.....)

  # calulate your OF (optimization function) that the model aims to MINIMIZE
  # I use the metric d, index of agreement but you could choose whatever you like and adjust the OF accordingly
  d_act_yield<-hydroGOF::d(sim_obs_yield$sim_yield, sim_obs_yield$obs_yield)
  # save the intermediate results to izualize later
  d_yield <<- c(d_act_yield_FAT_values, d_act_yield_FAT)## with <<- i can append values vectors defined outside the function!
  
  metric_optim<- 1- d_act_yield # In this case, d should be ideally =1, so OF= 1-d is to be minimized!
  crop_opt = metric_optim 
  print(crop_opt)
  return(crop_opt)  
}

#_______________________________________________________________________________
######    run DEoptim and define parameter ranges 
#_______________________________________________________________________________

set.seed(1234) # important, because otherwise the results will differ for each run
tic() # measure time # If I want to try smth quickly, set VTR=0.8 (stop when d=1-0.8) and itermax=2 e.g.

# These are the min and max ranges I chose, but they can be adjusted accordingly. 
# Note, that in the WOFOST reference manual somtimes restrictions on the range in some parameters can be varied are listed
low<-0.85
up<-1.15# varying by +/- 15%

outDEoptim <- DEoptim(crop_opt, 
                      lower=c(150*low, # TSUM
                              1.05*low, #FLTBb
                              1*low,  # FOTBc
                              0.869*low, #ALPHAcrit 
                              0.003*low, #SLATB
                              30*low, #AMAXTB
                              0.83*low, #FLTBc
                              1*low, # CF
                              0.1*low, #ADCRH
                              50*low), #RDC
                      upper=c(150*up, 1.05*up, 1*up, 0.869*up, 0.003*up, 30*up, 0.83*up,1*up, 0.1*up, 50*up),
                      control = DEoptim.control(trace = FALSE,NP=NP, itermax=iterations)) #, initialpop=initial_matrix))#, storepopfrom=1, storepopfreq = 1,)) 
toc()

#_______________________________________________________________________________
######    visualize and save results
#_______________________________________________________________________________


devAskNewPage(FALSE) # so R does not require me anymore to press return to see the next plot!!!

# When the simulation is done, I can save and visualize the output
plot(outDEoptim, type="b")
summary(outDEoptim)
summary(outDEoptim)
summary(outDEoptim)

#dev.new()
plot(outDEoptim, plot.type="bestvalit", type="b", col="blue")
summary(outDEoptim)
# dev.new()
# plot(outDEoptim, plot.type="storepop")

outDEoptim$member$pop # last popultaion generated
last_pop<-as.data.frame(outDEoptim$member$pop)
write.csv(last_pop, file = paste0(result_plots, "last_population_",iterations,".csv"),row.names = F)

par(mfrow=c(1,1))
plot.fn<-plot(outDEoptim$member$bestvalit, type="l")
plot.path = paste0(result_plots, "best_function_value_at_each_iteration",iterations,".png")
png(plot.path, width = 1200, height = 800, res = 180)
print(plot.fn)
dev.off()

OF_values<-outDEoptim$member$bestvalit
write.csv(OF_values, file = paste0(result_plots, "OF_values_",iterations,".csv"),row.names = F)

# plot evolution
devAskNewPage(FALSE)
par(mfrow=c(1,1))
plot(d_yield, type="l", col=alpha("blue", 0.6), ylim=c(0,1), main="Evolution of d_yield (black)")
lines(d_yield, col=alpha("black", 0.6),lty=2)
write.csv(d_yield, file = paste0(result_plots, "d_yield_",iterations,".csv"),row.names = F)


#Save results
TSUMEA= as.numeric(outDEoptim$optim$bestmem[1])
FLTBb= as.numeric(outDEoptim$optim$bestmem[2]) 
FOTBc=  as.numeric(outDEoptim$optim$bestmem[3])
ALPHACRIT= as.numeric(outDEoptim$optim$bestmem[4])
SLATB= as.numeric(outDEoptim$optim$bestmem[5])
AMAXTB= as.numeric(outDEoptim$optim$bestmem[6])
FLTBc= as.numeric(outDEoptim$optim$bestmem[7])
CF= as.numeric(outDEoptim$optim$bestmem[8])
ADCRL= as.numeric(outDEoptim$optim$bestmem[9])
RDC= as.numeric(outDEoptim$optim$bestmem[10])


final = data.frame(TSUMEA=TSUMEA, FLTBb=FLTBb, FOTBc=FOTBc, ALPHACRIT=ALPHACRIT, SLATB=SLATB, AMAXTB=AMAXTB, FLTBc=FLTBc, CF=CF, ADCRL=ADCRL, RDC=RDC)
write.csv(final, file = paste0(result_plots,'TSUMEA_FLTBb_FOTBc_ALPHACRIT_SLATB_AMAXTB_FLTBc_CF_ADCRL_RDC.csv'),row.names = F)
#write.csv(allmetrics, file = paste0(result_plots,"d_irrHAFL_and_d_meanFATyield"))

#_______________________________________________________________________________
# You can now save you newly calibrated crop file and compare the results for running SWAP with and wihtout it!
### END
#_______________________________________________________________________________


