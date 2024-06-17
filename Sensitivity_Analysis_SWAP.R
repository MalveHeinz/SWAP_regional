#####  Sensitivity Analysis (SA) using latin hypercube sampling (LHS) and Sobol indices

# Author: Malve HEinz
# Date: 02.04.2024

# Since SWAP has many parameters, optimizing all of them would take forever and would lead to overparameterization.
# Therefore, a SA is needed to identify those parameters that are most important to the model outputs you are interested in!

# we also require this script = crop_files/write_cropfile_for_SA.R

#_______________________________________________________________________________
### Packages needed
#_______________________________________________________________________________

library(lhs)
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
library(TDPanalysis)
library(sensobol)
library(tictoc)

#_______________________________________________________________________________
## Prior settings and notes
#_______________________________________________________________________________

# common error:

# running swap ....
# ERROR in rootextraction: 1 Max iterations for solver oxygen stress reached.
##### The trick is to increase the ISUBLAY 3 from 1x15cm to 2x7.5cm!

# If I want to know how long that takes
tic()

## irrigated or not? If yes, the irirgation module within the cropfile is activated and the model irrigates as defined there!
irr<-1 # 0=no, 1= yes

homedir<-"C:/Users/giub/Malve//AgroHydro/LHS_sensitivity_combined/KOBO/"

N<-3000 # number of samples, you can start with a smaller one, but a larger number gives more robsut results
samplesize<-paste0("sample",N)
R<-200 # number of bootstrap replicas, default = 200

# define path to cropfile which is overwritten at each iteration
pathy <-paste0(homedir, "samples/crop_files/secondbatch/PotatoD.crp") 

# create mat_lists to store the results in later
mat_listyield17 <- list()
mat_listyield18 <- list()
mat_listyield19 <- list()
mat_listyield_mean <- list()
mat_listirr_totaldepth <- list()

#_______________________________________________________________________________
### define the ranges to vary the parameters in question in
#_______________________________________________________________________________

#---------   parameters to look at with a +/-15 % range 
low = 0.85
up = 1.15

#----------------  A  ####  potential production =
CF<-seq(1*low,1*up, 0.01) # DVS 0 # instead CF2 will be CF1 * 1.1 !!!
RSC<-seq(100*low,100*up,0.1)
TDWI<-seq(75*low,75*up, 0.1)
RGRLAI<-seq(0.012*low,0.012*up, 0.0001)
SPAN<-seq(37*low,37*up,0.1)
SLATB<-seq(0.003*low,0.003*up,0.00001) # DVS: 0 und 1.1 (same)-> instead SLATB2 witll be SLATB1/2!!
AMAXTB<-seq(30*low,30*up,0.1) # DVS: 0 -> same DVS 0 and 1.5, (2=0)
Q10<-seq(2*low,2*up,0.01)
#-----------------  B  ####  actual production (limiting parameters)   =
TMPFTB<- seq(0.75*low,0.75*up,0.001) ### OBACHT can only be between 0 and 1!!! # here I only change the value at temp=10 and 26, because the shape should remain the same
Q10_microbial<-seq(2.8*low,2.8*up,0.001)
resp_humus<-seq(0.0016*low,0.0016*up,0.00001)
SRL<-seq(151375*low,151375*up,1)
root_radius<-seq(0.00015*low,0.00015*up,0.000001)
HLIM3H<- seq(-300*up,-300*low, 1) ### OBACH bei nagative values i need to switch up and low
HLIM3L<- seq(-500*up,-500*low, 1)
HLIM4<- seq(-10000*up,-10000*low, 1)
ADCRH<- seq(0.50*low,0.50*up,0.001)
ADCRL<- seq(0.1*low,0.1*up,0.001)
ALPHACRIT<- seq(0.84*low,0.84*up,0.001) # actually the default is 1.0, but that is the upper limit already
RRI<-seq(1.2*low,1.2*up, 0.001)
RDC<-seq(50*low,50*up, 0.1)
RDCTB<- seq(-0.38,-0.4155, -0.0001) ### This is less then 10 %, but it is the possible max to alter this value! # equation: (RDCTB*ln(x) + 0.9966)

# 10% deviation from default variables for modern potato cultivars from Den et al 2022
FRTB<-seq(0.2*low,0.2*up,0.001)# same for both DVS
FOTBa <- seq(7*low,7*up,0.01) 
FOTBb <- seq(1.05*low,1.05*up, 0.01) 
FOTBc <- seq(1*low,1*up, 0.01)
FLTBa <- seq(20*low,20*up,0.01)
FLTBb <- seq(1.05*low,1.05*up, 0.01)
FLTBc <- seq(0.83*low,0.83*up,0.01)

params=c("CF","RSC","TDWI","RGRLAI","SPAN","SLATB","AMAXTB","Q10","TMPFTB","Q10_microbial","resp_humus",
         "SRL","root_radius", "HLIM3H","HLIM3L","HLIM4","ADCRH","ADCRL","ALPHACRIT","RRI","RDC","RDCTB",
         "FRTB", "FOTBa" ,"FOTBb","FOTBc","FLTBa","FLTBb","FLTBc") 


#_______________________________________________________________________________
#### take a LHS sample from my hypercube
#_______________________________________________________________________________
# which type of matrix (in this case A) depends on the nature of my sampling= do i have many and interconnected paramters? -> A might not be the best choice
lhs_test<-sobol_matrices(matrics="A", N=N, params=params ,type = "LHS")
lhs_test<-sobol_matrices(N=N, params=params ,type = "LHS")


#_______________________________________________________________________________
### Now get the actual values for the sample taken!
#_______________________________________________________________________________

###------------ 03. from normalized sample to parameter values  -------------###

# assuming "parameter" is the vector of parameter values assigned in the ebginning of this^script (e.g. CF)
# This is out syntax --> parameter [lhs_test[ ,column of parameter] * length(parmeter) + 1]
# first we extract the column within our latin hypercube sample where the normalized latin hypercube sample (LHS) value of 
# the parameter of interest are located (between 0 and 1) --> lhs_test[ ,column of parameter]
# Then we multiply each element in this column by the length of our parameter --> lhs_test[, 1] * length(CF)
# As lenght(parameter) give us the number of elements in this sequence, this scaling maps the normalized lHS values to indices
# within the range of the parameter vector. The + 1 ensures, that the indices are in the range of 1 to length(CF)

# e.g. lhs_test[1, 1] = 0.7218774 ;  0.7218774 * length(CF) = 67134.6 ; 67134.6 + 1 = 67135.6; CF[67135.6]  = 1.03

CF <- CF[lhs_test[, 1] * length(CF) + 1]
RSC <- RSC[lhs_test[, 2] * length(RSC) + 1]
TDWI <- TDWI[lhs_test[, 3] * length(TDWI) + 1]
RGRLAI <- RGRLAI[lhs_test[, 4] * length(RGRLAI) + 1]
SPAN <- SPAN[lhs_test[, 5] * length(SPAN) + 1]
#SPAN <- paste0(as.character(SPAN), ".0")
SLATB <- SLATB[lhs_test[, 6] * length(SLATB) + 1]
AMAXTB <- AMAXTB[lhs_test[, 7] * length(AMAXTB) + 1]
Q10<- format(as.numeric(Q10[lhs_test[,8] * length(Q10) + 1]),nsmall=1) 
TMPFTB <- TMPFTB[lhs_test[, 9] * length(TMPFTB) + 1]
Q10_microbial<- Q10_microbial[lhs_test[, 10] * length(Q10_microbial) + 1]
resp_humus<- resp_humus[lhs_test[, 11] * length(resp_humus) + 1]
SRL<- SRL[lhs_test[, 12] * length(SRL) + 1]
#SRL<- paste0(as.character(SRL), ".0")
root_radius<- root_radius[lhs_test[, 13] * length(root_radius) + 1]
HLIM3H <- HLIM3H[lhs_test[, 14] * length(HLIM3H) + 1]
#HLIM3H<- paste0(as.character(HLIM3H), ".0")
HLIM3L <- HLIM3L[lhs_test[, 15] * length(HLIM3L) + 1]
#HLIM3L<- paste0(as.character(HLIM3L), ".0")
HLIM4 <- HLIM4[lhs_test[, 16] * length(HLIM4) + 1]
#HLIM4<- paste0(as.character(HLIM4), ".0")
ADCRH <- ADCRH[lhs_test[, 17] * length(ADCRH) + 1]
ADCRL<- ADCRL[lhs_test[, 18] * length(ADCRL) + 1]
ALPHACRIT <- ALPHACRIT[lhs_test[, 19] * length(ALPHACRIT) + 1]
RRI<-RRI[lhs_test[, 20] * length(RRI) + 1]
RDC<-RDC[lhs_test[, 21] * length(RDC) + 1]
#RDC <- paste0(as.character(RDC), ".0")
RDCTB<-RDCTB[lhs_test[, 22] * length(RDCTB) + 1]
FRTB<-FRTB[lhs_test[, 23] * length(FRTB) + 1]
FOTBa<-FOTBa[lhs_test[, 24] * length(FOTBa) + 1]
FOTBb<-FOTBb[lhs_test[, 25] * length(FOTBb) + 1]
FOTBc<-FOTBc[lhs_test[, 26] * length(FOTBc) + 1]
FLTBa<-FLTBa[lhs_test[, 27] * length(FLTBa) + 1]
FLTBb<-FLTBb[lhs_test[, 28] * length(FLTBb) + 1]
FLTBc<-FLTBc[lhs_test[, 29] * length(FLTBc) + 1]

lhs_paras=data.frame(CF,RSC,TDWI,RGRLAI,SPAN,SLATB,AMAXTB,Q10,TMPFTB,Q10_microbial,resp_humus,SRL,root_radius, HLIM3H,HLIM3L,HLIM4,ADCRH,ADCRL,
                     ALPHACRIT,RRI,RDC,RDCTB, FRTB, FOTBa, FOTBb, FOTBc, FLTBa,FLTBb,FLTBc) 

#_______________________________________________________________________________
### loop through my samples
#_______________________________________________________________________________

for (i in 1:nrow(lhs_paras)){
  
  # ----write crop text file --------------------------#
  
  source(paste0(homedir,"samples/crop_files/write_cropfile_for_SA.R"), local=TRUE)
  
  # ----execute swap file --------------------------# 
  
  setwd(paste0(homedir,"samples/",samplesize,"/" ))
  system('DeleteOldResults.cmd')
  system("runswap.bat")
  
  # ----save and post-process results --------------------------# 
  
  sim_y<-matrix(NA, nrow = nrow(lhs_paras), ncol=3)
  sim_y_pot<-matrix(NA, nrow = nrow(lhs_paras), ncol=3)
  
  ## Loop through each year (2017-2019)
  for (year in 2017:2019) {
    # Read the text file
    data<-read.table("result.crp", skip=7, header=TRUE, sep=",")
    # Convert the "date" column to Date format
    data$Date <- as.Date(data$Date)
    # Extract just the year from the "date" column
    data$years <- year(data$Date)
    # Filter data for the current year
    data_year <- subset(data, data$years == year)
    # Sort the data by date in ascending order
    data_year <- data_year[order(data_year$Date), ]
    # Get the last day's yield value for the current year
    last_actual_yield <- tail(data_year$CWSO, 1)
    last_potential_yield <- tail(data_year$CPWSO, 1)
    sim_y[i,year - 2016]<- last_actual_yield
    sim_y_pot[i,year - 2016]<- last_potential_yield
  }
  ## Print the resulting matrix
  
  # save yield for each year!
  mat_listyield17[[i]]<-sim_y[i,1]
  mat_listyield18[[i]]<-sim_y[i,2]
  mat_listyield19[[i]]<-sim_y[i,3]
  mat_listyield_mean[[i]]<-mean(sim_y[i,]) # save the matrices for each year and the mean of all years
  
  sim_irr<-matrix(NA, nrow = nrow(lhs_paras), ncol=1)
  irr_sim<-read.csv("result.irg" , skip=9, header=TRUE, sep=",")
  sim_irr[i]<-sum(irr_sim$Irrigation)*10
  mat_listirr_totaldepth[[i]]<- sim_irr[i]# depth in mm
  
}

toc()


####   first for actual yield
yield17<- do.call(rbind, lapply(mat_listyield17, as.vector))
yield17<-data.frame(yield17)
Y17<-as.numeric(yield17$yield)

yield18<- do.call(rbind, lapply(mat_listyield18, as.vector))
yield18<-data.frame(yield18)
Y18<-as.numeric(yield18$yield)

yield19<- do.call(rbind, lapply(mat_listyield19, as.vector))
yield19<-data.frame(yield19)
Y19<-as.numeric(yield19$yield)

yield_mean<- do.call(rbind, lapply(mat_listyield_mean, as.vector))
yield_mean<-data.frame(yield_mean)
Y_mean<-as.numeric(yield_mean$yield)

####  for irriation
irrdepth<- do.call(rbind, lapply(mat_listirr_totaldepth, as.vector))
irrdepth<-data.frame(irrdepth)
Y_irr<-as.numeric(irrdepth$irrdepth)

# save results
write.csv(Y17, paste0(homedir,"results/sample=",N,"/Y17.csv"))
write.csv(Y18, paste0(homedir,"results/sample=",N,"/Y18.csv"))
write.csv(Y19, paste0(homedir,"results/sample=",N,"/Y19.csv"))
write.csv(Y_mean, paste0(homedir,"results/sample=",N,"/Y_mean.csv"))
write.csv(Y_irr, paste0(homedir,"results/sample=",N,"/Y_irr.csv"))

#_______________________________________________________________________________
### apply and evaluate Sobol indices
#_______________________________________________________________________________


Y17<-read.csv(paste0(homedir,"results/sample=",N,"/Y17.csv"))
Y17<-Y17[,2]
Y18<-read.csv(paste0(homedir,"results/sample=",N,"/Y18.csv"))
Y18<-Y18[,2]
Y19<-read.csv(paste0(homedir,"results/sample=",N,"/Y19.csv"))
Y19<-Y19[,2]
Y_mean<-read.csv(paste0(homedir,"results/sample=",N,"/Y_mean.csv"))
Y_mean<-Y_mean[,2]
Y_irr<-read.csv(paste0(homedir,"results/sample=",N,"/Y_irr.csv"))
Y_irr<-Y_irr[,2]

############### actual yield
######### 2017
plot_uncertainty(Y = Y17, N = N) + labs(y = "Counts", x = "$y$")
S1<-plot_scatter(data = lhs_test, N = N, Y = Y17, params = params)
sobol_yield17<-sensobol::sobol_indices(Y=Y17, N=N, params=params, boot= TRUE, R=R, type= "norm", conf=0.95)
M1<-plot_multiscatter(data = lhs_test, N = N, Y = Y17, params = c(params[c(28,26,25,7,6,29)]))
cols<-colnames(sobol_yield17$results)[1:5]
#sobol_yield17$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
sobol_yield17
ind.dummy <- sobol_dummy(Y = Y17, N = N, params = params, boot = TRUE, R = R)
P1<-plot(sobol_yield17, dummy = ind.dummy)

######### 2018
plot_uncertainty(Y = Y18, N = N) + labs(y = "Counts", x = "$y$")
S2<-plot_scatter(data = lhs_test, N = N, Y = Y18, params = params)
sobol_yield18<-sensobol::sobol_indices(Y=Y18, N=N, params=params, boot= TRUE, R=R, type= "norm", conf=0.95)
M2<-plot_multiscatter(data = lhs_test, N = N, Y = Y18, params = c(params[c(28,26,25,7,6,29)]))
cols<-colnames(sobol_yield18$results)[1:5]
#sobol_yield18$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
sobol_yield18
ind.dummy <- sobol_dummy(Y = Y18, N = N, params = params, boot = TRUE, R = R)
P2<-plot(sobol_yield18, dummy = ind.dummy)

######### 2019
plot_uncertainty(Y = Y19, N = N) + labs(y = "Counts", x = "$y$")
S3<-plot_scatter(data = lhs_test, N = N, Y = Y19, params = params)
sobol_yield19<-sensobol::sobol_indices(Y=Y19, N=N, params=params, boot= TRUE, R=R, type= "norm", conf=0.95)
M3<-plot_multiscatter(data = lhs_test, N = N, Y = Y19, params = c(params[c(28,26,25,7,6,29)]))
cols<-colnames(sobol_yield19$results)[1:5]
#sobol_yield19$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
sobol_yield19
ind.dummy <- sobol_dummy(Y = Y19, N = N, params = params, boot = TRUE, R = R)
P3<-plot(sobol_yield19, dummy = ind.dummy)

######### mean yield
plot_uncertainty(Y = Y_mean, N = N) + labs(y = "Counts", x = "$y$")
S4<-plot_scatter(data = lhs_test, N = N, Y = Y_mean, params = params)
sobol_yield_mean<-sensobol::sobol_indices(Y=Y_mean, N=N, params=params, boot= TRUE, R=R, type= "norm", conf=0.95)
M4<-plot_multiscatter(data = lhs_test, N = N, Y = Y_mean, params = c(params[c(28,26,25,7,6,29)]))
cols<-colnames(sobol_yield_mean$results)[1:5]
#sobol_yield_mean$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
sobol_yield_mean
ind.dummy <- sobol_dummy(Y = Y_mean, N = N, params = params, boot = TRUE, R = R)
P4<-plot(sobol_yield_mean, dummy = ind.dummy)

####################### Irrigation

plot_uncertainty(Y = Y_irr, N = N) + labs(y = "Counts", x = "$y$")
S9<-plot_scatter(data = lhs_test, N = N, Y = Y_irr, params = params)
M9<-plot_multiscatter(data = lhs_test, N = N, Y = Y_irr, params = c(params[c(1,28, 6, 20,29, 26, 7)]))
sobol_irrdepth<-sensobol::sobol_indices(Y=Y_irr, N=N, params=params, boot= TRUE, R=R, type= "norm", conf=0.95)
cols<-colnames(sobol_irrdepth$results)[1:5]
#sobol_irrdepth$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
sobol_irrdepth
ind.dummy <- sobol_dummy(Y = Y_irr, N = N, params = params, boot = TRUE, R = R)
P9<-plot(sobol_irrdepth, dummy = ind.dummy)


#### save plots

pl_width = 3100
pl_heigh = 1300
pl_resol =200

plotfile = paste0(homedir,"results/sample=",N,"/Y17scatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(S1)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Y17multiscatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(M1)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Y17sobol.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(P1)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Y18scatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(S2)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Y18multiscatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(M2)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Y18sobol.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(P2)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Y19scatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(S3)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Y19multiscatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(M3)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Y19sobol.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(P3)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Ymeanscatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(S4)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Y1meanmultiscatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(M4)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Ymeansobol.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(P4)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Yirrscatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(S9)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Yirrmultiscatter.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(M9)
dev.off()

plotfile = paste0(homedir,"results/sample=",N,"/Yirrsobol.png")
png(plotfile, width = pl_width, height = pl_heigh, res = pl_resol)
print(P9)
dev.off()


#####---------------------------  End of script ----------------------------#####

