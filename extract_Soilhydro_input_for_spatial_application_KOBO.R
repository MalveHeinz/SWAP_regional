####   script to extract soil hydraulic data data for each cell/site i want to run the model for

# the goal is to create the csv file:soil_hydro_tex_your_area.csv" that is needed for input in my regional_application_SWAP.R script

# in this script we assume, that we already have extracted  soil textural information for each coordinate/cell!

# Ideally we would have as many properties as possible, but minimum requirements are: sand, silt clay and humus content.
# here the object soil_df contains exactly these parameters plus the soil depth for three depth, 0-30cm, 30-60cm, 60-120cm.
# These would be the columns in my case:
# "clay.0.30","clay.30.60","clay.60.120","sand.0.30","sand.30.60","sand.60.120","silt.0.30","silt.30.60","silt.60.120","soc.0.30","soc.30.60","soc.60.120","soil_depth","x","y"

# soil depth in my case is either 3= 90cm or 4 =120 cm, please adapt.

# We will use the soil textural information as inpüut for a pedotransferfunction to derive the mualem van Genuchten Parameters for each soillayer, as requested by SWAP

# required packages

library(raster)
library(ncdf4)
library(rasterVis)
library(maptools)
library(maps)
library(data.table)
library(lubridate) 
library(terra)
library(ggplot2)
library(tidyr)
library(rgdal)
library(sf)
library(dplyr)
library(exactextractr)
library(bfsMaps)
library(abind)
library(tictoc)
library(euptf2)
library(gdata)
library(zoo)

# define your homedir
path_r<-"C:/your_path"

# here is the raster with the gridcell i want to run my model for. in my case 30m x 30 m
raster_area<-raster(paste0(path_r, "your_path/raster_area.tif")) 
area<-raster_area
areaname<- "raster_area" 

# for each of the grid cell i before extracted the center coordinates. These coordinates are used in turn to extract climate and soil data
coords<-read.csv(paste0("C:path_to_your_coordinates/coords_",areaname,".csv"))
raster_areapoints<-coords
raster_area_ids<-data.frame(x=raster_areapoints[,1], y=raster_areapoints[,2])

###_____________________________ prepare data ________________________________

# read in pre-processed soil input file from your source ) as defined in the beginning of the script
soil_df<- read.csv(paste0("C:your_path/soil_input_",areaname,".csv"))
kobo<-soil_df # kobo a swiss soil dataset that I used (Boden"), KOBO (2023). Hinweiskarten für Bodeneigenschaften Landesweit modellierte Karten für Bodeneigenschaften für drei Tiefenstufen. )

# prepare the dataframe we will fill with the  hydraulic parameters later
hydro.tex.all<-matrix(ncol=49, nrow=nrow(kobo))

#_______________________________________________________________________________

tic()
for (x in 1:nrow(kobo)){  
  #---------------------- 04: prepare soil data input ----------------------------
  
  G<-kobo[x,13] # This is the soil depth, if unknown, put a default value
  soilmatrix<-as.data.frame(matrix(data=NA, nrow=5, ncol=5))
  colnames(soilmatrix)<-c("ISUBLAY","ISOILLAY","HSUBLAY","HCOMP","NCOMP")
  soilmatrix$ISUBLAY<-c(1:5)
  soilmatrix$ISOILLAY<-c(1,1,1,2,3)
  
  if (G==4){ # depending on 3 or 4, my entire soil profile is 90 or 120cm deep!
  soilmatrix$HSUBLAY<-format(c(5,10,15,30,60), nsmall=1)
  soilmatrix$NCOMP<-c(5,5,3,3,2)
  soilmatrix$HCOMP<-format(c(1,2,5,10,30), nsmall=1)
  } else {
  soilmatrix$HSUBLAY<-format(c(5,10,15,30,30), nsmall=1)
  soilmatrix$NCOMP<-c(5,5,3,3,3)
  soilmatrix$HCOMP<-format(c(1,2,5,10,10), nsmall=1)
  }
  
  if (G==4){depth="120.0"} else {depth="90.0"}
  
  #_______________________________________________________________________________
  #### soilhydro data
  
  # in my case i already have the data on soil organic carbon which i can extract here for my different soil layers
  # please note, that in my soilmatrix i defined 5 ISUBLAYs, but the first 3 ISUBLAYs are between 0 and 30cm.
  # Since I only have one value for 0-30 cm, i fill the first 3 ISUBLAYs with the same value.
  OC<-c(as.numeric(kobo[x,10]), as.numeric(kobo[x,10]), as.numeric(kobo[x,10]), as.numeric(kobo[x,11]) , as.numeric(kobo[x,12])) 
  
  # Since I dont have bulk density, I also derive it with another ptf
  BD = (1.66 - 0.318 *OC ^0.5)*1000 # De Vos, B., M. Van Meirvenne, P. Quataert, J. Deckers and B. Muys (2005). "Predictive Quality of Pedotransfer Functions for Estimating Bulk Density of Forest Soils." Soil Science Society of America Journal 69(2): 500-510.
  
  # now i prepare the inputdata in the required input format for the ptf I use.
  # The choice of ptf depends on your data availability!
  mydata = data.frame(SAMPLE_ID = c(1,2,3), DEPTH_M = c(30,60,120), USSAND = c(kobo[x,4], kobo[x,5], kobo[x,6]), 
                      USSILT = c(kobo[x,7], kobo[x,8], kobo[x,9]), USCLAY = c(kobo[x,1], kobo[x,2], kobo[x,3]), OC = OC[3:5])
  
  pred_MVG<- euptfFun(ptf = "PTF02", predictor = mydata, target = "MVG", query = "predictions")
  pred_MVG<-round(pred_MVG, digits=3)
  
  
  # Here I fill my before created dataframe with the newly derived mualem van genuchten parameters (MVG)
  hydromatrix<-as.data.frame(matrix(data=NA, ncol=11, nrow=3))
  colnames(hydromatrix)<-c("ISOILLAY1","ORES","OSAT","ALFA","NPAR","KSATFIT","LEXP","ALFAW","H_ENPR","KSATEXM","BDENS")
  hydromatrix[1,]<-c(1,pred_MVG[1,3],pred_MVG[1,2],pred_MVG[1,4],pred_MVG[1,5],format(pred_MVG[1,6]+0.01,nsmall=2, digits=2),pred_MVG[1,7],
                     pred_MVG[1,4],"0.0", format(pred_MVG[1,6]+0.01, nsmall=2, digits=2), BD[3]+0.1)
  hydromatrix[2,]<-c(2,pred_MVG[2,3],pred_MVG[2,2],pred_MVG[2,4],pred_MVG[2,5], format(pred_MVG[2,6]+0.01, nsmall=2, digits=2),pred_MVG[2,7],
                     pred_MVG[2,4],"0.0",format(pred_MVG[2,6]+0.01, nsmall=2, diits=2), BD[4]+0.1)
  hydromatrix[3,]<-c(3,pred_MVG[3,3],pred_MVG[3,2],pred_MVG[3,4],pred_MVG[3,5], format(pred_MVG[3,6]+0.01, nsmall=2, digits=2),pred_MVG[3,7],
                     pred_MVG[3,4],"0.0", format(pred_MVG[3,6]+0.01, nsmall=2, digits=2), BD[5]+0.1)
  
  # I also need to write the Humus content, so I assume that soil organic carbon makes up 58% of the total humus content (common value)
  OCtoHumus<-function(z) (z/58)*100 
  Humus<-lapply(OC, FUN=OCtoHumus)
  Humus<-unlist(Humus)
  Humus<-round(Humus, 2)
  
  # Now I can create the soil texture matrix, that is needed later in the SWAP main input file
  mydata = data.frame(SAMPLE_ID = c(1,2,3), DEPTH_M = c(30,60,120), USSAND = c(kobo[x,6], kobo[x,7], kobo[x,8]), 
                      USSILT = c(kobo[x,9], kobo[x,10], kobo[x,11]), USCLAY = c(kobo[x,3], kobo[x,4], kobo[x,5]), OC = OC[3:5])
  
  soiltexmatrix<-as.data.frame(matrix(data=NA, ncol=5, nrow=3))
  soiltexmatrix[1,]<-c("1",unlist(kobo[x,4])/100, unlist(kobo[x,7])/100, unlist(kobo[x,1])/100, Humus[3]/100)
  soiltexmatrix[2,]<-c("2",unlist(kobo[x,5])/100, unlist(kobo[x,8])/100, unlist(kobo[x,2])/100, Humus[4]/100)
  soiltexmatrix[3,]<-c("3",unlist(kobo[x,6])/100, unlist(kobo[x,9])/100, unlist(kobo[x,3])/100, Humus[5]/100)
  colnames(soiltexmatrix)<-c("ISOILLAY5","PSAND","PSILT","PCLAY","ORGMAT")
  
  hyd<- as.matrix(hydromatrix)
  stex<- as.matrix(soiltexmatrix)
  
  #  ________________________________________
  hydro.tex<-cbind(hydromatrix[1,],soiltexmatrix[1,],
              hydromatrix[2,],soiltexmatrix[2,],
              hydromatrix[3,],soiltexmatrix[3,]) # here i had a mistake before and was printing soiltexmatrix[2,] twice
  
  colnames(hydro.tex)<-c("ISOILLAY1.1","ORES.1","OSAT.1","ALFA.1","NPAR.1","KSATFIT.1","LEXP.1","ALFAW.1","H_ENPR.1","KSATEXM.1","BDENS.1",
                    "ISOILLAY5.1","PSAND.1","PSILT.1","PCLAY.1","ORGMAT.1",
                    "ISOILLAY1.2","ORES.2","OSAT.2","ALFA.2","NPAR.2","KSATFIT.2","LEXP.2","ALFAW.2","H_ENPR.2","KSATEXM.2","BDENS.2",
                    "ISOILLAY5.2","PSAND.2","PSILT.2","PCLAY.2","ORGMAT.2",
                    "ISOILLAY1.3","ORES.3","OSAT.3","ALFA.3","NPAR.3","KSATFIT.3","LEXP.3","ALFAW.3","H_ENPR.3","KSATEXM.3","BDENS.3",
                    "ISOILLAY5.3","PSAND.3","PSILT.3","PCLAY.3","ORGMAT.3")
  
  # ______________________________________
  
  hydro.tex<-as.numeric(hydro.tex)
  hydro.tex.all[x,]<-c(hydro.tex, depth)
  
}

toc() # This might take a couple of hours depending o the nbumber of cells
  
head(hydro.tex.all)
safe<-hydro.tex.all # safety copy

hydro.tex.all.final<-cbind(hydro.tex.all, coords[,2], coords[,3]) # add again the coordinates for each row/cell

colnames(hydro.tex.all.final)<-c("ISOILLAY1.1","ORES.1","OSAT.1","ALFA.1","NPAR.1","KSATFIT.1","LEXP.1","ALFAW.1","H_ENPR.1","KSATEXM.1","BDENS.1",
  "ISOILLAY5.1","PSAND.1","PSILT.1","PCLAY.1","ORGMAT.1",
  "ISOILLAY1.2","ORES.2","OSAT.2","ALFA.2","NPAR.2","KSATFIT.2","LEXP.2","ALFAW.2","H_ENPR.2","KSATEXM.2","BDENS.2",
  "ISOILLAY5.2","PSAND.2","PSILT.2","PCLAY.2","ORGMAT.2",
  "ISOILLAY1.3","ORES.3","OSAT.3","ALFA.3","NPAR.3","KSATFIT.3","LEXP.3","ALFAW.3","H_ENPR.3","KSATEXM.3","BDENS.3",
  "ISOILLAY5.3","PSAND.3","PSILT.3","PCLAY.3","ORGMAT.3", "depth", "x", "y")

# write large dataframe as csv that can now be used as an input to my regional_application_SWAP.R script!
write.csv(hydro.tex.all.final, paste0("C:/your_path/soil_hydro_tex_",areaname,".csv"))

#_____________________________________ END ___________________________________________