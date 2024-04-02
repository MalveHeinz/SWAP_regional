########  extract climate data and create climate input files for each cell ##########

# totally depends on the cliamte input data you have available. 
# i used the gridded cliamte data by MeteoSwiss and calculated ET0 beforehand.
# Ideally we would use ET0 from penman-Monteith, but I did not have all the parameters, so I calculated ET0 with Priestly-Taylor beforehand and used a large grid as input.
# If you have all required data for penman monteith just adjust the script accordingly

# You will also need the template_climatedata.csv file from the repository

library(raster)
library(dplyr)
library(euptf2)
library(gdata)
library(tictoc)
library(ncdf4)
library(ggplot2)
library(tidyr)


year<- 2020     # the year I want to write the meteofiles for
ye<-   "020"
yr<-     20 
path.base<-"C:you_home_directory"
# the directory where you want to write the cliamte input files to
climDir<- paste0(path.base,"your_path/") # 

# ___________________________________________________

# again, read in center coordinates of the grid cells you want to apply the model to
your_areapoints<-read.csv() 
your_area_ids<-data.frame(x=your_areapoints[,1], y=your_areapoints[,2])

#_____________________________________________________

####  read in climate data, in my case I had netCDF data, as well as ET0 and solrad in raster format
# this nc file contains the variables for Tmin, Tmax and precipitation
climdata<-nc_open(paste0(path.base, "climate_data/clim",yr,"_cropped.nc"))
time <- ncvar_get(climdata, "time.doy")
E <- ncvar_get(climdata, "meters_east") # E = X =lon = zb. 2574500
N <- ncvar_get(climdata, "meters_north")# N = y= lat =zb. 1059500
Tmin <- ncvar_get(climdata, "Tmin")
Tmax <- ncvar_get(climdata, "Tmax")
Prec <- ncvar_get(climdata, "RhiresD")

ET0<-read.csv(paste0(your_path), sep=";")

solrad<-read.csv(paste0(your_path), sep=",")
solrad<-solrad*100  # In my case I had the wrong unit, it needs to be in kj/m2

# template for year/month/day info
clim_muster<-read.csv("C:/your_path/template_climatedata.csv")
suse<-subset(clim_muster, clim_muster$YYYY== year) 
#________________________________________________________

tic()
for (x in 1:nrow(crop_coords)){
  
  #####  Tmin/Tmax
  # find index of nc file that is closest to the cell i am looping through at that step
  x_target <- crop_coords[x,1]  # Extract target coordinates
  y_target <- crop_coords[x,2] 
  
  # Find the index of the closest point in the climate array, please manually check if they match!
  x_index <- which.min(abs(E - x_target))
  y_index <- which.min(abs(N - y_target))
  
  # Extract all 365 values at the closest point to the target coordinate
  Tmin_daily<-Tmin[x_index,y_index,] 
  Tmax_daily<-Tmax[x_index,y_index,]
  Prec_daily<-Prec[x_index,y_index,]
  
  ET0_value<-as.numeric(ET0[matching_indices[x],])
  # ET0_value<-as.numeric(ET0[6017,]) # for the station of payerne
  
  solrad_value<-as.numeric(solrad[matching_indices[x],])
  
  ### write climate input file for that cell 
  #___________________________________________________________________________________________

  outfile <-  paste0(climDir, x ,"_siteX_clim.",ye) # write the file to the climDir directory and name it with the chosen year and the ID of the cell!

  write("*************************************************************************************************************************   ", outfile, append=T)
  write(paste0("*PAY_clim.", ye                                                                                                      ), outfile, append=T)  # or whichever year
  write("* Contents: SWAP  4.0 - Daily meterorological data                                                                          ", outfile, append=T)
  write("* Comment area:                                        ", outfile, append=T)
  write("*                                                                               ", outfile, append=T)
  write("*************************************************************************************************************************    ", outfile, append=T)
  write(" Station      DD      MM    YYYY         RAD       Tmin      Tmax        HUM      WIND      RAIN     ETref            WET    ", outfile, append=T)
  write("*             nr      nr      nr       kJ/m2         ºC        ºC        kPa       m/s        mm        mm              d    ", outfile, append=T)
  write("*************************************************************************************************************************    ", outfile, append=T)
  out_df <- data.frame(
    format(paste0(" 'PAY'"), width = 10),
    format(suse$DD, digits = 1, width = 5),
    format(suse$MM, digits = 1, width = 7),
    format(suse$YYYY, digits = 1, width = 7),
    format(solrad_value, digits=1, nsmall = 1, width = 11),
    format(Tmin_daily, digits = 1, width = 9),
    format(Tmax_daily, digits = 1, width = 9),
    format(round(1, 4), nsmall = 1, width = 10), # SWAP does not accept -99.9, although it is written in the manual, that this is how missing data should be flagged..
    format(round(1, 4), nsmall = 1, width = 10),
    format(Prec_daily, digits = 2, nsmall=1, width = 11),
    format(round(ET0_value, digits = 1), nsmall = 1, width = 10), 
    format(round(-99.9, 4), nsmall = 4, width = 13)
  )
  
  write.fwf(x=out_df , file =outfile, colnames = FALSE, append = TRUE, eol="\n")
  
}

toc()

###### If I want to delete files form just one year e.g. =
# 
# folder_path <- climDir
# files <- list.files(folder_path, full.names = TRUE)
# # Filter files with the ending .018
# files_to_delete <- files[grep("\\.020$", files)]
# # Remove the selected files
# file.remove(files_to_delete)


#________________________________________________________END_____________________________________________________________

 


