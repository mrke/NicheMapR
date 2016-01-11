#' build.global.climate
#'
#' Downloads and builds the global climate, elevation and soil moisture datasets used in the function micro_global. The global
#' climate dataset global.nc is generated from "New, M., Lister, D., Hulme, M. and Makin, I., 2002: A high-resolution
#' data set of surface climate over global land areas. Climate Research 21:1-25", specically 10' estimates of monthly
#' mean (over 1961-1990) preceipitation, wet-days, air temperature, diurnal temperature range, relative humidity and
#' wind speed. Cloud cover comes from a bilinear interpolation of a lower resolution version of this dataset (New, M.,
#' M. Hulme, and P. D. Jones. 1999. Representing twentieth century space-time climate variability. Part 1: development
#' of a 1961-90 mean monthly terrestrial climatology. Journal of Climate 12:829-856.). The micro_global function
#' optionaly uses a  global monthly soil moisture estimate from NOAA CPC Soil Moisture http://140.172.38.100/psd/
#' thredds/catalog/Datasets/cpcsoil/catalog.html
#' @param folder Path to the folder you want to install the global climate data in
#' @usage get.global.climate(folder)
#' @export
get.global.climate <- function(folder="c:/globalclimate"){
ANSWER<-readline(prompt = "This function downloads and extracts 0.93 GB of data, type 'y' if you want to continue: ")
if(substr(ANSWER, 1, 1) == "y"){
  if(dir.exists(folder)==FALSE){
    dir.create(folder)
  }

# soil moisture
soilmoist.file<-"http://140.172.38.100/psd/thredds/catalog/Datasets/cpcsoil/catalog.html?dataset=Datasets/cpcsoil/soilw.mon.ltm.v2.nc"
destin<-paste(folder,"/soilw.mon.ltm.v2.nc",sep="")
download.file(soilmoist.file, destin, mode="wb")

# global_climate.nc construction

# global climate files to download:
# "https://crudata.uea.ac.uk/cru/data/hrg/tmc/grid_10min_elv.dat.gz" # Elevation
# "https://crudata.uea.ac.uk/cru/data/hrg/tmc/grid_10min_pre.dat.gz" # Precipitation
# "https://crudata.uea.ac.uk/cru/data/hrg/tmc/grid_10min_rd0.dat.gz" # wet-days
# "https://crudata.uea.ac.uk/cru/data/hrg/tmc/grid_10min_tmp.dat.gz" # Mean temperature
# "https://crudata.uea.ac.uk/cru/data/hrg/tmc/grid_10min_dtr.dat.gz" # Mean diurnal temperature range
# "https://crudata.uea.ac.uk/cru/data/hrg/tmc/grid_10min_reh.dat.gz" # Mean relative humidity
# "https://crudata.uea.ac.uk/cru/data/hrg/tmc/grid_10min_wnd.dat.gz" # Mean wind speed
# "http://ipcc-ddc.cru.uea.ac.uk/download_data/observed/ccld6190.zip" # Mean wind speed

sourcepath="https://crudata.uea.ac.uk/cru/data/hrg/tmc/"
destfile1="grid_10min_"
destfile2=".dat.gz"
vars=c('elv','pre','rd0','wnd','tmp','dtr','reh')

# download and unzip all the data, construct rasters
library(R.utils)
library(raster)
library(ncdf4)
gridout <- raster(ncol=2160, nrow=1080, xmn=-180, xmx=180, ymn=-90, ymx=90)
global_climate=stack(replicate(97,gridout))
cat('downloading and decompressing the data')
for(i in 1:length(vars)){
  destin<-paste(folder,destfile1,vars[i],destfile2,sep="")
  download.file(paste(sourcepath,destfile1,vars[i],destfile2,sep=""), destin, mode="wb")
  cat(paste(vars[i],destfile2,' downloaded \n',sep=""))
  R.utils::gunzip(destin)
  cat(paste(vars[i],destfile2,' decompressed \n',sep=""))
}
cat('reading in each variable, rasterising and storing in global_climate stack \n')
for(i in 1:length(vars)){
  data<-read.table(paste(folder,destfile1,vars[i],".dat",sep=""))
  x<-cbind(data[,2],data[,1]) # list of co-ordinates
  if(i==1){ # elevation, just one colum
    global_climate[[1]] <- rasterize(x, gridout, data[,3])
    cat(paste(vars[i],' done \n',sep=""))
  }else{
    for(j in 1:12){
      global_climate[[1+(i-2)*12+j]] <- raster::rasterize(x, gridout, data[,2+j])
      cat(paste(vars[i],' month ',j,' done \n',sep=""))
    }
  }
}

meantemps=global_climate[[38:49]] #save mean temps
meandtrs=global_climate[[50:61]] #save mean diurnal temp ranges
meanhums=global_climate[[62:73]] #save mean hums

# construct max/min temps
for(i in 1:12){
  global_climate[[37+i]] = meantemps[[i]] - meandtrs[[i]]/2
  global_climate[[49+i]] = meantemps[[i]] + meandtrs[[i]]/2
}

# construct max/min relative humidities at new tmin/tmax values from vapour pressure at mean temperature/humidity
library(NicheMapR) # need functions WETAIR.rh and VAPPRS
for(i in 1:12){
  e=WETAIR.rh(rh=meanhums[[i]],db=meantemps[[i]])$e
  global_climate[[61+i]] = (e/VAPPRS(global_climate[[49+i]]))*100 # minhum = (e_meantemp/esat_mintemp)*100
  global_climate[[73+i]] = (e/VAPPRS(global_climate[[37+i]]))*100 # maxhum = (e_meantemp/esat_maxtemp)*100
}
# adjust values outside 0-100 range
for(i in 1:12){
  values(global_climate[[61+i]]) <- ifelse(values(global_climate[[61+i]]) < 0, 0, values(global_climate[[61+i]]))
  values(global_climate[[61+i]]) <- ifelse(values(global_climate[[61+i]]) > 100, 100, values(global_climate[[61+i]]))
  values(global_climate[[73+i]]) <- ifelse(values(global_climate[[73+i]]) < 0, 0, values(global_climate[[73+i]]))
  values(global_climate[[73+i]]) <- ifelse(values(global_climate[[73+i]]) > 100, 100, values(global_climate[[73+i]]))
}

# cloud cover

# download 0.5 degree data
cld<-"http://ipcc-ddc.cru.uea.ac.uk/download_data/observed/ccld6190.zip" # Mean wind speed
destin<-paste(folder,"ccld6190.zip",sep="")
download.file(cld, destin, mode="wb")
cat(paste("ccld6190.zip",' downloaded \n',sep=""))
untar(zipfile=destin,exdir=folder)
file.remove(destin)

# read in the ascii data
Lines <- readLines(paste(folder,"ccld6190.dat",sep="")) #1. break into lines
k <- count.fields(paste(folder,"ccld6190.dat",sep="")) #2. get field sizes
Lines2 <- gsub("^ *| *$", "", Lines[1:2]) # 3. trim whitespace from beginning and end for first two lines
Lines3 <- gsub(" +", ",", Lines2) # replace each string of spaces with a ,
header=read.table(text = Lines3, sep = ",", as.is = TRUE, nrows = 1,header=TRUE) # now read in the header and data
data<-Lines[3:length(Lines)] # get the rest of the table (actual data)
write.table(data,paste(folder,"data.txt",sep=""),col.names = FALSE,row.names = FALSE,sep="",quote=FALSE) # write it out, now that header is gone
data2<-read.fwf(paste(folder,"data.txt",sep=""),rep(5,720)) # read it in, assuming 5 spaces per value and 720 values per line

# empty grid -need to make it from -180 to 180
gridout <- raster(ncol=header$n_cols, nrow=header$n_rows, xmn=header$xmin-180, xmx=header$xmax-180, ymn=header$ymin, ymx=header$ymax)
monthly=stack(replicate(header$n_months,gridout))

longs=seq(header$xmin,header$xmax,header$grd_sz) # sequence of longitudes (0-360)
lats=rev(seq(header$ymin,header$ymax,header$grd_sz)) # sequence of lats (need to reverse them)
longlats=expand.grid(longs,lats) # get all combinations
for(i in 1:header$n_months){ # for each month
  grid=as.matrix(data2[(360*(i-1)+1):(360*(i-1)+360),]) # get the month's data chunk
  grid=as.vector(c(t(grid))) # transpose the grid, then vectorise it
  longlats[longlats[,1]>180,1]<-(longlats[longlats[,1]>180,1]-360) # turn longitudes > 180 into negative values going from -180 down to 0
  x<-cbind(longlats[1],longlats[2],grid) # list of co-ordinates plus values
  x<-x[x[,3]!=-9999,] # get rid of no data
  monthly[[i]] <- rasterize(cbind(x[,1],x[,2]), gridout, x[,3]) # rasterize!
}

# now resample down to 10 min and put in global_climate stack
one=(global_climate[[1]]+10000)/(global_climate[[1]]+10000)
gridout <- raster(ncol=2160, nrow=1080, xmn=-180, xmx=180, ymn=-90, ymx=90)
for(i in 1:header$n_months){
  global_climate[[85+i]]=resample(monthly[[i]], gridout)*one
}

writeRaster(global_climate,paste(folder,"global_climate.nc",sep=""))
file.remove(paste(folder,"data.txt",sep=""))

# cat("extracting climate data", '\n')
# CLIMATE <- raster::extract(global_climate,cbind(141,-34))
# ALTT<-as.numeric(CLIMATE[,1]*1000) # convert from km to m
# RAINFALL <- CLIMATE[,2:13]
# RAINYDAYS <- CLIMATE[,14:25]
# WNMAXX <- CLIMATE[,26:37]
# WNMINN<-WNMAXX*0.1 # impose diurnal cycle
# TMINN <- CLIMATE[,38:49]
# TMAXX <- CLIMATE[,50:61]
# ALLMINTEMPS<-TMINN
# ALLMAXTEMPS<-TMAXX
# ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)
# RHMINN <- CLIMATE[,62:73]
# RHMAXX <- CLIMATE[,74:85]
# CCMINN <- CLIMATE[,86:97]
# CCMAXX<-CCMINN
#
# global_climate_old<-raster::brick(paste("c:/global climate/","/global_climate.nc",sep=""))
# for(i in 1:12){
#   plot((100-global_climate[[48+i]])-monthly2[[i]]*one,zlim=c(-75,75),main=i)
# }
# plot(100-global_climate_old[[49]],zlim=c(0,100),main="Jan original")
# plot(global_climate[[86]]*one,zlim=c(0,100),main="Jan 0.5 res, resampled")
# plot(monthly[[1]],zlim=c(0,100),main="Jan 0.5 res")

} # end check for user confirmation
}

