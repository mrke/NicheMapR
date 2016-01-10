#' crudat2raster
#'
#' Converts CRU integer ascii file to raster with longitude range -180 to 180
#' @param folder Path to the folder you want to install the global climate data in
#' @return cru.raster CRU data in raster format
#' @return header Header data for the layer
#' @usage get.global.climate(folder)
#' @export
crudat2raster <- function(file=file,loc=loc){
  # read in the ascii data
  Lines <- readLines(paste(loc,file,sep="")) #1. break into lines
  k <- count.fields(paste(loc,file,sep="")) #2. get field sizes
  Lines2 <- gsub("^ *| *$", "", Lines[1:2]) # 3. trim whitespace from beginning and end for first two lines
  Lines3 <- gsub(" +", ",", Lines2) # replace each string of spaces with a ,
  header=read.table(text = Lines3, sep = ",", as.is = TRUE, nrows = 1,header=TRUE) # now read in the header and data
  data<-Lines[3:length(Lines)] # get the rest of the table (actual data)
  write.table(data,paste(loc,"data.txt",sep=""),col.names = FALSE,row.names = FALSE,sep="",quote=FALSE) # write it out, now that header is gone
  library(readr)
  data2<-read_fwf(paste(loc,"data.txt",sep=""),fwf_widths(rep(5,720)),progress = interactive())
  file.remove(paste(loc,"data.txt",sep=""))

  # empty grid -need to make it from -180 to 180
  library(raster)
  gridout <- raster(ncol=header$n_cols, nrow=header$n_rows, xmn=header$xmin-180, xmx=header$xmax-180, ymn=header$ymin, ymx=header$ymax)
  cru.raster=stack(replicate(header$n_months,gridout))

  longs=seq(header$xmin,header$xmax,header$grd_sz) # sequence of longitudes (0-360)
  lats=rev(seq(header$ymin,header$ymax,header$grd_sz)) # sequence of lats (need to reverse them)
  longlats=expand.grid(longs,lats) # get all combinations
  for(i in 1:header$n_months){ # for each month
    grid=as.matrix(data2[(360*(i-1)+1):(360*(i-1)+360),]) # get the month's data chunk
    grid=as.vector(c(t(grid))) # transpose the grid, then vectorise it
    longlats[longlats[,1]>180,1]<-(longlats[longlats[,1]>180,1]-360) # turn longitudes > 180 into negative values going from -180 down to 0
    x<-cbind(longlats[1],longlats[2],grid) # list of co-ordinates plus values
    x<-x[x[,3]!=-9999,] # get rid of no data
    cru.raster[[i]] <- rasterize(cbind(x[,1],x[,2]), gridout, x[,3]) # rasterize!
  }
  rm(data2)
  return(list(cru.raster=cru.raster,header=header))
}
