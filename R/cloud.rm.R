cloud.rm<-
function(BLUE, TEMP,nClus=30,nsample=100,layers=c(1,6),maskmin=5,maskmax=NA,
	 outname="clouds-masked.tif", outpath = getwd(), write=F, ...){
    
  if(class(BLUE)[[1]]!="RasterLayer"){
    BLUE = raster(BLUE)
  }
  if(class(TEMP)[[1]]!="RasterLayer"){
    TEMP = raster(TEMP)
  }
  lstack = stack(BLUE,TEMP)	
  lDF=as.data.frame(lstack[[layers]])
  lDF=na.omit(lDF)
  l.cl=clara(lDF, stand=T, metric="euclidean"
             ,k=nClus,trace=T,medoids.x=F,keep.data=F, samples=nsample)
  lDF$clus=l.cl$clustering
  
  cloud.mask=lstack[[1]]
  cloud.mask[][as.numeric(rownames(lDF))]=lDF$clus
  if(is.null(maskmax)) {
    cloud.mask[][which(cloud.mask[]<=maskmin)] = NA    			
  } else{
    cloud.mask[][which(cloud.mask[]<=maskmin | cloud.mask[]>=maskmax)] = NA    		
  }
  cloud.mask[][which(!is.na(cloud.mask[]))] = 1
  
  setwd(outpath)
  if(write==T){
   writeRaster(cloud.mask , filename = outname , "GTiff", overwrite=T, NAflag = -9999)  
  }
  return(cloud.mask)
}	


