ndvi<-
function(NIR,RED, write = F, outname= "ndvi.tif", outpath = getwd(), sensor = "tm",
	 Date = "1990-01-01", cal.reflectance = FALSE,
 	 Lmin.red = - 1.17*10^6,Lmax.red = 264*10^6, Lmin.nir = -1.51*10^6,
	 Lmax.nir = 221*10^6,  ESUN.band3 = 1536, ESUN.band4 = 1031, ...){

if(sensor == "etm+"){
	ESUN.band3 = 1533
	ESUN.band4 = 1039
	LminLow.red = - 5.0*10^6 ## can be found in calibration parameter file for Landsat scene. (im micrometers -> convert to m
	LmaxLow.red = 234.4*10^6 ##  "" values here for Lmin+Lmax for TM data
	LminHigh.red = -5.0*10^6 #can be found in calibration parameter file for Landsat scene. (im micrometers -> convert to m
	LmaxHigh.red = 152.9*10^6 ###  "" values here for Lmin+Lmax for TM data

	LminLow.nir = -5.1*10^6 ### can be found in calibration parameter file for Landsat scene. (im micrometers -> convert to m
	LmaxLow.nir = 241.1*10^6 ###  "" values here for Lmin+Lmax for TM data
	LminHigh.nir = -5.1*10^6 ### can be found in calibration parameter file for Landsat scene. (im micrometers -> convert to m
	LmaxHigh.nir = 157.4*10^6  ###  "" values here for Lmin+Lmax for TM data
}

if(class(NIR)[[1]]!="RasterLayer"){
NIR = raster(NIR)
}
if(class(RED)[[1]]!="RasterLayer"){
RED = raster(RED)
}
nirsplit = strsplit(names(NIR),"")
redsplit = strsplit(names(RED),"")

if(cal.reflectance == TRUE){
	zen = sun.zen(NIR,paste(Date," 10:00:00",sep=""))

	earth.sun.dist  <- read.csv(system.file("extdata", "EarthSunDist.txt", package="TriangleMethod"),sep=",")
	d = earth.sun.dist[as.numeric(strftime(Date,format = "%j")),2] 	# Earth - sun - distance [AU]

	### Step I: from DN to radiance
	gc()
	QcalMax=255
	QcalMin=1

	### Implementing the following formula in R:
	### L = ((Lmax-Lmin)/(QcalMAX-QcalMIN)) * (Qcal[i] - QcalMin)+Lmin
	### Qcal[i]: DN for each Pixel in the calculated scene

	### Step I: from DN to radiance

	if(sensor == "etm+"){
		RadLow_NIR=((LmaxLow.nir-LminLow.nir)/(QcalMax-QcalMin)) * (NIR  - QcalMin)+LminLow.nir
		RadHigh_NIR= ((LmaxHigh.nir-LminHigh.nir)/(QcalMax-QcalMin)) * (NIR  - QcalMin)+LminHigh.nir
		Radiance_NIR = (RadLow_NIR+RadHigh_NIR)/2 # Average of the Low Gain and High Gain 

		RadLow_RED= ((LmaxLow.red-LminLow.red)/(QcalMax-QcalMin)) * (RED  - QcalMin)+LminLow.red
		RadHigh_RED= ((LmaxHigh.red-LminHigh.red)/(QcalMax-QcalMin)) * (RED  - QcalMin)+LminHigh.red
		Radiance_RED = (RadLow_RED+RadHigh_RED)/2 # Average of the Low Gain and High Gain 
	}else{
		Radiance_NIR =((Lmax.red-Lmin.red)/(QcalMax-QcalMin)) * (NIR - QcalMin)+Lmin.nir
		Radiance_RED =((Lmax.nir-Lmin.nir)/(QcalMax-QcalMin)) * (RED - QcalMin)+Lmin.red
	}
	rm(RadLow_NIR, RadHigh_NIR, RadLow_RED, RadHigh_RED) 
	gc()
	
	### Convert to reflectance:

	RED = ((pi*Radiance_RED*round(round(d,5)^2,5))/(ESUN.band3*cos((zen/180*pi))))   # Equation according to CHANDER & MARKHAM 2003.
	# Band4:
	NIR = ((pi*Radiance_NIR*round(round(d,5)^2,5))/(ESUN.band4*cos((zen/180*pi))))
	rm(Radiance_RED, Radiance_NIR)
}

NDVI = (NIR -RED)/(NIR+RED)

setwd(outpath)

if(write == T){
	writeRaster(NDVI , filename=outname , "GTiff", overwrite=T, NAflag = -9999)  
}
return(NDVI)

invisible(y)
}
