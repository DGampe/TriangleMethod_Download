lstetm <-
  function(TEMPL, TEMPH, outpath = getwd(), outname = "ts.tif",
           LminLow=0.0*10^6, LmaxLow=17.04*10^6, LminHigh=3.2*10^6, LmaxHigh=12.65*10^6,
           QcalMax=255, QcalMin=1, unit="C", write=F, ...){

    lambda=11.478*10^-6#
    c1= 1.19104356*10^-16###
    c2= 1.43876869*10^-2

    if(class(TEMPL)[[1]]!="RasterLayer"){
      TEMPL = raster(TEMPL)
    }
    if(class(TEMPH)[[1]]!="RasterLayer"){
      TEMPH = raster(TEMPH)
    }

    TEMPL[which(TEMPL[]==0)] = NA
    TEMPH[which(TEMPH[]==0)] = NA
    ### Step I: from DN to radiance

    ### Implementing the following formula in R:
    ### L = ((Lmax-Lmin)/(QcalMAX-QcalMIN)) * (Qcal[i] - QcalMin)+Lmin
    ### Qcal[i]: DN for each Pixel in the calculated scene

    ### Step II: from radiance to brightness T in K

    k1= c1/lambda^5
    k2= c2/lambda

    ### Implementing the following formula in R:
    ### LST = k2/(ln(k1/L[i] + 1)
    ### L[i]: radiance for each Pixel in the calculated scene

    ### Step I: from DN to radiance

    RadLow= ((LmaxLow-LminLow)/(QcalMax-QcalMin)) * (TEMPL  - QcalMin)+LminLow
    RadHigh= ((LmaxHigh-LminHigh)/(QcalMax-QcalMin)) * (TEMPH  - QcalMin)+LminHigh
    RadAvg = (RadLow+RadHigh)/2 # Average of the Low Gain and High Gain thermal bands

    ### Step II: from radiance to LST in K
    LST = max(k2/(log(k1/RadAvg + 1)),1)

    if(unit=="C"){
      ### Step III: converting LST K in °C
      LST =LST  - 273.15
    }

    setwd(outpath)
    if(write ==  T){
      writeRaster(LST , filename=outname , "GTiff", overwrite=T, NAflag = -9999)
    }

    return(LST)
    invisible(y)
  }
