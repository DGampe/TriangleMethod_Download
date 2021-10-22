lsttm <-
  function(TEMP,outpath = getwd(), outname = "lst.tif", Lmin=1.2378*10^6, Lmax=15.3030*10^6,
           QcalMax=255, QcalMin=1, unit="C", write=F, ...){

    lambda=11.478*10^-6
    c1= 1.19104356*10^-16###
    c2= 1.43876869*10^-2

    if(class(TEMP)[[1]]!="RasterLayer"){
      TEMP = raster(TEMP)
    }

    TEMP[which(TEMP[]==0)] = NA

    ### Implementing the following formula in R:
    ### L = ((Lmax-Lmin)/(QcalMAX-QcalMIN)) * (Qcal[i] - QcalMin)+Lmin
    ### Qcal[i]: DN for each Pixel in the calculated scene

    ### Step II: from radiance to LST in K

    k1= c1/lambda^5
    k1 = k1
    k2= c2/lambda

    ### Implementing the following formula in R:
    ### LST = k2/(ln(k1/L[i] + 1)
    ### L[i]: radiance for each Pixel in the calculated scene

    ### Step I: from DN to radiance
    Radiance = max(((Lmax-Lmin)/(QcalMax-QcalMin)) * (TEMP - QcalMin)+Lmin,1)

    ### Step II: from radiance to LST in K
    LST = max(k2/(log(k1/Radiance + 1)),1)

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

