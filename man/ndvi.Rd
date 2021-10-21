\name{ndvi}
\alias{ndvi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate NDVI from raster.}
\description{
This function will calculate the Normalized Difference Vegetation Index from raster data sets - either in memory, or as files (which will be rastered), defined as: 

NDVI = (NIR - RED) / (NIR + RED)

The NDVI will be calculated from pure DN as default, however, an option to caculate the NDVI based on reflectance is implemented according to the procedure by Chander et al. (2009). See argument 'calc.ref' for more information.
}
\usage{
ndvi(NIR,RED, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NIR}{
File containing the NIR information for calculation. 
} 
\item{RED}{
File containing the RED information for calculation. 
}
\item{write}{Specify, if the resulting NDVI grid should be written to file.Default = F (T/F)}
\item{outpath}{User defined outpath to write the final NDVI file.}
\item{outname}{User defined name for the resulting file.(Default name / format: 'ndvi.tif'}
\item{calc.ref}{Should the NDVI be calculated from reflectance? Default = FALSE, NDVI will be calculated from DN only. If set TRUE, make sure to provide the date of the scene. See date for further information.}
\item{sensor}{Sensor for NDVI calculation. Default "tm", other option: "etm+". Only important if calc.ref = TRUE.}
\item{Date}{Date of the scene (overpassing date). Format POSIXct, e.g. "1990-01-10". Default is "1990-01-01". Only important if calc.ref = TRUE.}
\item{Lmin.red,...}{Calibration parameter files for the sensor. Landsat TM and ETM+ are implemented, if images originate from a different sensor, these parameters should be adjusted. Only relevant if calc.ref =T.}
}
\details{
Note: calculation of the NDVI for a larger scene (e.g. one entire Landsat image) will require a lot of time and memory! If only a smaller part of the image is needed, consider to crop the image first.
Conversion to reflectance according to Chander et al. (2009).
}
\value{
RasterLayer 
.tiff file
}

\seealso{
 
}

\examples{
NIR<- raster(ncol=30,nrow=20)
NIR [] <- runif(30*20,min = 1,max = 255)
RED <- raster(ncol=30,nrow=20)
RED [] <- runif(30*20,min = 1,max = 255)

NDVI = ndvi(NIR,RED)
plot(NDVI)

red <- system.file("extdata", "L4175038_03819880324_B30.TIF", package="TriangleMethod")
red <- raster(red)

nir <- system.file("extdata", "L4175038_03819880324_B40.TIF", package="TriangleMethod")
nir <- raster(nir)

NDVI = ndvi(nir,red)
plot(NDVI)


red <- system.file("extdata", "B30.TIF", package="TriangleMethod")
red = raster(red)
nir <- system.file("extdata", "B40.TIF", package="TriangleMethod")
nir <- raster(nir)
ndvi = ndvi(red,nir,cal.reflectance = TRUE, date = "1980-11-14")
plot(ndvi)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
