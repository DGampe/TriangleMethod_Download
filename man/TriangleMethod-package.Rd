\name{TriangleMethodPackage-package}
\alias{TriangleMethodPackage-package}
\alias{TriangleMethodPackage}
\docType{package}
\title{
Triangle Method
}
\description{
The Triangle Method (Jiang & Islam, 1999) approach provides estimations for actual evapotranspiration using remotely sensed information as the Normalized Differnce Vegetation Index (NDVI) and Land Surface Temperature (LST). The approach is based on the cooling effect of highly evaporating areas. Maximum evapotranspiration hence occurs over areas with high NDVI values and low LST. With help of the NDVI-LST-plot Tmax and Tmin are estimated to scale the Priestley-Taylor-coefficient Phi (often and originally referred to as Alpha, see Priestley & Taylor (1972) for further information). Maximum Phi for equilibrium conditions was determined as 1.26 (Eichinger et al. (1996)). This package provides functions to derive those information from common LANDSAT TM and ETM+ scenes (might be modified by the user to work also for other sensors). The function 'triangle' then reads in these two grid-files (NDVI & LST) and allows to estimate actual evpotranspiration. Furthermore, the function grid.plot provides a way to plot gridded data (as the final ETact, or NDVI etc.), based on the spplot-function and the package sp.
}
\details{
\tabular{ll}{
Package: \tab TriangleMethod\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2021-10-14\cr
License: \tab GNU General Public License v3.0\cr
}
Functions are:
ndvi, lstetm, lsttm, ndvilst.plot, cloud.rm, sun.zen
}
\author{
David Gampe (david.gampe@googlemail.com)
Verena Huber Garcia
Benjamin Mueller
Philip Marzahn
}
\references{

Allen, R., Pereira, L., Raes, D., and Smith, M. (1998). Crop evapotranspiration - guidelines for
computing crop water requirements. FAO Irrigation and drainage paper, 56.

Batra, N., Islam,S., Venturini,V.,Bisht,G.,Liang,L. (2006). Estimation and comparison of evapotranspiration from MODIS and AVHRR sensors for clear sky days over the Southern Great Plains, Remote Sensing of the Environment, 103, 1-15.

Chander,G., Markham, B., Helder,D. (2009). Summary of current radiometric calibration coefficients for Landsat MSS, TM, ETM+,and EO-1 ALI sensors, 113, 893-903.

Crago, R. and Brutsaert, W. (1996). Daytime evaporation and the self-preservation of the evaporative
fraction and the Bowen ratio. Journal of Hydrology , 178:241-255.

Davies, J. (1967). A note on the relationship between net radiation and solar radiation. Quarterly
Journal if the Royal Meteorological Society, 93:109-115.

Jiang,L. & Islam,S. (1999). A methodology for estimation of surface evapotranspiration over large areas using remote sensing observations.

Eichinger, E., Parlange,M.,Stricker,H. (1996). On the concept of equilibrium evaporation and the value of the Pristley-Taylor coefficient, water resources research, 32, 161-164. 

Priestley, C. & Taylor, R. (1972). On the assessment of surface heat flux and evaporation using
large-scale parameters. Monthly Weather Review, 100(2):81-92.

Stisen, S., Sandholt,I.,Norgaard,A., Fensholt,R., Hogh Jensen, K. (2008). Combining the triangle method with thermal inertia to estimate regional evapotranspiration - Applied to MSG-SEVIRI data in the Senegal River basin, Remote Sensing of the Environment, 112, 1242-1255.

Tang, R., Zhao-Liang, L., Bohui, T. (2010). An application of the Ts-VI triangle method with enhanced edges determination for evapotranspiration estimation from MODIS data in arid and semi-arid regions: Implementation and validation, Remote Sensing of Environment, 114, 540 - 551.

Wang, K., Li,Z., Cribb, M. (2006). Estimation of evaporative fraction from a combination of day and night land surface temperatures and NDVI: A new method to determine the Priestley-Taylor parameter, Remote Sensing of the Environment, 102, 293-305. 

}

\keyword{ package triangle evapotranspiration landsat }

\examples{

# calculate NDVI from raster in memory (within R):
NIR<- raster(ncol=30,nrow=20)
NIR [] <- runif(30*20,min = 1,max = 255)
RED <- raster(ncol=30,nrow=20)
RED [] <- runif(30*20,min = 1,max = 255)

NDVI = ndvi(NIR,RED)
plot(NDVI)


# calculate NDVI from file based on reflectance not DN:
# function will search for the band containing the string 'B30' & 'B40'.
path <- system.file("extdata", package="TriangleMethod")
NDVI = ndvi.file(path,calc.ref = TRUE, date = "1988-11-19")

# calculate LST for TM sensor from file
# function will search for the band containing the string 'B60'.
path <- system.file("extdata", package="TriangleMethod")
lst = lsttm(path)


# estimation of actual evapotranspiration with NDVI and LST:
path <- system.file("extdata", package="TriangleMethod")
lst = lsttm(path)
ndvi = ndvi.file(path)
tact <- 25  	# fictive homogeneous air temperature of 25 deg C
rad <- 225 	# fictive homogeneous global radiation of 225 W/m^2
et <- triangle(ndvi, lst, tact, rad,ti.method="tang")
plot(et)




}
