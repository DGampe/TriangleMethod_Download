\name{triangle}
\alias{triangle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimation of actual evapotranspiration with the Triangle Method.
}
\description{
  Function to estimate actual evapotranspiration according to Jiang and Islam (1999) using the Normalized Difference Vegetation Index (NDVI) and Land Surface Temperature (LST). The approach is based on the cooling effect of highly evaporating areas. Maximum evapotranspiration hence occurs over areas with high NDVI values and low LST. With help of the NDVI-LST-plot Tmax and Tmin are estimated to scale the Priestley-Taylor-coefficient Phi.
  A log-file will be written to the defined (or default) outpath containing all input parameters and calculation methods, as well as mean, max and min evapotranspiration from the scene.
  For details it is refered to the paper related to this package by Gampe et al. (2018).
}
\usage{
  triangle(NDVI, LST, Tact, Rad, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NDVI}{
    Can either be the filename of the raster containing the NDVI-information, or the allready raster variable in the workspace.
  }
  \item{LST}{
    Can either be the filename of the raster containing the LST-information, or the allready raster variable in the workspace.

  }
  \item{Tact}{
    Observed temperature at overpassing time. Can either be numeric, raster or filename of a grid (raster format, will then be rastered)[either deg C or K].
  }
  \item{Rad}{
    Observed global radiation at overpassing time. Can either be numeric, raster or filename of a grid (in raster format, will then be rastered) [in W/m^2].
  }
  \item{t.min}{Tmin from NDVI-LST-plot. If no user input is provided, it will be estimated from the automatic plot. }
  \item{t.max}{Tmax from NDVI-LST-plot. If no user input is provided, it will be estimated from the automatic plot. }
  \item{outpath.ras}{User defined outpath to write the final LST file.}
  \item{outname.ras}{User defined name for the resulting file.}
  \item{outname}{Output filename for automatically calculated NDVI-LST-plot.}
  \item{write.plot}{Should the ndvi-lst plot be written to file? Default: F)
\item{dem}{Input elevation, can either be numeric, raster or filename of a grid (in raster format) containing elevation. Default is 100m.}
\item{int.method}{Method for resampling of Tact, Rad and DEM to the input grid resolution. Default:  "ngb", see ?resample for further options.}
\item{ti.method}{Scaling of Ti in NDVI-LST plot. Default: "linear", further options are "tang"}
\item{delta.method}{Method to calculate the slope of the saturated vapor pressure as defined by Murray (1967). Default: "batra", further options are "fao" and "wang".}
\item{pressure.method}{Method to calculate pressure change with altitude. Default: "barometric" (approach according to Roedel (1994), further options are "wang".}
\item{gamma.method}{Method to estimate the psychrometric constant gamma. Default: "wang", further options are "fao"}
\item{ndvi.classify}{Should the NDVI be classified (in steps of 0.05). Default = T}
\item{net.radiation}{Net radiation (Rn) [W/m^2]. If existent can be either numeric, or file name. Otherwise Rn will be estimated as a share of the global radiation.}
\item{soil.heat}{Soil heat flux (G) in [W/m^2]. If existent can be either numeric, or file name. Otherwise G will be estimated as a share of the global radiation.}
\item{rn.perc}{If no input grid for Rn is provided, Rn will be calculated as a share of the global radiation. Default: 0.55}
\item{g.perc}{If no input grid for G is provided, Rn will be calculated as a share of the global radiation. Default: 0.07}
\item{fao.par1}{Parameter to convert ETact from [W/m^2] to [mm/d] according to the approach by FAO. Default: 0.0864}
\item{fao.par2}{Parameter to convert ETact from [W/m^2] to [mm/d] according to the approach by FAO. Default: 0.408}
\item{out.unit}{Output unit for final ETact. Default is "mm" for [mm/d], other option is "watt" for [W/m^2]"]}
  \item{ta.conv}{If Tact is read in in K, but calculation should be in deg C. Default: "F"}
}
\details{
\itemize{
\item{Methods:}{    For further description of the different methods see papers below.}

\item{Tmin and Tmax:}{Estimation from NDVI-LST plot might not be very accurate, check automatic plot carefully and re-plot using the function ndvilst.plot for better results.}

\item{Rn:}{   Default share of 55 percent of global radiation based on Davies (1967).}

\item{G:}{    Default share of 7 percent of global radiation based on Crago and Brutsaert (1996).}

\item{Fao pars:}{   see Allen et al. (1998) for details.}


}
}

\references{
Allen, R., Pereira, L., Raes, D., and Smith, M. (1998). Crop evapotranspiration - guidelines for
computing crop water requirements. FAO Irrigation and drainage paper, 56.

Batra, N., Islam,S., Venturini,V.,Bisht,G.,Liang,L. (2006): Estimation and comparison of evapotranspiration from MODIS and AVHRR sensors for clear sky days over the Southern Great Plains, Remote Sensing of the Environment, 103, 1-15.

Crago, R. and Brutsaert, W. (1996). Daytime evaporation and the self-preservation of the evaporative
fraction and the Bowen ratio. Journal of Hydrology , 178:241-255.

Davies, J. (1967). A note on the relationship between net radiation and solar radiation. Quarterly
Journal if the Royal Meteorological Society, 93:109-115.

Jiang,L. and Islam,S. (1999): A methodology for estimation of surface evapotranspiration over large areas using remote sensing observations.

Eichinger, E., Parlange,M.,Stricker,H. (1996): On the concept of equilibrium evaporation and the value of the Pristley-Taylor coefficient, water resources research, 32, 161-164.

Priestley, C. and Taylor, R. (1972). On the assessment of surface heat flux and evaporation using
large-scale parameters. Monthly Weather Review, 100(2):81-92.

Stisen, S., Sandholt,I.,Norgaard,A., Fensholt,R., Hogh Jensen, K. (2008): Combining the triangle method with thermal inertia to estimate regional evapotranspiration - Applied to MSG-SEVIRI data in the Senegal River basin, Remote Sensing of the Environment, 112, 1242-1255.

Tang, R., Zhao-Liang, L., Bohui, T. (2010): An application of the Ts-VI triangle method with enhanced edges determination for evapotranspiration estimation from MODIS data in arid and semi-arid regions: Implementation and validation, Remote Sensing of Environment, 114, 540 - 551.

Wang, K., Li,Z., Cribb, M. (2006): Estimation of evaporative fraction from a combination of day and night land surface temperatures and NDVI: A new method to determine the Priestley-Taylor parameter, Remote Sensing of the Environment, 102, 293-305.

}


\note{
The function will automatically write a log-file in *.txt format to the output folder. This file is supposed to provide an overwiew on the approaches used in the calculation process and make the results more understandable and transparent. Hence, the file provides information on the inputs used, for gridded information, the filename is logged, else the numeric value of the variable (e.g. obs. temperature: 25 or "temperature_grid_xy.tif"). Furthermore, it sums up the methods used for Delta, Gamma and the pressure estimation, as well as the resulting t.min and t.max from the NDVi - LST - plot. The filenames for this plot and the final ETact grid are also logged, as well as the output unit of the final evapotranspiration. Finally, also the mean, min and max act. evapotranspiration - for a first plausability check - are logged. The last line of the file contains the date of the calculation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link{raster} \link{lsttm} \link{lstetm} \link{ndvi}}
}

\examples{
ndvi <- raster(ncol=30,nrow=20)
ndvi [] <- runif(30*20,min = -1,max = 1)
lst <- raster(ncol=30,nrow=20)
lst [] <- runif(30*20,min = 15,max = 35)
tact <- 25  	# fictive homogeneous air temperature of 25 deg C
rad <- 225 	# fictive homogeneous global radiation of 225 W/m^2
et <- triangle(ndvi, lst, tact, rad,ti.method="tang")
plot(et)

# now input dem and gridded temperature following the approach of Wang et al.
ndvi <- raster(ncol=30,nrow=20)
ndvi [] <- runif(30*20,min = -1,max = 1)
lst <- raster(ncol=30,nrow=20)
lst [] <- runif(30*20,min = 15,max = 35)
dem <- raster(ncol=30,nrow=20)
dem [] <- runif(30*20,min = 10,max = 250)
tact <- raster(ncol=30,nrow=20)
tact [] <- runif(30*20,min = 5,max = 25)
et <- triangle(ndvi, lst, tact, rad, delta.method = "wang",pressure.method = "wang")
plot(et)


# same but output now in W/m^2 instaead of mm
et <- triangle(ndvi, lst, tact, rad, delta.method = "wang",pressure.method = "wang",out.unit="watt" )
plot(et)


# loading data sets from file:
path <- system.file("extdata", package="TriangleMethod")
lst = lsttm(path)
ndvi = ndvi.file(path)
elev = system.file("extdata","dem.TIF", package="TriangleMethod")
elev= raster(elev)
tact <- 25  	# fictive homogeneous air temperature of 25 deg C
rad <- 225 	# fictive homogeneous global radiation of 225 W/m^2
et_with <- triangle(ndvi, lst, tact, rad,ti.method="tang",dem = elev)
et_without <- triangle(ndvi, lst, tact, rad,ti.method="tang")
split.screen(c(1,2))
screen(1)
plot(et_with)
screen(2)
plot(et_without)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
