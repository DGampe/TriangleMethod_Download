\name{ndvilst.plot}
\alias{ndvilst.plot}
\title{
  Plot the NDVI-LST plot for estimation of Tmax and Tmin for further processing the triangle method.
}
\description{
  The Priestley Taylor coefficient PHI is scaled through the triangle shaped plot created from NDVI and LST. This function creates this type of plot, allows manual fitting and finally the more precise estimation of Tmin and Tmax for the estimation of actual evapotranspiration with the triangle method. The resulting plot combines the original NDVI-LST plot which can be found in Jiang & Islam (1999) with some modifications (e.g. adding the true dry edge) inspired by Stisen et al. (2008). However, the liniar relationship proposed by Jiang & Islam is kept here.
}
\usage{
  ndvilst.plot(NDVI, LST, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NDVI}{
    Can either be the filename of the raster containing the NDVI-information, or the allready raster variable in the workspace.
  }
  \item{LST}{
    Can either be the filename of the raster containing the LST/Ts-information, or the allready raster variable in the workspace.
  }
  \item{write}{specify, if the plot output should be writte to file. Default: T (alternative input F)}
  \item{Tmin}{Tmin from NDVI-LST plot. Lower end of vegetation cluster in plot. Default: c() (will be calculated)}
  \item{Tmax}{Tmax from NDVI-LST plot. Maximum temperature for vegetation cluster in plot. Default: c() (will be calculated)c()}
  \item{outpath}{Path to write the resulting plot file, if write =T. Default is current working directory.}
  \item{outname}{Name of the resulting .png file containg the plot. Default: c() (will be named automatically.)}
  \item{plot.unit}{Unit of the y-axis (LST). Default is depended on values of LST and will change automatically. Possible values are 'C' or 'K'.}
  \item{classify}{ Should NDVI be classified (classes of 0.05 -steps) ? Default = T}
  \item{plot.ylim}{Limit of the y - axis.}
  \item{plot.xlim}{Limit of the x - axis.}
  \item{plot.main}{Title for the plot. Default: 'NDVI-LST-plot'}
  \item{plot.xlab}{Label for x-axis. Default: 'NDVI'}
  \item{plot.ylab}{Label for y-axis. Default: LST deg C / LST K, dependend on scale of input LST.}
  \item{plot.wetcol}{Color for line wet edge. Default: 'blue'}
  \item{plot.drycol}{Color for observed dry edge. Default: 'red'}
  \item{plot.ndvi0}{Color for NDVI = 0 vertical line. Default: 'forestgreen'}
  \item{textwet.x}{X position for label of the wet edge.}
  \item{textwet.lab}{Label for wet edge line. Default: 'wet edge'}
  \item{textwet.col}{Color for wet edge line. Default: 'blue'}
  \item{textwet.font}{Font for wet edge line. Default: 2}
  \item{textwet.cex}{Siye of text. Default 0.9}
  \item{textobs.x}{X position for observed dry edge.}
  \item{textobs.lab}{see textwet. Default: 'observed dry edge'}
  \item{textobs.col}{see textwet. Default: 'red'}
  \item{textobs.font}{see textwet. Default: 2}
  \item{textobs.cex }{see textwet. Default: 0.9}
  \item{textdry.x}{X position for label of the true dry edge.}
  \item{textdry.lab}{see textwet. Default: 'true dry edge'}
  \item{textdry.col}{see textwet. Default: 'black'}
  \item{textdry.font}{see textwet. Default: 2}
  \item{textdry.cex}{see textwet. Default: 0.9}
  \item{textwater.x }{X position for label of the water-pixels}
  \item{textwater.lab }{see textwet. Default: 'water'}
  \item{textwater.col}{see textwet. Default: 'black'}
  \item{textwater.font}{see textwet. Default: 2}
  \item{textwater.cex}{see textwet. Default: 0.9}
  \item{textndvi.x}{X position for label of NDVI = 0.}
  \item{textndvi.lab}{see textwet.}
  \item{textndvi.col}{see textwet.}
  \item{textndvi.font}{see textwet.}
  \item{textndvi.cex}{see textwet.}
  \item{arrow.x0}{x 0 for arrow. see ?arrows for details}
  \item{arrow.y0}{y 0 for arrow. see ?arrows for details}
  \item{arrow.x1}{x 1 for arrow. see ?arrows for details}
  \item{arrow.angle}{Angle for arrow. see ?arrows for details}
  \item{arrow.code}{Default 2. see ?arrows for details}
  \item{arrow.length}{Default 0.15 see ?arrows for details}
  \item{arrow.lwd}{Lwd for arrow. Default 2. see ?arrows for details}
  \item{arrow.col}{Color for arrow. Default 'black'. see ?arrows for details}
  \item{pointmin.x}{X position for point for PHI-min.}
  \item{pointmin.col}{Color for point for PHI-min. Default: 'black'}
  \item{pointmin.type}{Type for point for PHI-min. Default: 'p'}
  \item{pointmin.lwd}{Lwd for point for PHI-min.}
  \item{pointmin.pch}{Symbol for point for PHI-min. Default: 16}
  \item{textmin.x}{X position for label of PHI-min.}
  \item{textmin.lab}{Label for label of PHI-min.}
  \item{textmin.col}{Color for lable for label of PHI-min.}
  \item{textmin.font}{Font for for label of PHI-min.}
  \item{textmin.cex}{Size of label of PHI-min.}
  \item{pointmax.x}{X position for point for PHI-max.}
  \item{pointmax.col}{see pointmin.}
  \item{pointmax.type }{see pointmin.}
  \item{pointmax.lwd}{see pointmin.}
  \item{pointmax.pch}{see pointmin.}
  \item{textmax.x}{see textmin}
  \item{textmax.lab}{see textmin}
  \item{textmax.col}{see textmin}
  \item{textmax.font}{see textmin}
  \item{textmax.cex}{see textmin}
  \item{textwet.y}{Y-position for textwet. Default will be determined by the function.}
  \item{textobs.y}{Y-position for textobs. Default will be determined by the function.}
  \item{textdry.y}{Y-position for textdry. Default will be determined by the function.}
  \item{textndvi.y}{Y-position for textndvi. Default will be determined by the function.}
  \item{textwater.y}{Y-position for textwater. Default will be determined by the function.}
  \item{arrow.y1}{Y-position for arrow y 1. Default will be determined by the function.}
  \item{pointmin.y }{Y-position for pointmin. Default will be determined by the function.}
  \item{textmin.y}{Y-position for textmin. Default will be determined by the function.}
  \item{pointmax.y }{Y-position for pointmax y. Default will be determined by the function.}
  \item{textmax.y}{Y-position for textmax. Default will be determined by the function.}
  \item{n.highest}{The nth highest temperature to automatically define the dry edge temperature. Default is n.highest = 2. Numeric.}
  \item{...}{Plot is based on smoothScatter function. Additional arguments can be passed.}
}

\references{
  Jiang and Islam (1999): A methodology for estimation of surface evapotranspiration over large areas using remote sensing observatios

  Batra et al. (2006): Estimation and comparison of evapotranspiration from MODIS and AVHRR sensors for clear sky days over the Southern Great Plains, Remote Sensing of the Environment, 103, 1-15

  Stisen et al. (2008): Combining the triangle method with thermal inertia to estimate regional evapotranspiration - Applied to MSG-SEVIRI data in the Senegal River basin, Remote Sensing of the Environment, 112, 1242-1255

  Wang et al. (2006): Estimation of evaporative fraction from a combination of day and night land surface temperatures and NDVI: A new method to determine the Priestley-Taylor parameter, Remote Sensing of the Environment, 102, 293-305
}


\seealso{
  \code{\link{smoothScatter} \link{ndvi} \link{lsttm} \link{lstetm}}
}
\examples{
  # ndvilst.plot('Your-NDVI-File-Here','Your-LST-File-Here')
  ndvi <- raster(ncol=30,nrow=20)
  ndvi [] <- runif(30*20,min = -1,max = 1)
  lst <- raster(ncol=30,nrow=20)
  lst [] <- runif(30*20,min = 15,max = 35)
  ndvilst.plot(ndvi,lst,write = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

