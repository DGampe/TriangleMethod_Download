\name{lsttm}
\alias{lsttm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Brightness Temperature as proxy for Surface Temperature (Ts) from Landsat TM scenes.
}
\description{
Brightness Temperature as proxy for Surface Temperature (Ts) will be calculated in degC or K from Landsat TM thermal band here. Default parameters to calibrate the temperature at the sensor etc. are implemented for the TM sensor. If ETM+ scenes should be processed, see function lstetm. If LST should be calculated from other remote sensing data adjust parameters for calibration. NOTE: the parameters here match most of the scenes, however to be sure, check the corresponding calibration file for lmax and lmin! }
\usage{
lsttm(TEMP, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TEMP}{
Thermal band for Landsat TM or other sensors image, as raster data set - in memory, or as files (which will be rastered).
}
 \item{write}{Specify, if the resulting Ts grid should be written to file.Default = F (T/F)}
 \item{outpath}{User defined outpath to write the final Ts file.}
 \item{outname}{User defined name for the resulting file.(Default name / format: 'ts.tif'}
 \item{unit}{Define, if Ts should be calculated in degC or K.}
 \item{LminLow,...}{Calibration parameter files for the sensor. Landsat TM is the default parameterization, if images originate from a different sensor, these parameters should be adjusted.}

}
\details{
\itemize{
Default calibration parameter values are:
 
 \item{lmin:}{  1.2378*10^6	(can be found in calibration parameter file for Landsat scene.) }
 \item{lmax:}{  15.3030*10^6	}

\item{Calculation of radiation:}{

RadLow= ((LmaxLow-LminLow)/(QcalMax-QcalMin)) * (TEMPL  - QcalMin)+LminLow 

RadHigh= ((LmaxHigh-LminHigh)/(QcalMax-QcalMin)) * (TEMPH  - QcalMin)+LminHigh  

RadAvg = (RadLow+RadHigh)/2 # Average of the Low Gain and High Gain thermal bands 

see lstetm for additional information.
}

}

#\value{
#RasterLayer, .tif file
#}
#\examples{
#
#}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
