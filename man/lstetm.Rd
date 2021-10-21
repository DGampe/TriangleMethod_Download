\name{lstetm}
\alias{lstetm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Brightness Temperature as proxy for Surface Temperature (Ts) from Landsat ETM+ scenes.
}
\description{
Brightness Temperature as proxy for Surface Temperature (Ts) will be calculated in degC or K from Landsat ETM+ thermal bands (as ETM+ has upper and lower thermal bands, the average radiation will be used.) here. Default parameters to calibrate the temperature at the sensor etc. are implemented for the ETM+ sensor. If TM scenes should be processed, see function lsttm. If LST should be calculated from other remote sensing data adjust parameters for calibration. NOTE: the parameters here match most of the scenes, however to be sure, check the corresponding calibration file for lmax and lmin! 
}
\usage{
lstetm(TEMPL,TEMPH, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{TEMPL}{
Lower thermal band (B61) for Landsat ETM+ or other sensors image, as raster data set - in memory, or as files (which will be rastered).
}
 \item{TEMPH}{
Upper thermal band (B62) for Landsat ETM+, as raster data set - in memory, or as files (which will be rastered).
}
 \item{write}{Specify, if the resulting Ts grid should be written to file.Default = F (T/F)}
 \item{outpath}{User defined outpath to write the final Ts file.}
 \item{outname}{User defined name for the resulting file.(Default name / format: 'ts.tif'}
 \item{unit}{Define, if Ts should be calculated in degC or K.}
 \item{LminLow,...}{Calibration parameter files for the sensor. Landsat ETM+ is the default parameterization, if images originate from a different sensor, these parameters should be adjusted.}
}
\details{
\itemize{
Calibration for ETM+ Low gain thermal band (band B61)
\item{LminLow:}{    0.0*10^6   (can be found in calibration parameter file for Landsat scene)}
\item{LmaxLow:}{    17.04*10^6  }

# Calibration for ETM+ High gain thermal band (band B62)
\item{LminHigh:}{   3.2*10^6}
\item{LmaxHigh:}{   12.65*10^6 # values here for Lmin+Lmax for ETM+ B5 Low Gain data}

\item{QcalMin:}{Minimum digital number, default 1}
\item{QcalMax:}{Maximum digital number, default 255}

\item{Calculation of radiation:}{

RadLow= ((LmaxLow-LminLow)/(QcalMax-QcalMin)) * (TEMPL  - QcalMin)+LminLow 

RadHigh= ((LmaxHigh-LminHigh)/(QcalMax-QcalMin)) * (TEMPH  - QcalMin)+LminHigh  

RadAvg = (RadLow+RadHigh)/2 # Average of the Low Gain and High Gain thermal bands 

see lsttm for additional information.
}

}

#\value{
#RasterLayer, .tif file
#}

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
