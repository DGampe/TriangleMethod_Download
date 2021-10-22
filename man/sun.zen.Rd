\name{sun.zen}
\alias{sun.zen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the solar zenit angle from raster.}
\description{
  This function will return the mean solar zenith angle for a given input raster with a given projection. Based on the zenith function (package GeoLight). Function will convert the coordinates to WGS 84 Lat Lon and use the mean longitude and latitude for the calculation of the solar zenith angle.
}
\usage{
  sun.zen(ras,time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ras}{
    Raster containing a spatial reference system. Function will convert the coordinates to WGS 84 Lat Lon.
  }
  \item{time}{
    Overpassing time / time when the scene was aquired. Time format 'POSIXct', e.g. "1990-01-01 10:00:00"
  }
  \item{crs.code}{
    CRS code for the respective image. Default: crs.code = CRS("+init=epsg:4326").
  }
}
\details{
  Most likely, an occuring error stems from a to large grid -> crop grid before the calculation of the solar zenit angle!
}
\value{
  scalar
}

\examples{

  red <- system.file("extdata", "B30.TIF", package="TriangleMethod")
  red <- raster(red)
  zen = sun.zen(red,"1980-11-14 10:00:00")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
