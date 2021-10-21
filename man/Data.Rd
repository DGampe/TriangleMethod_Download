\name{Data}
\alias{Data}
\docType{data}
\title{
External data sets.
}
\description{
Small subset of Landsat TM scene (Path/Row 174/38) with bands 3,4 and 6 as .tif file. 
}
\usage{
nir <- system.file("extdata", "B40.TIF", package="TriangleMethod")
}
\format{
  The format is:
*.tif
}
\details{
These grids can be loaded by using the code below and will be used in the examples within the package to examplify the calculation with external data (e.g. from a certain path containg the Landsat scenes).
}
\source{
The Landsat TM scene (Path/Row 174/38) was originally downloaded from http://glovis.usgs.gov/.
}
\references{
http://glovis.usgs.gov/
}

\examples{
nir <- system.file("extdata", "B40.TIF", package="TriangleMethod")
nir = raster(nir)
plot(nir)
}
\keyword{datasets data landsat}
